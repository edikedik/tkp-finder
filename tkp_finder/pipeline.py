from __future__ import annotations

import multiprocessing
from collections import abc

import lXtractor.chain as lxc
from easydict import EasyDict
from lXtractor.core import Segment
from lXtractor.core.exceptions import LengthMismatch
from lXtractor.core.segment import resolve_overlaps
from lXtractor.util import apply
from loguru import logger

from tkp_finder.annotation import Annotator, get_score
from tkp_finder.constants import PSKD_NAME, HIT_SEP
from tkp_finder.discovery import PKDiscoverer


class Pipeline:
    def __init__(self, cfg: EasyDict):
        self.cfg = cfg
        self.discoverer = PKDiscoverer(cfg.discovery, cfg.get("PKD_HMM", None))

        match self.cfg.get("num_proc"):
            case None:
                self.cfg.num_proc = 1
            case 0:
                self.cfg.num_proc = multiprocessing.cpu_count()

        self.annotator = Annotator(cfg.annotation)
        self.plant_hmms = None

    def _apply_wrapper(self, fn, it, task):
        yield from apply(
            fn, it, True, task, self.cfg.num_proc, use_joblib=self.cfg.use_joblib
        )

    def annotate(self, chains: abc.Sequence[lxc.ChainSequence], is_plant: bool):
        keys = ["pfam_domain", "pfam_family", "pfam_motif"]
        if is_plant:
            keys.append("plants_pkfam")
        num_hits = sum(1 for _ in self.annotator.annotate(keys, chains))
        logger.info(f"Total hits={num_hits}.")

    def discover_one(self, cs: lxc.ChainSequence) -> lxc.ChainSequence:
        try:
            hits = list(self.discoverer.annotate([cs]))
            if len(hits) > 0:
                self.discoverer.join_hits(cs)
        except LengthMismatch as e:
            logger.error(
                f"Sequence {cs}. Length mismatch: {e}. Inspect and debug the "
                f"exception below. Skipping domain annotation."
            )
            logger.exception(e)
            cs.children = lxc.ChainList([])
        return cs

    @staticmethod
    def resolve_overlaps(
        cs: lxc.ChainSequence,
        value_fn: abc.Callable[[Segment], float] = get_score,
    ):
        groups = dict(cs.children.groupby(lambda x: tuple(x.categories)))
        for g, gg in groups.items():
            gg_no_ov = lxc.ChainList(
                resolve_overlaps(gg, value_fn=value_fn, max_it=1000)
            )
            if len(gg_no_ov) == len(gg):
                continue
            logger.debug(
                f"Group={g}, Overlapping_seqs={gg}, Non-overlapping: {gg_no_ov}"
            )
            logger.debug(f"Children before: {cs.children}")
            cs.children = cs.children.filter(
                lambda x: tuple(x.categories) != g or x in gg_no_ov
            ).sort(key=lambda x: x.start)
            logger.debug(f"Children after: {cs.children}")
        return cs

    def fmt_hit_composition(self, cs: lxc.ChainSequence) -> lxc.ChainSequence:
        def fmt(hit: lxc.ChainSequence):
            return f"<{hit.start}-{hit.name}-{hit.end}>"

        for cat in cs.children.categories:
            hits = cs.children[cat].sort(key=lambda x: x.start)
            cs.meta[f"{cat}_hits"] = HIT_SEP.join(map(fmt, hits))
        cs.meta["hits"] = HIT_SEP.join(
            map(fmt, cs.children.sort(key=lambda x: x.start))
        )

        return cs

    def run(
        self, chains: abc.Iterable[lxc.ChainSequence], is_plant: bool
    ) -> lxc.ChainList[lxc.ChainSequence]:
        chains = map(self.discover_one, chains)
        chains = lxc.ChainList(filter(lambda x: len(x.children[PSKD_NAME]) > 0, chains))
        logger.info(f"Discovered {len(chains)} chains with PKD hits.")
        chains = self.discoverer.post_join_filter(chains)
        logger.info(f"After post-join filter: {len(chains)} chains with PKD hits.")
        if len(chains) > 0:
            self.annotate(chains, is_plant)
        chains = chains.apply(
            self.resolve_overlaps,
            verbose=True,
            desc="Overlaps resolved",
            num_proc=self.cfg.num_proc,
        )
        chains = chains.apply(self.fmt_hit_composition)
        return chains


if __name__ == "__main__":
    raise RuntimeError
