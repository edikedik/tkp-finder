from collections import abc
from pathlib import Path

import lXtractor.chain as lxc
import lXtractor.ext.hmm as lxh
from lXtractor.core.segment import resolve_overlaps, Segment
from lXtractor.util import fetch_to_file
from loguru import logger
from more_itertools import consume, one
from pyhmmer.plan7 import HMM, HMMFile
from tqdm import tqdm

from tkp_finder.constants import PLANT_HMM_URL, RESOURCES_PATH, PLANT_HMMS_PATH


def split_hmm(
    path: Path, get_path: abc.Callable[[HMM], Path], verbose: bool = False
):
    with HMMFile(path) as hmms:
        if verbose:
            hmms = tqdm(hmms, desc="Splitting HMM")
        for hmm in hmms:
            hmm_path = get_path(hmm)
            hmm_path.parent.mkdir(exist_ok=True, parents=True)
            with hmm_path.open("wb") as f:
                hmm.write(f)


def load_pfam() -> lxh.Pfam:
    pfam = lxh.Pfam()
    try:
        pfam.read(hmm=False)
    except FileNotFoundError:
        logger.warning(
            "Pfam was not initialized. "
            "Will fetch and parse Pfam from scratch. "
            "This might take a while..."
        )
        pfam.fetch()
        pfam.parse()
    return pfam


def fetch_plant_hmms():
    path_plants = fetch_to_file(PLANT_HMM_URL, root_dir=RESOURCES_PATH)
    return path_plants


def split_plant_hmms(path: Path, verbose: bool = False):
    def make_path(hmm: HMM):
        return (
            PLANT_HMMS_PATH
            / f"{hmm.name.decode('utf-8').replace(' ', '_').replace('-', '_')}.hmm"
        )

    split_hmm(path, make_path, verbose)


def get_score(cs: Segment) -> float:
    return one(v for k, v in cs.meta.items() if k.endswith("_score"))


class Annotator:
    """
    A simple wrapper class helping to annotate a sequence by HMM models.
    """

    def __init__(self, chain_seq: lxc.ChainSequence):
        self.chain_seq = chain_seq

    def annotate(self, hmms: abc.Sequence[lxh.PyHMMer], **kwargs):
        for hmm in hmms:
            consume(hmm.annotate([self.chain_seq], **kwargs))

    def annotate_pfam(self, categories: tuple[str, ...], **kwargs):
        pfam = load_pfam()
        df = pfam.read(categories=categories, hmm=True)
        for cat, gg in df.groupby["category"]:
            self.annotate(gg["PyHMMer"], category=cat, **kwargs)

    def resolve_overlaps(self, value_fn: abc.Callable[[Segment], float] = get_score):
        for g, gg in self.chain_seq.children.groupby(lambda x: x.categories):
            gg_no_ov = lxc.ChainList(
                resolve_overlaps(gg, value_fn=value_fn, max_it=1000)
            )
            if len(gg_no_ov) == len(gg):
                continue
            other_children_ids = self.chain_seq.children.ids
            self.chain_seq.children = self.chain_seq.children.filter(
                lambda x: x.id in other_children_ids or x.id in gg_no_ov
            )


if __name__ == "__main__":
    raise RuntimeError
