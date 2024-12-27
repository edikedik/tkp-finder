from __future__ import annotations

import operator as op
import typing as t
from collections import abc
from itertools import pairwise, chain, repeat, islice, starmap
from pathlib import Path

import lXtractor.chain as lxc
import lXtractor.ext.hmm as lxh
from easydict import EasyDict
from lXtractor.ext import PyHMMer
from lXtractor.util import read_fasta
from loguru import logger
from more_itertools import split_when, mark_ends, one
from more_itertools.more import always_iterable

from tkp_finder.constants import MSA_PATH, PSKD_NAME

_T = t.TypeVar("_T")


def get_pskd_hmm(msa: Path = MSA_PATH, name: str = PSKD_NAME):
    return lxh.PyHMMer.from_msa(read_fasta(msa), name, "amino")


def take_n_filtered(
    pred: abc.Callable[[_T], bool], it: abc.Iterable[_T], n
) -> abc.Iterator[_T]:
    return islice(filter(pred, it), n)


def join_non_overlapping_hits(
    parent: lxc.ChainSequence, hits: abc.Sequence[lxc.ChainSequence], map_name: str
) -> lxc.ChainSequence:
    """
    Join hits within a parent sequence that don't overlap in their `map_name`
    numbering.

    :param parent: Parent chain sequence.
    :param hits: Hits to join produced by the same hmm annotation act and
        existing within ``parent.children``.
    :param map_name: Map name pointing to HMM numbering within each hit.
    """
    seg_start, seg_end = hits[0].start, hits[-1].end
    cats = hits[0].categories
    cat = cats[0] if cats else None
    hmm_numbering = always_iterable(None)
    for is_first, is_last, (h1, h2) in mark_ends(pairwise(hits)):
        h1_num, h2_num = h1[map_name], h2[map_name]
        insert_size = (h2.start - 1) - (h1.end + 1) + 1
        insert = repeat(None, insert_size)
        if not is_first and is_last:
            append = chain(insert, h2_num)
        else:
            append = chain(h1_num, insert, h2_num)
        hmm_numbering = chain(hmm_numbering, append)

    loc_h1 = one(i for i, c in enumerate(parent.children) if c.id == hits[0].id)
    child = parent.spawn_child(
        seg_start,
        seg_end,
        map_name,
        category=cat,
        keep=False,
    )
    child[map_name] = list(hmm_numbering)
    parent.children.insert(loc_h1, child)
    child.parent = parent
    hit_ids = tuple(h.id for h in hits)
    parent.children = parent.children.filter(lambda x: x.id not in hit_ids)
    return child


def join_overlapping_hits(
    parent: lxc.ChainSequence,
    hits: abc.Sequence[lxc.ChainSequence],
    map_name: str,
    pyhmm: PyHMMer,
    boundary_offset: int = 5,
) -> lxc.ChainSequence | None:
    """
    Join hits within a parent sequence that overlap in their `map_name`
    numberings. In contrast to :func:`join_non_overlapping_hits`, here,
    hits are assumed to have overlaps and are therefore can't be merged
    directly. In this case, hit sequences will be joined and the resulting
    sequence annotated anew to infer correct non-overlapping HMM numbering.

    :param parent: Parent chain sequence.
    :param hits: Hits to join produced by the same hmm annotation act and
        existing within ``parent.children``.
    :param map_name: Map name pointing to HMM numbering within each hit.
    :param pyhmm: PyHMMer object holding the reference HMM.
    """

    def iter_hit_seqs():
        parent_seq = parent.seq1
        for is_first, is_last, h in mark_ends(hits):
            start = h.start - boundary_offset if is_first else h.start
            end = h.end + boundary_offset if is_last else h.end
            yield parent_seq[start:end]

    merged_seq = "".join(iter_hit_seqs())
    merged_cs = lxc.ChainSequence.from_string(merged_seq)
    try:
        hit = one(
            pyhmm.annotate([merged_cs], new_map_name=map_name, keep=False),
            too_long=RuntimeError("More than one hit"),
        )
    except RuntimeError:
        return None

    # creates "S" parent numbering within the hit
    numbering = hit.map_numbering(parent)
    start, end = numbering[0], numbering[-1]

    # Infer category
    cats = hits[0].categories
    cat = cats[0] if cats else None

    # Create a new merged hit child subsequence
    new_child = parent.spawn_child(start, end, map_name, cat, keep=False)
    hit.relate(new_child, map_name, "i", "S")
    new_child.meta.update(filter(lambda it: map_name in it[0], hit.meta.items()))
    new_child.meta["merged"] = True
    new_child.meta["merge_kind"] = "overlapping"

    # Handle parent-child relationship and remove old hits
    loc_h1 = one(i for i, c in enumerate(parent.children) if c.id == hits[0].id)
    parent.children.insert(loc_h1, new_child)
    new_child.parent = parent
    hit_ids = tuple(h.id for h in hits)
    parent.children = parent.children.filter(lambda x: x.id not in hit_ids)
    return new_child


class PKDiscoverer:
    """
    A class for PK domain discovery.
    """

    def __init__(
        self,
        cfg: EasyDict,
        pyhmm: lxh.PyHMMer | lxh._HmmInpT | None,
        assign_name: str = PSKD_NAME,
    ):
        """
        :param pyhmm: PyHMMer object with HMM of a PK domain.
        :param assign_name: Internal domain name.
        """
        self.cfg = cfg
        if pyhmm is None:
            logger.info("Initializing default PSKD HMM.")
            pyhmm = get_pskd_hmm(name=assign_name)
        if not isinstance(pyhmm, lxh.PyHMMer):
            pyhmm = lxh.PyHMMer(pyhmm)
        self.pyhmm = pyhmm
        self.assign_name = assign_name

    def annotate(
        self, chains: abc.Iterable[lxc.ChainSequence], strip_hmm_ord: bool = True
    ) -> abc.Iterator[lxc.ChainSequence]:
        """
        Run initial domain annotation. A basic wrapper around :meth:`PyHMMer.annotate`.

        :param chains: Candidate chain sequences.
        :param strip_hmm_ord: Strip hit numeration if multiple hits are found
            in the same sequence. Eg, ``"PSKD_1" -> "PSKD"``.
        :return: An iterator over chain sequence hits.
        """
        hits = self.pyhmm.annotate(
            chains, new_map_name=self.assign_name, **self.cfg.annotation
        )
        if strip_hmm_ord:
            for hit in hits:
                hit.name = hit.name.split("_")[0]
                yield hit
        else:
            yield from hits

    def join_hits(
        self,
        parent_chain: lxc.ChainSequence,
    ) -> None:
        """
        Oftentimes large inserts cause fragmentary hits. This method attempts
        to join such fragments. There are two round of joining:

            #. For directly mergeable hits having continuous HMM numbering.
            #. For indirectly mergeable hits having HMM numbering overlaps.

        For the former, hits are interspersed with inserts and concatenated
        into a single hit. For the latter, an additional annotation round
        is executed on concatenated hit sequences (without inserts) and HMM
        numbering is then transferred to a hit variant that has the original
        inserts.

        :param parent_chain: Chain whose hits should be joined. Hits are found
            via :attr:`assign_name`.
        :param max_insert: Maximum allowed insert size between fragments.
        :param min_size: Minimum fragment size to consider joining.
        :param max_overlap: Maximum overlap size for a pair of hits to qualify
            as indirectly mergeable.
        :return: Modifies the provided parent sequence and returns nothing.
        """

        def get_ends(
            c1: lxc.ChainSequence, c2: lxc.ChainSequence
        ) -> tuple[int, int, int] | None:
            c1_end_item = c1.get_closest(
                self.assign_name, self.pyhmm.hmm.M, reverse=True
            )
            c2_start_item = c2.get_closest(self.assign_name, 1)

            if c1_end_item is None or c2_start_item is None:
                return None

            c1_end, c2_start = map(
                lambda x: x._asdict()[self.assign_name], [c1_end_item, c2_start_item]
            )
            insert_size = c2_start_item.i - c1_end_item.i + 1
            return c1_end, c2_start, insert_size

        def not_directly_mergeable(
            c1: lxc.ChainSequence, c2: lxc.ChainSequence
        ) -> bool:
            ends = get_ends(c1, c2)
            if ends is None:
                # throw a warning?
                return True
            c1_end, c2_start, insert_size = ends
            return (
                len(c1) < self.cfg.min_size
                or len(c2) < self.cfg.min_size
                or c1_end >= c2_start
                or insert_size > self.cfg.max_insert
            )

        def not_indirectly_mergeable(
            c1: lxc.ChainSequence, c2: lxc.ChainSequence
        ) -> bool:
            ends = get_ends(c1, c2)
            if ends is None:
                # throw a warning?
                return True
            insert_size = ends[-1]
            items_c1_end = take_n_filtered(
                lambda x: x is not None, c1[self.assign_name], self.cfg.max_overlap
            )
            items_c2_start = take_n_filtered(
                lambda x: x is not None,
                c2[self.assign_name][::-1],
                self.cfg.max_overlap,
            )
            overlap = set(items_c1_end) & set(items_c2_start)

            return (
                len(c1) < self.cfg.min_size
                or len(c2) < self.cfg.min_size
                or insert_size > self.cfg.max_insert
                or len(overlap) > self.cfg.max_overlap
            )

        def populate_meta(
            merged_hit: lxc.ChainSequence, hits: abc.Sequence[lxc.ChainSequence]
        ):
            num_nodes_covered = sum(x is not None for x in merged_hit[self.assign_name])
            cov_hmm_name = f"{self.assign_name}_cov_hmm"
            cov_seq_name = f"{self.assign_name}_cov_seq"
            score_name = f"{self.assign_name}_score"
            bias_name = f"{self.assign_name}_bias"
            pval_name = f"{self.assign_name}_pvalue"
            score_sum = sum(h.meta[score_name] for h in hits)
            pval_max = max(h.meta[pval_name] for h in hits)
            meta_upd = {
                cov_hmm_name: num_nodes_covered / self.pyhmm.hmm.M,
                cov_seq_name: num_nodes_covered / len(merged_hit),
                score_name: score_sum,
                pval_name: pval_max,
                bias_name: None,
                "merged": True,
                "merge_kind": "non-overlapping",
            }
            merged_hit.meta.update(meta_upd)

        # First passage for directly mergeable hits
        hits = parent_chain.children.filter(lambda x: x.name == self.assign_name)
        splits = filter(lambda x: len(x) > 1, split_when(hits, not_directly_mergeable))
        splits = list(splits)
        for split in splits:
            logger.debug(f"Merging {split} for {parent_chain}")
            merged = join_non_overlapping_hits(parent_chain, split, self.assign_name)
            populate_meta(merged, split)

        # Second passage for indirectly mergeable hits
        hits = parent_chain.children.filter(lambda x: x.name == self.assign_name)
        splits = filter(
            lambda x: len(x) > 1, split_when(hits, not_indirectly_mergeable)
        )
        for split in splits:
            join_overlapping_hits(parent_chain, split, self.assign_name, self.pyhmm)

    def post_join_filter(self, chains: abc.Sequence[lxc.ChainSequence]):
        def get_meta(c: lxc.ChainSequence, postfix: str):
            key = f"{self.assign_name}_{postfix}"
            try:
                return float(c.meta[key])
            except KeyError as e:
                raise KeyError(f"No key {key} in {c}'s meta {c.meta}")

        def accept_hit(hit: lxc.ChainSequence):
            values = (
                get_meta(hit, "cov_hmm"),
                get_meta(hit, "cov_seq"),
                get_meta(hit, "score"),
                len(hit),
                get_meta(hit, "pvalue"),
            )
            criteria = starmap(
                lambda k, v, o: self.cfg.get(k, None) is None or o(v, self.cfg[k]),
                zip(keys, values, ops),
            )
            accepted = all(criteria)
            criteria_fmt = ",".join(
                "{}={:.2f}".format(k, v) for k, v in zip(keys, values)
            )
            logger.debug(f"Criteria={criteria_fmt}. Accepted={accepted}.")
            return accepted

        def iter_accepted():
            for c in chains:
                accepted_hits = c.children[self.assign_name].filter(accept_hit)
                if accepted_hits:
                    c.children = c.children.filter(
                        lambda x: x.name != self.assign_name or x in accepted_hits
                    )
                    yield c

        keys = ("min_cov_hmm", "min_cov_seq", "min_score", "min_size", "max_pvalue")
        ops = (op.ge, op.ge, op.ge, op.ge, op.lt)

        return lxc.ChainList(iter_accepted())


if __name__ == "__main__":
    raise RuntimeError
