import operator as op
from collections import abc
from itertools import starmap
from pathlib import Path

import lXtractor.chain as lxc
import lXtractor.ext.hmm as lxh
import pyhmmer
from easydict import EasyDict
from lXtractor.core.segment import Segment
from lXtractor.util import fetch_to_file
from loguru import logger
from more_itertools import one
from pyhmmer.plan7 import HMM, HMMFile, Domain, OptimizedProfile, Alignment
from toolz import valmap
from tqdm import tqdm

from tkp_finder.constants import PLANT_HMM_URL, RESOURCES_PATH, PLANT_HMMS_PATH

_HMM_BASE = Path(__file__).parent / "resources" / "hmm"
HMM_PATHS = EasyDict(
    {
        "pfam_domain": _HMM_BASE / "Domain.h3m",
        "pfam_family": _HMM_BASE / "Family.h3m",
        "pfam_motif": _HMM_BASE / "Motif.h3m",
        "plants_pkfam": _HMM_BASE / "PKFamilyPlants.h3m",
    }
)


def split_hmm(path: Path, get_path: abc.Callable[[HMM], Path], verbose: bool = False):
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


def make_map_name(hmm: HMM) -> str:
    return hmm.name.decode("utf-8").replace(" ", "_").replace("-", "_")


def load_profiles() -> EasyDict[str, dict[tuple[str, str], OptimizedProfile]]:
    def load(p: Path):
        with HMMFile(p) as hmms:
            try:
                hmms = hmms.optimized_profiles()
            except ValueError:
                pass
            return {(p.accession, p.name): p for p in hmms}

    return valmap(load, HMM_PATHS)


class Annotator:
    """
    A simple wrapper class helping to annotate a sequence by HMM models.
    """

    def __init__(self, cfg: EasyDict):
        self.cfg = cfg
        self.profiles: EasyDict = load_profiles()
        report = "; ".join(f"{k} {len(v)}" for k, v in self.profiles.items())
        logger.info(f"Loaded profiles. {report}. Config: {cfg}")

    def annotate(
        self,
        profile_keys: abc.Iterable[str],
        chains: abc.Sequence[lxc.ChainSequence],
        hmm2name: abc.Callable[[HMM], str] = make_map_name,
        **kwargs,
    ) -> abc.Iterator[lxc.ChainSequence]:
        """
        Run `hmmscan` against HMM models.

        :param hmms: A sequence of HMM models to query against.
        :return: Iterator over hits sorted by bitscore from highest to lowest.
        """

        def accept_domain(d: Domain, cov_hmm: float, cov_seq: float) -> bool:
            size = d.alignment.target_to - d.alignment.target_from

            keys = ("min_cov_hmm", "min_cov_seq", "min_score", "min_size", "max_pvalue")
            values = (cov_hmm, cov_seq, d.score, size, d.pvalue)
            ops = (op.ge, op.ge, op.ge, op.ge, op.lt)

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

        def calculate_coverage(
            aln: Alignment, hmm: HMM
        ) -> tuple[list[int | None], float, float]:
            num = [hmm_i for seq_i, hmm_i in lxh._enumerate_numbering(aln) if seq_i]

            # n = the number of valid HMM nodes covered
            # => 2 3 4 5 => 4
            n = sum(1 for x in num if x is not None)

            # SEQ coverage 4 / 6
            cov_seq = n / len(num)
            # HMM coverage 4 / M
            cov_hmm = n / hmm.M
            return num, cov_hmm, cov_seq

        logger.info(f"HMM annotation. Sequences={len(chains)}, profiles={profile_keys}")

        dseq = list(map(lxh.digitize_seq, chains))
        for pk in profile_keys:
            logger.info(f"Starting annotation for {pk}.")
            hmms_map = self.profiles[pk]
            hmms = list(hmms_map.values())
            top_hits_iter = pyhmmer.hmmer.hmmscan(
                dseq, hmms, cpus=self.cfg.get("num_proc", 0)
            )
            for top_hits, chain_seq in zip(top_hits_iter, chains):
                for hit in top_hits:
                    if hit.included:
                        dom = hit.best_domain
                        aln = dom.alignment
                        hmm = hmms_map[(hit.accession, hit.name)]
                        num, cov_hmm, cov_seq = calculate_coverage(aln, hmm)

                        if accept_domain(dom, cov_hmm, cov_seq):
                            map_name = hmm2name(hmm)
                            offset = chain_seq.start - 1
                            child = chain_seq.spawn_child(
                                aln.target_from + offset,
                                aln.target_to + offset,
                                map_name,
                                category=pk,
                                **kwargs,
                            )
                            child.add_seq(map_name, num)
                            child.meta[f"{map_name}_pvalue"] = dom.pvalue
                            child.meta[f"{map_name}_score"] = dom.score
                            child.meta[f"{map_name}_bias"] = dom.bias
                            child.meta[f"{map_name}_cov_seq"] = cov_seq
                            child.meta[f"{map_name}_cov_hmm"] = cov_hmm

                            yield child


if __name__ == "__main__":
    raise RuntimeError
