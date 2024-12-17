from __future__ import annotations

import gzip
import re
import typing as t
from abc import abstractmethod
from collections import abc
from io import StringIO
from pathlib import Path

import lXtractor.chain as lxc
from lXtractor.util import read_fasta


def read_fasta_gz(p: Path):
    with gzip.open(p, "rb") as f:
        text = f.read().decode("utf-8")
        yield from read_fasta(StringIO(text), strip_id=False)


def select_transcripts_by_meta(
    chains: abc.Iterable[lxc.ChainSequence],
) -> abc.Iterable[lxc.ChainSequence]:
    def get_transcript(c: lxc.ChainSequence) -> str:
        return c.meta["transcript"]

    chains = lxc.ChainList(chains)
    for g, gg in chains.sort(get_transcript).groupby(get_transcript):
        sel_transcript = max(gg, key=len)
        sel_transcript.meta["num_transcripts"] = len(gg)
        yield sel_transcript


class Proteome:
    def __init__(self, fasta_path: Path, category: str):
        self.fasta_path = fasta_path
        self.file_meta: dict[str, t.Any] = self.parse_filename(fasta_path.name)
        self.file_meta["category"] = category
        self.file_meta["filename"] = fasta_path.name

    @abstractmethod
    def init_chain_sequence(self, fasta_item: tuple[str, str]) -> lxc.ChainSequence:
        ...

    @staticmethod
    @abstractmethod
    def parse_filename(filename: str) -> dict[str, t.Any]:
        ...

    @staticmethod
    def select_transcripts(
        chains: abc.Iterable[lxc.ChainSequence],
    ) -> abc.Iterator[lxc.ChainSequence]:
        return iter(chains)

    def preprocess(self) -> abc.Iterator[lxc.ChainSequence]:
        if self.fasta_path.name.endswith(".gz"):
            fasta_iter = read_fasta_gz(self.fasta_path)
        else:
            fasta_iter = read_fasta(self.fasta_path, strip_id=False)
        chains = map(self.init_chain_sequence, fasta_iter)
        chains = self.select_transcripts(chains)
        yield from chains


class EnsemblProteome(Proteome):
    @staticmethod
    def parse_filename(filename: str) -> dict[str, t.Any]:
        organism = filename.split(".")[0]
        return {"organism": organism}

    @staticmethod
    def select_transcripts(
        chains: abc.Iterable[lxc.ChainSequence],
    ) -> abc.Iterator[lxc.ChainSequence]:
        yield from select_transcripts_by_meta(chains)

    def init_chain_sequence(self, fasta_item: tuple[str, str]) -> lxc.ChainSequence:
        header, seq = fasta_item

        try:
            gene = re.findall(r"gene:([^\s]+)", header)[0]
        except IndexError:
            gene = None

        try:
            transcript = re.findall(r"transcript:([^\s]+)", header)[0]
        except IndexError:
            transcript = None

        cs = lxc.ChainSequence.from_string(seq, name=header.split()[0])
        cs.meta.update(self.file_meta)
        cs.meta.update(dict(gene=gene, transcript=transcript, header=header))
        return cs


class NCBIProteome(Proteome):
    @staticmethod
    def parse_filename(filename: str) -> dict[str, t.Any]:
        organism = filename.split(".")[0]
        return {"organism": organism}

    def init_chain_sequence(self, fasta_item: tuple[str, str]) -> lxc.ChainSequence:
        header, seq = fasta_item
        cs = lxc.ChainSequence.from_string(seq, name=header.split()[0])
        cs.meta.update(self.file_meta)
        cs.meta.update(dict(header=header))
        return cs

    def preprocess(self):
        chains = super().preprocess()
        chains = filter(lambda c: "partial" in c.meta["header"], chains)
        yield from chains


class Phytozome_proteome(Proteome):
    @staticmethod
    def parse_filename(filename: str) -> dict[str, t.Any]:
        organism = filename.split("_")[0]
        return {"organism": organism}

    @staticmethod
    def select_transcripts(
        chains: abc.Iterable[lxc.ChainSequence],
    ) -> abc.Iterator[lxc.ChainSequence]:
        yield from select_transcripts_by_meta(chains)

    def init_chain_sequence(self, fasta_item: tuple[str, str]) -> lxc.ChainSequence:
        header, seq = fasta_item
        seq_id = header.split()[0]
        seq = seq.removesuffix("*")

        meta = {}

        transcript = None
        if "transcript" in header:
            transcript = re.findall(r"transcript=([^\s]+)", header)[0]
        if transcript is None:
            transcript = seq_id

        meta["transcript"] = transcript

        if "locus" in header:
            locus = re.findall(r"locus=([^\s]+)", header)[0]
            meta["locus"] = locus

        if "pacid" in header:
            pacid = re.findall(r"pacid=([^\s]+)", header)[0]
            meta["pacid"] = pacid

        meta["header"] = header
        cs = lxc.ChainSequence.from_string(seq, name=seq_id)
        cs.meta.update(self.file_meta)
        cs.meta.update(meta)
        return cs


if __name__ == "__main__":
    raise RuntimeError
