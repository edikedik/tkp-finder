from collections import abc
from itertools import filterfalse
from pathlib import Path
from tempfile import NamedTemporaryFile

import biolib
from lXtractor.core import ChainList, ChainSequence
from lXtractor.core.segment import Segment
from lXtractor.util.seq import write_fasta
from more_itertools import peekable, split_at


class DeepTMHMM:
    def __init__(self):
        self.interface = biolib.load('DTU/DeepTMHMM')

    def run(
            self, seqs: ChainList | abc.Iterable[ChainSequence] | abc.Iterable[tuple[str, str]] | Path
    ) ->  list[tuple[str, list[Segment]]]:
        match seqs:
            case abc.Iterable():
                peek = peekable(seqs)
                first = peek.peek()
                match first:
                    case [str(), str()]:
                        pass
                    case ChainSequence():
                        seqs = [(c.id, c.seq1) for c in seqs]
                    case other:
                        raise ValueError(
                            f'Expected to find a (header, seq) pair or `ChainSequence`, '
                            f'but found {other}')

                with NamedTemporaryFile('w') as f:
                    write_fasta(seqs, f)
                    f.seek(0)
                    res = self.run_cli(f.name)

            case Path():
                res = self.run_cli(seqs)

            case _:
                raise TypeError(f'Invalid input type {type(seqs)}')

        return list(self.parse_output_gff(res.get_output_file('/TMRs.gff3')))

    def run_cli(self, inp_fasta: Path | str):
        return self.interface.cli(args=f'--fasta {inp_fasta}')

    @staticmethod
    def parse_output_gff(inp: Path | str | bytes) -> abc.Iterator[tuple[str, list[Segment]]]:
        def parse_chunk(xs: abc.Iterable[str]) -> tuple[str, list[Segment]]:
            xs = map(lambda x: x.split('\t'), xs)
            peek = peekable(xs)
            first = peek.peek(None)
            if first is None:
                raise ValueError('Empty GFF3 output')
            chain_id = first[0]
            return chain_id, [Segment(int(x[2]), int(x[3]), x[1]) for x in xs]

        if isinstance(inp, Path):
            inp = inp.read_text()
        if isinstance(inp, bytes):
            inp = inp.decode('utf-8')

        lines = filterfalse(
            lambda x: not x or x.startswith('#'), inp.split('\n'))

        chunks = split_at(lines, lambda x: x.startswith('\\'))
        return map(parse_chunk, chunks)
