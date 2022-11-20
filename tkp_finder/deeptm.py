import logging
from collections import abc
from itertools import filterfalse
from pathlib import Path
from tempfile import NamedTemporaryFile

from lXtractor.core import ChainSequence, ChainList
from lXtractor.core.segment import Segment
from lXtractor.util.seq import write_fasta
from more_itertools import peekable, zip_equal


class DeepTMHMM:
    def __init__(self):
        import biolib
        biolib.biolib_logging.logger.setLevel(logging.CRITICAL)
        self.interface = biolib.load('DTU/DeepTMHMM')

    def run(
            self, seqs: abc.Iterable[ChainSequence] | abc.Iterable[tuple[str, str]] | Path
    ) -> abc.Iterator[tuple[str, list[Segment]]]:
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

        return self.parse_output_gff(res.get_output_file('/TMRs.gff3').get_data())

    def run_cli(self, inp_fasta: Path | str):
        return self.interface.cli(args=f'--fasta {inp_fasta}')

    @staticmethod
    def parse_output_gff(inp: Path | str | bytes) -> abc.Iterator[tuple[str, list[Segment]]]:
        def parse_chunk(c: str) -> tuple[str, list[Segment]]:
            lines = map(
                lambda x: x.split('\t'),
                filterfalse(lambda x: not x or x.startswith('#'), c.split('\n')))
            lines = peekable(lines)
            fst = lines.peek(None)
            assert fst, 'Non-empty chunk of data'
            chain_id = fst[0]
            return chain_id, [Segment(int(x[2]), int(x[3]), x[1]) for x in lines]

        if isinstance(inp, Path):
            inp = inp.read_text()
        if isinstance(inp, bytes):
            inp = inp.decode('utf-8')

        chunks = inp.split('//')
        return map(parse_chunk, chunks)

    def annotate(
            self, chains: abc.Iterable[ChainSequence], category: str = 'DeepTMHMM', **kwargs
    ) -> abc.Generator[ChainSequence]:
        if not isinstance(chains, ChainList):
            chains: ChainList[ChainSequence] = ChainList(chains)
        results = self.run(chains)
        for (c_id, segments), c in zip_equal(results, chains):
            assert c_id == c.id, "IDs and order are preserved"
            for s in segments:
                yield c.spawn_child(s.start, s.end, f'{category}_{s.name}', **kwargs)


if __name__ == '__main__':
    raise RuntimeError
