import logging
import operator as op
from collections import abc
from itertools import filterfalse, chain
from pathlib import Path
from tempfile import NamedTemporaryFile
from time import sleep

from lXtractor.core import ChainSequence, ChainList
from lXtractor.core.segment import Segment
from lXtractor.util.seq import write_fasta
from more_itertools import peekable, zip_equal, chunked
from tqdm.auto import tqdm

LOGGER = logging.getLogger(__name__)


class DeepTMHMM:
    def __init__(self, biolib_account_login: bool = False):
        import biolib
        # biolib.set_api_token('m4YbOukiwOilinBDYLlNnS6VIxuayGiq2DlejEdikVE')
        if biolib_account_login:
            try:
                print('Send KeyboardInterrupt to omit')
                biolib.login()
            except KeyboardInterrupt:
                pass

        biolib.biolib_logging.logger.setLevel(logging.CRITICAL)
        self.biolib = biolib
        self.interface = biolib.load('DTU/DeepTMHMM')

    def run(
            self, seqs: abc.Iterable[ChainSequence] | abc.Iterable[tuple[str, str]] | Path,
            blocking: bool = True, parse: bool = True
    ):
        # proper output annotation requires global import of biolib, causing loggin artifacts
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
                    res = self.run_cli(f.name, blocking=blocking)

            case Path():
                res = self.run_cli(seqs, blocking=blocking)

            case _:
                raise TypeError(f'Invalid input type {type(seqs)}')

        if not parse or not blocking:
            return res

        return self.parse_output_gff(res.get_output_file('/TMRs.gff3').get_data())

    def run_cli(self, inp_fasta: Path | str, blocking: bool = True):
        return self.interface.cli(args=f'--fasta {inp_fasta}', blocking=blocking)

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
            self, chains: abc.Iterable[ChainSequence], category: str = 'DeepTMHMM',
            chunk_size: int | None = None, **kwargs
    ) -> abc.Generator[ChainSequence]:
        if not isinstance(chains, ChainList):
            chains: ChainList[ChainSequence] = ChainList(chains)
        if chunk_size is None or chunk_size >= len(chains):
            LOGGER.info('Running a single job')
            results = self.run(chains)
        else:
            LOGGER.info('Submitting parallel jobs')
            jobs = []
            with self.biolib.Experiment('tkp-finder'):
                for i, chunk in enumerate(chunked(chains, chunk_size), start=1):
                    jobs.append((i, self.run(chunk, blocking=False, parse=False)))
            experiment = self.biolib.get_experiment('tkp-finder')
            LOGGER.info(f'Submitted {i} jobs: {experiment.show_jobs()}')
            outputs = []
            bar = tqdm(total=len(jobs), desc='Waiting for jobs')
            while jobs:
                for x in jobs:
                    i, j = x
                    if j.is_finished():
                        try:
                            outputs.append(
                                (i, self.parse_output_gff(j.get_output_file('/TMRs.gff3').get_data()))
                            )
                        except Exception as e:
                            LOGGER.exception(e)
                            LOGGER.warning(f'Failed to complete a job {j.id} due to {e}')
                        jobs.remove(x)
                        bar.update(1)
                if jobs:
                    sleep(2)

            bar.close()
            results = list(chain.from_iterable(
                map(op.itemgetter(1), sorted(outputs, key=op.itemgetter(0)))))

        for (c_id, segments), c in zip_equal(results, chains):
            assert c_id == c.id, "IDs and order are preserved"
            for s in segments:
                yield c.spawn_child(s.start, s.end, f'{category}_{s.name}', **kwargs)


if __name__ == '__main__':
    raise RuntimeError
