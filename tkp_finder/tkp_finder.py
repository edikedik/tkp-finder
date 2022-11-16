import gzip
import logging
import operator as op
import shutil
from collections import abc
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from warnings import warn

import click
import pandas as pd
from lXtractor.core import ChainSequence, ChainList, ChainInitializer, ChainIO
from lXtractor.core.segment import resolve_overlaps
from lXtractor.ext import PyHMMer
from lXtractor.util.io import download_to_file, get_files, get_dirs
from lXtractor.util.seq import read_fasta
from lXtractor.variables.base import SequenceVariable
from lXtractor.variables.manager import Manager
from lXtractor.variables.sequential import SeqEl
from more_itertools import consume, unique_everseen, split_at, zip_equal
from pyhmmer.plan7 import HMMFile
from toolz import keyfilter, valfilter, compose_left, pipe, curry
from tqdm.auto import tqdm

PFAM_A_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_DAT_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"
PFAM_A_NAME = 'Pfam-A.hmm'
PFAM_DAT_NAME = 'Pfam-A.hmm.dat'
PFAM_PK_NAME = 'PF00069'
PK_NAME = 'PK'
PPK_NAME = 'PPK'
UNK_HMM = 'unknown'
VARIABLES = (
    SeqEl(30),  # Beta-3 Lys
    SeqEl(48),  # aC Glu
    SeqEl(121), SeqEl(122), SeqEl(123),  # HRD
    SeqEl(141), SeqEl(142), SeqEl(143),  # DFG
)
MOTIF = 'KEXXDDXX'

LOGGER = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
LOGGER.addHandler(handler)


@click.group(
    'tkp_finder',
    context_settings=dict(
        help_option_names=['-h', '--help'],
        ignore_unknown_options=True
    ),
    invoke_without_command=True)
def tkp_finder(): ...


@tkp_finder.command('setup')
@click.option(
    '-H', '--hmm_dir', required=True,
    type=click.Path(dir_okay=True, file_okay=False, writable=True),
    help='Path to a directory to store hmm-related data.'
)
@click.option(
    '-d', '--download', is_flag=True, default=False, show_default=True,
    help='If True, download the Pfam data from interpro.'
)
@click.option(
    '-q', '--quiet', is_flag=True, default=False, show_default=True,
    help='Disable verbose output.'
)
@click.option(
    '--path_pfam_a', type=click.Path(dir_okay=False, file_okay=True, exists=True),
    help='A path to downloaded Pfam-A HMM profiles. By default, if `download` is ``False``,'
         'will try to find it within the `hmm_dir`.'
)
@click.option(
    '--path_pfam_dat', type=click.Path(dir_okay=False, file_okay=True, exists=True),
    help='A path to downloaded Pfam-A (meta)data file. By default, if `download` is ``False``,'
         'will try to find it within the `hmm_dir`.'
)
def setup(hmm_dir, download, quiet, path_pfam_a, path_pfam_dat):
    """
    Command to initialize the HMM data needed for TKPs' annotation.
    """
    if not quiet:
        LOGGER.setLevel(logging.INFO)

    hmm_dir = Path(hmm_dir)
    hmm_dir.mkdir(exist_ok=True, parents=True)

    if download:
        LOGGER.info('Fetching and unpacking Pfam data')
        path_pfam_a, path_pfam_dat = fetch_pfam(hmm_dir)
    else:
        files = get_files(hmm_dir)
        if path_pfam_a is None:
            if PFAM_A_NAME not in files:
                raise ValueError(
                    f'If `download` is false, {hmm_dir} must contain {PFAM_A_NAME}'
                )
            path_pfam_a = files[PFAM_A_NAME]
        if path_pfam_dat is None:
            if PFAM_DAT_NAME not in files:
                raise ValueError(
                    f'If `download` is false, {hmm_dir} must contain {PFAM_DAT_NAME}'
                )
            path_pfam_dat = files[PFAM_DAT_NAME]
        LOGGER.info(f'Using existing {path_pfam_a, path_pfam_dat}')

    df = parse_pfam_dat(path_pfam_dat)
    df.to_csv(hmm_dir / 'pfam_entries.tsv', sep='\t', index=False)
    LOGGER.info(f'Obtained {len(df)} metadata entries from {path_pfam_dat}')

    with HMMFile(path_pfam_a) as f:
        hmms = list(f)
    LOGGER.info(f'Obtained {len(hmms)} profiles from {path_pfam_a}')

    acc2type = {
        dom_id: dom_type for dom_id, dom_type in
        df[['Accession', 'Type']].itertuples(index=False)
    }

    profiles_dir = hmm_dir / 'profiles'
    type_dirs = {x: (profiles_dir / x) for x in df['Type'].unique()}
    type_dirs[UNK_HMM] = profiles_dir / UNK_HMM
    for d in type_dirs.values():
        d.mkdir(exist_ok=True, parents=True)

    if not quiet:
        hmms = tqdm(hmms, desc='Writing hmms')

    for hmm in hmms:
        acc = hmm.accession.decode('utf-8').split('.')[0]
        hmm_type = acc2type.get(acc, UNK_HMM)
        hmm_path = type_dirs[hmm_type] / f'{acc}.hmm'
        with hmm_path.open('wb') as f:
            hmm.write(f)

    domains_dir = profiles_dir / 'Domain'
    pk_path = domains_dir / f'{PFAM_PK_NAME}.hmm'
    if not pk_path.exists():
        raise ValueError(f'Expected to find {PFAM_PK_NAME} in {domains_dir}')
    shutil.copy(pk_path, hmm_dir / pk_path.name)
    LOGGER.info(f'Copied PK profile {pk_path.name} to {hmm_dir / pk_path.name}')
    LOGGER.info('Finished HMM setup')


def gunzip(path_in: Path, path_out: Path | None = None, rm: bool = True) -> Path:
    if path_out is None:
        name_out = path_in.name.removesuffix('.gz')
        if path_in.name == name_out:
            name_out = f'{name_out}_unpacked'
        path_out = path_in.parent / name_out
    if path_out.exists():
        warn(f'Overwriting existing {path_out}')
    with gzip.open(path_in, 'rb') as f_in:
        with path_out.open('wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    if rm:
        path_in.unlink()
    return path_out


def fetch_pfam(base: Path):
    return (gunzip(download_to_file(PFAM_A_URL, root_dir=base, text=False)),
            gunzip(download_to_file(PFAM_DAT_URL, root_dir=base, text=False)))


def parse_pfam_dat(path: Path):
    def wrap_chunk(xs: list[str]):
        _id, _acc, _desc, _, _type = map(lambda x: x.split('   ')[-1], xs[:5])
        _acc = _acc.split('.')[0]
        return _id, _acc, _desc, _type

    with path.open() as f:
        lines = filter(bool, map(lambda x: x.rstrip(), f))
        chunks = filter(bool, split_at(lines, lambda x: x.startswith('# ')))
        return pd.DataFrame(
            map(wrap_chunk, chunks),
            columns=['ID', 'Accession', 'Description', 'Type'])


@tkp_finder.command('find', context_settings={"ignore_unknown_options": True})
@click.argument('fasta', nargs=-1, type=click.Path())
@click.option(
    '-H', '--hmm_dir', type=click.Path(exists=True, dir_okay=True, file_okay=False),
    help='Directory with HMM profiles. Expected to contain `profiles` dir and target '
         'PK profile (PF00069.hmm). See `tkp-finder setup` on how to prepare this dir.'
)
@click.option(
    '-t', '--hmm_type', multiple=True, default=['Family', 'Domain', 'Motif'], show_default=True,
    help='Which HMM types to use for annotating the discovered TKPs. The names must correspond to '
         'folders within he `hmm_dir`.'
)
@click.option(
    '-p', '--pk_profile', type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='A path to the PK HMM profile. By default, will try to find it within the `hmm_dir`.'
)
@click.option(
    '-m', '--motif', default=MOTIF, show_default=True,
    help='A motif to discriminate between PKs and pseudo PKs. This corresponds to the following '
         'conserved elements:: \n'
         '(1) b3-Lys'
         '(2) aC-helix Glu'
         '(3-4-5) HRD motif'
         '(6-7-8) DFG motif'
)
@click.option(
    '-o', '--output', type=click.Path(file_okay=False, dir_okay=True, writable=True),
    help='Output directory to store the results. Be default, will store within `./tkp-finder`.'
)
@click.option(
    '-n', '--num_proc', type=int, default=None,
    help='The number of cpus for data parallelism: each input fasta will be annotated within '
         'separate process. HINT: one may split large fasta files for faster processing.'
)
@click.option(
    '-q', '--quiet', is_flag=True, default=False,
    help='Disable logging and progress bar'
)
@click.option(
    '--pk_map_name', default=PK_NAME, show_default=True,
    help='Use this name for the protein kinase domain.'
)
@click.option(
    '--ppk_map_name', default=PPK_NAME, show_default=True,
    help='Use this name for pseudo protein kinases.'
)
@click.option(
    '--min_domain_size', type=int, default=150, show_default=True,
    help='The minimum number of amino acid residues within a PK domain.'
)
@click.option(
    '--min_domains', type=int, default=2,
    help='The number of domains to classify a protein as TKP.'
)
def find(
        fasta, hmm_dir, hmm_type, pk_profile, motif, output, num_proc, quiet,
        pk_map_name, ppk_map_name, min_domain_size, min_domains
):
    if not quiet:
        LOGGER.setLevel(logging.INFO)
    fasta = [Path(f) for f in fasta]
    if not fasta:
        raise ValueError('No inputs provided. Use -h or --help to invoke help.')
    for f in fasta:
        if not f.exists():
            raise ValueError(f'File {f} does not exist')
    if hmm_dir is None:
        hmm_dir = Path.cwd() / 'hmm'
        LOGGER.info(f'Assuming hmm dir to be {hmm_dir}')
        if not hmm_dir.exists():
            raise ValueError(f'HMM dir {hmm_dir} does not exist.')
        dirs = get_dirs(hmm_dir)
        if 'profiles' not in dirs:
            raise ValueError(f'Expected to find `profiles` dir in {hmm_dir}')
    if output is None:
        output = Path.cwd() / 'tkp-finder'
        output.mkdir(exist_ok=True)
        LOGGER.info(f'Setting output dir to {output}')
    else:
        output = Path(output)
    if pk_profile is None:
        pk_profile = hmm_dir / f'{PFAM_PK_NAME}.hmm'
        if not pk_profile.exists():
            raise ValueError(f'Expected to find profile {PFAM_PK_NAME} within {hmm_dir}')

    use_parallel = num_proc is not None and num_proc > 1 and len(fasta) > 1

    pipe_one = discover_and_annotate(
        pk_profile=pk_profile,
        hmm_base_dir=hmm_dir / 'profiles',
        hmm_types=hmm_type,
        min_size=min_domain_size,
        min_domains=min_domains,
        motif=motif,
        pk_map_name=pk_map_name,
        ppk_name=ppk_map_name,
        seq_variables=VARIABLES,
        quiet=True if use_parallel else quiet
    )

    if use_parallel:
        LOGGER.info(f'Processing {len(fasta)} files in parallel')
        results = yield_parallel(pipe_one, num_proc, fasta)
    else:
        results = yield_sequentially(pipe_one, fasta)

    if not quiet and len(fasta) > 1:
        results = tqdm(results, desc='Processing inputs', total=len(fasta))

    io = ChainIO(verbose=False if use_parallel else not quiet)
    for f, chains in zip_equal(fasta, results):
        if chains is None or len(chains) == 0:
            continue
        f_dir = output / f.name
        consume(io.write(chains, f_dir, write_children=True))
        df = aggregate_annotations(
            chains.collapse_children()
        ).sort_values(['ObjectID', 'HMM_type', 'Start'])
        df.to_csv(f_dir / 'summary.tsv', sep='\t', index=False)


@curry
def find_tkps(
        path: Path, profile: Path, min_size: int = 150,
        min_domains: int = 2, map_name: str = PK_NAME, quiet: bool = True
) -> ChainList:
    def wrap_pbar(it):
        if quiet:
            return it
        return tqdm(it, desc='Initializing chains')

    initializer = ChainInitializer()
    annotator = PyHMMer(profile, bit_cutoffs='trusted')

    chains = pipe(
        path,
        curry(read_fasta)(strip_id=True),
        curry(unique_everseen)(key=op.itemgetter(1)),
        initializer.from_iterable,
        # wrap_pbar,
        ChainList
    )

    consume(annotator.annotate(chains, min_size=min_size, new_map_name=map_name))

    return pipe(
        chains,
        lambda x: x.filter(lambda c: len(c.children) >= min_domains),
        filter_child_overlaps,
        lambda x: x.filter(lambda c: len(c.children) >= min_domains),
    )


@curry
def calculate_variables(
        chains: ChainList,
        vs: abc.Sequence[SequenceVariable] = VARIABLES,
        map_name: str = PK_NAME,
) -> pd.DataFrame:
    manager = Manager()
    df = manager.aggregate_from_it(
        manager.calculate(chains, vs, map_name=map_name))
    return df


@curry
def annotate_by_hmms(
        chains: ChainList, hmm_paths: abc.Iterable[Path], hmm_type: str,
        quiet: bool = True
) -> ChainList:
    if not quiet:
        hmm_paths = tqdm(hmm_paths, desc=f'Annotating by HMM {hmm_type}')
    for path in hmm_paths:
        map_name = f'{hmm_type}_{path.stem}'
        annotator = PyHMMer(path, bit_cutoffs='trusted')
        consume(annotator.annotate(chains, new_map_name=map_name))
    return chains


@curry
def filter_child_overlaps(
        chains: ChainList,
        filt_fn: abc.Callable[[ChainSequence], bool] = lambda _: True,
        val_fn: abc.Callable[[ChainSequence], float] = len,
) -> ChainList:
    for c in chains:
        target_children = filter(filt_fn, next(c.iter_children()))
        non_overlapping = resolve_overlaps(target_children, value_fn=val_fn)
        non_overlapping_ids = [x.id for x in non_overlapping]
        c.children = valfilter(
            lambda x: not filt_fn(x) or x.id in non_overlapping_ids, c.children)
    return chains


def annotate_ppks(
        chains: ChainList, vs_df: pd.DataFrame,
        pk_name: str = PK_NAME, ppk_name: str = PPK_NAME, motif: str = MOTIF,
):
    def is_ppk(motif_observed):
        assert len(motif) == len(motif_observed), 'motif sizes must match'
        for c1, c2 in zip(motif, motif_observed):
            if c1 != 'X' and c1 != c2:
                return True
        return False

    vs_df = vs_df.copy().fillna('X')
    vs_df['Motif'] = [''.join(x[1:]) for x in vs_df.itertuples(index=False)]
    id2motif = dict(x for x in vs_df[['ObjectID', 'Motif']].itertuples(index=False))

    for c in chains:
        try:
            chain_motif = id2motif[c.id]
            is_pseudo = is_ppk(chain_motif)
            if is_pseudo:
                c.name = c.name.replace(pk_name, ppk_name)
            c.meta['motif'] = chain_motif
        except KeyError:
            c.meta['motif'] = '-'

    return chains, vs_df


def aggregate_annotations(
        chains: abc.Iterable[ChainSequence],
        pk_name: str = PK_NAME, ppk_name: str = PPK_NAME
) -> pd.DataFrame:
    def agg_one(c):
        if pk_name in c.name or ppk_name in c.name:
            hmm_type = 'Target'
            hmm_name = c.id.split('_')[0]
        else:
            hmm_type, hmm_name = c.id.split('_')[:2]
        score = next(filter(lambda x: x[0].endswith('score'), c.meta.items()))[1]
        parent_name = c.parent.name
        return hmm_type, hmm_name, c.id, parent_name, c.start, c.end, score

    return pd.DataFrame(
        map(agg_one, chains),
        columns=['HMM_type', 'HMM', 'ObjectID', 'ParentName', 'Start', 'End', 'BitScore']
    )


@curry
def discover_and_annotate(
        path: Path, pk_profile: Path, hmm_base_dir: Path,
        hmm_types: abc.Iterable[str] = ('Family', 'Domain', 'Motif'),
        min_size: int = 150, min_domains: int = 2, motif=MOTIF,
        pk_map_name: str = PK_NAME, ppk_name: str = PPK_NAME,
        seq_variables: abc.Sequence[SequenceVariable] = VARIABLES,
        quiet: bool = True,
) -> ChainList[ChainSequence]:
    @curry
    def value_fn(c: ChainSequence, hmm_type: str):
        scores = keyfilter(
            lambda x: x.startswith(hmm_type) and x.endswith('score'),
            c.meta)
        if len(scores) > 1:
            raise ValueError(
                f'Expected exactly one score for hmm type {hmm_type}, '
                f'got {len(scores)}: {scores}')
        return scores.popitem()[1]

    @curry
    def annotate_and_filter(chains, hmm_type):
        hmms = list((hmm_base_dir / hmm_type).glob('*hmm'))
        return pipe(
            chains,
            annotate_by_hmms(hmm_paths=hmms, hmm_type=hmm_type),
            filter_child_overlaps(
                filt_fn=lambda c: any(
                    (x.startswith(hmm_type) for x in c.fields)
                ),
                val_fn=value_fn(hmm_type=hmm_type)
            )
        )

    # if not pk_map_name.startswith('Domain'):
    #     pk_map_name = f'Domain_{pk_map_name}'

    chains = find_tkps(
        path,
        min_size=min_size, min_domains=min_domains,
        map_name=pk_map_name, profile=pk_profile, quiet=quiet
    )
    if len(chains) == 0:
        LOGGER.info(f'Found no TKPs in {path}')
        return chains

    chains = compose_left(*(annotate_and_filter(hmm_type=x) for x in hmm_types))(chains)

    pk_children = filter(lambda x: pk_map_name in x, chains.collapse_children())
    vs_df = calculate_variables(pk_children, seq_variables, pk_map_name)

    annotate_ppks(
        chains.collapse_children(), vs_df,
        pk_name=pk_map_name, ppk_name=ppk_name, motif=motif
    )

    return chains


def yield_sequentially(fn, *args):
    yield from map(fn, *args)


def yield_parallel(fn, num_proc, objs):
    with ProcessPoolExecutor(num_proc) as executor:
        futures = [(o, executor.submit(fn, o)) for o in objs]
        for o, f in futures:
            try:
                yield f.result()
            except Exception as e:
                LOGGER.error(f'Failed on input {o} with {e}; stacktrace below')
                LOGGER.exception(e)
                yield None


if __name__ == '__main__':
    tkp_finder()
