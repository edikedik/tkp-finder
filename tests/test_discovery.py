from copy import deepcopy
from pathlib import Path

import lXtractor.chain as lxc
import pytest
from lXtractor.util import read_fasta

from tkp_finder.discovery import PKDiscoverer, get_pskd_hmm


def large_insert_pk_chains() -> lxc.ChainList[lxc.ChainSequence]:
    seqs = list(read_fasta(Path("data/large_insert_cases.fasta")))
    seqs = ((s[0].split("|")[0], s[1]) for s in seqs)
    return lxc.ChainList(map(lxc.ChainSequence.from_tuple, seqs))


@pytest.fixture
def large_insert_pks() -> list[tuple[str, str]]:
    return list(read_fasta(Path("data/large_insert_cases.fasta")))


@pytest.fixture
def sbmsa_hmm():
    # Returns a protein kinase HMM derived from structure-based MSA
    return get_pskd_hmm()


@pytest.fixture
def one_domain_pk():
    # F4I4F2|1-499
    # PK boundaries 115-436
    return (
        "MLLKPGNKLVSPETSHHRDSASNSSNHKCQQQKPRKDKQKQVEQNTKKIEEHQIKSESTLLISNHNVNMS"
        "SQSNNSESTSTNNSSKPHTGGDIRWDAVNSLKSRGIKLGISDFRVLKRLGYGDIGSVYLVELKGANPTTY"
        "FAMKVMDKASLVSRNKLLRAQTEREILSQLDHPFLPTLYSHFETDKFYCLVMEFCSGGNLYSLRQKQPNK"
        "CFTEDAARFFASEVLLALEYLHMLGIVYRDLKPENVLVRDDGHIMLSDFDLSLRCSVNPTLVKSFNGGGT"
        "TGIIDDNAAVQGCYQPSAFFPRMLQSSKKNRKSKSDFDGSLPELMAEPTNVKSMSFVGTHEYLAPEIIKN"
        "EGHGSAVDWWTFGIFIYELLHGATPFKGQGNKATLYNVIGQPLRFPEYSQVSSTAKDLIKGLLVKEPQNR"
        "IAYKRGATEIKQHPFFEGVNWALIRGETPPHLPEPVDFSCYVKKEKESLPPAATEKKSKMFDEANKSGSD"
        "PDYIVFEYF"
    )


@pytest.fixture
def make_pk_with_insert():
    def _make(pos: int, insert_size: int, ref_seq: str) -> str:
        """Insert sequence of N's at specified position"""
        return ref_seq[:pos] + "X" * insert_size + ref_seq[pos:]

    return _make


@pytest.fixture
def make_pk_multiple_inserts():
    def _make(positions: list[int], sizes: list[int], ref_seq: str) -> str:
        result = ref_seq
        for pos, size in sorted(zip(positions, sizes), reverse=True):
            result = result[:pos] + "X" * size + result[pos:]
        return result

    return _make


def test_cases_single_insert():
    return [
        # Core functionality tests
        (
            "simple_insert",
            {
                "pos": 250,  # Middle of PK domain
                "size": 100,
                "expected_fragments": 2,
                "should_merge": True,
            },
        ),
        (
            "simple_insert 3",
            {
                "pos": 170,  # Towards the start
                "size": 100,
                "expected_fragments": 2,
                "should_merge": True,
            },
        ),
        (
            "simple_insert 3",
            {
                "pos": 350,  # Towards the end
                "size": 100,
                "expected_fragments": 2,
                "should_merge": True,
            },
        ),
        # Edge cases
        (
            "insert_at_boundary",
            {
                "pos": 110,  # At domain start
                "size": 50,
                "expected_fragments": 1,
                "should_merge": False,
            },
        ),
        # Failure cases
        (
            "large_insert",
            {
                "pos": 300,
                "size": 500,  # Too large to merge
                "expected_fragments": 2,
                "should_merge": False,
            },
        ),
    ]


@pytest.mark.parametrize("case_name,params", test_cases_single_insert())
def test_pk_discovery_with_inserts(
    case_name: str, params: dict, make_pk_with_insert, one_domain_pk: str, sbmsa_hmm
):
    # Setup
    seq_with_insert = make_pk_with_insert(params["pos"], params["size"], one_domain_pk)
    cs = lxc.ChainSequence.from_string(seq_with_insert)

    # Run discovery
    discoverer = PKDiscoverer(sbmsa_hmm)

    # Initial annotation
    hits = list(discoverer.annotate([cs]))
    assert len(hits) == params["expected_fragments"]

    # Test merging
    discoverer.join_hits(cs, max_insert=250, min_size=30, max_overlap=5)
    merged_hits = [h for h in cs.children if h.name == discoverer.assign_name]

    if params["should_merge"]:
        assert len(merged_hits) == 1
        # Validate merged hit properties
        merged = merged_hits[0]
        assert merged.start < merged.end
        assert discoverer.assign_name == merged.name
        assert discoverer.assign_name in merged.fields
    else:
        assert len(merged_hits) == params["expected_fragments"]


def test_cases_multiple_inserts():
    return [
        {
            "pos": [200, 300],
            "size": [110, 110],
            "expected_fragments": 3,
            "should_merge": True,
            "max_insert": 200,
            "max_overlap": 10,
        },
        {
            # Should "eat" the starting five positions
            "pos": [120, 250],
            "size": [100, 100],
            "expected_fragments": 2,
            "should_merge": True,
            "max_insert": 120,
            "max_overlap": 10,
        },
        {
            # Should "eat" the ending dozen positions
            "pos": [250, 420],
            "size": [100, 100],
            "expected_fragments": 2,
            "should_merge": True,
            "max_insert": 120,
            "max_overlap": 10,
        },
    ]


@pytest.mark.parametrize("params", test_cases_multiple_inserts())
def test_pk_discovery_with_multiple_inserts(
    params, make_pk_multiple_inserts, one_domain_pk, sbmsa_hmm
):
    seq_with_insert = make_pk_multiple_inserts(
        params["pos"], params["size"], one_domain_pk
    )
    cs = lxc.ChainSequence.from_string(seq_with_insert)
    discoverer = PKDiscoverer(sbmsa_hmm)
    hits = list(discoverer.annotate([cs]))
    assert len(hits) == params["expected_fragments"]
    discoverer.join_hits(
        cs,
        max_insert=params["max_insert"],
        min_size=30,
        max_overlap=params["max_overlap"],
    )
    merged_hits = [h for h in cs.children if h.name == discoverer.assign_name]
    if params["should_merge"]:
        assert len(merged_hits) == 1
    else:
        assert len(merged_hits) == params["expected_fragments"]


@pytest.mark.parametrize("chain_seq", large_insert_pk_chains())
@pytest.mark.parametrize("max_insert", [100, 200, 300, 500])
def test_fragment_joining(chain_seq, sbmsa_hmm, max_insert):
    chain_seq = deepcopy(chain_seq)
    pkd = PKDiscoverer(sbmsa_hmm)

    hits = list(pkd.annotate([chain_seq]))
    assert len(hits) == 2
    insert_size = hits[1].start - hits[0].end + 1

    pkd.join_hits(chain_seq, max_insert=max_insert, min_size=30, max_overlap=10)
    expected = 1 if insert_size < max_insert else 2
    assert len(chain_seq.children) == expected
