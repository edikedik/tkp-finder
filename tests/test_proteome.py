from pathlib import Path

import pytest

from tkp_finder.proteome import EnsemblProteome, NCBIProteome, Phytozome_proteome


def ensemble_animals_path() -> Path:
    return (
        Path(__file__).parent
        / "data"
        / "test_proteomes"
        / "Ensembl_animals"
        / "Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz"
    )


def ensemble_plants_path() -> Path:
    return (
        Path(__file__).parent
        / "data"
        / "test_proteomes"
        / "Ensembl_plants"
        / "Chondrus_crispus.ASM35022v2.pep.all.fa.gz"
    )


def ncbi_path() -> Path:
    return (
        Path(__file__).parent
        / "data"
        / "test_proteomes"
        / "NCBI_plants"
        / "Ostreobium_quekettii.faa.tar.gz"
    )


def phytozome_path() -> Path:
    return (
        Path(__file__).parent
        / "data"
        / "test_proteomes"
        / "Phytozome_plants"
        / "Ppatens_318_v3.3.protein_primaryTranscriptOnly.fa.gz"
    )


@pytest.mark.parametrize(
    "path,size", [(ensemble_animals_path(), 6600), (ensemble_plants_path(), 9807)]
)
def test_ensemble(path, size):
    proteome = EnsemblProteome(path, "")
    chains = list(proteome.preprocess())
    assert len(chains) == size
    for c in chains:
        assert c.meta["filename"] == path.name
        assert "transcript" in c.meta
        assert "gene" in c.meta


def test_ncbi_proteome():
    path = ncbi_path()
    proteome = NCBIProteome(path, "")
    chains = list(proteome.preprocess())
    for c in chains:
        assert c.meta["filename"] == path.name
        assert c.meta["organism"] == path.name.split(".")[0]


def test_phytozome():
    path = phytozome_path()
    proteome = Phytozome_proteome(path, "")
    chains = list(proteome.preprocess())
    for c in chains:
        assert c.meta["filename"] == path.name
        assert "transcript" in c.meta
