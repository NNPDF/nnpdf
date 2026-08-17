"""
Configuration for nnpdf-data tests

Provides fixtures using real datasets from the commondata package
and mock bibtex entries for tests that exercise the INSPIRE API layer.
"""

from pathlib import Path

import pytest

from nnpdf_data.cli.index import build_index
from nnpdf_data.cli.inspire import BibliographyEntry

# Real dataset names from the commondata package used across tests
DATASET_1 = "ATLAS_1JET_13TEV_DIF_PT-Y"
DATASET_2 = "LHCB_DY_7TEV_MUON_Y"
DATASET_NO_INSPIRE = "ATLAS_DY_13TEV_TOT"

# Mock bibtex entries for the INSPIRE record IDs of DATASET_1 and DATASET_2
MOCK_BIBTEX = {
    "1634970": BibliographyEntry(
        key="ATLAS:2020xxx",
        raw_bibtex="@article{ATLAS:2020xxx,\n  author = {ATLAS Collab},\n  year = {2020}\n}",
    ),
    "1373300": BibliographyEntry(
        key="LHCB:2015abc",
        raw_bibtex="@article{LHCB:2015abc,\n  author = {LHCb Collab},\n  year = {2015}\n}",
    ),
}


def _mock_fetch_bibtex(record_id):
    """Return a mock bibliography entry for the given INSPIRE record ID."""
    return MOCK_BIBTEX.get(record_id)


@pytest.fixture(scope="session")
def real_index():
    """Build the full dataset index once for the entire test session."""
    return build_index()


@pytest.fixture
def dataset_1(real_index):
    """Return the first real dataset (ATLAS_1JET_13TEV_DIF_PT-Y)."""
    return real_index[DATASET_1]


@pytest.fixture
def dataset_2(real_index):
    """Return the second real dataset (LHCB_DY_7TEV_MUON_Y)."""
    return real_index[DATASET_2]


@pytest.fixture
def dataset_no_inspire(real_index):
    """Return a real dataset that has no INSPIRE URL."""
    return real_index[DATASET_NO_INSPIRE]


@pytest.fixture
def sample_bib_entry():
    """Return a simple BibliographyEntry fixture for use in tests."""
    return BibliographyEntry(
        key="Test:2024abc",
        raw_bibtex="@article{Test:2024abc,\n  author = {Test},\n  year = {2024}\n}",
    )


@pytest.fixture
def tmp_runcard(tmp_path):
    """Create a minimal runcard YAML file on a temp directory."""
    runcard = tmp_path / "test_runcard.yml"
    runcard.write_text(
        "dataset_inputs:\n" f"  - {{dataset: {DATASET_1}}}\n" f"  - {{dataset: {DATASET_2}}}\n"
    )
    return runcard
