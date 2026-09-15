"""
Tests for the LaTeX generation module (nnpdf_data.cli.latex).

Verifies the LaTeX character escaper, row builder, and table generator
across a range of scenarios using real datasets from the commondata package.
"""

from unittest.mock import MagicMock, patch

from nnpdf_data.cli.inspire import BibliographyEntry
from nnpdf_data.cli.latex import LatexDatasetRow, _escape_latex, build_latex_rows, generate_table
from nnpdf_data.tests.conftest import _mock_fetch_bibtex


def test_escape_latex():
    """Test that special characters are escaped but LaTeX math mode is preserved."""
    assert _escape_latex("foo_bar") == r"foo\_bar"
    assert _escape_latex("100%") == r"100\%"
    # Math mode substrings should be left intact
    assert r"$\sigma$" in _escape_latex(r"Cross section $\sigma$")
    # Existing LaTeX commands should pass through
    assert _escape_latex(r"foo\bar") == r"foo\bar"


def test_build_latex_rows(dataset_1):
    """Test row building with a real dataset that has a valid INSPIRE URL."""
    with patch("nnpdf_data.cli.latex.fetch_bibtex", side_effect=_mock_fetch_bibtex):
        result = build_latex_rows({dataset_1.name: dataset_1})

    # Should have one row in "unsorted"
    assert "unsorted" in result
    assert len(result["unsorted"]) == 1
    row = result["unsorted"][0]
    assert row.experiment == dataset_1.experiment
    assert row.process == dataset_1.process
    # bib_entry should have been resolved
    assert row.bib_entry is not None
    assert row.bib_key is not None


def test_build_latex_rows_grouping(real_index):
    """Test grouping real datasets by experiment."""
    # Grab a few ATLAS datasets with inspire URLs
    atlas_datasets = {
        name: ds for name, ds in real_index.items() if ds.experiment == "ATLAS" and ds.inspire_url
    }
    atlas_datasets = dict(list(atlas_datasets.items())[:3])

    with patch("nnpdf_data.cli.latex.fetch_bibtex", side_effect=_mock_fetch_bibtex):
        grouped = build_latex_rows(atlas_datasets, group_by="experiment")

    assert "ATLAS" in grouped
    assert len(grouped["ATLAS"]) == len(atlas_datasets)


def test_build_latex_rows_missing_inspire(dataset_no_inspire):
    """Test that real datasets without an INSPIRE URL do not crash and have bib_entry=None."""
    with patch("nnpdf_data.cli.latex.fetch_bibtex", side_effect=_mock_fetch_bibtex):
        result = build_latex_rows({dataset_no_inspire.name: dataset_no_inspire})

    assert result["unsorted"][0].bib_entry is None


def test_build_latex_rows_bad_inspire_url():
    """Test that a non-standard INSPIRE URL is handled gracefully (bib_entry=None)."""
    ds = MagicMock()
    ds.experiment = "TEST"
    ds.inspire_url = "https://example.com/not-inspire"
    ds.plotting.dataset_label = "Test_Label"
    ds.description = "Test"
    ds.process = "P"
    ds.nnpdf31_process = "P31"

    with patch("nnpdf_data.cli.latex.fetch_bibtex", side_effect=_mock_fetch_bibtex):
        result = build_latex_rows({"TEST_DS": ds})

    assert result["unsorted"][0].bib_entry is None


def test_build_latex_rows_empty_inspire_url():
    """Test that an empty INSPIRE URL is handled gracefully (bib_entry=None)."""
    ds = MagicMock()
    ds.experiment = "TEST"
    ds.inspire_url = ""
    ds.plotting.dataset_label = "Test_Label"
    ds.description = "Test"
    ds.process = "P"
    ds.nnpdf31_process = "P31"

    with patch("nnpdf_data.cli.latex.fetch_bibtex", side_effect=_mock_fetch_bibtex):
        result = build_latex_rows({"TEST_DS": ds})

    assert result["unsorted"][0].bib_entry is None


def test_generate_table_with_reference():
    """Test table generation when bibliography entries are available."""
    bib = BibliographyEntry(key="Key1", raw_bibtex="@article{Key1}")
    row = LatexDatasetRow("DATA", "Label", "Obs", "TEST", "P", "P31", bib)
    table, bibtex = generate_table([row])

    assert "\\begin{table}" in table
    assert "\\cite{Key1}" in table
    assert len(bibtex) == 1


def test_generate_table_without_reference():
    """Test table generation when bib_entry is None (shows 'Reference unavailable')."""
    row_none = LatexDatasetRow("DATA", "Label", "Obs", "TEST", "P", "P31", None)
    table2, bibtex2 = generate_table([row_none])

    assert r"\textit{Reference unavailable}" in table2
    assert len(bibtex2) == 0
