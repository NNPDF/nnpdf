"""
Tests for the nnpdf-data CLI (nnpdf_data.cli.script) list and latex commands.
"""

from pathlib import Path
import sys
from unittest.mock import patch

import pytest

from nnpdf_data.cli import script
from nnpdf_data.tests.conftest import _mock_fetch_bibtex

# -- Tests for the ``list`` command --


def test_cmd_list(capsys, real_index):
    """Test the list command prints datasets and respects filter patterns."""
    # Unfiltered list should contain entries
    script._cmd_list(real_index)
    assert "ATLAS" in capsys.readouterr().out

    # Filtered list should match only ATLAS 1JET datasets
    script._cmd_list(real_index, filter_pattern="ATLAS_1JET*")
    out = capsys.readouterr().out
    assert "ATLAS" in out and "CMS" not in out

    # Nonexistent filter should show "No entries found"
    script._cmd_list(real_index, filter_pattern="NONEXISTENT_12345*")
    assert "No entries found" in capsys.readouterr().out


# -- Tests for the ``latex`` command --


def test_cmd_latex_stdout(dataset_1, dataset_2, tmp_runcard, capsys):
    """Test that the latex command prints LaTeX and BibTeX to stdout."""
    entries = {dataset_1.name: dataset_1, dataset_2.name: dataset_2}

    with patch("nnpdf_data.cli.latex.fetch_bibtex", side_effect=_mock_fetch_bibtex):
        script._cmd_latex(entries, tmp_runcard)

    captured = capsys.readouterr()
    assert "\\cite{" in captured.out


def test_cmd_latex_file_output(dataset_1, dataset_2, tmp_runcard, tmp_path, capsys):
    """Test that the latex command writes a .bib file when --output-bib is given."""
    entries = {dataset_1.name: dataset_1, dataset_2.name: dataset_2}
    bib_path = tmp_path / "out.bib"

    with patch("nnpdf_data.cli.latex.fetch_bibtex", side_effect=_mock_fetch_bibtex):
        script._cmd_latex(entries, tmp_runcard, output_bib=bib_path)

    assert bib_path.is_file()
    content = bib_path.read_text()
    assert "@article{" in content


def test_cmd_latex_missing_runcard(real_index):
    """Test that a missing runcard path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        script._cmd_latex(real_index, Path("/nonexistent.yml"))


def test_cmd_latex_missing_dataset_inputs(real_index, tmp_path):
    """Test that a runcard without dataset_inputs raises KeyError."""
    runcard = tmp_path / "_bad_runcard.yml"
    runcard.write_text("theory: {theoryid: 123}\n")
    try:
        with pytest.raises(KeyError):
            script._cmd_latex(real_index, runcard)
    finally:
        runcard.unlink(missing_ok=True)


# -- Tests for the main entry point --


def test_main_no_command(capsys):
    """Test that running without a subcommand prints help and returns 0."""
    with patch.object(sys, "argv", ["nnpdf-data"]):
        with patch("nnpdf_data.cli.script.index.build_index"):
            assert script.main() == 0

    assert "usage" in capsys.readouterr().out.lower()


def test_main_group_by_requires_sort_mode():
    """Test that --group-by without a sort key causes a SystemExit."""
    with patch.object(sys, "argv", ["nnpdf-data", "latex", "dummy.yml", "--group-by"]):
        with pytest.raises(SystemExit):
            script.main()
