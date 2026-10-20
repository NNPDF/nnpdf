"""
Tests for the INSPIRE API client module (nnpdf_data.cli.inspire).

Exercises the BibTeX key extractor, the network fetcher (with mocked
urlopen), and the retry/error-handling logic so that no real HTTP
requests are made during CI.
"""

from unittest.mock import MagicMock, patch
from urllib.error import HTTPError

from nnpdf_data.cli.inspire import BibliographyEntry, _extract_key, fetch_bibtex


def test_extract_key():
    """Test that the BibTeX key is correctly extracted from raw entries."""
    # Valid entries with different whitespace
    assert _extract_key("@article{ATLAS:2017rso,\n}") == "ATLAS:2017rso"
    assert _extract_key("  @article{Key123,\n}") == "Key123"
    # Invalid entries should return None
    assert _extract_key("") is None
    assert _extract_key("invalid bibtex") is None


def test_fetch_success():
    """Test that fetch_bibtex returns a BibliographyEntry on a successful response."""
    mock_resp = MagicMock()
    mock_resp.read.return_value = b"@article{Test:2021xyz,\n}"
    mock_resp.__enter__ = MagicMock(return_value=mock_resp)
    mock_resp.__exit__ = MagicMock(return_value=False)

    with patch("nnpdf_data.cli.inspire.urlopen", return_value=mock_resp):
        result = fetch_bibtex(1234567)

    # Verify the key was extracted correctly
    assert result.key == "Test:2021xyz"


def test_fetch_http_404():
    """Test that a 404 error returns None without retrying."""
    err404 = HTTPError("url", 404, "Not Found", {}, None)
    with patch("nnpdf_data.cli.inspire.urlopen", side_effect=err404):
        assert fetch_bibtex(1) is None


def test_fetch_http_500_retry():
    """Test that a 500 error triggers a retry and succeeds on the second attempt."""
    err500 = HTTPError("url", 500, "Server Error", {}, None)
    mock_resp = MagicMock()
    mock_resp.read.return_value = b"@article{Retry,\n}"
    mock_resp.__enter__ = MagicMock(return_value=mock_resp)
    mock_resp.__exit__ = MagicMock(return_value=False)

    with patch("nnpdf_data.cli.inspire.urlopen", side_effect=[err500, mock_resp]):
        result = fetch_bibtex(2)

    assert result.key == "Retry"
