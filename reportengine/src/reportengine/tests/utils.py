"""
utils.py

Common utilities for testing.
"""
from pathlib import Path

import pytest

@pytest.fixture
def tmp(tmpdir):
    """A fixture that returns a pathlib object representing
    a newly created temporary path."""
    return Path(tmpdir)

