"""
    Add markers for pytest
"""
import sys
import pytest
import pathlib


@pytest.fixture
def tmp(tmpdir):
    """A tempdir that is manipulated like pathlib Paths"""
    return pathlib.Path(tmpdir)


def pytest_runtest_setup(item):
    ALL = {"darwin", "linux"}
    supported_platforms = ALL.intersection(mark.name for mark in item.iter_markers())
    plat = sys.platform
    if supported_platforms and plat not in supported_platforms:
        pytest.skip("cannot run on platform {}".format(plat))


def pytest_configure(config):
    config.addinivalue_line("markers", "linux: mark test to run only on linux")
