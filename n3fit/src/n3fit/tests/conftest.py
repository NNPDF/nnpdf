"""
Add markers for pytest
"""

import sys

import pytest

from validphys.loader import FallbackLoader

THEORYID = 40_000_000


@pytest.fixture(scope='module')
def nnpdf_theory_card():
    """Return a theory card already loaded as a dictionary"""
    th = FallbackLoader().check_theoryID(THEORYID)
    return th.get_description()


def pytest_runtest_setup(item):
    ALL = {"darwin", "linux"}
    supported_platforms = ALL.intersection(mark.name for mark in item.iter_markers())
    plat = sys.platform
    if supported_platforms and plat not in supported_platforms:
        pytest.skip("cannot run on platform {}".format(plat))


def pytest_configure(config):
    config.addinivalue_line("markers", "linux: mark test to run only on linux")
