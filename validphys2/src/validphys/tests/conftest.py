"""
conftest.py

Pytest fixtures.
"""
import pathlib
import sys

import pytest
from hypothesis import settings

# Adding this here to change the time of deadline from default (200ms) to 1500ms
settings.register_profile("extratime", deadline=1500)
settings.load_profile("extratime")

#Fortunately py.test works much like reportengine and providers are
#connected by argument names.
@pytest.fixture
def tmp(tmpdir):
    """A tempdir that is manipulated like pathlib Paths"""
    return pathlib.Path(tmpdir)

# Here define the default config items like the PDF, theory and experiment specs

DATA = [
    {'dataset': 'NMC'},
    {'dataset': 'ATLASTTBARTOT', 'cfac':['QCD']},
    {'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10}
]


SINGLE_EXP = [
    {
        'experiment': 'pseudo experiment',
        'datasets': [
            {'dataset': 'NMC'},
            {'dataset': 'ATLASTTBARTOT', 'cfac':['QCD']},
            {'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10}]}]

WEIGHTED_DATA = [
    {'dataset': 'NMC'},
    {'dataset': 'NMC', 'weight': 100},
]

PDF = "NNPDF40_nnlo_as_01180"
HESSIAN_PDF = "NNPDF40_nnlo_as_01180_hessian"
THEORYID = 162
FIT = "NNPDF40_nnlo_lowprecision"
FIT_ITERATED = "NNPDF40_nnlo_lowprecision_iterated"

base_config = dict(
        pdf=PDF,
        use_cuts='nocuts',
        dataset_inputs=DATA,
        theoryid=THEORYID,
        use_fitthcovmat=False
    )

@pytest.fixture(scope='module')
def data_config():
    return base_config

@pytest.fixture(scope='module')
def data_internal_cuts_config(data_config):
    config_dict = dict(data_config)
    config_dict.update(use_cuts='internal')
    return config_dict

@pytest.fixture(scope='module')
def single_data_internal_cuts_config(data_internal_cuts_config):
    """Like data_internal_cuts_config but for a single dataset"""
    config_dict = dict(data_internal_cuts_config)
    config_dict.pop("dataset_inputs")
    config_dict.update(dataset_input=DATA[0])
    return config_dict

@pytest.fixture(scope='module')
def data_witht0_config():
    config_dict = dict(
        **base_config,
        use_t0=True,
        t0pdfset=PDF)
    return config_dict

@pytest.fixture(scope='module')
def data_witht0_internal_cuts_config(data_witht0_config):
    config_dict = dict(data_witht0_config)
    config_dict.update(use_cuts='internal')
    return config_dict

@pytest.fixture(scope='module')
def data_singleexp_witht0_config(data_witht0_config):
    config_dict = dict(data_witht0_config)
    config_dict.pop("dataset_inputs")
    config_dict.update({'experiments': SINGLE_EXP})
    config_dict.update(use_cuts='internal')
    return config_dict

@pytest.fixture(scope='module')
def weighted_data_witht0_config(data_witht0_config):
    config_dict = dict(data_witht0_config)
    config_dict.update({'dataset_inputs': WEIGHTED_DATA})
    return config_dict

@pytest.fixture(scope='module')
def weighted_data_witht0_internal_cuts_config(data_witht0_internal_cuts_config):
    config_dict = dict(data_witht0_internal_cuts_config)
    config_dict.update({'dataset_inputs': WEIGHTED_DATA})
    return config_dict

def pytest_runtest_setup(item):
    ALL = {"darwin", "linux"}
    supported_platforms = ALL.intersection(mark.name for mark in item.iter_markers())
    plat = sys.platform
    if supported_platforms and plat not in supported_platforms:
        pytest.skip("cannot run on platform {}".format(plat))

def pytest_configure(config):
    config.addinivalue_line(
        "markers", "linux: mark test to run only on linux"
    )

@pytest.fixture(scope='module')
def flavour_basis_initial_scale_config():
    """Set the basis to ``flavour`` and fix ``Q`` to be 1.651, which is
    about the initial scale at the point of creating this fixture.

    In practice these values don't matter, but we will fix them here for the
    sake of having stable tests.

    """
    return {"basis": "flavour", "Q": 1.651}

@pytest.fixture(scope='module')
def mc_pdf_config(flavour_basis_initial_scale_config):
    """``flavour_basis_initial_scale_config`` with pdf set to be
    a MC pdf

    """
    return {"pdf": PDF, **flavour_basis_initial_scale_config}

@pytest.fixture(scope='module')
def hessian_pdf_config(flavour_basis_initial_scale_config):
    """``flavour_basis_initial_scale_config`` with pdf set to be
    a hessian pdf

    """
    return {"pdf": HESSIAN_PDF, **flavour_basis_initial_scale_config}

@pytest.fixture(scope='module')
def hessian_data_config(data_config):
    """Same as data config but with hessian PDF"""
    new_config = dict(data_config)
    new_config["pdf"] = HESSIAN_PDF
    return new_config

@pytest.fixture(scope='module')
def hessian_data_internal_cuts_config(data_internal_cuts_config):
    """Same as data config but with hessian PDF"""
    new_config = dict(data_internal_cuts_config)
    new_config["pdf"] = HESSIAN_PDF
    return new_config

@pytest.fixture(scope='module')
def hessian_single_data_internal_cuts_config(single_data_internal_cuts_config):
    """Same as single data config but with hessian PDF"""
    new_config = dict(single_data_internal_cuts_config)
    new_config["pdf"] = HESSIAN_PDF
    return new_config
