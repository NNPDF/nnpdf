"""
conftest.py

Pytest fixtures.
"""

import contextlib
import pathlib
import sys

from hypothesis import settings
import lhapdf
import pytest

# Adding this here to change the time of deadline from default (200ms) to 1500ms
settings.register_profile("extratime", deadline=1500)
settings.load_profile("extratime")


# Fortunately py.test works much like reportengine and providers are
# connected by argument names.
@pytest.fixture
def tmp(tmpdir):
    """A tempdir that is manipulated like pathlib Paths"""
    return pathlib.Path(tmpdir)


# Here define the default config items like the PDF, theory and experiment specs
SINGLE_DATAPOINT = {'dataset': 'ATLAS_TTBAR_8TEV_TOT_X-SEC', 'variant': 'legacy'}

SINGLE_DATASET = {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy'}

SINGLE_CATEGORICAL = {"dataset": "ATLAS_DY_13TEV_TOT", 'variant': 'legacy'}

DATA = [
    {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy'},
    {'dataset': 'ATLAS_TTBAR_7TEV_TOT_X-SEC', 'variant': 'legacy'},
    {'dataset': 'CMS_Z0J_8TEV_PT-Y', 'cfac': ['NRM'], 'variant': 'legacy_10'},
    # Explicitly put a CMS dataset between the two ATLAS
    SINGLE_DATAPOINT,
]


SINGLE_EXP = [
    {
        'experiment': 'pseudo experiment',
        'datasets': [
            {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy'},
            {'dataset': 'ATLAS_TTBAR_7TEV_TOT_X-SEC', 'variant': 'legacy'},
            {'dataset': 'CMS_Z0J_8TEV_PT-Y', 'cfac': ['NRM'], 'variant': 'legacy_10'},
        ],
    }
]

WEIGHTED_DATA = [
    {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy'},
    {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy', 'weight': 100},
]

DATA_THCOVMAT = [
    {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy'},
    {'dataset': 'CHORUS_CC_NOTFIXED_PB_DW_NU-SIGMARED', 'variant': 'legacy'},
    {'dataset': 'CMS_Z0J_8TEV_PT-Y', 'cfac': ['NRM'], 'variant': 'legacy_10'},
    {'dataset': 'ATLAS_WJ_8TEV_WP-PT', 'variant': 'legacy'},
    {'dataset': 'LHCB_Z0_8TEV_MUON_Y', 'cfac': ['NRM']},
]

POSITIVITIES = ["NNPDF_POS_2P24GEV_DYU", "NNPDF_POS_2P24GEV_F2S"]

PDF = "NNPDF40_nnlo_as_01180"
HESSIAN_PDF = "NNPDF40_nnlo_as_01180_hessian"
THEORYID = 162
THEORYID_NEW = 399
THEORY_QED = 398
FIT = "NNPDF40_nnlo_low_precision_240916"
FIT_3REPLICAS = "Basic_runcard_3replicas_lowprec_221130"
FIT_3REPLICAS_DCUTS = "Basic_runcard_3replicas_diffcuts_230221"
FIT_ITERATED = "NNPDF40_nnlo_low_precision_240916_iterated"
PSEUDODATA_FIT = "pseudodata_test_fit_n3fit_240916"


base_config = dict(pdf=PDF, use_cuts='nocuts', dataset_inputs=DATA, theoryid=THEORYID_NEW, Q=10)


@pytest.fixture(scope='module')
def data_config():
    return base_config


@pytest.fixture(scope='module')
def thcovmat_config(data_config):
    """Same as data_config but with additional info for the thcovmat production."""
    new_config = dict(data_config)
    new_config["point_prescription"] = "3 point"
    new_config["use_theorycovmat"] = "true"
    new_config["use_cuts"] = "internal"
    new_config.update(theoryid=708)
    new_config["theoryids"] = {"from_": "scale_variation_theories"}
    new_config.update(dataset_inputs=DATA_THCOVMAT)
    return new_config


@pytest.fixture(scope='module')
def data_internal_cuts_config(data_config):
    config_dict = dict(data_config)
    config_dict.update(use_cuts='internal')
    return config_dict


@pytest.fixture(scope='module')
def data_internal_cuts_new_theory_config(data_internal_cuts_config):
    config = dict(data_internal_cuts_config)
    config["theoryid"] = THEORYID_NEW
    return config


@pytest.fixture(scope='module')
def data_fromfit_cuts_config(data_internal_cuts_new_theory_config):
    config = dict(data_internal_cuts_new_theory_config)
    config.update(use_cuts="fromfit")
    return config


@pytest.fixture(scope='module')
def single_data_internal_cuts_config(data_internal_cuts_config):
    """Like data_internal_cuts_config but for a single dataset"""
    config_dict = dict(data_internal_cuts_config)
    config_dict.pop("dataset_inputs")
    config_dict.update(dataset_input=DATA[0])
    return config_dict


@pytest.fixture(scope='module')
def single_data_categorical_internal_cuts_config(data_internal_cuts_config):
    """Test dataset with categorical plotting variables"""
    return {
        **data_internal_cuts_config,
        'dataset_input': SINGLE_CATEGORICAL,
        'theoryid': THEORYID_NEW,
    }


@pytest.fixture(scope='module')
def single_data_single_point_internal_cuts_config(single_data_internal_cuts_config):
    config_dict = dict(single_data_internal_cuts_config)
    config_dict.update(dataset_input=SINGLE_DATAPOINT)
    return config_dict


@pytest.fixture(scope='module')
def data_witht0_config():
    config_dict = dict(**base_config, use_t0=True, t0pdfset=PDF)
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


@pytest.fixture(scope='module')
def fromfit_closure_config():
    """A configuration useful for closure test where everything is
    read from the fit"""
    config = {
        "dataset_inputs": {"from_": "fit"},
        "datacuts": {"from_": "fit"},
        "use_cuts": "fromfit",
        "fakepdf": {"from_": "closuretest"},
        "theory": {"from_": "fit"},
        "theoryid": {"from_": "theory"},
        "pdf": {"from_": "fit"},
        "closuretest": {"from_": "fit"},
        "filterseed": {"from_": "closuretest"},
        "use_fitcommondata": True,
        "use_t0": True,
        "t0pdfset": {"from_": "datacuts"},
    }
    return config


def pytest_runtest_setup(item):
    ALL = {"darwin", "linux"}
    supported_platforms = ALL.intersection(mark.name for mark in item.iter_markers())
    plat = sys.platform
    if supported_platforms and plat not in supported_platforms:
        pytest.skip(f"cannot run on platform {plat}")


def pytest_configure(config):
    config.addinivalue_line("markers", "linux: mark test to run only on linux")


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


@contextlib.contextmanager
def temp_lhapdf_path(folder):
    """Modify the data path for LHAPDF sets"""
    oldpaths = lhapdf.paths()
    lhapdf.setPaths([str(folder)])
    try:
        yield
    finally:
        lhapdf.setPaths(oldpaths)
