"""
test_covmats.py

Tests related to the computation of the covariance matrix and its derivatives
"""
import numpy as np
import pytest

from validphys.api import API
from validphys.commondataparser import load_commondata
from validphys.covmats import dataset_t0_predictions, sqrt_covmat
from validphys.tests.conftest import DATA, HESSIAN_PDF, PDF, THEORYID

# Experiments which have non trivial correlations between their datasets
CORR_DATA = [
    {'dataset': 'ATLASWZRAP36PB', 'cfac': ['QCD']},
    {'dataset': 'ATLASZHIGHMASS49FB', 'cfac': ['QCD']},
    {'dataset': 'ATLASLOMASSDY11EXT', 'cfac': ['QCD']},
    {'dataset': 'ATLASWZRAP11', 'frac': 0.5, 'cfac': ['QCD']},
    {'dataset': 'CMSZDIFF12', 'cfac': ('QCD', 'NRM'), 'sys': 10},
    {'dataset': 'CMSJETS11', 'frac': 0.5, 'sys': 10},
]


def test_self_consistent_covmat_from_systematics(data_internal_cuts_config):
    """Test which checks that the single dataset implementation of
    ``covmat_from_systematics`` matches ``dataset_inputs_covmat_from_systematics``
    when the latter is given a list containing a single dataset.

    """
    base_config = dict(data_internal_cuts_config)
    dataset_inputs = base_config.pop("dataset_inputs")

    for dsinp in dataset_inputs:
        covmat_a = API.covmat_from_systematics(**base_config, dataset_input=dsinp)
        covmat_b = API.dataset_inputs_covmat_from_systematics(**base_config, dataset_inputs=[dsinp])
        np.testing.assert_allclose(covmat_a, covmat_b)


@pytest.mark.parametrize("use_cuts", ["nocuts", "internal"])
@pytest.mark.parametrize("dataset_inputs", [DATA, CORR_DATA])
def test_covmat_from_systematics(data_config, use_cuts, dataset_inputs):
    """Test which checks the python computation of the covmat relating to a
    collection of datasets from dataset_inputs matches the direct call to groups_covmat

    Tests all combinations of internal/no cuts and correlated/uncorrelated data.
    """
    config = dict(data_config)
    config["use_cuts"] = use_cuts
    config["dataset_inputs"] = dataset_inputs

    covmat = API.dataset_inputs_covmat_from_systematics(**config)
    another_covmat = API.groups_covmat(**config)

    np.testing.assert_allclose(another_covmat, covmat)


def test_covmat_with_one_systematic():
    """Test that a dataset with 1 systematic successfully builds covmat.
    This special case can break the covmat construction in python because of pandas indexing.
    """
    dsinput = {"dataset": "D0ZRAP", "frac": 1.0, "cfac": ["QCD"]}
    config = dict(dataset_input=dsinput, theoryid=THEORYID, use_cuts="nocuts")

    ds = API.dataset(**config)
    # double check that the dataset does indeed only have 1 systematic.
    assert ds.commondata.nsys == 1

    # Test the covmat can be constructed
    _ = API.covmat_from_systematics(**config)


def test_sqrt_covmat(data_config):
    """Tests the python implementation of the sqrt of the covariance matrix, namely
    :py:func:`validphys.covmats.sqrt_covmat`.
    """
    rectangular_covmat = np.random.randint(10, size=(4, 5))

    with pytest.raises(ValueError):
        # Check whether ValueError is raised for a
        # rectangular covmat matrix
        sqrt_covmat(rectangular_covmat)

    # Test for empty covmat
    assert sqrt_covmat(np.array([])).size == 0

    config = dict(data_config)
    config["dataset_inputs"] = CORR_DATA
    covmat = API.dataset_inputs_covmat_from_systematics(**config)

    cholesky_cov = sqrt_covmat(covmat)
    np.testing.assert_allclose(cholesky_cov @ cholesky_cov.T, covmat)


@pytest.mark.parametrize("t0pdfset", [PDF, HESSIAN_PDF])
@pytest.mark.parametrize("dataset_inputs", [DATA, CORR_DATA])
def test_python_t0_covmat_matches_variations(data_internal_cuts_config, t0pdfset, dataset_inputs):
    """Test which checks the python computation of the t0 covmat relating to a
    collection of datasets

    Tests all combinations of hessian/MC t0pdfset and correlated/uncorrelated
    data.

    """
    config = dict(data_internal_cuts_config)
    config["dataset_inputs"] = dataset_inputs
    config["t0pdfset"] = t0pdfset
    config["use_t0"] = True
    covmat = API.dataset_inputs_t0_covmat_from_systematics(**config)
    another_covmat = API.groups_covmat(**config)
    # use allclose defaults or it fails
    np.testing.assert_allclose(another_covmat, covmat, rtol=1e-05, atol=1e-08)
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(covmat, API.dataset_inputs_covmat_from_systematics(**config))


@pytest.mark.parametrize("use_cuts", ["nocuts", "internal"])
@pytest.mark.parametrize("dataset_input", DATA)
def test_systematic_matrix(data_config, use_cuts, dataset_input):
    """Test which checks the python computation of the t0 covmat relating to a
    collection of datasets is equivalent using different functions

    Tests all combinations of hessian/MC t0pdfset and correlated/uncorrelated
    data.

    """
    config = dict(data_config)
    config["dataset_input"] = dataset_input
    config["use_cuts"] = use_cuts
    covmat = API.covmat_from_systematics(**config)
    sys_mat = API.systematics_matrix_from_commondata(**config)
    covmat_from_sys_mat = sys_mat @ sys_mat.T
    np.testing.assert_allclose(covmat_from_sys_mat, covmat)


def test_single_datapoint(single_data_single_point_internal_cuts_config):
    # Make the t0 predictions
    ds = API.dataset(**single_data_single_point_internal_cuts_config)
    t0set = API.pdf(**single_data_single_point_internal_cuts_config)
    t0_predictions = dataset_t0_predictions(ds, t0set)

    cd = API.commondata(**single_data_single_point_internal_cuts_config)
    ld = load_commondata(cd)
    # Ensure the dataset is only a single datapoint
    assert ld.ndata == 1
    ld.systematic_errors(t0_predictions)
