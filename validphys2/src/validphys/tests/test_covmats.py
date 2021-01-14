"""
test_covmats.py

Tests related to the computation of the covariance matrix and its derivatives
"""
import random

import pytest

import numpy as np


from validphys.api import API
from validphys.covmats import (
    sqrt_covmat,
)
from validphys.loader import Loader
from validphys.tests.conftest import THEORYID


def test_covmat_from_systematics_correlated(data_with_correlations_config):
    """Test the covariance matrix generation from a set of correlated datasets
    given their systematic errors
    """
    covmat = API.dataset_inputs_experimental_covmat(**data_with_correlations_config)
    cpp_covmat = API.groups_covmat(**data_with_correlations_config)
    np.testing.assert_allclose(cpp_covmat, covmat)


def test_self_consistent_covmat_from_systematics(data_internal_cuts_config):
    """Test which checks that the single dataset implementation of
    ``covmat_from_systematics`` matches ``datasets_covmat_from_systematics``
    when the latter is given a list containing a single dataset.

    """
    base_config = dict(data_internal_cuts_config)
    #TODO: update tests to use new input.
    exps = base_config.pop("experiments")
    dataset_inputs = [dsinp for exp in exps for dsinp in exp["datasets"]]

    for dsinp in dataset_inputs:
        covmat_a = API.experimental_covmat(
            **base_config, dataset_input=dsinp)
        covmat_b = API.dataset_inputs_experimental_covmat(
            **base_config, dataset_inputs=[dsinp])
        np.testing.assert_allclose(covmat_a, covmat_b)


def test_covmat_from_systematics(data_internal_cuts_config):
    """Test which checks the python computation of the covmat relating to a
    collection of datasets matches that of the C++ computation. Note that the
    datasets are cut using the internal rules, but the datasets are not correlated.
    """
    covmat = API.dataset_inputs_experimental_covmat(**data_internal_cuts_config)
    cpp_covmat = API.groups_covmat(**data_internal_cuts_config)

    np.testing.assert_allclose(cpp_covmat, covmat)

def test_covmat_with_one_systematic():
    """Test that a dataset with 1 systematic successfully builds covmat, and
    that it agrees with cpp code. This special case can break the covmat
    construction in python because of pandas indexing.

    """
    dsinput = {"dataset": "D0ZRAP", "frac": 1.0, "cfac": ["QCD"]}
    config = dict(dataset_input=dsinput, theoryid=THEORYID, use_cuts="nocuts")

    covmat = API.experimental_covmat(**config)
    ds = API.dataset(**config)
    cpp_covmat = ds.load().get_covmat()

    np.testing.assert_allclose(cpp_covmat, covmat)


def test_cpp_sqrtcovmat():
    """Test that the sqrt of the covariance matrix is computed correctly for a
    random sample of 10 datasets. This uses the get methods of a loaded dataset
    which currently call the C++ code. This therefore currently tests the
    computation of the sqrt of the covariance matrix in the C++ code. In time
    the get_sqrtcovmat method should call the python code, in which case this
    test can be merged with :py:func:`test_sqrt_covmat`.

    """
    l = Loader()
    # Only test 10 datasets to avoid test taking too long
    datasets = random.sample(l.available_datasets, 10)
    cuts = (None, "internal")

    for ds_name in datasets:
        try:
            for cut in cuts:
                ds = l.check_dataset(ds_name, theoryid=THEORYID, cuts=cut)
                ds_ld = ds.load()
                sqrt_cov = ds_ld.get_sqrtcovmat()
                assert np.allclose(sqrt_cov @ sqrt_cov.T, ds_ld.get_covmat())
        except FileNotFoundError:
            continue


def test_sqrt_covmat(data_config):
    """In contrast to :py:func:`test_cpp_sqrtcovmat` this tests the python
    implementation of the sqrt of the covariance matrix, namely
    :py:func:`validphys.covmats.sqrt_covmat`.

    """
    rectangular_covmat = np.random.randint(10, size=(4, 5))

    with pytest.raises(ValueError):
        # Check whether ValueError is raised for a
        # rectangular covmat matrix
        sqrt_covmat(rectangular_covmat)

        # Check whether an empty covmat input raises
        # a ValueError
        sqrt_covmat(np.array([]))

    exps = API.experiments_data(**data_config)

    for exp in exps:
        ld_exp = exp.load()
        covmat = ld_exp.get_covmat()
        cholesky_cov = sqrt_covmat(covmat)
        np.testing.assert_allclose(cholesky_cov @ cholesky_cov.T, covmat)

def test_t0_covmat(data_witht0_config):
    #TODO: expand t0 tests
    base_config = dict(data_witht0_config)
    # just compare cut data
    base_config["use_cuts"] = "internal"
    covmat = API.dataset_inputs_t0_covmat(**base_config)
    cpp_covmat = API.groups_covmat(**base_config)
    # use allclose defaults or it fails
    np.testing.assert_allclose(cpp_covmat, covmat, rtol=1e-05, atol=1e-08)
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(
            covmat, API.dataset_inputs_experimental_covmat(**base_config)
        )
