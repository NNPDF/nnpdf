"""
test_covmats.py

Tests related to the computation of the covariance matrix and its derivatives
"""
import random

import pytest

import numpy as np


from validphys.api import API
from validphys.commondataparser import load_commondata
from validphys.covmats import (
    sqrt_covmat,
    datasets_covmat_from_systematics,
    covmat_from_systematics
)
from validphys.loader import Loader
from validphys.tests.conftest import THEORYID


def test_covmat_from_systematics_correlated(data_with_correlations_config):
    """Test the covariance matrix generation from a set of correlated datasets
    given their systematic errors
    """
    data = API.data(**data_with_correlations_config)
    cds = [ds.commondata for ds in data.datasets]

    ld_cds = list(map(load_commondata, cds))

    covmat = datasets_covmat_from_systematics(ld_cds)

    cpp_covmat = API.groups_covmat(**data_with_correlations_config)

    np.testing.assert_allclose(cpp_covmat, covmat)


def test_self_consistent_covmat_from_systematics(data_internal_cuts_config):
    """Test which checks that the single dataset implementation of
    ``covmat_from_systematics`` matches ``datasets_covmat_from_systematics``
    when the latter is given a list containing a single dataset.

    """
    data = API.data(**data_internal_cuts_config)
    cds = [ds.commondata for ds in data.datasets]

    ld_cds = list(map(load_commondata, cds))

    internal_cuts = [ds.cuts for ds in data.datasets]
    cut_ld_cds = list(map(lambda x: x[0].with_cuts(x[1]), zip(ld_cds, internal_cuts)))

    for cut_ld_cd in cut_ld_cds:
        covmat_a = covmat_from_systematics(cut_ld_cd)
        covmat_b = datasets_covmat_from_systematics([cut_ld_cd])
        np.testing.assert_allclose(covmat_a, covmat_b)


def test_covmat_from_systematics(data_internal_cuts_config):
    """Test which checks the python computation of the covmat relating to a
    collection of datasets matches that of the C++ computation. Note that the
    datasets are cut using the internal rules, but the datasets are not correlated.
    """
    data = API.data(**data_internal_cuts_config)
    cds = [ds.commondata for ds in data.datasets]

    ld_cds = list(map(load_commondata, cds))

    internal_cuts = [ds.cuts for ds in data.datasets]
    cut_ld_cds = list(map(lambda x: x[0].with_cuts(x[1]), zip(ld_cds, internal_cuts)))

    covmat = datasets_covmat_from_systematics(cut_ld_cds)

    cpp_covmat = API.groups_covmat(**data_internal_cuts_config)

    np.testing.assert_allclose(cpp_covmat, covmat)

def test_covmat_with_one_systematic():
    """Test that a dataset with 1 systematic successfully builds covmat, and
    that it agrees with cpp code. This special case can break the covmat
    construction in python because of pandas indexing.

    """
    dsinput = {"dataset": "D0ZRAP", "frac": 1.0, "cfac": ["QCD"]}
    cd = API.commondata(dataset_input=dsinput)
    l_cd = load_commondata(cd)
    covmat = covmat_from_systematics(l_cd)

    ds = API.dataset(dataset_input=dsinput, theoryid=THEORYID, use_cuts="nocuts")
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
