"""
test_covmats.py

Tests related to the computation of the covariance matrix and its derivatives
"""
import pathlib

import pytest

import numpy as np
import random


from validphys.api import API
from validphys.commondataparser import load_commondata
from validphys.covmats import (
    sqrt_covmat,
    datasets_covmat_from_systematics
)
from validphys.loader import Loader
from validphys.tests.conftest import THEORYID

REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regressions")


def test_covmat_from_systematics_correlated():
    """Test the covariance matrix generation from a set of correlated datasets
    given their systematic errors
    """
    l = Loader()
    correlated_datasets = [
        "ATLASWZRAP36PB",
        "ATLASZHIGHMASS49FB",
        "ATLASLOMASSDY11EXT",
        "ATLASWZRAP11",
    ]

    cds = list(map(l.check_commondata, correlated_datasets))
    ld_cds = list(map(load_commondata, cds))

    covmat = datasets_covmat_from_systematics(ld_cds)

    dss = [
        l.check_dataset(i, theoryid=THEORYID, cuts=None) for i in correlated_datasets
    ]
    exp = l.check_experiment("null", dss)
    ld_exp = exp.load()

    np.testing.assert_allclose(ld_exp.get_covmat(), covmat)


def test_covmat_from_systematics(data_config):
    """Test which checks the python computation of the covmat relating to a
    collection of datasets matches that of the C++ computation. Note that the
    datasets are cut using the internal rules, but the datasets are not correlated.
    """
    l = Loader()
    exps = API.experiments(**data_config)

    cds = [l.check_commondata(i.name) for exp in exps for i in exp.datasets]
    ld_cds = list(map(load_commondata, cds))

    internal_cuts = [
        l.check_internal_cuts(cd, API.rules(theoryid=THEORYID, use_cuts="internal"))
        for cd in cds
    ]
    cut_ld_cds = list(map(lambda x: x[0].with_cuts(x[1]), zip(ld_cds, internal_cuts)))

    covmat = datasets_covmat_from_systematics(cut_ld_cds)

    dss = [
        l.check_dataset(i.name, theoryid=THEORYID, cuts="internal")
        for exp in exps
        for i in exp.datasets
    ]
    exp = l.check_experiment("null", dss)
    ld = exp.load()
    cpp_covmat = ld.get_covmat()

    assert np.allclose(cpp_covmat, covmat)


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

    exps = API.experiments(**data_config)

    for exp in exps:
        ld_exp = exp.load()
        covmat = ld_exp.get_covmat()
        cholesky_cov = sqrt_covmat(covmat)
        np.testing.assert_allclose(cholesky_cov @ cholesky_cov.T, covmat)
