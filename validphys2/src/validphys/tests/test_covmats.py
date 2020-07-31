"""
test_covmats.py

Tests related to the computation of the covariance matrix and its derivatives
"""
import pytest

import numpy as np
import random

from validphys.loader import Loader
from validphys.tests.conftest import THEORYID
from validphys.api import API
from validphys.covmats import sqrt_covmat


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
