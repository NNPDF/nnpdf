"""
test_sqrtcovmat.py

Tests that the sqrt of the covariance matrix is computed correctly
"""
import numpy as np
import random

from validphys.loader import Loader
from validphys.tests.conftest import THEORYID


def test_sqrtcovmat():
    """Test that the sqrt of the covariance matrix is computed correctly for a
    random sample of 10 datasets"""
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
