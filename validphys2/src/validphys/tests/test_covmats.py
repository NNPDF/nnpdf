import pytest

import numpy as np

from validphys.api import API
from validphys.covmats import sqrt_covmat


def test_sqrt_covmat(data_config):
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
