from hypothesis import given
from hypothesis.extra.numpy import arrays
from hypothesis.strategies import floats
import numpy as np
import scipy.linalg as la

from validphys import calcutils

sane_floats = floats(min_value=-1, max_value=1, allow_nan=False, allow_infinity=False)
diffs = arrays(dtype=float, shape=10, elements=sane_floats)
sqrtcov = arrays(dtype=float, shape=(10, 10), elements=sane_floats)


@given(sqrtcov, diffs)
def test_calc_chi2(s, d):
    cov = s @ s.T
    # Handle zero matrices and so on
    np.fill_diagonal(cov, np.diag(cov) + 1)
    chi2 = d @ la.inv(cov) @ d
    chol = la.cholesky(cov, lower=True)
    calc = calcutils.calc_chi2(chol, d)
    assert np.allclose(chi2, calc)
    dd = np.repeat(d, 5).reshape(len(d), 5)
    calcdd = calcutils.calc_chi2(chol, dd)
    assert np.allclose(chi2, calcdd)
