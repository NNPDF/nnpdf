"""
calcutils.py

Low level utilities to calculate χ² and such. These are used to implement the
higher level functions in results.py
"""
import logging
from typing import Callable

import numpy as np
import pandas as pd
import scipy.linalg as la

log = logging.getLogger(__name__)


def calc_chi2(sqrtcov, diffs):
    """Elementary function to compute the chi², given a Cholesky decomposed
    lower triangular part and a vector of differences.

    Parameters
    ----------
    sqrtcov : matrix
        A lower tringular matrix corresponding to the lower part of
        the Cholesky decomposition of the covariance matrix.
    diffs : array
        A vector of differences (e.g. between data and theory).
        The first dimenssion must match the shape of `sqrtcov`.
        The computation will be broadcast over the other dimensions.

    Returns
    -------
    chi2 : array
        The result of the χ² for each vector of differences.
        Will have the same shape as ``diffs.shape[1:]``.

    Notes
    -----
    This function computes the χ² more efficiently and accurately than
    following the direct definition of inverting the covariance matrix,
    :math:`\chi^2 = d\Sigma^{-1}d`,
    by solving the triangular linear system instead.

    Examples
    --------

    >>> from validphys.calcutils import calc_chi2
    >>> import numpy as np
    >>> import scipy.linalg as la
    >>> np.random.seed(0)
    >>> diffs = np.random.rand(10)
    >>> s = np.random.rand(10,10)
    >>> cov = s@s.T
    >>> calc_chi2(la.cholesky(cov, lower=True), diffs)
    44.64401691354948
    >>> diffs@la.inv(cov)@diffs
    44.64401691354948

    """
    # handle empty data
    if not diffs.size:
        return np.full(diffs.shape[1:], np.nan)
    # Note la.cho_solve doesn't really improve things here
    # NOTE: Do not enable check_finite. The upper triangular part is not
    # guaranteed to make any sense.
    vec = la.solve_triangular(sqrtcov, diffs, lower=True, check_finite=False)
    # This sums up the result for the chi² for any input shape.
    # Sum the squares over the first dimension and leave the others alone
    return np.einsum('i...,i...->...', vec, vec)


def all_chi2(results):
    """Return the chi² for all elements in the result, regardless of the Stats class
    Note that the interpretation of the result will depend on the PDF error type"""
    data_result, th_result = results
    diffs = th_result.rawdata - data_result.central_value[:, np.newaxis]
    return calc_chi2(sqrtcov=data_result.sqrtcovmat, diffs=diffs)


def central_chi2(results):
    """Calculate the chi² from the central value of the theory prediction to
    the data"""
    data_result, th_result = results
    central_diff = th_result.central_value - data_result.central_value
    return calc_chi2(data_result.sqrtcovmat, central_diff)


def all_chi2_theory(results, totcov):
    """Like all_chi2 but here the chi² are calculated using a covariance matrix
    that is the sum of the experimental covmat and the theory covmat."""
    data_result, th_result = results
    diffs = th_result.rawdata - data_result.central_value[:, np.newaxis]
    total_covmat = np.array(totcov)
    return calc_chi2(sqrtcov=la.cholesky(total_covmat, lower=True), diffs=diffs)


def central_chi2_theory(results, totcov):
    """Like central_chi2 but here the chi² is calculated using a covariance matrix
    that is the sum of the experimental covmat and the theory covmat."""
    data_result, th_result = results
    central_diff = th_result.central_value - data_result.central_value
    total_covmat = np.array(totcov)
    return calc_chi2(la.cholesky(total_covmat, lower=True), central_diff)


def calc_phi(sqrtcov, diffs):
    """Low level function which calculates phi given a Cholesky decomposed
    lower triangular part and a vector of differences. Primarily used when phi
    is to be calculated independently from chi2.

    The vector of differences `diffs` is expected to have N_bins on the first
    axis
    """
    diffs = np.array(diffs)
    return np.sqrt(
        (np.mean(calc_chi2(sqrtcov, diffs), axis=0) - calc_chi2(sqrtcov, diffs.mean(axis=1)))
        / diffs.shape[0]
    )


def bootstrap_values(
    data, nresamples, *, boot_seed: int = None, apply_func: Callable = None, args=None
):
    """General bootstrap sample

    `data` is the data which is to be sampled, replicas is assumed to
    be on the final axis e.g N_bins*N_replicas

    `boot_seed` can be specified if the user wishes to be able to
    take exact same bootstrap samples multiple times, as default it is
    set as None, in which case a random seed is used.

    If just `data` and `nresamples` is provided, then `bootstrap_values`
    creates N resamples of the data, where each resample is a Monte Carlo
    selection of the data across replicas. The mean of each resample is
    returned

    Alternatively, the user can specify a function to be sampled `apply_func`
    plus any additional arguments required by that function.
    `bootstrap_values` then returns `apply_func(bootstrap_data, *args)`
    where `bootstrap_data.shape = (data.shape, nresamples)`. It is
    critical that `apply_func` can handle data input in this format.
    """
    data = np.atleast_2d(data)
    N_reps = data.shape[-1]
    bootstrap_data = data[
        ..., np.random.RandomState(boot_seed).randint(N_reps, size=(N_reps, nresamples))
    ]
    if apply_func is None:
        return np.mean(bootstrap_data, axis=-2)
    else:
        return apply_func(bootstrap_data, *args)


def get_df_block(matrix: pd.DataFrame, key: str, level):
    """Given a pandas dataframe whose index and column keys match, and data represents a symmetric
    matrix returns a diagonal block of this matrix corresponding to `matrix`[key`, key`] as a numpy
    array

    addtitionally, the user can specify the `level` of the key for which the cross section is being
    taken, by default it is set to 1 which corresponds to the dataset level of a theory covariance
    matrix
    """
    block = matrix.xs(key, level=level, axis=0).xs(key, level=level, axis=1).values
    return block


def regularize_covmat(covmat: np.array, norm_threshold=4):
    """Given a covariance matrix, performs a regularization which is equivalent
    to performing `regularize_l2` on the sqrt of `covmat`: the l2 norm of
    the inverse of the correlation matrix calculated from `covmat` is set to be
    less than or equal to `norm_threshold`. If the input covmat already fulfills
    this criterion it is returned.

    Parameters
    ----------
    covmat : array
        a covariance matrix which is to be regularized.
    norm_threshold : float
        The acceptable l2 norm of the sqrt correlation matrix, by default
        set to 4.

    Returns
    -------
    new_covmat : array
        A new covariance matrix which has been regularized according to
        prescription above.

    """
    # square up threshold since we have cov not sqrtcov
    sqr_threshold = norm_threshold**2
    d = np.sqrt(np.diag(covmat))[:, np.newaxis]
    corr = covmat / d / d.T
    # eigh gives eigenvals in ascending order
    e_val, e_vec = la.eigh(corr)
    # if eigenvalues are close to zero, can be negative
    if e_val[0] < 0:
        log.warning(
            "Negative eigenvalue encountered in correlation matrix: %s. "
            "Assuming eigenvalue should be zero and is negative due to numerical "
            "precision.",
            e_val[0],
        )
    if e_val[0] > 1 / sqr_threshold:
        return covmat
    new_e_val = np.clip(e_val, a_min=1 / sqr_threshold, a_max=None)
    return ((e_vec * new_e_val) @ e_vec.T) * d * d.T


def regularize_l2(sqrtcov, norm_threshold=4):
    r"""Return a regularized version of `sqrtcov`.

    Given `sqrtcov` an (N, nsys) matrix, such that it's
    gram matrix is the covariance matrix (`covmat = sqrtcov@sqrtcov.T`), first
    decompose it like ``sqrtcov = D@A``, where `D` is a positive diagonal matrix
    of standard deviations and `A` is the "square root" of the correlation
    matrix, ``corrmat = A@A.T``. Then produce a new version of `A` which removes
    the unstable behaviour and assemble a new square root covariance matrix,
    which is returned.

    The stability condition is controlled by `norm_threshold`. It is

    .. math::

        \left\Vert A^+ \right\Vert_{L2}
        \leq \frac{1}{\text{norm_threshold}}

    A+ is the pseudoinverse of A, `norm_threshold` roughly corresponds to the
    sqrt of the maximimum relative uncertainty in any systematic.

    Parameters
    ----------

    sqrtcov : 2d array
        An (N, nsys) matrix specifying the uncertainties.
    norm_threshold : float
        The tolerance for the regularization.

    Returns
    -------

    newsqrtcov : 2d array
        A regularized version of `sqrtcov`.
    """

    d = np.sqrt(np.sum(sqrtcov**2, axis=1))[:, np.newaxis]
    sqrtcorr = sqrtcov / d
    u, s, vt = la.svd(sqrtcorr, full_matrices=False)
    if 1 / s[-1] <= norm_threshold:
        return sqrtcov
    snew = np.clip(s, a_min=1 / norm_threshold, a_max=None)
    return u * (snew * d) @ vt
