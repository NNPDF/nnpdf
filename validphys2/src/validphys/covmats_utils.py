"""
covmat_utils.py

Utils functions for constructing covariance matrices from systematics.
Leveraged by :py:mod:`validphys.covmats` which contains relevant
actions/providers.

"""
import numpy as np
import pandas as pd
import scipy.linalg as la

def systematics_matrix(stat_errors: np.array, sys_errors: pd.DataFrame):
    """Basic function to create a systematics matrix , :math:`A`, such that:

    .. math::

        C = A A^T

    Where :math:`C` is the covariance matrix. This is achieved by creating a
    block diagonal matrix by adding the uncorrelated systematics in quadrature
    then taking the square-root and concatenating the correlated systematics,
    schematically:

    .. code::python

        A = concat([diag(sqrt(A_uncorr.sum(axis=1))), A_corr])

    Parameters
    ----------
    stat_errors: np.array
        a 1-D array of statistical uncertainties
    sys_errors: pd.DataFrame
        a dataframe with shape (N_data * N_sys) and systematic name as the
        column headers. The uncertainties should be in the same units as the
        data.

    Notes
    -----
    This function doesn't contain any logic to ignore certain contributions to
    the covmat, if you wanted to not include a particular systematic/set of
    systematics i.e all uncertainties with MULT errors, then filter those out
    of ``sys_errors`` before passing that to this function.

    """
    diagonal = stat_errors ** 2

    is_uncorr = sys_errors.columns.isin(("UNCORR", "THEORYUNCORR"))
    diagonal += (sys_errors.loc[:, is_uncorr].to_numpy() ** 2).sum(axis=1)

    corr_sys_mat = sys_errors.loc[:, ~is_uncorr].to_numpy()
    return np.concatenate((np.diag(np.sqrt(diagonal)), corr_sys_mat), axis=1)


def construct_covmat(stat_errors: np.array, sys_errors: pd.DataFrame):
    """Basic function to construct a covariance matrix (covmat), given the
    statistical error and a dataframe of systematics.

    Errors with name UNCORR or THEORYUNCORR are added in quadrature with
    the statistical error to the diagonal of the covmat.

    Other systematics are treated as correlated; their covmat contribution is
    found by multiplying them by their transpose.

    Parameters
    ----------
    stat_errors: np.array
        a 1-D array of statistical uncertainties
    sys_errors: pd.DataFrame
        a dataframe with shape (N_data * N_sys) and systematic name as the
        column headers. The uncertainties should be in the same units as the
        data.

    Notes
    -----
    This function doesn't contain any logic to ignore certain contributions to
    the covmat, if you wanted to not include a particular systematic/set of
    systematics i.e all uncertainties with MULT errors, then filter those out
    of ``sys_errors`` before passing that to this function.

    """
    diagonal = stat_errors ** 2

    is_uncorr = sys_errors.columns.isin(("UNCORR", "THEORYUNCORR"))
    diagonal += (sys_errors.loc[:, is_uncorr].to_numpy() ** 2).sum(axis=1)

    corr_sys_mat = sys_errors.loc[:, ~is_uncorr].to_numpy()
    return np.diag(diagonal) + corr_sys_mat @ corr_sys_mat.T

def diagonalize_sampling_covmat(
    dataset_inputs_sampling_covmat
):
    return la.eig(dataset_inputs_sampling_covmat)