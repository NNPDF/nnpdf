"""
Module contains general utility functions for general use in implementing datasets.
"""

import numpy as np


def correlation_to_covariance(correlation, uncertainties):
    """
    Converts a correlation matrix into a covariance matrix
    using a list of uncertainties.

    Parameters:
    -----------
    correlation : np.ndarray
        A square matrix of correlations.
    uncertainties : np.ndarray
        A 1D array of uncertainties.

    Returns:
    --------
    np.ndarray
        The corresponding covariance matrix.
    """
    covariance = np.outer(uncertainties, uncertainties) * correlation
    return covariance


def decompose_covmat(covmat):
    """Given a covmat it return an array sys with shape (ndat,ndat)
    giving ndat correlated systematics for each of the ndat point.
    The original covmat is obtained by doing sys@sys.T"""

    lamb, mat = np.linalg.eig(covmat)
    sys = np.multiply(np.sqrt(lamb), mat)
    return sys
