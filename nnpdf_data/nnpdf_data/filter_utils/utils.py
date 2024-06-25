"""
Module contains general utility functions for general use in implementing datasets.
"""

import numpy as np
import yaml


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


def prettify_float(dumper, value):
    """
    Override the default yaml representer:
    https://github.com/yaml/pyyaml/blob/48838a3c768e3d1bcab44197d800145cfd0719d6/lib/yaml/representer.py#L189

    This function is used to prettify the float representation in the yaml file.
    If the float has more than 8 digits, it will be represented in scientific notation with 8 digits.

    Note:
    -----
    When importing yaml in a module,

    yaml.add_representer(float, prettify_float)

    must be called to use this function.
    """

    ret = dumper.represent_float(value)
    if len(ret.value) > 8:
        ret_str = f"{value:.8e}"
        ret = dumper.represent_scalar('tag:yaml.org,2002:float', ret_str)
    return ret
