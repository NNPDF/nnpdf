import numpy as np
from numpy.linalg import eig
from nnpdf_data.filter_utils.utils import covmat_to_artunc


def upper_triangular_to_symmetric(ut, dim):
    """Build a symmetric matrix from the upper diagonal"""
    corr = np.zeros((dim, dim))
    last = dim
    first = 0
    for i in range(dim):
        corr[i, i:] = ut[first:last]
        last += dim - i - 1
        first += dim - i
    return corr


def compute_covmat(corrmat: np.ndarray, unc: np.ndarray, ndata: int) -> list:
    """Compute the covariance matrix with the artificial stat uncertainties."""
    # multiply by stat err
    cov_mat = np.einsum("i,ij,j->ij", unc, corrmat, unc)
    return covmat_to_artunc(ndata, cov_mat.flatten().tolist())
