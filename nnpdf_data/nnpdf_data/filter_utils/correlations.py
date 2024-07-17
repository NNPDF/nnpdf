import numpy as np
from numpy.linalg import eig


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


def covmat_to_artunc(ndata, covmat_list, no_of_norm_mat=0):
    r"""Convert the covariance matrix to a matrix of
    artificial uncertainties.

    NOTE: This function has been taken from validphys.newcommondata_utils.
    If those utils get merged in the future, we can replace this.

    Parameters
    ----------
    ndata : integer
        Number of data points
    covmat_list : list
        A one dimensional list which contains the elements of
        the covariance matrix row by row. Since experimental
        datasets provide these matrices in a list form, this
        simplifies the implementation for the user.
    no_of_norm_mat : int
        Normalized covariance matrices may have an eigenvalue
        of 0 due to the last data point not being linearly
        independent. To allow for this, the user should input
        the number of normalized matrices that are being treated
        in an instance. For example, if a single covariance matrix
        of a normalized distribution is being processed, the input
        would be 1. If a covariance matrix contains pertains to
        3 normalized datasets (i.e. cross covmat for 3
        distributions), the input would be 3. The default value is
        0 for when the covariance matrix pertains to an absolute
        distribution.

    Returns
    -------
    artunc : list
        A two dimensional matrix (given as a list of lists)
        which contains artificial uncertainties to be added
        to the commondata. i^th row (or list) contains the
        artificial uncertainties of the i^th data point.

    """
    epsilon = -0.0000000001
    neg_eval_count = 0
    psd_check = True
    covmat = np.zeros((ndata, ndata))
    artunc = np.zeros((ndata, ndata))
    for i in range(len(covmat_list)):
        a = i // ndata
        b = i % ndata
        covmat[a][b] = covmat_list[i]
    eigval, eigvec = eig(covmat)
    for j in range(len(eigval)):
        if eigval[j] < epsilon:
            psd_check = False
        elif eigval[j] > epsilon and eigval[j] <= 0:
            neg_eval_count = neg_eval_count + 1
            if neg_eval_count == (no_of_norm_mat + 1):
                psd_check = False
        elif eigval[j] > 0:
            continue
    if psd_check == False:
        raise ValueError("The covariance matrix is not positive-semidefinite")
    else:
        for i in range(ndata):
            for j in range(ndata):
                if eigval[j] < 0:
                    continue
                else:
                    artunc[i][j] = eigvec[i][j] * np.sqrt(eigval[j])
    return artunc.tolist()
