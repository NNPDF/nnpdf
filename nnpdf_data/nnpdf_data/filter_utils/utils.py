"""
General Python utilities for commondata implementation.

This module provides helpful functions that automate a few
tasks that are regularly needed for the implementation of
experimental data to the commondata format. If there are
any additional functions that could be added here as they
could simplify some repetitve tasks, please do suggest.

Before the usage of any functions, it is recommended to read
the docstrings of the function to understand the inputs and
outputs.
"""

import copy
from math import sqrt
import os
import shutil

import numpy as np
from numpy.linalg import eig
import yaml


def symmetrize_errors(delta_plus, delta_minus):
    r"""Compute the symmterized uncertainty and the shift in data point.

    Parameters
    ----------
    delta_plus : float
        The top/plus uncertainty with sign
    delta_minus : float
        The bottom/minus uncertainty with sign

    Returns
    -------
    se_delta : float
        The value to be added to the data point
    se_sigma : float
        The symmetrized uncertainty to be used in commondata

    """
    semi_diff = (delta_plus + delta_minus) / 2
    average = (delta_plus - delta_minus) / 2
    se_delta = semi_diff
    se_sigma = sqrt(average * average + 2 * semi_diff * semi_diff)
    return se_delta, se_sigma


def percentage_to_absolute(percentage, value):
    r"""Compute the absolute value of uncertainty from percentage.

    Parameters
    ----------
    percentage : string/float
        Experimental datasets can provide the percentage
        uncertainties with a % sign or without one.
        The function will autostrip % sign and convert to
        a float type in case the percentage uncertainty
        comes with a % sign. Else, it will directly perform
        the computation.
    value : float
        The data point

    Returns
    -------
    absolute : float
        The absolute value of the uncertainty

    """
    if type(percentage) is str:
        percentage = float(percentage.replace("%", ""))
        absolute = percentage * value * 0.01
        return absolute
    else:
        absolute = percentage * value * 0.01
        return absolute


def cormat_to_covmat(err_list, cormat_list):
    r"""Convert correlation matrix elements to covariance
    matrix elements.

    Parameters
    ----------
    err_list : list
        A one dimensional list which contains the uncertainty
        associated to each data point in order.
    cormat_list : list
        A one dimensional list which contains the elements of
        the correlation matrix row by row. Since experimental
        datasets provide these matrices in a list form, this
        simplifies the implementation for the user.

    Returns
    -------
    covmat_list : list
        A one dimensional list which contains the elements of
        the covariance matrix row by row.

    """
    covmat_list = []
    for i in range(len(cormat_list)):
        a = i // len(err_list)
        b = i % len(err_list)
        covmat_list.append(cormat_list[i] * err_list[a] * err_list[b])
    return covmat_list


def covmat_to_artunc(ndata, covmat_list, no_of_norm_mat=0):
    r"""Convert the covariance matrix to a matrix of
    artificial uncertainties.

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
    if not psd_check:
        raise ValueError("The covariance matrix is not positive-semidefinite")
    else:
        for i in range(ndata):
            for j in range(ndata):
                if eigval[j] < 0:
                    continue
                else:
                    artunc[i][j] = eigvec[i][j] * sqrt(eigval[j])
    return artunc.tolist()


def cross_cormat_to_covmat(row_err_list, col_err_list, cormat_list):
    r"""Convert cross correlation matrix elements
    (i.e. those between different different variables or
    observables) to covariance matrix elements.

    Parameters
    ----------
    row_err_list : list
        A one dimensional list which contains the uncertainty
        associated to each data point of the variable that is
        given on the vertical axis.
    col_err_list : list
        A one dimensional list which contains the uncertainty
        associated to each data point of the variable that is
        given on the horizontal axis.
    cormat_list : list
        A one dimensional list which contains the elements of
        the correlation matrix row by row. Since experimental
        datasets provide these matrices in a list form, this
        simplifies the implementation for the user.

    Returns
    -------
    covmat_list : list
        A one dimensional list which contains the elements of
        the covariance matrix row by row.

    """
    covmat_list = []
    for i in range(len(cormat_list)):
        a = i // len(col_err_list)
        b = i % len(col_err_list)
        covmat_list.append(cormat_list[i] * row_err_list[a] * col_err_list[b])
    return covmat_list


def matlist_to_matrix(rows, columns, mat_list):
    r"""Convert a 1d list to a 2d matrix.

    Note: This utils function is not strictly needed for
    data implementation, however, it is provided for
    the aid of the user due to how matrices are treated
    throughout all the other functions. This function
    allows the user to convert a list that contains the
    elemnets of matrix row by row to a proper matrix, if
    need be for any reason.

    Parameters
    ----------
    rows : int
        No. of rows in the matrix
    columns : int
        No. of columns in the matrix
    mat_list : list
        A one dimensional list which contains the elements of
        the matrix row by row.

    Returns
    -------
    matrix : numpy.ndarray
        The matrix as a numpy 2d array.

    """
    if rows * columns == len(mat_list):
        matrix = np.zeros((rows, columns))
        for i in range(rows):
            for j in range(columns):
                matrix[i][j] = mat_list[j + i * columns]
        matrix = np.array(matrix)
        return matrix
    else:
        raise Exception("rows * columns != len(mat_list)")


def concat_matrices(rows, columns, list_of_matrices):
    r"""Join smaller matrices into a large matrix.

    This function aims to simplify the process of joining multiple
    smaller matrices into one large matrix. Such a need could arise,
    for instance, when cross variable covariance matrices are provided
    but the user needs to join all the matrices to generate the full
    covariance matrix corresponding to the entire dataset.

    Parameters
    ----------
    rows : int
        No. of rows of matrices to be concatenated. E.g., if 6
        matrices: A, B, C, D, E, F need to be joined as
        [[A, B, C],
        [D, E, F]],
        the number of rows would be 2.
    columns : int
        No. of columns of matrices to be concatenated. In the
        above example, this would be 3.
    list_of_matrices : list
        A list of the matrices that have to concatenated row by
        row. In the above example, this would be [A, B, C, D, E, F].
        The matrices themselves need to be provided as a list of lists,
        or a numpy 2d array. If the user has the matrix in a 1d row by
        row form, use matList_to_matrix() to convert it. It is assumed
        the user verifies that all the input matrices have the correct
        dimensions. Matrices with incompatible dimensions will lead to
        undesired behavior.

    Returns
    -------
    final_mat_list : list
        A one dimensional list which contains the elements of
        the final, fully concatenated matrix row by row.

    """
    for i in range(len(list_of_matrices)):
        list_of_matrices[i] = np.array(list_of_matrices[i])
    col_list = []
    for i in range(rows):
        row_list = []
        for j in range(columns):
            row_list.append(list_of_matrices[j + i * columns])
        col_list.append(np.concatenate(tuple(row_list), axis=1))
    final_mat = np.concatenate(tuple(col_list), axis=0)
    final_mat_list = []
    for i in range(len(final_mat)):
        for j in range(len(final_mat[i])):
            final_mat_list.append(final_mat[i][j])
    return final_mat_list


def trimat_to_fullmat(mode, tri_mat_list):
    r"""Convert a list of values of a triangular matrix
    to a symmetric matrix.

    Experimental datasets can provide the entries of
    correlation or covariance matrices as a triangular
    matrix, as these matrices are symmetric by their
    very nature. This function can convert these list to
    a complete symmetric matrix, that can be used for the
    dataset implementation.

    mode : bool
        Enter 0 or 1 based on the following scenarios:
        Use mode 0 if matrix entries are given row by
        row such as:
        0 1 2 3
          4 5 6
            7 8
              9
        Use mode 1 if the matrix entries are given column
        by column such as:
        0 1 3 6
          2 4 7
            5 8
              9
        Please note that the numbers above (0-9) are not
        entries of the matrix but rather the index of the
        entries of the list which contains the elements of
        the triangular matrix.
    tri_mat_list : list
        A list containing the elements of the triangular matrix,
        for example, for a 4*4 matrix, the list of
        triangular matrix entries could be:
        [a, b, c, d, e, f, g, h, i, j]

    Returns
    -------
    mat_list : list
        A one dimensional list which contains the elements of
        the fully populated, symmetric matrix row by row.

    """
    dim = int((np.sqrt(1 + 8 * len(tri_mat_list)) - 1) / 2)
    matrix = np.zeros((dim, dim))
    if mode == 0:
        for i in range(dim):
            for j in range(i + 1):
                list_el = len(tri_mat_list) - 1 - ((i * (i + 1)) // 2 + j)
                if i == j:
                    matrix[dim - 1 - i][dim - 1 - j] = tri_mat_list[list_el]
                else:
                    matrix[dim - 1 - i][dim - 1 - j] = tri_mat_list[list_el]
                    matrix[dim - 1 - j][dim - 1 - i] = tri_mat_list[list_el]
    elif mode == 1:
        for i in range(dim):
            for j in range(i + 1):
                list_el = (i * (i + 1)) // 2 + j
                if i == j:
                    matrix[i][j] = tri_mat_list[list_el]
                else:
                    matrix[i][j] = tri_mat_list[list_el]
                    matrix[j][i] = tri_mat_list[list_el]
    else:
        raise Exception("Mode should be 0 or 1, refer to the function for usage")
    mat_list = []
    for i in range(dim):
        for j in range(dim):
            mat_list.append(matrix[i][j])
    return mat_list


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
        ret = dumper.represent_scalar("tag:yaml.org,2002:float", ret_str)
    return ret


def check_xq2_degenearcy(Q2, x):
    """Check is the pair of (x,Q2) is unique."""
    size = len(x)
    unique_pairs = np.unique(np.array([Q2, x]), axis=1)
    try:
        assert unique_pairs.shape[1] == size
    except AssertionError as exc:
        raise ValueError(
            f"""(x,Q2) kinematic is degenerate need to store 3rd kinematic variable as well.
            unique kinematics are: {unique_pairs.shape[1]}, original size: {size}"""
        ) from exc

def uncert_skip_variant(source_file, skip_file, uncert_file, uncert_name, remove_source=True):
    r"""
    Create two new uncertainty files, one where the specified uncertainty
    is skipped, and one with only the specified uncertainty.

    This function should necessarily be used in the filter.py file only
    after the main uncertainty file has been created and saved. To ensure
    that the filter.py file works for everyone, the source_file and
    destination_file should be just the file names and not the full path.

    Parameters
    ----------
    source_file : str
        The name of the original uncertainty file
    skip_file : str
        The name of the uncertainty file where the uncertainty is skipped
    uncert_file : str
        The name of the uncertainty file where only the specified uncertainty is present
    uncert_name : str
        The name of the uncertainty to be skipped
    remove_source : bool
        If True, the source_file will be removed after the operation is complete
    """
    shutil.copy(source_file, skip_file)
    with open(skip_file, 'r') as file:
        content = yaml.safe_load(file)

    yaml.add_representer(float, prettify_float)
    content_uncert = {}

    if 'definitions' in content and uncert_name in content['definitions']:
        content_uncert['definitions'] = {uncert_name: copy.deepcopy(content['definitions'][uncert_name])}
        content_uncert['bins'] = {}
        bins = []
        for i in range(len(content['bins'])):
            bins.append({uncert_name: copy.deepcopy(content['bins'][i][uncert_name])})
        content_uncert['bins'] = bins

        del content['definitions'][uncert_name]
        for i in range(len(content['bins'])):
            del content['bins'][i][uncert_name]

    with open(skip_file, 'w') as file:
        yaml.dump(content, file, sort_keys=False)
    with open(uncert_file, 'w') as file:
        yaml.dump(content_uncert, file, sort_keys=False)
    if remove_source:
        os.remove(source_file)
