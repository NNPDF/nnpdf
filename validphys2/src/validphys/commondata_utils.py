"""
Python utilities for commondata implementation.

This module provides helpful functions that automate a few
tasks that are regularly needed for the implementation of
experimental data to the commondata format. If there are 
any additional functions that should be added here as they
could simply some repetitve tasks, please do suggest.

Before the usage of any functions, it is recommended to read
the docstrings of the function to understand the inputs and
outputs.

@author: Tanishq Sharma
"""

import numpy as np

from math import sqrt
from numpy.linalg import eig

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
    semi_diff = (delta_plus + delta_minus)/2
    average = (delta_plus - delta_minus)/2
    se_delta = semi_diff
    se_sigma = sqrt(average*average + 2*semi_diff*semi_diff)
    return se_delta, se_sigma

def percentage_to_absolute(percentage, value):
    r"""Compute the absolute value of uncertainty from percentage.

    Parameters
    ----------
    percentage : string
        Experimental datasets can provide the percentage 
        uncertainties with a % sign or without one. 
        In case a % sign is used, it will be automatically 
        stripped and the remaining string converted to a 
        float. If a % sign is not present, just input the 
        percentage uncertainty as a string, i.e. str(value).
        This is important as this functions uses replace()
        function which requires a string input.
    value : float
        The data point

    Returns
    -------
    absolute : float
        The absolute value of the uncertainty
    
    """
    percentage = float(percentage.replace("%", ""))
    absolute = percentage * value * 0.01
    return absolute 

def corMat_to_covMat(errList, corMatList):
    r"""Convert correlation matrix elements to covariance
    matrix elements.

    Parameters
    ----------
    errList : list
        A one dimensional list which contains the uncertainty
        associated to each data point in order.
    corMatList : list
        A one dimensional list which contains the elements of 
        the correlation matrix row by row. Since experimental
        datasets provide these matrices in a list form, this 
        simplifies the implementation for the user.

    Returns
    -------
    covMatList : list
        A one dimensional list which contains the elements of
        the covariance matrix row by row.
    
    """
    covMatList = []
    for i in range(len(corMatList)):
        a = i // len(errList)
        b = i % len(errList)
        covMatList.append(corMatList[i] * errList[a] * errList[b])
    return covMatList

def covMat_to_artUnc(ndata, covMatList, is_normalized):
    r"""Convert the covariance matrix to a matrix of 
    artificial uncertainties.

    Parameters
    ----------
    ndata : integer
        Number of data points
    covMatList : list
        A one dimensional list which contains the elements of
        the covariance matrix row by row. Since experimental
        datasets provide these matrices in a list form, this 
        simplifies the implementation for the user.
    is_normalized : boolean
        True if the dataset contains normalized values, False 
        if the dataset contains absolute values. Needed for 
        proper determination of whether the matrix is postive-
        semidefinite.

    Returns
    -------
    artUnc : list
        A two dimensional matrix (given as a list of lists)
        which contains artificial uncertainties to be added 
        to the commondata. i^th row (or list) contains the 
        artificial uncertainties of the i^th data point.
            
    """
    epsilon = -0.0000000001
    negEValCount = 0
    psdCheck = True
    covMat = np.zeros((ndata, ndata))
    artUnc = np.zeros((ndata, ndata))
    for i in range(len(covMatList)):
        a = i // ndata
        b = i % ndata
        covMat[a][b] = covMatList[i]
    eigVal, eigVec = eig(covMat)
    if is_normalized == False:
        for i in range(ndata):
            for j in range(ndata):
                if eigVal[j] <= 0:
                    raise ValueError('The covariance matrix is not positive-semidefinite')
                else:
                    artUnc[i][j] = eigVec[i][j] * sqrt(eigVal[j])
    elif is_normalized == True:
        for j in range(len(eigVal)):
            if eigVal[j] < epsilon:
                psdCheck = False
            elif eigVal[j] > epsilon and eigVal[j] <= 0:
                negEValCount = negEValCount + 1
                if negEValCount == 2:
                    psdCheck = False
            elif eigVal[j] > 0:
                continue
        if psdCheck == False:
            raise ValueError('The covariance matrix is not positive-semidefinite')
        else:
            for i in range(ndata):
                for j in range(ndata):
                    if eigVal[j] < 0:
                        continue
                    else:
                        artUnc[i][j] = eigVec[i][j] * sqrt(eigVal[j]) 
    return artUnc.tolist()

def cross_corMat_to_covMat(rowErrList, colErrList, corMatList):
    r"""Convert cross correlation matrix elements 
    (i.e. those between different different variables or 
    observables) to covariance matrix elements.
    
    Parameters
    ----------
    rowErrList : list
        A one dimensional list which contains the uncertainty
        associated to each data point of the variable that is
        given on the vertical axis.
    colErrList : list
        A one dimensional list which contains the uncertainty
        associated to each data point of the variable that is
        given on the horizontal axis.
    corMatList : list
        A one dimensional list which contains the elements of 
        the correlation matrix row by row. Since experimental
        datasets provide these matrices in a list form, this 
        simplifies the implementation for the user.
        
    Returns
    -------
    covMatList : list
        A one dimensional list which contains the elements of
        the covariance matrix row by row.
    
    """
    covMatList = []
    for i in range(len(corMatList)):
        a = i // len(colErrList)
        b = i % len(colErrList)
        covMatList.append(corMatList[i] * rowErrList[a] * colErrList[b])
    return covMatList

def matList_to_matrix(rows, columns, matList):
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
    matList : list
        A one dimensional list which contains the elements of
        the matrix row by row.

    Returns
    -------
    matrix : numpy.ndarray
        The matrix as a numpy 2d array.
    
    """
    if rows * columns == len(matList):
        matrix = np.zeros((rows, columns))
        for i in range(rows):
            for j in range(columns):
                matrix[i][j] = matList[j + i * columns]
        matrix = np.array(matrix)
        return matrix
    else:
        raise Exception('rows * columns != len(matList)')
    
def concatMatrices(rows, columns, listOfMatrices):
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
    listOfMatrices : list
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
    finalMatList : list
        A one dimensional list which contains the elements of
        the final, fully concatenated matrix row by row.
        
    """
    for i in range(len(listOfMatrices)):
        listOfMatrices[i] = np.array(listOfMatrices[i])
    colList = []
    for i in range(rows):
        rowList = []
        for j in range(columns):
            rowList.append(listOfMatrices[j + i * columns])
        colList.append(np.concatenate(tuple(rowList), axis=1))
    finalMat = np.concatenate(tuple(colList), axis=0)
    finalMatList = []
    for i in range(len(finalMat)):
        for j in range(len(finalMat[i])):
            finalMatList.append(finalMat[i][j])
    return finalMatList

def triMat_to_fullMat(mode, triMatList):
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
    triMatList : list
        A list containing the elements of the triangular matrix,
        for example, for a 4*4 matrix, the list of 
        triangular matrix entries could be: 
        [a, b, c, d, e, f, g, h, i, j] 

    Returns
    -------
    matList : list
        A one dimensional list which contains the elements of
        the fully populated, symmetric matrix row by row.    
    
    """
    dim = int((np.sqrt(1 + 8*len(triMatList)) - 1)/2)
    matrix = np.zeros((dim, dim))
    if mode == 0:
        for i in range(dim):
            for j in range(i + 1):
                listEl = len(triMatList) - 1 - ((i*(i + 1))//2 + j)
                if i == j:
                    matrix[dim - 1 - i][dim - 1 - j] = triMatList[listEl]
                else:
                    matrix[dim - 1 - i][dim - 1 - j] = triMatList[listEl]
                    matrix[dim - 1 - j][dim - 1 - i] = triMatList[listEl]
    elif mode == 1:
        for i in range(dim):
            for j in range(i + 1):
                listEl = (i*(i + 1))//2 + j
                if i == j:
                    matrix[i][j] = triMatList[listEl]
                else:
                    matrix[i][j] = triMatList[listEl]
                    matrix[j][i] = triMatList[listEl]
    else:
        raise Exception('Mode should be 0 or 1, refer to the function for usage')
    matList = []
    for i in range(dim):
        for j in range(dim):
            matList.append(matrix[i][j])
    return matList
