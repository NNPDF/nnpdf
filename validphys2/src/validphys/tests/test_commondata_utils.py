import pytest

import numpy as np

from validphys import commondata_utils as cdu


def test_symmetrize_errors():
    right = [1, 2, -1, -2]
    left = [-2, -1, -1, 1]
    delta = [-1/2, 1/2, -1, -1/2]
    sigma = [1.65831239518, 1.65831239518,
              1.41421356237, 1.65831239518]
    for i in range(4):
        se_del, se_sig = cdu.symmetrize_errors(right[i], left[i])
        np.testing.assert_equal(delta[i], se_del)
        np.testing.assert_almost_equal(sigma[i], se_sig)

def test_percentage_to_absolute():
    perc = [5, '5%', -5, '-5%']
    value = 50
    res = [2.5, 2.5, -2.5, -2.5]
    for i in range(4):
        np.testing.assert_equal(cdu.percentage_to_absolute(perc[i], value), res[i])

def test_cormat_to_covmat():
    cormat = [1, 0.625, 0.25, 0.625, 1, 0.5, 0.25, 0.5, 1]
    errlist = [4, 2, 2]
    res = [16, 5, 2, 5, 4, 2, 2, 2, 4]
    covmat = cdu.cormat_to_covmat(errlist, cormat)
    np.testing.assert_allclose(res, covmat)

def test_covmat_to_artunc_and_matlist_to_matrix():
    covmat_list = [16, 5, 2, 5, 4, 2, 2, 2, 4]
    artunc = np.array(cdu.covmat_to_artunc(3, covmat_list))
    artuncT = artunc.transpose()
    covmat = cdu.matlist_to_matrix(3, 3, covmat_list)
    np.testing.assert_array_almost_equal(covmat, np.matmul(artunc, artuncT))

def test_cross_cormat_to_covmat():
    rowerrlist = [5, 4, 3, 5]
    colerrlist = [1, 2, 2]
    cormat = [0.2, 0.3, 0.4, 0.25, 0.35, 0.45, 
              0.15, 0.2, 0.25, 0.6, 0.4, 0.2]
    res = [1.0, 3.0, 4.0, 1.0, 2.8, 3.6, 
              0.45, 1.2, 1.5, 3.0, 4.0, 2.0]
    covmat = cdu.cross_cormat_to_covmat(rowerrlist, colerrlist, cormat)
    np.testing.assert_allclose(res, covmat, atol=10e-10)

def test_concat_matrices():
    arrA = np.array([[1, 2], 
                     [3, 4]])
    arrB = np.array([[5, 6], 
                     [7, 8]])
    arrC = np.array([[9, 10], 
                     [11, 12]])
    arrD = np.array([[13, 14], 
                     [15, 16]])
    arrE = np.array([[17, 18], 
                     [19, 20]])
    arrF = np.array([[21, 22], 
                     [23, 24]])
    res1 = [1, 2, 5, 6, 3, 4, 7, 8, 9, 10, 13, 14, 11, 
            12, 15, 16, 17, 18, 21, 22, 19, 20, 23, 24]
    res2 = [1, 2, 5, 6, 9, 10, 3, 4, 7, 8, 11, 12, 13, 
            14, 17, 18, 21, 22, 15, 16, 19, 20, 23, 24]
    np.testing.assert_allclose(res1, cdu.concat_matrices(3, 2, [arrA, arrB, arrC, arrD, arrE, arrF]))
    np.testing.assert_allclose(res2, cdu.concat_matrices(2, 3, [arrA, arrB, arrC, arrD, arrE, arrF]))

def test_trimat_to_fullmat():
    trimat = [1, 2, 3, 4, 5, 6]
    mode0 = [1, 2, 3, 2, 4, 5, 3, 5, 6]
    mode1 = [1, 2, 4, 2, 3, 5, 4, 5, 6]
    np.testing.assert_allclose(mode0, cdu.trimat_to_fullmat(0, trimat))
    np.testing.assert_allclose(mode1, cdu.trimat_to_fullmat(1, trimat))
