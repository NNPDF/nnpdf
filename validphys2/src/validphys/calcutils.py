"""
calcutils.py

Low level utilities to calculate χ² and such. These are used to implement the
higher level functions in results.py
"""
import numpy as np
import scipy.linalg as la

def calc_chi2(sqrtcov, diffs):
    """Elementary function to compute the chi², given a Cholesky decomposed
    lower triangular part and a vector of differences"""
    #Note la.cho_solve doesn't really improve things here
    #NOTE: Do not enable check_finite. The upper triangular part is not
    #guaranteed to make any sense. If this causes a problem, it is a bug in
    #ibnnpdf.
    vec = la.solve_triangular(sqrtcov, diffs, lower=True, check_finite=False)
    #This sums up the result for the chi² for any input shape.
    #Sum the squares over the first dimension and leave the others alone
    return np.einsum('i...,i...->...', vec,vec)

def all_chi2(results):
    """Return the chi² for all elements in the result. Note that the
    interpretation of the result will depend on the PDF error type"""
    data_result, th_result = results
    diffs = th_result._rawdata - data_result.central_value[:,np.newaxis]
    return calc_chi2(sqrtcov=data_result.sqrtcovmat, diffs=diffs)

def central_chi2(results):
    """Calculate the chi² from the central value of the theory prediction to
    the data"""
    data_result, th_result = results
    central_diff = th_result.central_value - data_result.central_value
    return calc_chi2(data_result.sqrtcovmat, central_diff)

def all_chi2_theory(results, theory_covmat_dataset_3pt):
    data_result, th_result = results
    diffs = th_result._rawdata - data_result.central_value[:,np.newaxis]
    total_covmat = np.array(theory_covmat_dataset_3pt)
    return calc_chi2(sqrtcov=la.cholesky(total_covmat), diffs=diffs)

def central_chi2_theory(results, theory_covmat_dataset_3pt):
    data_result, th_result = results
    central_diff = th_result.central_value - data_result.central_value
    total_covmat = np.array(theory_covmat_dataset_3pt)
    return calc_chi2(la.cholesky(total_covmat), central_diff)
