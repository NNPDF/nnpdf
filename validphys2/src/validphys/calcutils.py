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

def bootstrap_error(Stats_Object, nresamples, 
                    apply_func:Callable=None, Data_Result_Object=None):
    """Performs bootstrap sample on either the input Stats_Object.data, a
    function applied to the Stats_Object.data or a function applied to a 
    tuple of (Data_Result_Object, Stats_Object.data) and returns error
    according to bootstrap sample
    """
    resample_data = np.empty(nresamples)
    if ~apply_func:
        for i in range(nresamples):
            resample_data[i] = Stats_Object.bootstrap_values().central_value()
    elif apply_func && ~Data_Result_Object:
        for i in range(nresamples):
            resample_data[i] = apply_func(Stats_Object.bootstrap_values().data())
    elif apply_func && Data_Result_Object:
        for i in range(nresamples):
            results = Data_Result_Object, Stats_Object.bootstrap_values().data()
            resample_data[i] = apply_func(results)
    return np.std(resample_data)
