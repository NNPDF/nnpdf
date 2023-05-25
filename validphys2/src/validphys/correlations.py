# -*- coding: utf-8 -*-
"""
Utilities for computing correlations in batch.

@author: Zahari Kassabov
"""
import numpy as np
import numpy.linalg as la

from reportengine import collect

from validphys.core import Stats
from validphys.checks import check_pdf_is_montecarlo
from validphys.results import ThPredictionsResult

#This would be a good candidate to be optimized to calculate everything in one
#pass over x,
def _basic_obs_pdf_correlation(pdf_val, obs_val):
    """Calculate the correlation between pdfs and observables.
    The expected format is two arrays

    obs_val: (nbin x nreplicas) as returned from thresults.error_members
    pdf_val: (nreplicas x nf x nf) as returned from xplotting_grid.grid_values.error_members

    The returned array contains the PDF correlation between
    the value of the obsevable and the PDF at the corresponding point in (fl,x)
    space. The dimensions are:
    (nbins x nf x nx), compatible with grid_values, upon
    changing replicas->bins.
    """
    x = pdf_val - np.mean(pdf_val, axis=0)
    y = (obs_val - np.mean(obs_val, axis=-1, keepdims=True)).T

    #We want to compute:
    #sum(x*y)/(norm(x)*norm(y))
    #broadcast to the appropriate dimensions

    num = np.einsum('ij,ikm->jkm',y,x)

    xnorm = la.norm(x, axis=0)
    ynorm = la.norm(y, axis=0)
    #like np.outer, but keeping the right shape
    den = np.einsum('i,jk->ijk',ynorm,  xnorm)

    return num/den

def _basic_obs_obs_correlation(obs1, obs2):
    """Calculate the correlation between two observables. The expected format is
    arrays instances of:
    
    obs1: (nbins, nreplicas)
    obs2: (nbins, nreplicas)

    The result is (nbins1 , nbins2), a mareix containing the correlation
    coefficients between the two sets.
    """
    x = obs1 - np.mean(obs1, axis=1, keepdims=True)
    y = (obs2 - np.mean(obs2, axis=1, keepdims=True)).T

    return x@y/np.outer(la.norm(x,axis=1),la.norm(y,axis=0))

@check_pdf_is_montecarlo
def obs_pdf_correlations(pdf, results, xplotting_grid):
    """Return the correlations between each point in a dataset and the PDF
    values on a grid of (x,f) points in a format similar to `xplotting_grid`."""
    _, th = results
    # Wrap the result in a standard Stats class
    # since the result is (npoints, flavours, ndata) and has nothing to do with the PDF replicas
    pdf_val = xplotting_grid.grid_values.error_members()
    obs_val = th.error_members
    corrs = Stats(_basic_obs_pdf_correlation(pdf_val, obs_val))
    return xplotting_grid.copy_grid(grid_values=corrs)


corrpair_results = collect("results", ["corrpair"])
corrpair_datasets = collect("dataset", ["corrpair"])

@check_pdf_is_montecarlo
def obs_obs_correlations(pdf, corrpair_results):
    """Return the theoretical correlation matrix between a pair of observables."""
    (_, th1), (_, th2) = corrpair_results
    return _basic_obs_obs_correlation(th1.error_members, th2.error_members)


@check_pdf_is_montecarlo
def mc_dataset_correlation(dataset, pdf):
    """
    Returns the theoretical correlation between the datapoints of one
    dataset
    """
    th = ThPredictionsResult.from_convolution(pdf, dataset)
    pdf_corr = np.corrcoef(th.error_members, rowvar=True)
    return pdf_corr

@check_pdf_is_montecarlo
def mc_dataset_covariance(dataset, pdf):
    """
    Returns the theoretical correlation between the datapoints of one
    dataset
    """
    th = ThPredictionsResult.from_convolution(pdf, dataset)
    pdf_cov = np.cov(th.error_members, rowvar=True)
    return pdf_cov

def hessian_dataset_correlation(dataset, pdf):
    """
    Returns the theoretical correlation between the datapoints of one
    dataset
    """
    th = ThPredictionsResult.from_convolution(pdf, dataset)
    hessian_eigenvectors = th.error_members
    central_predictions = th.central_value
    
    # Need to subtract the central set, which is not the average of the hessian eigenvectors
    # hence we cannot use np.cov(hessian_eigenvectors)
    X = hessian_eigenvectors - central_predictions.reshape((central_predictions.shape[0],1))

    # Covariance Matrix
    cov_matrix = np.einsum("ij,kj->ik", X, X)

    # compute the correlation matrix
    std_dev = np.sqrt(np.diag(cov_matrix))
    corr_matrix = cov_matrix / np.outer(std_dev, std_dev)

    return corr_matrix
    



