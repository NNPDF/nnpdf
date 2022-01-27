# -*- coding: utf-8 -*-
"""
Utilities for computing correlations in batch.

@author: Zahari Kassabov
"""
import numpy as np
import numpy.linalg as la

from reportengine import collect

from validphys.checks import check_pdf_is_montecarlo


#This would be a good candidate to be optimized to calculate everything in one
#pass over x,
def _basic_obs_pdf_correlation(pdfarr, obsarr):
    """Calculate the correlation between pdfs and observables.
    The expected format is:

    obsarr: (nbins x nreplicas), as returned from thresult.
    pdfarr: (nreplicas x nf x nx), as returned from xplotting_grid.grid_values

    The returned array contains the PDF correlation between
    the value of the obsevable and the PDF at the corresponding point in (fl,x)
    space. The dimensions are:
    (nbins x nf x nx), compatible with grid_values, upon
    changing replicas->bins.
    """

    #Remove mean
    #TODO: This should be done at the Result level
    x = pdfarr  - np.mean(pdfarr, axis=0)
    y = obsarr.T - np.mean(obsarr, axis=1)

    #We want to compute:
    #sum(x*y)/(norm(x)*norm(y))
    #broadcast to the appropriate dimensions

    num = np.einsum('ij,ikm->jkm',y,x)

    xnorm = la.norm(x, axis=0)
    ynorm = la.norm(y, axis=0)
    #like np.outer, but keeping the right shape
    den = np.einsum('i,jk->ijk',ynorm,  xnorm)

    return num/den

def _basic_obs_obs_correlation(obsarr1, obsarr2):
    """Calculate the correlation between two observables. The expected format is
    obsarr1: (nbins1, nreplicas)
    obsarr2: (nbins2, nreplicas)

    The result is (nbins1 , nbins2), a mareix containing the correlation
    coefficients between the two sets.
    """
    #TODO: Do this at Result level taking into account error type
    x = (obsarr1.T - np.mean(obsarr1, axis=1)).T
    y = (obsarr2.T - np.mean(obsarr2, axis=1))

    return x@y/np.outer(la.norm(x,axis=1),la.norm(y,axis=0))

#TODO: Implement for other error types.
@check_pdf_is_montecarlo
def obs_pdf_correlations(pdf, results, xplotting_grid):
    """Return the correlations between each point in a dataset and the PDF
    values on a grid of (x,f) points in a format similar to `xplotting_grid`."""
    _ , th = results
    corrs = _basic_obs_pdf_correlation(xplotting_grid.grid_values, th.error_members)
    return xplotting_grid._replace(grid_values=corrs)


corrpair_results = collect('results', ['corrpair'])
corrpair_datasets = collect('dataset', ['corrpair'])

@check_pdf_is_montecarlo
def obs_obs_correlations(pdf, corrpair_results):
    """Return the theoretical correlation matrix between a pair of observables."""
    (_,th1), (_,th2) = corrpair_results
    return _basic_obs_obs_correlation(th1.error_members, th2.error_members)
