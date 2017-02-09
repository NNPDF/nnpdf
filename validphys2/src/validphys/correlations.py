# -*- coding: utf-8 -*-
"""
Utilities for computing correlations in batch.

@author: Zahari Kassabov
"""
import numpy as np
import numpy.linalg as la

from validphys.checks import check_pdf_is_montecarlo


#This would be a good candidate to be optimized to calculate everything in one
#pass over x,
def _basic_obs_pdf_correlation(pdfarr, obsarr):
    """Calculate the correlation between pdfs and observables.
    The expected format is:

    obsarr: (nbins x nreplicas), as returned from thresult.
    pdfarr: (nreplicas x nf x nx), as returned from xplotting_grid.grid_values

    The returned array has dimensions and contains the PDF correlation between
    the value of the obsevable and the PDF at the corresponding point in (fl,x)
    space:
    (nbins x nf x nx), compatible with grid_values, upon
    changing replicas->bins.
    """

    #Remove mean
    x = pdfarr  - np.mean(pdfarr, axis=0)
    y = obsarr.T - np.mean(obsarr, axis=1)

    #We want to compute:
    #sum(x*y)/(norm(x)*norm(y))
    #broatcast to the appropriate dimensions

    num = np.einsum('ij,ikm->jkm',y,x)

    xnorm = la.norm(x, axis=0)
    ynorm = la.norm(y, axis=0)
    #like np.outer, but keeping the right shape
    den = np.einsum('i,jk->ijk',ynorm,  xnorm)

    return num/den


#TODO: Implement for other error types. Do not use the _rawdata.
@check_pdf_is_montecarlo
def obs_pdf_correlations(pdf, results, xplotting_grid):
    _ , th = results
    corrs = _basic_obs_pdf_correlation(xplotting_grid.grid_values, th._rawdata)
    return xplotting_grid._replace(grid_values=corrs)
