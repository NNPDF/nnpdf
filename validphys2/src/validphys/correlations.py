# -*- coding: utf-8 -*-
"""
Utilities for computing correlations in batch.

@author: Zahari Kassabov
"""
import numpy as np
import numpy.linalg as la

from reportengine import collect

from validphys.core import Stats

#This would be a good candidate to be optimized to calculate everything in one
#pass over x,
def _basic_obs_pdf_correlation(pdf_stats, obs_stats):
    """Calculate the correlation between pdfs and observables.
    The expected format is are Stats instances of:

    obs_stats: (nbins x nreplicas), as returned from thresult.
    pdf_stats: (nreplicas x nf x nx), as returned from xplotting_grid.grid_values

    The returned array contains the PDF correlation between
    the value of the obsevable and the PDF at the corresponding point in (fl,x)
    space. The dimensions are:
    (nbins x nf x nx), compatible with grid_values, upon
    changing replicas->bins.
    """
    x = pdf_stats.error_members() - pdf_stats.central_value()
    y = obs_stats.error_members() - obs_stats.central_value()

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
    Stats instances of:

    obs1: (nreplicas, nbins1)
    obs2: (nreplicas, nbins2)

    The result is (nbins1 , nbins2), a mareix containing the correlation
    coefficients between the two sets.
    """
    x = (obs1.error_members() - obs1.central_value()).T
    y = obs2.error_members() - obs2.central_value()

    return x@y/np.outer(la.norm(x,axis=1),la.norm(y,axis=0))

def obs_pdf_correlations(pdf, results, xplotting_grid):
    """Return the correlations between each point in a dataset and the PDF
    values on a grid of (x,f) points in a format similar to `xplotting_grid`."""
    _, th = results
    # Wrap the result in a standard Stats class
    # since the result is (npoints, flavours, ndata) and has nothing to do with the PDF replicas
    corrs = Stats(_basic_obs_pdf_correlation(xplotting_grid.grid_values, th.stats))
    return xplotting_grid.copy_grid(grid_values=corrs)


corrpair_results = collect("results", ["corrpair"])
corrpair_datasets = collect("dataset", ["corrpair"])


def obs_obs_correlations(pdf, corrpair_results):
    """Return the theoretical correlation matrix between a pair of observables."""
    (_, th1), (_, th2) = corrpair_results
    return _basic_obs_obs_correlation(th1.stats, th2.stats)
