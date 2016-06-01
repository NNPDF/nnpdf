# -*- coding: utf-8 -*-
"""
Utilities for reweighting studies

Implements utilities for calculating the NNPDF weights and

Created on Mon May 30 12:50:16 2016

@author: Zahari Kassabov
"""
import numpy as np

from reportengine import table

from validphys.results import chi2_data, results


#TODO: implement this using reportengine expand
#use_t0 is to request that parameter to be set explicitly
def chi2_data_for_reweighting_experiments(reweighting_experiments, pdf, use_t0, t0set=None):
    return [chi2_data(results(exp,pdf,t0set)) for exp in reweighting_experiments]

def nnpdf_weights(chi2_data):
    ndata, nreplicas = chi2_data.shape
    chi2s = chi2_data #??
    logw = ((ndata - 1)/2)*np.log(chi2s) - 0.5*chi2s
    logw -= np.max(logw)
    w = np.exp(logw)
    w /= sum(w)
    return w

def reweighting_stats():...

def unweighted_index(nnpdf_weights, nsamples:int=100):
    return np.random.choice(len(nnpdf_weights), size=nsamples, p=nnpdf_weights)

def make_unweighted_index(pdf, unweight_mask):
    ...