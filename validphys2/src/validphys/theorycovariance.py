# -*- coding: utf-8 -*-
"""
theorycovariance.py

Tools for constructing and studying theory covariance matrices.
"""
from __future__ import generator_stop

from collections import OrderedDict, namedtuple, Sequence
import itertools
import logging

import numpy as np
import scipy.linalg as la
import pandas as pd

from NNPDF import ThPredictions, CommonData, Experiment
from reportengine.checks import require_one, remove_outer, check_not_empty, make_argcheck, CheckError
from reportengine.table import table
from reportengine import collect

from validphys.checks import assert_use_cuts_true
from validphys.core import DataSetSpec, PDF, ExperimentSpec
from validphys.results import results, experiment_results, experiments_central_values

log = logging.getLogger(__name__)



theoryids_experiments_central_values = collect(experiments_central_values, ('theoryids',))

@make_argcheck
def check_have_three_theories(theoryids):
    l = len(theoryids)
    if l!=3:
        raise CheckError(f"Expecting exactly 3 theories, but got {l}.")

@table
@check_have_three_theories
def theory_covmat_3pt(theoryids_experiments_central_values, experiments, experiments_index):
    """Calculates the theory covariance matrix for 3-point scale variations."""
    number_theories = len(theoryids_experiments_central_values)
    central, low, high = np.array(theoryids_experiments_central_values)
    lowdiff  = low - central
    highdiff = high - central
    s = np.zeros((len(central),len(central)))
    s = 0.5*(np.outer(lowdiff,lowdiff) + np.outer(highdiff,highdiff))
    df = pd.DataFrame(s, index=experiments_index, columns=experiments_index)
    return df

theoryids_results = collect(results, ('theoryids',))

@check_have_three_theories
def theory_covmat_datasets_3pt(theoryids_experiments_central_values, each_dataset_results_theory):
    number_theories = len(theoryids_experiments_central_values)
    print(np.shape(each_dataset_results_theory))
    for dataset in each_dataset_results_theory:
        data_centrals = [x[0].central_value for x in dataset]
        theory_centrals = [x[1].central_value for x in dataset]
        central, low, high = theory_centrals
        print(low)
        print("////////////////////")
        lowdiff = low - central
        highdiff = high - central
        s = np.zeros((len(central),len(central)))
        s = 0.5*(np.outer(lowdiff,lowdiff) + np.outer(highdiff,highdiff))
        sigmas = [x[0].covmat for x in dataset]
        sigma = sigmas[0]
        cov = s + sigma
        dataset_cent_th = dataset[0]
        for x in dataset_cent_th:
            x.total_covmat = cov
    dataset_cent = [dataset[0] for dataset in each_dataset_results_theory]
    print(np.shape(dataset_cent))
    dataset_covmats = [x[0].total_covmat for x in dataset_cent]
    return dataset_covmats

def theory_covmat_experiments_3pt(theoryids_experiments_central_values, experiments_results_theory):
    number_theories = len(theoryids_experiments_central_values)
    experiments_results_theory = np.swapaxes(experiments_results_theory, 0, 1)
    for experiment in experiments_results_theory:
        data_centrals = [x[0].central_value for x in experiment]
        theory_centrals = [x[1].central_value for x in experiment]
        central, low, high = theory_centrals
        lowdiff = low - central
        highdiff = high - central
        s = np.zeros((len(central),len(central)))
        s = 0.5*(np.outer(lowdiff,lowdiff) + np.outer(highdiff,highdiff))
        sigmas = [x[0].covmat for x in experiment]
        sigma = sigmas[0]
        cov = s + sigma
        experiment_cent_th = experiment[0]
        for x in experiment_cent_th:
            x.total_covmat = cov
    experiment_cent = [experiment[0] for experiment in experiments_results_theory]
    experiment_covmats = [x[0].total_covmat for x in experiment_cent]
    return experiment_covmats

@table
def theory_corrmat_3pt(theory_covmat_3pt):
    """Calculates the theory correlation matrix for 3-point scale variations."""
    df = theory_covmat_3pt
    covmat = df.as_matrix()
    diag_minus_half = (np.diagonal(covmat))**(-0.5)
    mat = diag_minus_half[:,np.newaxis]*df*diag_minus_half
    return mat

@table
def theory_normcovmat_3pt(theory_covmat_3pt, experiments_data):
    """Calculates the theory covariance matrix for 3-point scale variations normalised to data."""
    df = theory_covmat_3pt
    experiments_data_array = np.array(experiments_data)
    mat = df/np.outer(experiments_data_array, experiments_data_array)
    return mat

@table
def experimentsplustheory_normcovmat_3pt(experiments_covmat, theory_covmat_3pt, experiments_data):
    """Calculates the experiment + theory covariance matrix for 3-point scale variations normalised to data."""
    df = experiments_covmat + theory_covmat_3pt
    experiments_data_array = np.array(experiments_data)
    mat = df/np.outer(experiments_data_array, experiments_data_array)
    return mat

@table
def experimentsplustheory_corrmat_3pt(experiments_covmat, theory_covmat_3pt):
    """Calculates the correlation matrix for the experimental  plus theory (3 pt) covariance matrices."""
    exp_df = experiments_covmat
    theory_df = theory_covmat_3pt
    total_df = experiments_covmat + theory_covmat_3pt
    exp_covmat = exp_df.as_matrix()
    theory_covmat = theory_df.as_matrix()
    total_covmat = exp_covmat + theory_covmat
    diag_minus_half = (np.diagonal(total_covmat))**(-0.5)
    corrmat = diag_minus_half[:,np.newaxis]*total_df*diag_minus_half
    return corrmat

def chi2_impact(theory_covmat_3pt, experiments_covmat, experiments_results):
    dataresults = [ x[0] for x in experiments_results ]
    theoryresults = [ x[1] for x in experiments_results ]
    dat_central_list = [x.central_value for x in dataresults]
    th_central_list = [x.central_value for x in theoryresults]
    dat_central = np.concatenate([x for x in dat_central_list])
    th_central  = np.concatenate([x for x in th_central_list])
    central_diff = dat_central - th_central
    cov = theory_covmat_3pt.as_matrix() + experiments_covmat.as_matrix()
    elements = np.dot(central_diff.T,np.dot(la.inv(cov),central_diff))
    chi2 = (1/len(central_diff))*np.sum(elements)
    return chi2

def test(experiments_results_theory, experiments_results_theory2):
 #   print(np.shape(experiment_results))
#    print(np.shape(experiment_results_theoryids))
    print(np.shape(experiments_results_theory))
    print(np.shape(experiments_results_theory2))

experiments_results = collect(experiment_results, ('experiments',))
theoryids_experiments_results = collect('experiments_results', ('theoryids',))
each_dataset_results = collect(results, ('experiments', 'experiment'))
results_theoryids = collect(results,('theoryids',))
experiment_results_theoryids = collect(experiment_results, ('theoryids',))
each_dataset_results_theory = collect('results_theoryids', ('experiments', 'experiment'))
experiments_results_theory2 = collect('experiment_results_theoryids', ('experiments', 'experiment'))
#experiments_results_theory = collect('experiment_results_theoryids', ('experiments',))
#experiments_results_theory = collect('results_theoryids', ('experiments',))

experiments_results_theory = collect('experiments_results', ('theoryids',))
