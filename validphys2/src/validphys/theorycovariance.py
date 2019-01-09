# -*- coding: utf-8 -*-
"""
theorycovariance.py
Tools for constructing and studying theory covariance matrices.
"""
from __future__ import generator_stop

import logging

from collections import defaultdict, namedtuple
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors
import pandas as pd

from reportengine.figure import figure
from reportengine.checks import make_argcheck, check
from reportengine.table import table
from reportengine import collect

from validphys.results import experiments_central_values, results
from validphys.results import Chi2Data, experiments_chi2_table
from validphys.calcutils import calc_chi2, all_chi2_theory, central_chi2_theory
from validphys.plotoptions import get_info
from validphys import plotutils
from validphys.checks import check_two_dataspecs

from IPython import embed

log = logging.getLogger(__name__)

theoryids_experiments_central_values = collect(experiments_central_values,
                                               ('theoryids',))

@make_argcheck
def _check_correct_theory_combination(theoryids, fivetheories):
    """Checks that a valid theory combination corresponding to an existing
    prescription has been inputted"""
    l = len(theoryids)
    check(l in {3, 5, 7, 9},
          "Expecting exactly 3, 5, 7 or 9 theories, but got {l}.")
    opts = {'bar', 'nobar'}
    xifs = [theoryid.get_description()['XIF'] for theoryid in theoryids]
    xirs = [theoryid.get_description()['XIR'] for theoryid in theoryids]
    if l == 3:
        correct_xifs = [1.0, 2.0, 0.5]
        correct_xirs = [1.0, 2.0, 0.5]
    elif l == 5:
        check(
            fivetheories is not None,
            "For five input theories a prescription bar or nobar for "
            "the flag fivetheories must be specified.")
        check(fivetheories in opts,
              "Invalid choice of prescription for 5 points", fivetheories,
              opts)
        if fivetheories == "nobar":
            correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0]
            correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5]
        else:
            correct_xifs = [1.0, 2.0, 0.5, 2.0, 0.5]
            correct_xirs = [1.0, 2.0, 0.5, 0.5, 2.0]
    elif l == 7:
        correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0, 2.0, 0.5]
        correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5]
    else:
        correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5]
        correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5, 0.5, 2.0]
    check(
        xifs == correct_xifs and xirs == correct_xirs,
        "Choice of input theories does not correspond to a valid "
        "prescription for theory covariance matrix calculation")

@make_argcheck
def _check_valid_shift_matrix_threshold_method(shift_threshold:(int, float, None) = None,
                                               method:(int, None) = None):
    """Checks that a valid method 1 or 2 is chosen where a threshold for
    removing elements of the shift correlation matrix has been specified"""
    opts = {1,2}
    if shift_threshold != None:
        check(method is not None, "A threshold for removing elements of the "
               "shift correlation matrix has been specified but no choice of "
               "method (1 or 2) was provided")
        check(method in opts,
              "Invalid choice of method for removing shift correlation matrix "
              "elements. Please choose 1 or 2.")

def make_scale_var_covmat(predictions):
    """Takes N theory predictions at different scales and applies N-pt scale
    variations to produce a covariance matrix."""
    l = len(predictions)
    central, *others = predictions
    deltas = (other - central for other in others)
    if l==3:
        norm = 0.5
    elif l==5:
        norm = 0.5
    elif l==7:
        norm = 1/3
    elif l==9:
        norm = 0.25
    s = norm*sum(np.outer(d, d) for d in deltas)
    return s

@table
@_check_correct_theory_combination
def theory_covmat(theoryids_experiments_central_values, experiments_index, theoryids):
    """Calculates the theory covariance matrix for scale variations.
    The matrix is a dataframe indexed by experiments_index."""
    s = make_scale_var_covmat(theoryids_experiments_central_values)
    df = pd.DataFrame(s, index=experiments_index, columns=experiments_index)
    return df

results_bytheoryids = collect(results,('theoryids',))
each_dataset_results_bytheory = collect('results_bytheoryids',
                                        ('experiments', 'experiment'))

@_check_correct_theory_combination
def theory_covmat_datasets(each_dataset_results_bytheory):
    """Produces an array of theory covariance matrices. Each matrix corresponds
    to a different dataset, which must be specified in the runcard. """
    dataset_covmats=[]
    for dataset in each_dataset_results_bytheory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        dataset_covmats.append(s)
    return dataset_covmats

@_check_correct_theory_combination
def total_covmat_datasets(each_dataset_results_bytheory):
    """Produces an array of total covariance matrices; the sum of experimental
    and scale-varied theory covariance matrices. Each matrix corresponds
    to a different dataset, which must be specified in the runcard.
    These are needed for calculation of chi2 per dataset. """
    dataset_covmats=[]
    for dataset in each_dataset_results_bytheory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        sigma = dataset[0][0].covmat
        cov = s + sigma
        dataset_covmats.append(cov)
    return dataset_covmats

@_check_correct_theory_combination
def total_covmat_diagtheory_datasets(each_dataset_results_bytheory):
    """Same as total_covmat_theory_datasets but for diagonal theory only"""
    dataset_covmats=[]
    for dataset in each_dataset_results_bytheory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        s_diag = np.zeros((len(s),len(s)))
        np.fill_diagonal(s_diag, np.diag(s))
        sigma = dataset[0][0].covmat
        cov = s_diag + sigma
        dataset_covmats.append(cov)
    return dataset_covmats


def theory_block_diag_covmat(theory_covmat_datasets, experiments_index):
    """Takes the theory covariance matrices for individual datasets and
    returns a data frame with a block diagonal theory covariance matrix
    by dataset"""
    s  = la.block_diag(*theory_covmat_datasets)
    df = pd.DataFrame(s, index=experiments_index, columns=experiments_index)
    return df

experiments_results_theory = collect('experiments_results', ('theoryids',))

@_check_correct_theory_combination
def total_covmat_experiments(experiments_results_theory):
    """Same as total_covmat_datasets but per experiment rather than
    per dataset. Needed for calculation of chi2 per experiment."""
    exp_result_covmats = []
    for exp_result in zip(*experiments_results_theory):
        theory_centrals = [x[1].central_value for x in exp_result]
        s = make_scale_var_covmat(theory_centrals)
        sigma = exp_result[0][0].covmat
        cov = s + sigma
        exp_result_covmats.append(cov)
    return exp_result_covmats

commondata_experiments = collect('commondata', ['experiments', 'experiment'])

# TODO: Improve how processes are assigned. Currently we group manually into
# Drell-Yan, Heavy Quarks and Jets but adding more processes could break
# this assignment
def process_lookup(commondata_experiments):
    """Produces a dictionary with keys corresponding to dataset names
    and values corresponding to process types. Process types are
    regrouped into the five categories 'Drell-Yan', 'Heavy Quarks', Jets',
    'DIS NC' and 'DIS CC'."""
    d = {commondata.name: get_info(commondata).process_description
         for commondata in commondata_experiments}
    for key, value in d.items():
        if "Deep Inelastic Scattering" in value:
            if ("CHORUS" in key) or ("NTV" in key) or ("HERACOMBCC" in key):
                d[key] = "DIS CC"
            else:
                d[key] = "DIS NC"
        elif "Drell-Yan" in value:
            d[key] = "Drell-Yan"
        elif "Heavy Quarks" in value:
            d[key] = "Heavy Quarks"
        elif "Jet" in value:
             d[key] = "Jets"
        else:
            pass
    return d

def dataset_names(commondata_experiments):
    """Returns a list of the names of the datasets, in the same order as
    they are inputted in the runcard"""
    names = [commondata.name for commondata in commondata_experiments]
    return names

ProcessInfo = namedtuple("ProcessInfo", ('theory', 'namelist', 'sizes'))


def combine_by_type(process_lookup,
                    each_dataset_results_bytheory, dataset_names):
    """Groups the datasets according to processes and returns three objects:
    theories_by_process: the relevant theories grouped by process type
    ordered_names: dictionary with keys of process type and values being the
                   corresponding list of names of datasets, in the order they
                   are appended to theories_by_process
    dataset_size:  dictionary with keys of dataset name and values being the
                   number of points in that dataset"""
    dataset_size = defaultdict(list)
    theories_by_process = defaultdict(list)
    ordered_names = defaultdict(list)
    for dataset, name in zip(each_dataset_results_bytheory, dataset_names):
        theory_centrals = [x[1].central_value for x in dataset]
        dataset_size[name] = len(theory_centrals[0])
        proc_type = process_lookup[name]
        ordered_names[proc_type].append(name)
        theories_by_process[proc_type].append(theory_centrals)
    for key, item in theories_by_process.items():
        theories_by_process[key] = np.concatenate(item, axis=1)
    process_info = ProcessInfo(theory = theories_by_process,
                               namelist = ordered_names,
                               sizes = dataset_size)
    return process_info


def process_starting_points(combine_by_type):
    """Returns a dictionary of indices in the covariance matrix corresponding
    to the starting point of each process."""
    process_info = combine_by_type
    running_index = 0
    start_proc = defaultdict(list)
    for name in process_info.theory:
        size = len(process_info.theory[name][0])
        start_proc[name] = running_index
        running_index += size
    return start_proc

def covmap(combine_by_type, dataset_names):
    """Creates a map between the covmat indices from matrices ordered by
    process to matrices ordered by experiment as listed in the runcard"""
    mapping = defaultdict(list)
    start_exp = defaultdict(list)
    process_info = combine_by_type
    running_index = 0
    for dataset in dataset_names:
        size = process_info.sizes[dataset]
        start_exp[dataset] = running_index
        running_index += size
    start = 0
    names_by_proc_list = [item for sublist in process_info.namelist.values() for item in sublist]
    for dataset in names_by_proc_list:
        for i in range(process_info.sizes[dataset]):
            mapping[start+i] = start_exp[dataset] + i
        start += process_info.sizes[dataset]
    return mapping

@_check_correct_theory_combination
def covs_pt_prescrip(combine_by_type, process_starting_points, theoryids,
                     fivetheories:(str, type(None)) = None):
    """Produces the sub-matrices of the theory covariance matrix according
    to a point prescription which matches the number of input theories.
    If 5 theories are provided, a scheme 'bar' or 'nobar' must be
    chosen in the runcard in order to specify the prescription. Sub-matrices
    correspond to applying the scale variation prescription to each pair of
    processes in turn, using a different procedure for the case where the
    processes are the same relative to when they are different."""
    l = len(theoryids)
    start_proc = process_starting_points
    covmats = defaultdict(list)
    process_info = combine_by_type
    for name1 in process_info.theory:
        for name2 in process_info.theory:
            central1, *others1 = process_info.theory[name1]
            deltas1 = list((other - central1 for other in others1))
            central2, *others2 = process_info.theory[name2]
            deltas2 = list((other - central2 for other in others2))
            if l==3:
                if name1 == name2:
                    s = 0.5*sum(np.outer(d, d) for d in deltas1)
                else:
                    s = 0.25*(np.outer((deltas1[0]+deltas1[1]),
                                       (deltas2[0]+deltas2[1])))
                start_locs = (start_proc[name1], start_proc[name2])
                covmats[start_locs] = s
            elif l==5:
                if name1 == name2:
                    s = 0.5*sum(np.outer(d, d) for d in deltas1)
            # 5 point --------------------------------------------------------------------
                elif fivetheories=='nobar':
                    s = 0.5*(np.outer(deltas1[0], deltas2[0]) + np.outer(
                                             deltas1[1], deltas2[1])) + 0.25*(
                        np.outer((deltas1[2] + deltas1[3]),
                                 (deltas2[2] + deltas2[3])))
             # 5bar-point -----------------------------------------------------------------
                else:
                    s = 0.25*(np.outer((deltas1[0]+deltas1[2]),
                                        (deltas2[0]+deltas2[2]))
                               + np.outer((deltas1[1]+deltas1[3]),
                                          (deltas2[1]+deltas2[3])))
             #  -----------------------------------------------------------------
                start_locs = (start_proc[name1], start_proc[name2])
                covmats[start_locs] = s
            elif l==7:
                if name1 == name2:
                    s = (1/3)*sum(np.outer(d, d) for d in deltas1)
                else:
                    s = (1/6)*(np.outer((deltas1[0]+ deltas1[4]),
                               (deltas2[0] + deltas2[4]))
                               + np.outer((deltas1[1]+ deltas1[5]),
                                          (deltas2[1] + deltas2[5]))
                               + np.outer((deltas1[2]+deltas1[3]), (
                                       deltas2[2]+ deltas2[3])))
                start_locs = (start_proc[name1], start_proc[name2])
                covmats[start_locs] = s
            elif l==9:
                if name1 == name2:
                    s = 0.25*sum(np.outer(d, d) for d in deltas1)
                else:
                    s = (1/12)*(np.outer((deltas1[0]+deltas1[4]+deltas1[6]),
                                         (deltas2[0]+deltas2[4]+deltas2[6]))
                                + np.outer((deltas1[1]+deltas1[5]+deltas1[7]),
                                           (deltas2[1]+deltas2[5]+deltas2[7]))) + (1/8)*(
                                           np.outer((deltas1[2]+deltas1[3]),
                                            (deltas2[2]+deltas2[3])))
                start_locs = (start_proc[name1], start_proc[name2])
                covmats[start_locs] = s
    return covmats

def theory_covmat_custom(covs_pt_prescrip, covmap, experiments_index):
    """Takes the individual sub-covmats between each two processes and assembles
    them into a full covmat. Then reshuffles the order from ordering by process
    to ordering by experiment as listed in the runcard"""
    matlength = int(sum([len(covmat) for covmat in covs_pt_prescrip.values()]
                        )/int(np.sqrt(len(covs_pt_prescrip))))
    mat = np.zeros((matlength,matlength))
    cov_by_exp = np.zeros((matlength,matlength))
    for locs in covs_pt_prescrip:
        cov = covs_pt_prescrip[locs]
        mat[locs[0]:(len(cov) + locs[0]),locs[1]:(len(cov.T)+locs[1])] = cov
    for i in range(matlength):
        for j in range(matlength):
            cov_by_exp[covmap[i]][covmap[j]] = mat[i][j]
    df = pd.DataFrame(cov_by_exp, index=experiments_index,
                      columns=experiments_index)
    return df

@_check_correct_theory_combination
def total_covmat_diagtheory_experiments(experiments_results_theory):
    """Same as total_covmat_datasets but per experiment rather than
    per dataset. Needed for calculation of chi2 per experiment."""
    exp_result_covmats = []
    for exp_result in zip(*experiments_results_theory):
        theory_centrals = [x[1].central_value for x in exp_result]
        s = make_scale_var_covmat(theory_centrals)
        s_diag = np.zeros((len(s),len(s)))
        np.fill_diagonal(s_diag, np.diag(s))
        sigma = exp_result[0][0].covmat
        cov = s_diag + sigma
        exp_result_covmats.append(cov)
    return exp_result_covmats

@table
def theory_corrmat(theory_covmat):
    """Calculates the theory correlation matrix for scale variations."""
    df = theory_covmat
    covmat = df.values
    diag_minus_half = (np.diagonal(covmat))**(-0.5)
    mat = diag_minus_half[:,np.newaxis]*df*diag_minus_half
    return mat

@table
def theory_blockcorrmat(theory_block_diag_covmat):
    """Calculates the theory correlation matrix for scale variations
    with block diagonal entries by dataset only"""
    mat = theory_corrmat(theory_block_diag_covmat)
    return mat

@table
def theory_corrmat_custom(theory_covmat_custom):
    """Calculates the theory correlation matrix for scale variations
    with variations by process type"""
    mat = theory_corrmat(theory_covmat_custom)
    return mat

@table
def theory_normcovmat(theory_covmat, experiments_data):
    """Calculates the theory covariance matrix for scale variations normalised
    to data."""
    df = theory_covmat
    experiments_data_array = np.array(experiments_data)
    mat = df/np.outer(experiments_data_array, experiments_data_array)
    return mat

@table
def theory_normblockcovmat(theory_block_diag_covmat, experiments_data):
    """Calculates the theory covariance matrix for scale variations
    normalised to data, block diagonal by dataset."""
    df = theory_block_diag_covmat
    experiments_data_array = np.array(experiments_data)
    mat = df/np.outer(experiments_data_array, experiments_data_array)
    return mat

@table
def theory_normcovmat_custom(theory_covmat_custom, experiments_data):
    """Calculates the theory covariance matrix for scale variations normalised
    to data, with variations according to the relevant prescription."""
    df = theory_covmat_custom
    experiments_data_array = np.array(experiments_data)
    mat = df/np.outer(experiments_data_array, experiments_data_array)
    return mat

@table
def experimentsplustheory_covmat(experiments_covmat, theory_covmat):
    """Calculates the experiment + theory covariance matrix for
    scale variations."""
    df = experiments_covmat + theory_covmat
    return df

@table
def experimentsplusblocktheory_covmat(experiments_covmat,
                                      theory_block_diag_covmat):
    """Calculates the experiment + theory covariance
    matrix for scale variations."""
    df = experiments_covmat + theory_block_diag_covmat
    return df

@table
def experimentsplustheory_covmat_custom(experiments_covmat,
                                        theory_covmat_custom):
    """Calculates the experiment + theory covariance matrix for
    scale variations correlated according to the relevant prescription."""
    df = experiments_covmat + theory_covmat_custom
    return df

@table
def experimentsplustheory_normcovmat(experiments_covmat, theory_covmat,
                                     experiments_data):
    """Calculates the experiment + theory covariance matrix for scale
       variations normalised to data."""
    df = experiments_covmat + theory_covmat
    experiments_data_array = np.array(experiments_data)
    mat = df/np.outer(experiments_data_array, experiments_data_array)
    return mat

@table
def experimentsplusblocktheory_normcovmat(experiments_covmat,
                                          theory_block_diag_covmat,
                                          experiments_data,
                                          experimentsplustheory_normcovmat):
    """Calculates the experiment + theory covariance matrix for scale
       variations normalised to data, block diagonal by data set."""
    mat = experimentsplustheory_normcovmat(experiments_covmat,
                                           theory_block_diag_covmat,
                                           experiments_data)
    return mat

@table
def experimentsplustheory_normcovmat_custom(experiments_covmat,
                                            theory_covmat_custom,
                                            experiments_data,
                                            experimentsplustheory_normcovmat):
    """Calculates the experiment + theory covariance matrix for scale
       variations normalised to data, correlations by process type."""
    mat = experimentsplustheory_normcovmat(experiments_covmat,
                                           theory_covmat_custom,
                                           experiments_data)

    return mat
@table
def experimentsplustheory_corrmat(experiments_covmat, theory_covmat):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices."""
    total_df = experiments_covmat + theory_covmat
    total_cov = (experiments_covmat + theory_covmat).values
    diag_minus_half = (np.diagonal(total_cov))**(-0.5)
    corrmat = diag_minus_half[:,np.newaxis]*total_df*diag_minus_half
    return corrmat

@table
def experimentsplusblocktheory_corrmat(experiments_covmat,
                                       theory_block_diag_covmat):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices, block diagonal by dataset."""
    corrmat = experimentsplustheory_corrmat(experiments_covmat,
                                            theory_block_diag_covmat)
    return corrmat

@table
def experimentsplustheory_corrmat_custom(experiments_covmat,
                                         theory_covmat_custom):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices, correlations by prescription."""
    corrmat = experimentsplustheory_corrmat(experiments_covmat,
                                            theory_covmat_custom)
    return corrmat

def chi2_impact(theory_covmat, experiments_covmat, experiments_results):
    """Returns total chi2 including theory cov mat"""
    dataresults, theoryresults = zip(*experiments_results)
    dat_central_list = [x.central_value for x in dataresults]
    th_central_list = [x.central_value for x in theoryresults]
    dat_central = np.concatenate(dat_central_list)
    th_central  = np.concatenate([x for x in th_central_list])
    central_diff = dat_central - th_central
    cov = theory_covmat.values + experiments_covmat.values
    return calc_chi2(la.cholesky(cov, lower=True), central_diff)/len(central_diff)

def data_theory_diff(experiments_results):
    """Returns (D-T) for central theory, for use in chi2 calculations"""
    dataresults, theoryresults = zip(*experiments_results)
    dat_central_list = [x.central_value for x in dataresults]
    th_central_list = [x.central_value for x in theoryresults]
    dat_central = np.concatenate(dat_central_list)
    th_central  = np.concatenate(th_central_list)
    central_diff = dat_central - th_central
    return central_diff

def chi2_block_impact(theory_block_diag_covmat, experiments_covmat,
                      experiments_results):
    """ Returns total chi2 including theory cov mat """
    chi2 = chi2_impact(theory_block_diag_covmat, experiments_covmat,
                       experiments_results)
    return chi2


def chi2_impact_custom(theory_covmat_custom, experiments_covmat,
                       experiments_results):
    """ Returns total chi2 including theory cov mat """
    chi2 = chi2_impact(theory_covmat_custom, experiments_covmat,
                       experiments_results)
    return chi2

def theory_diagcovmat(theory_covmat):
    """Returns theory covmat with only diagonal values"""
    s = theory_covmat.values
    s_diag = np.zeros((len(s),len(s)))
    np.fill_diagonal(s_diag, np.diag(s))
    return s_diag

def chi2_diag_only(theory_diagcovmat, experiments_covmat, data_theory_diff):
    """ Returns total chi2 including only diags of theory cov mat """
    cov = theory_diagcovmat + experiments_covmat.values
    elements = np.dot(data_theory_diff.T,np.dot(la.inv(cov),data_theory_diff))
    chi2 = (1/len(data_theory_diff))*np.sum(elements)
    return chi2

each_dataset_results = collect(results, ('experiments', 'experiment'))

def abs_chi2_data_theory_dataset(each_dataset_results, total_covmat_datasets):
    """ Returns an array of tuples (member_chi², central_chi², numpoints)
    corresponding to each data set, where theory errors are included"""
    chi2data_array = []
    for datresults, covmat in zip(each_dataset_results, total_covmat_datasets):
        data_result, th_result = datresults
        chi2s = all_chi2_theory(datresults,covmat)
        central_result = central_chi2_theory(datresults, covmat)
        chi2data_array.append(Chi2Data(th_result.stats_class(chi2s[:,np.newaxis]),
                                   central_result, len(data_result)))
    return chi2data_array

def abs_chi2_data_theory_experiment(experiments_results, total_covmat_experiments):
    """ Like abs_chi2_data_theory_dataset but for experiments not datasets"""
    chi2data_array = []
    for expresults, covmat in zip(experiments_results, total_covmat_experiments):
        data_result, th_result = expresults
        chi2s = all_chi2_theory(expresults, covmat)
        central_result = central_chi2_theory(expresults, covmat)
        chi2data_array.append(Chi2Data(th_result.stats_class(chi2s[:,np.newaxis]),
                              central_result, len(data_result)))
    return chi2data_array

def abs_chi2_data_diagtheory_experiment(experiments_results,
                                        total_covmat_diagtheory_experiments):
    """ For a diagonal theory covmat """
    return abs_chi2_data_theory_experiment(experiments_results,
                                           total_covmat_diagtheory_experiments)

def abs_chi2_data_diagtheory_dataset(each_dataset_results,
                                     total_covmat_diagtheory_datasets):
    """ For a diagonal theory covmat """
    return abs_chi2_data_theory_dataset(each_dataset_results,
                                        total_covmat_diagtheory_datasets)

@table
def experiments_chi2_table_theory(experiments, pdf,
                                  abs_chi2_data_theory_experiment,
                                  abs_chi2_data_theory_dataset):
    """Same as experiments_chi2_table but including theory covariance matrix"""
    return experiments_chi2_table(experiments, pdf,
                                  abs_chi2_data_theory_experiment,
                                  abs_chi2_data_theory_dataset)
@table
def experiments_chi2_table_diagtheory(experiments, pdf,
                                      abs_chi2_data_diagtheory_experiment,
                                      abs_chi2_data_diagtheory_dataset):
    """Same as experiments_chi2_table but including
    diagonal theory covariance matrix"""
    return experiments_chi2_table(experiments, pdf,
                                  abs_chi2_data_diagtheory_experiment,
                                  abs_chi2_data_diagtheory_dataset)

def matrix_plot_labels(df):
    explabels = [x[0] for x in df.columns]
    points = [x[2] for x in df.columns]
    unique_exp = []
    unique_exp.append([explabels[0],points[0]])
    for x in range(len(explabels)-1):
        if explabels[x+1] != explabels[x]:
            unique_exp.append([explabels[x+1],x+1])
    ticklabels = [unique_exp[x][0] for x in range(len(unique_exp))]
    startlocs = [unique_exp[x][1] for x in range(len(unique_exp))]
    startlocs += [len(explabels)]
    ticklocs = [0 for x in range(len(startlocs)-1)]
    for i in range(len(startlocs)-1):
        ticklocs[i] = 0.5*(startlocs[i+1]+startlocs[i])
    return ticklocs, ticklabels

@figure
def plot_covmat_heatmap(covmat, title):
    """Matrix plot of a covariance matrix"""
    df = covmat
    matrix = df.values
    fig,ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(100*matrix,
                            cmap=cm.Spectral_r,
                            norm=mcolors.SymLogNorm(linthresh=0.01,
                            linscale=10,
                            vmin=-100*matrix.max(),
                            vmax=100*matrix.max()))
    fig.colorbar(matrixplot, label="% of data")
    ax.set_title(title)
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks(ticklocs, ticklabels, rotation=30, ha="right")
    plt.gca().xaxis.tick_bottom()
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_corrmat_heatmap(corrmat, title):
    """Matrix plot of a correlation matrix"""
    df = corrmat
    matrix = df.values
    fig, ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    fig.colorbar(matrixplot)
    ax.set_title(title)
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks(ticklocs, ticklabels, rotation=30, ha="right")
    plt.gca().xaxis.tick_bottom()
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_normexpcovmat_heatmap(experiments_normcovmat):
    """Matrix plot of the experiment covariance matrix normalised to data."""
    fig = plot_covmat_heatmap(experiments_normcovmat,
                              "Experiment covariance matrix")
    return fig

@figure
def plot_expcorrmat_heatmap(experiments_corrmat):
    """Matrix plot of the experiment correlation matrix"""
    fig = plot_corrmat_heatmap(experiments_corrmat,
                               "Experiment correlation matrix")
    return fig

@figure
def plot_normthblockcovmat_heatmap(theory_normblockcovmat):
    """Matrix plot for block diagonal theory covariance matrix"""
    fig = plot_covmat_heatmap(theory_normblockcovmat,
                              "Block diagonal theory covariance matrix by dataset")
    return fig

@figure
def plot_normthcovmat_heatmap_custom(theory_normcovmat_custom, theoryids):
    """Matrix plot for block diagonal theory covariance matrix by process type"""
    l = len(theoryids)
    fig = plot_covmat_heatmap(theory_normcovmat_custom,
                              f"Theory covariance matrix for {l} points")
    return fig

@figure
def plot_thblockcorrmat_heatmap(theory_blockcorrmat):
    """Matrix plot of the theory correlation matrix"""
    fig = plot_corrmat_heatmap(theory_blockcorrmat,
                               "Theory correlation matrix block diagonal by dataset")
    return fig

@figure
def plot_thcorrmat_heatmap_custom(theory_corrmat_custom, theoryids):
    """Matrix plot of the theory correlation matrix, correlations by process type"""
    l = len(theoryids)
    fig = plot_corrmat_heatmap(theory_corrmat_custom,
                               f"Theory correlation matrix for {l} points")
    return fig

@figure
def plot_normexpplusblockthcovmat_heatmap(experimentsplusblocktheory_normcovmat):
    """Matrix plot of the exp + theory covariance matrix normalised to data"""
    fig = plot_covmat_heatmap(experimentsplusblocktheory_normcovmat,
                              "Experiment + theory (block diagonal by dataset) covariance matrix")
    return fig

@figure
def plot_normexpplusthcovmat_heatmap_custom(experimentsplustheory_normcovmat_custom, theoryids):
    """Matrix plot of the exp + theory covariance matrix normalised to data"""
    l = len(theoryids)
    fig = plot_covmat_heatmap(experimentsplustheory_normcovmat_custom,
                              f"Experiment + theory covariance matrix for {l} points")
    return fig

@figure
def plot_expplusblockthcorrmat_heatmap(experimentsplusblocktheory_corrmat):
    """Matrix plot of the exp + theory correlation matrix"""
    fig = plot_corrmat_heatmap(experimentsplusblocktheory_corrmat,
                               "Experiment + theory (block diagonal by dataset) correlation matrix")
    return fig

@figure
def plot_expplusthcorrmat_heatmap_custom(experimentsplustheory_corrmat_custom, theoryids):
    """Matrix plot of the exp + theory correlation matrix"""
    l = len(theoryids)
    fig = plot_corrmat_heatmap(experimentsplustheory_corrmat_custom,
                               f"Experiment + theory correlation matrix for {l} points")
    return fig

@figure
def plot_blockcovdiff_heatmap(theory_block_diag_covmat, experiments_covmat):
    """Matrix plot (thcov + expcov)/expcov"""
    df = (theory_block_diag_covmat.as_matrix()+experiments_covmat.values
          )/np.mean(experiments_covmat.values)
    fig = plot_covmat_heatmap(df,"(Theory + experiment)/mean(experiment)" +
                              "for block diagonal theory covmat by dataset")
    return fig

@figure
def plot_covdiff_heatmap_custom(theory_covmat_custom, experiments_covmat, theoryids):
    """Matrix plot (thcov + expcov)/expcov"""
    l = len(theoryids)
    df = (theory_covmat_custom+experiments_covmat
          )/np.mean(experiments_covmat.values)
    fig = plot_covmat_heatmap(df,
                              "(Theory + experiment)/mean(experiment)"
                              + f"covariance matrices for {l} points")
    return fig

@figure
def plot_diag_cov_comparison(theory_covmat_custom, experiments_covmat, experiments_data, theoryids):
    """Plot of sqrt(cov_ii)/|data_i| for cov = exp, theory, exp+theory"""
    l = len(theoryids)
    data = np.abs(experiments_data)
    df_theory = theory_covmat_custom
    df_experiment = experiments_covmat
    df_total = df_theory + df_experiment
    sqrtdiags1 = np.sqrt(np.diag(df_theory.values))
    sqrtdiags2 = np.sqrt(np.diag(df_experiment.values))
    sqrtdiags3 = np.sqrt(np.diag(df_total.values))
    fig,ax = plt.subplots(figsize=(20,10))
    ax.plot((sqrtdiags2/data).values, '.', label="Experiment", color="orange")
    ax.plot((sqrtdiags1/data).values, '.', label="Theory", color = "red")
    ax.plot((sqrtdiags3/data).values, '.', label="Total", color = "blue")
    ticklocs, ticklabels = matrix_plot_labels(df_experiment)
    plt.xticks(ticklocs, ticklabels, rotation=45, fontsize=6)
    ax.set_ylabel(r"$\frac{\sqrt{cov_{ii}}}{|D_i|}$")
    ax.set_ylim([0,0.5])
    ax.set_title(f"Square root of diagonal elements of covariances matrices for {l} points, "
                 + "normalised to absolute value of data")
    ax.legend()
    return fig

@figure
def plot_diag_cov_impact(theory_covmat_custom, experiments_covmat,
                         experiments_data, theoryids):
    """Plot ((expcov)^-1_ii)^-0.5 versus ((expcov + thcov)^-1_ii)^-0.5"""
    l = len(theoryids)
    data = experiments_data
    df_theory = theory_covmat_custom
    df_experiment = experiments_covmat
    matrix_theory = df_theory.values
    matrix_experiment = df_experiment.values
    a = (np.diag(la.inv(matrix_experiment)))**(-0.5)
    b = (np.diag(la.inv(matrix_theory+matrix_experiment)))**(-0.5)
    fig,ax = plt.subplots()
    ax.plot((a/data).values, '.', label="Experiment", color="orange")
    ax.plot((b/data).values, '.', label="Experiment + Theory", color="mediumseagreen")
    ticklocs, ticklabels = matrix_plot_labels(df_experiment)
    plt.xticks(ticklocs, ticklabels, rotation="vertical")
    ax.set_ylabel(r"$\frac{1}{D_i}\frac{1}{\sqrt{[cov^{-1}_]{ii}}}$")
    ax.set_title(f"Diagonal impact of adding theory covariance matrix for {l} points")
    ax.legend()
    return fig

@figure
def plot_datasets_chi2_theory(experiments,
                              each_dataset_chi2,
                              abs_chi2_data_theory_dataset):
    """Plot the chi² of all datasets, before and after adding theory errors, with bars."""
    ds = iter(each_dataset_chi2)
    dstheory = iter(abs_chi2_data_theory_dataset)
    dschi2 = []
    dschi2theory = []
    xticks = []
    for experiment in experiments:
        for dataset, dsres in zip(experiment, ds):
            dschi2.append(dsres.central_result/dsres.ndata)
            xticks.append(dataset.name)
    for experiment in experiments:
        for dataset, dsres in zip(experiment, dstheory):
            dschi2theory.append(dsres.central_result/dsres.ndata)
    plotvalues = np.stack((dschi2theory, dschi2))
    fig,ax = plotutils.barplot(plotvalues, collabels=xticks,
                               datalabels=["experiment + theory", "experiment"])
    ax.set_title(r"$\chi^2$ distribution for datasets")
    ax.legend(fontsize=14)
    return fig


matched_dataspecs_results = collect('results', ['dataspecs'])

LabeledShifts = namedtuple('LabeledShifts',
    ('experiment_name', 'dataset_name', 'shifts'))
@check_two_dataspecs
def dataspecs_dataset_prediction_shift(matched_dataspecs_results, experiment_name,
                                       dataset_name):
    """Compute the differnce in theory predictions between two dataspecs.
    This can be used in combination with `matched_datasets_from_dataspecs`

    It returns a ``LabeledShifts`` containing ``dataset_name``,
    ``experiment_name`` and``shifts``.
    """
    r1, r2 = matched_dataspecs_results
    res =  r1[1].central_value - r2[1].central_value
    return LabeledShifts(dataset_name=dataset_name,
                         experiment_name=experiment_name, shifts=res)

matched_dataspecs_dataset_prediction_shift = collect(
    'dataspecs_dataset_prediction_shift', ['dataspecs'])


#Not sure we want to export this, as it is 231 Mb...
#@table
@_check_valid_shift_matrix_threshold_method
def matched_datasets_shift_matrix(matched_dataspecs_dataset_prediction_shift,
                                  matched_dataspecs_dataset_theory,
                                  shift_threshold:(int, float, type(None)) = None,
                                  method:(int, type(None)) = None):
    """Produce a matrix out of the outer product of
    ``dataspecs_dataset_prediction_shift``. The matrix will be a
    pandas DataFrame, indexed similarly to ``experiments_index``.
    Note that this produces the normalised shift matrix, i.e. it is
    computed from shifts normalised to central theory."""
    all_shifts = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_prediction_shift])
    all_theory = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_theory])
    norm_shifts = all_shifts/all_theory
    mat = np.outer(norm_shifts, norm_shifts)
    for i, ival in enumerate(norm_shifts):
        for j, jval in enumerate(norm_shifts):
            if method == 1:
                if (np.abs(ival) < shift_threshold) or (np.abs(jval) < shift_threshold):
                    mat[i][j] = 0
                elif mat[i][j] > 0:
                    mat[i][j] = 1
                else:
                    mat[i][j] = -1
            elif method == 2:
                if (ival!=0) and (jval!=0):
                    if 1/shift_threshold <= np.abs(ival/jval) <= shift_threshold:
                        if mat[i][j] > 0:
                            mat[i][j] = 1
                        else:
                            mat[i][j] = -1
                    else:
                        mat[i][j] = 0
                else:
                    pass
    #build index
    expnames = np.concatenate([
        np.full(len(val.shifts), val.experiment_name, dtype=object)
        for val in matched_dataspecs_dataset_prediction_shift
    ])
    dsnames = np.concatenate([
        np.full(len(val.shifts), val.dataset_name, dtype=object)
        for val in matched_dataspecs_dataset_prediction_shift
    ])
    point_indexes = np.concatenate([
        np.arange(len(val.shifts))
        for val in matched_dataspecs_dataset_prediction_shift
    ])

    index = pd.MultiIndex.from_arrays(
        [expnames, dsnames, point_indexes],
        names=["Experiment name", "Dataset name", "Point"])

    return pd.DataFrame(mat, columns=index, index=index)

def shift_vector(matched_dataspecs_dataset_prediction_shift,
                 matched_dataspecs_dataset_theory):
    all_shifts = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_prediction_shift])
    all_theory = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_theory])
    norm_shifts = all_shifts/all_theory
     #build index
    expnames = np.concatenate([
        np.full(len(val.shifts), val.experiment_name, dtype=object)
        for val in matched_dataspecs_dataset_prediction_shift
    ])
    dsnames = np.concatenate([
        np.full(len(val.shifts), val.dataset_name, dtype=object)
        for val in matched_dataspecs_dataset_prediction_shift
    ])
    point_indexes = np.concatenate([
        np.arange(len(val.shifts))
        for val in matched_dataspecs_dataset_prediction_shift
    ])

    index = pd.MultiIndex.from_arrays(
        [expnames, dsnames, point_indexes],
        names=["Experiment name", "Dataset name", "Point"])
    return pd.DataFrame(norm_shifts, index=index)

def dataspecs_dataset_theory(matched_dataspecs_results, experiment_name, dataset_name):
    central, *others = matched_dataspecs_results
    res = central[1].central_value
    return LabeledShifts(dataset_name=dataset_name,
                         experiment_name=experiment_name, shifts=res)

matched_dataspecs_dataset_theory = collect('dataspecs_dataset_theory', ['dataspecs'])

def theory_vector(matched_dataspecs_dataset_theory):
    all_theory = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_theory])
     #build index
    expnames = np.concatenate([
        np.full(len(val.shifts), val.experiment_name, dtype=object)
        for val in matched_dataspecs_dataset_theory
    ])
    dsnames = np.concatenate([
        np.full(len(val.shifts), val.dataset_name, dtype=object)
        for val in matched_dataspecs_dataset_theory
    ])
    point_indexes = np.concatenate([
        np.arange(len(val.shifts))
        for val in matched_dataspecs_dataset_theory
    ])

    index = pd.MultiIndex.from_arrays(
        [expnames, dsnames, point_indexes],
        names=["Experiment name", "Dataset name", "Point"])
    return pd.DataFrame(all_theory, index=index)

def dataspecs_dataset_alltheory(matched_dataspecs_results, experiment_name, dataset_name):
    central, *others = matched_dataspecs_results
    res = [other[1].central_value for other in others]
    return LabeledShifts(dataset_name=dataset_name,
                         experiment_name=experiment_name, shifts=res)

matched_dataspecs_dataset_alltheory = collect('dataspecs_dataset_alltheory', ['dataspecs'])

def alltheory_vector(matched_dataspecs_dataset_alltheory, matched_dataspecs_dataset_theory):
    all_theory = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_alltheory], axis=1)
    expnames = np.concatenate([
        np.full(len(val.shifts),
        val.experiment_name, dtype=object)
        for val in matched_dataspecs_dataset_theory
    ])
    dsnames = np.concatenate([
        np.full(len(val.shifts),
        val.dataset_name, dtype=object)
        for val in matched_dataspecs_dataset_theory
    ])
    point_indexes = np.concatenate([
        np.arange(len(val.shifts))
        for val in matched_dataspecs_dataset_theory
    ])
    index = pd.MultiIndex.from_arrays(
        [expnames, dsnames, point_indexes],
        names=["Experiment name", "Dataset name", "Point"])
    theory_vectors = []
    for theoryvector in all_theory:
        theory_vectors.append(pd.DataFrame(theoryvector, index=index))
    return theory_vectors


@figure
def plot_matched_datasets_shift_matrix(matched_datasets_shift_matrix):
    """Heatmap plot of matched_datasets_shift_matrix"""
    return plot_covmat_heatmap(matched_datasets_shift_matrix,

    "Shift outer product matrix")

@table
def matched_datasets_shift_matrix_correlations(matched_datasets_shift_matrix):
    mat = matched_datasets_shift_matrix.values
    diag_minus_half = (np.diagonal(mat))**(-0.5)
    corrmat = diag_minus_half[:, np.newaxis] * mat * diag_minus_half
    corrmat = pd.DataFrame(
        corrmat,
        columns=matched_datasets_shift_matrix.columns,
        index=matched_datasets_shift_matrix.index)
    return corrmat

@figure
def plot_matched_datasets_shift_matrix_correlations(
        matched_datasets_shift_matrix):
    """Heatmap plot of the correlations of
    matched_datasets_shift_matrix. By construction these are
    zero or one."""
    corrmat = matched_datasets_shift_matrix_correlations
    return plot_corrmat_heatmap(
        corrmat, "Shift outer product normalized (correlation) matrix")

all_matched_results = collect('matched_dataspecs_results',
                              ['dataspecs'])

def combine_by_type_dataspecs(process_lookup, all_matched_results, matched_dataspecs_dataset_name):
    return combine_by_type(process_lookup, all_matched_results, matched_dataspecs_dataset_name)

datapsecs_theoryids = collect('theoryid', ['theoryconfig', 'original', 'dataspecs'])

def process_starting_points_dataspecs(combine_by_type_dataspecs):
    return process_starting_points(combine_by_type_dataspecs)

@make_argcheck
def _check_correct_theory_combination_dataspecs(datapsecs_theoryids,
                                                fivetheories):
    return _check_correct_theory_combination.__wrapped__(
        datapsecs_theoryids, fivetheories)

#@_check_correct_theory_combination_dataspecs
def covs_pt_prescrip_dataspecs(combine_by_type_dataspecs,
                      process_starting_points_dataspecs,
                      datapsecs_theoryids,
                      fivetheories: (str, type(None)) = None):
    return covs_pt_prescrip(combine_by_type_dataspecs, process_starting_points_dataspecs,
                            datapsecs_theoryids, fivetheories)

def covmap_dataspecs(combine_by_type_dataspecs, matched_dataspecs_dataset_name):
    return covmap(combine_by_type_dataspecs, matched_dataspecs_dataset_name)

matched_dataspecs_experiment_name = collect(
    'experiment_name', ['dataspecs'])
matched_dataspecs_dataset_name = collect(
    'dataset_name', ['dataspecs'])
matched_cuts_datasets = collect('dataset', ['dataspecs'])
all_matched_datasets = collect('matched_cuts_datasets',
                               ['dataspecs'])


def all_matched_data_lengths(all_matched_datasets):
    lens = []
    for rlist in all_matched_datasets:
        lens.append(rlist[0].load().GetNData())
    return lens

def matched_experiments_index(matched_dataspecs_experiment_name,
                              matched_dataspecs_dataset_name,
                              all_matched_data_lengths):

    enames = matched_dataspecs_experiment_name
    dsnames = matched_dataspecs_dataset_name
    lens = all_matched_data_lengths
    #build index
    expnames = np.concatenate([
        np.full(l, ename, dtype=object)
        for (l, ename) in zip(lens, enames)
    ])
    dsnames = np.concatenate([
        np.full(l, dsname, dtype=object)
        for (l, dsname) in zip(lens, dsnames)
    ])
    point_indexes = np.concatenate([
        np.arange(l)
        for l in lens
    ])

    index = pd.MultiIndex.from_arrays(
        [expnames, dsnames, point_indexes],
        names=["Experiment name", "Dataset name", "Point"])
    return index

@table
def theory_covmat_custom_dataspecs(covs_pt_prescrip_dataspecs, covmap_dataspecs,
                          matched_experiments_index):
    return theory_covmat_custom(covs_pt_prescrip_dataspecs, covmap_dataspecs,
                                matched_experiments_index)

thx_corrmat = collect('theory_corrmat_custom_dataspecs',
                      ['combined_shift_and_theory_dataspecs', 'theoryconfig'])

shx_corrmat = collect('matched_datasets_shift_matrix_correlations',
                      ['combined_shift_and_theory_dataspecs', 'shiftconfig'])

thx_covmat = collect('theory_covmat_custom_dataspecs',
                      ['combined_shift_and_theory_dataspecs', 'theoryconfig'])

combined_dataspecs_results = collect('all_matched_results',
                	['combined_shift_and_theory_dataspecs', 'theoryconfig'])

shx_vector = collect('shift_vector', ['combined_shift_and_theory_dataspecs', 'shiftconfig'])

thx_vector = collect('theory_vector', ['combined_shift_and_theory_dataspecs', 'theoryconfig'])

allthx_vector = collect('alltheory_vector', ['combined_shift_and_theory_dataspecs', 'theoryconfig'])

def theory_matrix_threshold(theory_threshold:(int, float) = 0):
    """Returns the threshold below which theory correlation elements are set to
    zero when comparing to shift correlation matrix"""
    return theory_threshold

@table
def shift_to_theory_ratio(thx_corrmat, shx_corrmat):
    ratio = (thx_corrmat[0]/shx_corrmat[0]).fillna(0)
    return ratio

@figure
def shift_to_theory_ratio_plot(shift_to_theory_ratio):
    matrix = shift_to_theory_ratio
    matrix[((matrix==np.inf) | (matrix==-np.inf))] = 0
#    bins = np.linspace(-1, 1, num=21)
#    digmatrix = np.digitize(matrix, bins)
#    symdigmatrix = np.zeros((len(digmatrix), len(digmatrix)))
#    for binnumber in range(len(bins)):
#        symdigmatrix[digmatrix==binnumber+1] = bins[binnumber]
    fig, ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r)
    fig.colorbar(matrixplot)
    ax.set_title("Ratio of theory to shift correlation matrices")
    ticklocs, ticklabels = matrix_plot_labels(matrix)
    plt.xticks(ticklocs, ticklabels, rotation=30, ha="right")
    plt.gca().xaxis.tick_bottom()
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def shift_corrmat_plot(shx_corrmat):
    # removing nans and setting them to 0
    fig = plot_corrmat_heatmap(shx_corrmat[0].fillna(0),
                               "Shift correlation matrix")
    return fig

@figure
def theory_corrmat_plot(thx_corrmat):
    fig = plot_corrmat_heatmap(thx_corrmat[0],
                               "Theory correlation matrix")
    return fig

@table
def shift_corrmat_value_fractions(shx_corrmat):
    mat = shx_corrmat[0].fillna(0).values
    matsize = np.size(mat)
    fracplus = 100*np.size(np.where(mat>0))/(2*matsize)
    fracminus = 100*np.size(np.where(mat<0))/(2*matsize)
    fraczero = 100*np.size(np.where(mat==0))/(2*matsize)
    table = pd.DataFrame([fracplus, fracminus, fraczero],
                         index=[r'$f_{\rho=+1}$',
                                r'$f_{\rho=-1}$',
                                r'$f_{\rho=0}$'],
                         columns = ["% of total entries"])
    return table

@table
def theory_corrmat_value_fractions(thx_corrmat,
                                   theory_threshold:(int, float) = 0):
    mat = thx_corrmat[0].values
    newmat = np.zeros((len(mat),len(mat)))
    #coarse graining for comparison
    newmat[mat>=theory_threshold]=1
    newmat[mat<=-theory_threshold]=-1
#    mask = ((mat!=-1) & (mat!=1))
#    mat[mask]=0
    matsize = np.size(mat)
    fracplus = 100*np.size(np.where(newmat>0))/(2*matsize)
    fracminus = 100*np.size(np.where(newmat<0))/(2*matsize)
    fraczero = 100*np.size(np.where(newmat==0))/(2*matsize)
    table = pd.DataFrame([fracplus, fracminus, fraczero],
                         index=[r'$\tilde{f}_{\rho=+1}$',
                                r'$\tilde{f}_{\rho=-1}$',
                                r'$\tilde{f}_{\rho=0}$'],
                         columns = ["% of total entries"])
    return table

@table
def shift_theory_element_comparison(shx_corrmat, thx_corrmat,
                                    theory_threshold:(int, float) = 0):
    #coarse graining for comparison
    thmat = thx_corrmat[0].values
    newthmat = np.zeros((len(thmat),len(thmat)))
    newthmat[thmat>=theory_threshold]=1
    newthmat[thmat<=-theory_threshold]=-1
#    mask = ((thmat!=-1) & (thmat!=1))
#    thmat[mask]=0
    shmat = shx_corrmat[0].fillna(0).values
    num_non_zero = np.size(np.where(shmat!=0))/2
    fracsame = 100*np.size(np.where((shmat==newthmat) & (shmat!=0)))/(2*num_non_zero)
    fracdiff = 100*np.size(np.where((shmat!=newthmat) & (shmat!=0)))/(2*num_non_zero)
    table = pd.DataFrame([fracsame, fracdiff],
                         index=['same sign',
                                'different signs'],
                         columns = ["% of total entries (where shift matrix is non-zero)"])
    return table


@table
def theory_corrmat_custom_dataspecs(theory_covmat_custom_dataspecs):
    """Calculates the theory correlation matrix for scale variations
    with variations by process type"""
    mat = theory_corrmat(theory_covmat_custom_dataspecs)
    return mat

@figure
def plot_thcorrmat_heatmap_custom_dataspecs(theory_corrmat_custom_dataspecs, theoryids):
    """Matrix plot of the theory correlation matrix, correlations by process type"""
    l = len(theoryids)
    fig = plot_corrmat_heatmap(theory_corrmat_custom_dataspecs,
                               f"Theory correlation matrix for {l} points")
    return fig

def evals_nonzero_basis(allthx_vector, thx_covmat, thx_vector):
    orig_matrix = thx_covmat[0]/(np.outer(thx_vector[0], thx_vector[0]))
    # constructing shift vectors
    scalevartheory_vectors = allthx_vector[0]
    deltas = [(thx_vector[0] - scalevarvector)/thx_vector[0] for scalevarvector in allthx_vector[0]]
    # iteratively orthogonalising deltas
    ys = [delta/np.linalg.norm(delta) for delta in deltas]
    xdash = deltas[1] - ys[0]*np.dot(ys[0].T, deltas[1])[0]
    ys[1] = xdash/np.linalg.norm(xdash)
    P = np.column_stack(ys)
    projected_matrix = np.dot(P.T, np.dot(orig_matrix, P))
    w, v_projected = la.eigh(projected_matrix)
    v = np.dot(P, v_projected)
    return w, v

def theory_shift_test(shx_vector, thx_vector, evals_nonzero_basis,
                     num_evals:(int, type(None)) = None,
		     evalue_cutoff:(float, type(None)) = None):
    w, v = evals_nonzero_basis
    # Sorting real part of eigenvalues
    w = np.real(w)
    v = np.real(v)
    sort_indices = np.argsort(w)
    w = w[sort_indices]
    v = v[:, sort_indices]
    w_max = w[np.argmax(w)]
    f = -shx_vector[0].values.T[0]
    all_projectors = np.sum(f*v.T, axis=1)
#    if num_evals is not None:
#        w_nonzero = w[-num_evals:]
#        nonzero_locs = range(len(w)-num_evals, len(w))
#    elif evalue_cutoff is not None:
#        w_nonzero = w[w>evalue_cutoff*w_max]
#        nonzero_locs = np.nonzero(w>evalue_cutoff*w_max)[0]
#    else:
#        mod_larg_neg_eval = np.abs(w[0])
#        nonzero_locs = np.nonzero(w>10*mod_larg_neg_eval)[0]
##        ratio = np.abs(all_projectors/np.sqrt(np.abs(w)))
#        ratio_nonzero = ratio[ratio<3]
#        w_nonzero = []
#        for loc in nonzero_locs:
#            if loc >=0:
#                w_nonzero.append(w[loc])
    # ^ taking 0th element to extract list from tuple
    nonzero_locs = range(len(w))
    w_nonzero = w[nonzero_locs]
    v_nonzero = []
    for loc in nonzero_locs:
        if loc >=0:
            v_nonzero.append(v[:,loc])
    projectors = np.sum(f*v_nonzero, axis=1)
    projected_evectors = np.zeros((len(projectors), (len(f))))
    for i in range(len(projectors)):
        projected_evectors[i] = projectors[i]*v_nonzero[i]
    fmiss = f - np.sum(projected_evectors, axis=0)
    embed()
    return w_nonzero, v_nonzero, projectors, f, fmiss, w_max, w, all_projectors

def cutoff(theory_shift_test, num_evals:(int, type(None)) = None,
           evalue_cutoff:(float, type(None)) = None):
    w_max = theory_shift_test[5]
    if num_evals is not None:
        cutoff = f"{num_evals} largest eigenvalues were selected"
    elif evalue_cutoff is not None:
        cutoff = evalue_cutoff*w_max
    else:
        cutoff = "10 times modulus of largest 'negative' eigenvalue"
    print(f"cutoff = {cutoff}")
    return cutoff

@table
def theory_covmat_eigenvalues(theory_shift_test):
    w_nonzero, v_nonzero, projectors = theory_shift_test[:3]
    s_scrambled = np.sqrt(np.abs(w_nonzero))
    projectors_scrambled = np.ndarray.tolist(projectors)
    ratio_scrambled = projectors_scrambled/s_scrambled
    table = pd.DataFrame([s_scrambled[::-1], projectors_scrambled[::-1], ratio_scrambled[::-1]],
         		index = [r'$s_a$', r'$\delta_a$', r'$\delta_a/s_a$'])
    return table

def efficiency(theory_shift_test):
    f = theory_shift_test[3]
    fmiss = theory_shift_test[4]
    fmod = np.sqrt(np.sum(f**2))
    fmiss_mod = np.sqrt(np.sum(fmiss**2))
    efficiency = 1 - fmiss_mod/fmod
    print(f"efficiency = {efficiency}")
    return efficiency

def maxrat(theory_shift_test):
    f = theory_shift_test[3]
    fmiss = theory_shift_test[4]
    maxrat = np.max(np.abs(fmiss))/np.max(np.abs(f))
#    print(f"maxrat = {maxrat}")
    return maxrat

def validation_theory_chi2(theory_shift_test):
    projectors = theory_shift_test[2]
    evals = theory_shift_test[0]
    ratio = projectors/np.sqrt(np.abs(evals))
    th_chi2 = 1/len(evals)*np.sum(ratio**2)
    print(f"Theory chi2 = {th_chi2}")
    return th_chi2

@figure
def projector_eigenvalue_ratio(theory_shift_test):
    surviving_evals = theory_shift_test[0][::-1]
    all_projectors = theory_shift_test[7][::-1]
    all_evals = theory_shift_test[6][::-1]
    ratio = np.abs(all_projectors)/np.sqrt(np.abs(all_evals))
    masked_evals = np.zeros((len(all_evals)))
    for loc, eval in enumerate(all_evals):
        if eval in surviving_evals:
            masked_evals[loc] = eval
  #  num_evals_ignored = len(all_evals)-len(surviving_evals)
  #  # Ordering eigenvalues and their projectors at random
  #  randomise = np.arange(len(evals))
  #  np.random.shuffle(randomise)
  #  evals = np.asarray(evals)[randomise]
  #  projectors = projectors[randomise]
  #  ratio = ratio[randomise]
    fig, (ax1, ax2) = plt.subplots(2, figsize=(5,5))
    ax1.plot(np.abs(all_projectors), 'o', label = r'|$\delta_a$|')
    ax1.plot(np.sqrt(np.abs(all_evals)), 'o', label = r'$|s_a|$')
    ax1.plot(np.sqrt(np.abs(masked_evals)), 'o', label = r'surviving $|s_a|$', color='k')
    ax2.plot(ratio, 'o', color="red")
    ax1.set_title(f"Number of surviving eigenvalues = {len(surviving_evals)}", fontsize=10)
    ax1.set_xlabel("a", fontsize=20)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
#    ax1.set_ylim([all_evals[np.argmin(np.sqrt(np.abs(all_evals)))],
#		all_evals[np.argmax(np.sqrt(np.abs(all_evals)))]])
    ax1.legend()
  #  ax1.axvline(x=num_evals_ignored, color='k')
    ax2.axhline(y=3, color='k', label=r'|$\delta_a$/$s_a$| = 3')
    ax2.legend()
    ax2.set_ylabel(r"|$\delta_a$/$s_a$|")
    print(f"Subspace dimension = {len(all_evals)}")
    return fig

@figure
def shift_diag_cov_comparison(shx_vector, thx_covmat, thx_vector):
    matrix = thx_covmat[0]/(np.outer(thx_vector[0], thx_vector[0]))
    fnorm = -shx_vector[0]
    sqrtdiags = np.sqrt(np.diag(matrix))
    fig, ax = plt.subplots(figsize=(20,10))
    ax.plot(sqrtdiags*100, '.-', label="Theory", color = "red")
    ax.plot(fnorm.values*100, '.-', label="NNLO-NLO Shift", color = "black")
    ticklocs, ticklabels = matrix_plot_labels(matrix)
    plt.xticks(ticklocs, ticklabels, rotation=45, fontsize=20)
    ax.set_ylabel("% of central theory", fontsize=20)
    ax.legend(fontsize=20)
    return fig

