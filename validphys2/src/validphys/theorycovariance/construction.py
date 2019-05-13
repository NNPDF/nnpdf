# -*- coding: utf-8 -*-
"""
construction.py
Tools for constructing theory covariance matrices and computing their chi2s.
"""
from __future__ import generator_stop

import logging

from collections import defaultdict, namedtuple
import numpy as np
import scipy.linalg as la
import pandas as pd

from reportengine.table import table
from reportengine import collect

from validphys.results import experiments_central_values, experiments_central_values_no_table
from validphys.results import Chi2Data, results
from validphys.calcutils import calc_chi2, all_chi2_theory, central_chi2_theory
from validphys.theorycovariance.theorycovarianceutils import process_lookup, check_correct_theory_combination


log = logging.getLogger(__name__)

theoryids_experiments_central_values = collect(experiments_central_values,
                                               ('theoryids',))

theoryids_experiments_central_values_no_table = collect(experiments_central_values_no_table,
                                               ('theoryids',))

collected_theoryids = collect('theoryids',
                              ['theoryconfig',])

def dataset_index_byprocess(experiments_index):
    """Return a multiindex with index
       per dataset per point, ordered by process"""
    dsnames = []
    ids = experiments_index.get_level_values("id")
    for dsname in experiments_index.get_level_values("dataset"):
        dsnames.append(dsname)
    processnames = [process_lookup(dsname) for dsname in dsnames]
    experiments_index.droplevel(level="experiment")
    newindex = pd.MultiIndex.from_arrays([processnames, dsnames, ids],
                                         names = ("process", "dataset", "id"))
    return newindex

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

@check_correct_theory_combination
def theory_covmat_no_table(theoryids_experiments_central_values_no_table,
                           experiments_index, theoryids,
                           fivetheories:(str, type(None)) = None):

    """Calculates the theory covariance matrix for scale variations.
    The matrix is a dataframe indexed by experiments_index."""
    s = make_scale_var_covmat(theoryids_experiments_central_values_no_table)
    df = pd.DataFrame(s, index=experiments_index, columns=experiments_index)
    return df

@table
@check_correct_theory_combination
def theory_covmat(theory_covmat_no_table, fivetheories:(str, type(None)) = None):
    """Duplicate of theory_covmat_no_table but with a table decorator."""
    return theory_covmat_no_table

results_bytheoryids = collect(results,('theoryids',))
each_dataset_results_bytheory = collect('results_bytheoryids',
                                        ('experiments', 'experiment'))

@check_correct_theory_combination
def theory_covmat_datasets(each_dataset_results_bytheory,
                           fivetheories:(str, type(None)) = None):
    """Produces an array of theory covariance matrices. Each matrix corresponds
    to a different dataset, which must be specified in the runcard."""
    dataset_covmats=[]
    for dataset in each_dataset_results_bytheory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        dataset_covmats.append(s)
    return dataset_covmats

@check_correct_theory_combination
def total_covmat_datasets(each_dataset_results_bytheory,
                          fivetheories:(str, type(None)) = None):
    """Produces an array of total covariance matrices; the sum of experimental
    and scale-varied theory covariance matrices. Each matrix corresponds
    to a different dataset, which must be specified in the runcard.
    These are needed for calculation of chi2 per dataset."""
    dataset_covmats=[]
    for dataset in each_dataset_results_bytheory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        sigma = dataset[0][0].covmat
        cov = s + sigma
        dataset_covmats.append(cov)
    return dataset_covmats

@check_correct_theory_combination
def total_covmat_diagtheory_datasets(each_dataset_results_bytheory,
                                     fivetheories:(str, type(None)) = None):
    """Same as total_covmat_theory_datasets but for diagonal theory only"""
    dataset_covmats=[]
    for dataset in each_dataset_results_bytheory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        # Initialise array of zeros and set precision to same as FK tables
        s_diag = np.zeros((len(s),len(s)), dtype=np.float32)
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

@check_correct_theory_combination
def total_covmat_experiments(experiments_results_theory,
                             fivetheories:(str, type(None)) = None):
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

def dataset_names(commondata_experiments):
    """Returns a list of the names of the datasets, in the same order as
    they are inputted in the runcard"""
    names = [commondata.name for commondata in commondata_experiments]
    return names

ProcessInfo = namedtuple("ProcessInfo", ('theory', 'namelist', 'sizes'))


def combine_by_type(each_dataset_results_bytheory, dataset_names):
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
        proc_type = process_lookup(name)
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

def covmat_3pt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 3pt prescription,
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = 0.5*sum(np.outer(d, d) for d in deltas1)
    else:
        s = 0.25*(np.outer((deltas1[0]+deltas1[1]),
                           (deltas2[0]+deltas2[1])))
    return s

def covmat_5pt_linear(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for Zahari's 5pt
    linear prescription, given two dataset names and collections
    of scale variation shifts"""
    if name1 == name2:
        s = 0.25*(np.outer(deltas1[0], deltas2[0])
                - np.outer(deltas1[0], deltas2[1])
                - np.outer(deltas1[1], deltas2[0])
                + np.outer(deltas1[1], deltas2[1])
                + np.outer(deltas1[2], deltas2[2])
                - np.outer(deltas1[2], deltas2[3])
                - np.outer(deltas1[3], deltas2[2])
                + np.outer(deltas1[3], deltas2[3]))
    else:
        s = 0.25*(np.outer(deltas1[0], deltas2[0])
                - np.outer(deltas1[0], deltas2[1])
                - np.outer(deltas1[1], deltas2[0])
                + np.outer(deltas1[1], deltas2[1]))
    return s

def covmat_5pt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 5pt prescription,
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = 0.5*sum(np.outer(d, d) for d in deltas1)
    else:
        s = 0.5*(np.outer(deltas1[0], deltas2[0])
                + np.outer(deltas1[1], deltas2[1])) + 0.25*(
                np.outer((deltas1[2] + deltas1[3]),
                (deltas2[2] + deltas2[3])))
    return s

def covmat_5barpt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 5barpt prescription,
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = 0.5*sum(np.outer(d, d) for d in deltas1)
    else:
        s = 0.25*(np.outer((deltas1[0]+deltas1[2]),
                           (deltas2[0]+deltas2[2]))
            + np.outer((deltas1[1]+deltas1[3]),
                       (deltas2[1]+deltas2[3])))
    return s

def covmat_7pt_orig(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for original 7pt prescription,
    now deprecated but kept for posterity,
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = (1/3)*sum(np.outer(d, d) for d in deltas1)
    else:
        s = (1/6)*(np.outer((deltas1[0]+ deltas1[4]),
                   (deltas2[0] + deltas2[4]))
            + np.outer((deltas1[1]+ deltas1[5]),
                       (deltas2[1] + deltas2[5]))
            + np.outer((deltas1[2]+deltas1[3]), (
                    deltas2[2]+ deltas2[3])))
    return s

def covmat_7pt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 7pt prescription (Gavin),
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = (1/3)*sum(np.outer(d, d) for d in deltas1)
    else:
        s = (1/6)*(2*(np.outer(deltas1[0], deltas2[0])
        + np.outer(deltas1[1], deltas2[1]))
        + (np.outer((deltas1[2] + deltas1[3]),
                    (deltas2[2] + deltas2[3]))
                + np.outer((deltas1[4] + deltas1[5]),
                           (deltas2[4] + deltas2[5]))))
    return s

def covmat_9pt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 9pt prescription,
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = 0.25*sum(np.outer(d, d) for d in deltas1)
    else:
        s = (1/12)*(np.outer((deltas1[0]+deltas1[4]+deltas1[6]),
                    (deltas2[0]+deltas2[4]+deltas2[6]))
                + np.outer((deltas1[1]+deltas1[5]+deltas1[7]),
                           (deltas2[1]+deltas2[5]+deltas2[7]))) + (1/8)*(
                           np.outer((deltas1[2]+deltas1[3]),
                                    (deltas2[2]+deltas2[3])))
    return s

@check_correct_theory_combination
def covs_pt_prescrip(combine_by_type, process_starting_points, theoryids,
                     fivetheories:(str, type(None)) = None,
                     seventheories:(str, type(None)) = None):
    """Produces the sub-matrices of the theory covariance matrix according
    to a point prescription which matches the number of input theories.
    If 5 theories are provided, a scheme 'bar' or 'nobar' must be
    chosen in the runcard in order to specify the prescription. Sub-matrices
    correspond to applying the scale variation prescription to each pair of
    processes in turn, using a different procedure for the case where the
    processes are the same relative to when they are different."""
    l = len(theoryids)
    start_proc = process_starting_points
    process_info = combine_by_type
    covmats = defaultdict(list)
    for name1 in process_info.theory:
        for name2 in process_info.theory:
            central1, *others1 = process_info.theory[name1]
            deltas1 = list((other - central1 for other in others1))
            central2, *others2 = process_info.theory[name2]
            deltas2 = list((other - central2 for other in others2))
            if l==3:
                s = covmat_3pt(name1, name2, deltas1, deltas2)
            elif l==5:
             # Zahari's proposal for the theoretical covariance matrix --------------
                if fivetheories=='linear':
                    s = covmat_5pt_linear(name1, name2, deltas1, deltas2)
             # 5 point --------------------------------------------------------------
                elif fivetheories=='nobar':
                    s = covmat_5pt(name1, name2, deltas1, deltas2)
             # 5bar-point -----------------------------------------------------------
                else:
                    s = covmat_5barpt(name1, name2, deltas1, deltas2)
             #  ---------------------------------------------------------------------
            elif l==7:
             # Outdated 7pts implementation: left for posterity ---------------------
                if seventheories=='original':
                    s = covmat_7pt_orig(name1, name2, deltas1, deltas2)
             # 7pt (Gavin) ----------------------------------------------------------
                else:
                    s = covmat_7pt(name1, name2, deltas1, deltas2)
            elif l==9:
                s = covmat_9pt(name1, name2, deltas1, deltas2)
            start_locs = (start_proc[name1], start_proc[name2])
            covmats[start_locs] = s
    return covmats

@table
def theory_covmat_custom(covs_pt_prescrip, covmap, experiments_index):
    """Takes the individual sub-covmats between each two processes and assembles
    them into a full covmat. Then reshuffles the order from ordering by process
    to ordering by experiment as listed in the runcard"""
    matlength = int(sum([len(covmat) for covmat in covs_pt_prescrip.values()]
                        )/int(np.sqrt(len(covs_pt_prescrip))))
    # Initialise arrays of zeros and set precision to same as FK tables
    mat = np.zeros((matlength,matlength), dtype=np.float32)
    cov_by_exp = np.zeros((matlength,matlength), dtype=np.float32)
    for locs in covs_pt_prescrip:
        cov = covs_pt_prescrip[locs]
        mat[locs[0]:(len(cov) + locs[0]),locs[1]:(len(cov.T)+locs[1])] = cov
    for i in range(matlength):
        for j in range(matlength):
            cov_by_exp[covmap[i]][covmap[j]] = mat[i][j]
    df = pd.DataFrame(cov_by_exp, index=experiments_index,
                      columns=experiments_index)
    return df

@check_correct_theory_combination
def total_covmat_diagtheory_experiments(experiments_results_theory,
                                        fivetheories:(str, type(None)) = None):
    """Same as total_covmat_datasets but per experiment rather than
    per dataset. Needed for calculation of chi2 per experiment."""
    exp_result_covmats = []
    for exp_result in zip(*experiments_results_theory):
        theory_centrals = [x[1].central_value for x in exp_result]
        s = make_scale_var_covmat(theory_centrals)
        # Initialise array of zeros and set precision to same as FK tables
        s_diag = np.zeros((len(s),len(s)), dtype=np.float32)
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
def experimentsplustheory_covmat(experiments_covmat_no_table, theory_covmat_no_table):
    """Calculates the experiment + theory covariance matrix for
    scale variations."""
    df = experiments_covmat_no_table + theory_covmat_no_table
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
    """Returns total chi2 including theory cov mat"""
    chi2 = chi2_impact(theory_block_diag_covmat, experiments_covmat,
                       experiments_results)
    return chi2


def chi2_impact_custom(theory_covmat_custom, experiments_covmat,
                       experiments_results):
    """Returns total chi2 including theory cov mat"""
    chi2 = chi2_impact(theory_covmat_custom, experiments_covmat,
                       experiments_results)
    return chi2

def theory_diagcovmat(theory_covmat):
    """Returns theory covmat with only diagonal values"""
    s = theory_covmat.values
    # Initialise array of zeros and set precision to same as FK tables
    s_diag = np.zeros((len(s),len(s)), dtype=np.float32)
    np.fill_diagonal(s_diag, np.diag(s))
    return s_diag

def chi2_diag_only(theory_diagcovmat, experiments_covmat, data_theory_diff):
    """Returns total chi2 including only diags of theory cov mat"""
    cov = theory_diagcovmat + experiments_covmat.values
    elements = np.dot(data_theory_diff.T,np.dot(la.inv(cov),data_theory_diff))
    chi2 = (1/len(data_theory_diff))*np.sum(elements)
    return chi2

each_dataset_results = collect(results, ('experiments', 'experiment'))

def abs_chi2_data_theory_dataset(each_dataset_results, total_covmat_datasets):
    """Returns an array of tuples (member_chi², central_chi², numpoints)
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
    """Like abs_chi2_data_theory_dataset but for experiments not datasets"""
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
    """For a diagonal theory covmat"""
    return abs_chi2_data_theory_experiment(experiments_results,
                                           total_covmat_diagtheory_experiments)

def abs_chi2_data_diagtheory_dataset(each_dataset_results,
                                     total_covmat_diagtheory_datasets):
    """For a diagonal theory covmat"""
    return abs_chi2_data_theory_dataset(each_dataset_results,
                                        total_covmat_diagtheory_datasets)
