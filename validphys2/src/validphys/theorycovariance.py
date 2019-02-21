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
from validphys import plotutils
from validphys.checks import check_two_dataspecs

from IPython import embed

log = logging.getLogger(__name__)

theoryids_experiments_central_values = collect(experiments_central_values,
                                               ('theoryids',))

@make_argcheck
def _check_correct_theory_combination(theoryids,
                                      fivetheories:(str, type(None)) = None):
    """Checks that a valid theory combination corresponding to an existing
    prescription has been inputted"""
    l = len(theoryids)
    check(l in {3, 5, 7, 9},
          "Expecting exactly 3, 5, 7 or 9 theories, but got {l}.")
    opts = {'bar', 'nobar', 'linear'}
    xifs = [theoryid.get_description()['XIF'] for theoryid in theoryids]
    xirs = [theoryid.get_description()['XIR'] for theoryid in theoryids]
    if l == 3:
        correct_xifs = [1.0, 2.0, 0.5]
        correct_xirs = [1.0, 2.0, 0.5]
    elif l == 5:
        check(
            fivetheories is not None,
            "For five input theories a prescription bar, nobar or linear "
            "for the flag fivetheories must be specified.")
        check(fivetheories in opts,
              "Invalid choice of prescription for 5 points", fivetheories,
              opts)
        if fivetheories == "nobar" or "linear":
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

def dataset_index_byprocess(experiments_index):
    """Return a multiindex with index
       per dataset per point, ordered by process"""
    dsnames = []
    ids = experiments_index.get_level_values("id")
    for dsname in experiments_index.get_level_values("dataset"):
        dsnames.append(dsname)
    processnames = [_process_lookup(dsname) for dsname in dsnames]
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

@table
@_check_correct_theory_combination
def theory_covmat(theoryids_experiments_central_values, experiments_index, theoryids,
                  fivetheories:(str, type(None)) = None):
    """Calculates the theory covariance matrix for scale variations.
    The matrix is a dataframe indexed by experiments_index."""
    s = make_scale_var_covmat(theoryids_experiments_central_values)
    df = pd.DataFrame(s, index=experiments_index, columns=experiments_index)
    return df

results_bytheoryids = collect(results,('theoryids',))
each_dataset_results_bytheory = collect('results_bytheoryids',
                                        ('experiments', 'experiment'))

@_check_correct_theory_combination
def theory_covmat_datasets(each_dataset_results_bytheory,
                           fivetheories:(str, type(None)) = None):
    """Produces an array of theory covariance matrices. Each matrix corresponds
    to a different dataset, which must be specified in the runcard. """
    dataset_covmats=[]
    for dataset in each_dataset_results_bytheory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        dataset_covmats.append(s)
    return dataset_covmats

@_check_correct_theory_combination
def total_covmat_datasets(each_dataset_results_bytheory,
                          fivetheories:(str, type(None)) = None):
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

@_check_correct_theory_combination
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

def _process_lookup(name):
    """Produces a dictionary with keys corresponding to dataset names
    and values corresponding to process types. Process types are
    regrouped into the five categories 'Drell-Yan', 'Top', Jets',
    'DIS NC' and 'DIS CC'."""
    process_dictionary = {	"ATLASZPT8TEVMDIST": 			"DY",
				"ATLASZPT8TEVYDIST":			"DY",
				"CMSZDIFF12":				"DY",
				"ATLAS1JET11":				"JETS",
				"CMSJETS11":				"JETS",
				"CDFR2KT":				"JETS",
				"CMSTOPDIFF8TEVTTRAPNORM":		"TOP",
				"ATLASTOPDIFF8TEVTRAPNORM":		"TOP",
				"ATLASTTBARTOT":			"TOP",
				"CMSTTBARTOT":				"TOP",
				"DYE605":				"DY",
				"DYE886P":				"DY",
				"DYE886R":				"DY",
				"ATLASWZRAP36PB":			"DY",
				"ATLASZHIGHMASS49FB":			"DY",
				"ATLASLOMASSDY11EXT":			"DY",
				"ATLASWZRAP11":				"DY",
				"CMSWEASY840PB":			"DY",
				"CMSWMASY47FB":				"DY",
				"CMSDY2D11":				"DY",
				"CMSWMU8TEV":				"DY",
				"CMSWCHARMRAT":				"DY",
				"LHCBZ940PB":				"DY",
				"LHCBZEE2FB":				"DY",
				"LHCBWZMU7TEV":				"DY",
				"LHCBWZMU8TEV":				"DY",
				"D0WEASY":				"DY",
				"D0WMASY":				"DY",
				"D0ZRAP":				"DY",
				"CDFZRAP":				"DY",
				"H1HERAF2B":				"DIS NC",
				"HERACOMBCCEM":				"DIS CC",
				"HERACOMBCCEP":				"DIS CC",
				"HERACOMBNCEM":				"DIS NC",
				"HERACOMBNCEP460":			"DIS NC",
				"HERACOMBNCEP575":			"DIS NC",
				"HERACOMBNCEP820":	 		"DIS NC",
				"HERACOMBNCEP920":			"DIS NC",
				"HERAF2CHARM":				"DIS NC",
				"ZEUSHERAF2B":				"DIS NC",
				"NMCPD":				"DIS NC",
				"NMC":					"DIS NC",
				"SLACP":				"DIS NC",
				"SLACD":				"DIS NC",
				"BCDMSP":				"DIS NC",
				"BCDMSD":				"DIS NC",
				"CHORUSNU":				"DIS CC",
				"CHORUSNB":				"DIS CC",
				"NTVNUDMN":				"DIS CC",
				"NTVNBDMN":				"DIS CC"	}
    proc = process_dictionary[name]
    return proc

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
        proc_type = _process_lookup(name)
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
             # Zahari's proposal for the theoretical covariance matrix --------------------
                if fivetheories=='linear':
                    if name1 == name2:
                        s = 0.25*(np.outer(deltas1[0], deltas2[0]) - np.outer(deltas1[0], deltas2[1])
                            - np.outer(deltas1[1], deltas2[0]) + np.outer(deltas1[1], deltas2[1])
                            + np.outer(deltas1[2], deltas2[2]) - np.outer(deltas1[2], deltas2[3])
                            - np.outer(deltas1[3], deltas2[2]) + np.outer(deltas1[3], deltas2[3]))
                    else:
                        s = 0.25*(np.outer(deltas1[0], deltas2[0]) - np.outer(deltas1[0], deltas2[1])
                            - np.outer(deltas1[1], deltas2[0]) + np.outer(deltas1[1], deltas2[1]))
                else:
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
                elif seventheories=='original':
                    s = (1/6)*(np.outer((deltas1[0]+ deltas1[4]),
                               (deltas2[0] + deltas2[4]))
                               + np.outer((deltas1[1]+ deltas1[5]),
                                          (deltas2[1] + deltas2[5]))
                               + np.outer((deltas1[2]+deltas1[3]), (
                                       deltas2[2]+ deltas2[3])))
                else:
                    s = (1/6)*(2*(np.outer(deltas1[0], deltas2[0])
                                  + np.outer(deltas1[1], deltas2[1]))
                             + (np.outer((deltas1[2] + deltas1[3]),
                                         (deltas2[2] + deltas2[3]))
                             + np.outer((deltas1[4] + deltas1[5]),
                                        (deltas2[4] + deltas2[5]))))
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

@_check_correct_theory_combination
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
    # Initialise array of zeros and set precision to same as FK tables
    s_diag = np.zeros((len(s),len(s)), dtype=np.float32)
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
    """Returns the tick locations and labels based on a dataframe
    to be plotted. The dataframe is assumed to be multiindexed by
    (process, dataset, points) or else (dataset, points). The tick
    location is in the centre of the dataset, and labelling is by
    the outermost index of the multiindex."""
    if len(df.index[0]) == 3:
        proclabels = [x[0] for x in df.index]
        dslabels = [x[1] for x in df.index]
        points = [x[2] for x in df.index]
        labels = proclabels
    elif len(df.index[0]) == 2:
        dslabels = [x[0] for x in df.index]
        points = [x[1] for x in df.index]
        labels = dslabels
    unique_ds = []
    unique_ds.append([labels[0],points[0]])
    for x in range(len(labels)-1):
        if labels[x+1] != labels[x]:
            unique_ds.append([labels[x+1],x+1])
    ticklabels = [unique_ds[x][0] for x in range(len(unique_ds))]
    startlocs = [unique_ds[x][1] for x in range(len(unique_ds))]
    startlocs += [len(labels)]
    ticklocs = [0 for x in range(len(startlocs)-1)]
    for i in range(len(startlocs)-1):
        ticklocs[i] = 0.5*(startlocs[i+1]+startlocs[i])
    return ticklocs, ticklabels

@figure
def plot_covmat_heatmap(covmat, title, dataset_index_byprocess):
    """Matrix plot of a covariance matrix"""
    df = pd.DataFrame(covmat.values, index=dataset_index_byprocess,
			columns=dataset_index_byprocess)
    df.sort_index(0, inplace=True)
    df.sort_index(1, inplace=True)
    procorder = ['DIS NC', 'DIS CC', 'DY', 'JETS', 'TOP']
    dsorder = ['BCDMSP', 'BCDMSD', 'SLACP', 'SLACD', 'NMC', 'NMCPD',
               'HERAF2CHARM', 'HERACOMBNCEP460', 'HERACOMBNCEP575',
               'HERACOMBNCEP820', 'HERACOMBNCEP920', 'HERACOMBNCEM',
               'CHORUSNU', 'CHORUSNB', 'NTVNUDMN', 'NTVNBDMN',
               'HERACOMBCCEP', 'HERACOMBCCEM', 'CDFZRAP', 'D0ZRAP',
               'D0WEASY', 'D0WMASY', 'ATLASWZRAP36PB', 'ATLASZHIGHMASS49FB',
               'ATLASLOMASSDY11EXT', 'ATLASWZRAP11', 'ATLASZPT8TEVMDIST',
               'ATLASZPT8TEVYDIST', 'CMSWEASY840PB', 'CMSWMASY47FB',
               'CMSWCHARMRAT', 'CMSDY2D11', 'CMSWMU8TEV', 'LHCBZ940PB',
               'LHCBZEE2FB', 'ATLAS1JET11', 'CMSJETS11', 'CDFR2KT',
               'ATLASTTBARTOT', 'ATLASTOPDIFF8TEVTRAPNORM', 'CMSTTBARTOT',
               'CMSTOPDIFF8TEVTTRAPNORM']
    oldindex = df.index.tolist()
    newindex = sorted(oldindex, key=lambda r: (procorder.index(r[0]), dsorder.index(r[1]), r[2]))
    # reindex index
    newdf = df.reindex(newindex)
    #reindex columns by transposing, reindexing, then transposing back
    newdf = (newdf.T.reindex(newindex)).T
    matrix = newdf.values
    fig,ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(100*matrix,
                            cmap=cm.Spectral_r,
                            norm=mcolors.SymLogNorm(linthresh=0.01,
                            linscale=10,
                            vmin=-100*matrix.max(),
                            vmax=100*matrix.max()))
    fig.colorbar(matrixplot, label="% of data")
    ax.set_title(title)
    ticklocs, ticklabels = matrix_plot_labels(newdf)
    plt.xticks(ticklocs, ticklabels, rotation=30, ha="right")
    plt.gca().xaxis.tick_bottom()
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_corrmat_heatmap(corrmat, title, dataset_index_byprocess):
    """Matrix plot of a correlation matrix"""
    df = pd.DataFrame(corrmat.values, index=dataset_index_byprocess,
		 	columns=dataset_index_byprocess)
    df.sort_index(0, inplace=True)
    df.sort_index(1, inplace=True)
    procorder = ['DIS NC', 'DIS CC', 'DY', 'JETS', 'TOP']
    dsorder = ['BCDMSP', 'BCDMSD', 'SLACP', 'SLACD', 'NMC', 'NMCPD',
               'HERAF2CHARM', 'HERACOMBNCEP460', 'HERACOMBNCEP575',
               'HERACOMBNCEP820', 'HERACOMBNCEP920', 'HERACOMBNCEM',
               'CHORUSNU', 'CHORUSNB', 'NTVNUDMN', 'NTVNBDMN',
               'HERACOMBCCEP', 'HERACOMBCCEM', 'CDFZRAP', 'D0ZRAP',
               'D0WEASY', 'D0WMASY', 'ATLASWZRAP36PB', 'ATLASZHIGHMASS49FB',
               'ATLASLOMASSDY11EXT', 'ATLASWZRAP11', 'ATLASZPT8TEVMDIST',
               'ATLASZPT8TEVYDIST', 'CMSWEASY840PB', 'CMSWMASY47FB',
               'CMSWCHARMRAT', 'CMSDY2D11', 'CMSWMU8TEV', 'LHCBZ940PB',
               'LHCBZEE2FB', 'ATLAS1JET11', 'CMSJETS11', 'CDFR2KT',
               'ATLASTTBARTOT', 'ATLASTOPDIFF8TEVTRAPNORM', 'CMSTTBARTOT',
               'CMSTOPDIFF8TEVTTRAPNORM']
    oldindex = df.index.tolist()
    newindex = sorted(oldindex, key=lambda r: (procorder.index(r[0]), dsorder.index(r[1]), r[2]))
    # reindex index
    newdf = df.reindex(newindex)
    #reindex columns by transposing, reindexing, then transposing back
    newdf = (newdf.T.reindex(newindex)).T
    matrix = newdf.values
    fig, ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    fig.colorbar(matrixplot)
    ax.set_title(title)
    ticklocs, ticklabels = matrix_plot_labels(newdf)
    plt.xticks(ticklocs, ticklabels, rotation=30, ha="right")
    plt.gca().xaxis.tick_bottom()
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_normexpcovmat_heatmap(experiments_normcovmat, dataset_index_byprocess):
    """Matrix plot of the experiment covariance matrix normalised to data."""
    fig = plot_covmat_heatmap(experiments_normcovmat,
                              "Experiment covariance matrix",
				dataset_index_byprocess)
    return fig

@figure
def plot_expcorrmat_heatmap(experiments_corrmat, dataset_index_byprocess):
    """Matrix plot of the experiment correlation matrix"""
    fig = plot_corrmat_heatmap(experiments_corrmat,
                               "Experiment correlation matrix",
				dataset_index_byprocess)
    return fig

@figure
def plot_normthblockcovmat_heatmap(theory_normblockcovmat, dataset_index_byprocess):
    """Matrix plot for block diagonal theory covariance matrix"""
    fig = plot_covmat_heatmap(theory_normblockcovmat,
                              "Block diagonal theory covariance matrix by dataset",
				dataset_index_byprocess)
    return fig

@figure
def plot_normthcovmat_heatmap_custom(theory_normcovmat_custom, theoryids,
					dataset_index_byprocess):
    """Matrix plot for block diagonal theory covariance matrix by process type"""
    l = len(theoryids)
    fig = plot_covmat_heatmap(theory_normcovmat_custom,
                              f"Theory covariance matrix for {l} points",
				dataset_index_byprocess)
    return fig

@figure
def plot_thblockcorrmat_heatmap(theory_blockcorrmat, dataset_index_byprocess):
    """Matrix plot of the theory correlation matrix"""
    fig = plot_corrmat_heatmap(theory_blockcorrmat,
                               "Theory correlation matrix block diagonal by dataset",
				dataset_index_byprocess)
    return fig

@figure
def plot_thcorrmat_heatmap_custom(theory_corrmat_custom, theoryids,
					dataset_index_byprocess):
    """Matrix plot of the theory correlation matrix, correlations by process type"""
    l = len(theoryids)
    fig = plot_corrmat_heatmap(theory_corrmat_custom,
                               f"Theory correlation matrix for {l} points",
				dataset_index_byprocess)
    return fig

@figure
def plot_normexpplusblockthcovmat_heatmap(experimentsplusblocktheory_normcovmat,
						dataset_index_byprocess):
    """Matrix plot of the exp + theory covariance matrix normalised to data"""
    fig = plot_covmat_heatmap(experimentsplusblocktheory_normcovmat,
                              "Experiment + theory (block diagonal by dataset) covariance matrix",
	dataset_index_byprocess)
    return fig

@figure
def plot_normexpplusthcovmat_heatmap_custom(experimentsplustheory_normcovmat_custom,
					theoryids, dataset_index_byprocess):
    """Matrix plot of the exp + theory covariance matrix normalised to data"""
    l = len(theoryids)
    fig = plot_covmat_heatmap(experimentsplustheory_normcovmat_custom,
                              f"Experiment + theory covariance matrix for {l} points",
				dataset_index_byprocess)
    return fig

@figure
def plot_expplusblockthcorrmat_heatmap(experimentsplusblocktheory_corrmat,
					dataset_index_byprocess):
    """Matrix plot of the exp + theory correlation matrix"""
    fig = plot_corrmat_heatmap(experimentsplusblocktheory_corrmat,
                               "Experiment + theory (block diagonal by dataset) correlation matrix",
	dataset_index_byprocess)
    return fig

@figure
def plot_expplusthcorrmat_heatmap_custom(experimentsplustheory_corrmat_custom,
					theoryids, dataset_index_byprocess):
    """Matrix plot of the exp + theory correlation matrix"""
    l = len(theoryids)
    fig = plot_corrmat_heatmap(experimentsplustheory_corrmat_custom,
                               f"Experiment + theory correlation matrix for {l} points",
                               dataset_index_byprocess)
    return fig

@figure
def plot_blockcovdiff_heatmap(theory_block_diag_covmat, experiments_covmat,
				dataset_index_byprocess):
    """Matrix plot (thcov + expcov)/expcov"""
    df = (theory_block_diag_covmat.as_matrix()+experiments_covmat.values
          )/np.mean(experiments_covmat.values)
    fig = plot_covmat_heatmap(df,"(Theory + experiment)/mean(experiment)" +
                              "for block diagonal theory covmat by dataset",
				dataset_index_byprocess)
    return fig

@figure
def plot_covdiff_heatmap_custom(theory_covmat_custom, experiments_covmat,
				theoryids, dataset_index_byprocess):
    """Matrix plot (thcov + expcov)/expcov"""
    l = len(theoryids)
    df = (theory_covmat_custom+experiments_covmat
          )/np.mean(experiments_covmat.values)
    fig = plot_covmat_heatmap(df,
                              "(Theory + experiment)/mean(experiment)"
                              + f"covariance matrices for {l} points",
				dataset_index_byprocess)
    return fig

@figure
def plot_diag_cov_comparison(theory_covmat_custom, experiments_covmat,
			experiments_data, theoryids, dataset_index_byprocess):
    """Plot of sqrt(cov_ii)/|data_i| for cov = exp, theory, exp+theory"""
    l = len(theoryids)
    data = np.abs(experiments_data)
    sqrtdiags_th = np.sqrt(np.diag(theory_covmat_custom))/data
    sqrtdiags_th = pd.DataFrame(sqrtdiags_th.values, index=dataset_index_byprocess)
    sqrtdiags_th.sort_index(0,inplace=True)
    sqrtdiags_exp = np.sqrt(np.diag(experiments_covmat))/data
    sqrtdiags_exp = pd.DataFrame(sqrtdiags_exp.values, index=dataset_index_byprocess)
    sqrtdiags_exp.sort_index(0,inplace=True)
    df_total = theory_covmat_custom + experiments_covmat
    sqrtdiags_tot = np.sqrt(np.diag(df_total))/data
    sqrtdiags_tot = pd.DataFrame(sqrtdiags_tot.values, index=dataset_index_byprocess)
    sqrtdiags_tot.sort_index(0,inplace=True)
    fig,ax = plt.subplots(figsize=(20,10))
    ax.plot(sqrtdiags_exp.values, '.', label="Experiment", color="orange")
    ax.plot(sqrtdiags_th.values, '.', label="Theory", color = "red")
    ax.plot(sqrtdiags_tot.values, '.', label="Total", color = "blue")
    ticklocs, ticklabels = matrix_plot_labels(sqrtdiags_th)
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
    matrix_theory = theory_covmat_custom.values
    matrix_experiment = experiments_covmat.values
    inv_exp = (np.diag(la.inv(matrix_experiment)))**(-0.5)/experiments_data
    inv_tot = (np.diag(la.inv(matrix_theory+matrix_experiment)))**(-0.5)/experiments_data
    df_inv_exp = pd.DataFrame(inv_exp, index=dataset_index_byprocess)
    df_inv_exp.sort_index(0,inplace=True)
    df_inv_tot = pd.DataFrame(inv_tot, index=dataset_index_byprocess)
    df_inv_tot.sort_index(0,inplace=True)
    fig,ax = plt.subplots()
    ax.plot(df_inv_exp.values, '.', label="Experiment", color="orange")
    ax.plot(df_inv_tot.values, '.', label="Experiment + Theory", color="mediumseagreen")
    ticklocs, ticklabels = matrix_plot_labels(df_inv_exp)
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


matched_dataspecs_results = collect('results', ['dataspecs_with_matched_cuts'])

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
    'dataspecs_dataset_prediction_shift', ['matched_datasets_from_dataspecs'])


#Not sure we want to export this, as it is 231 Mb...
#@table
def matched_datasets_shift_matrix(matched_dataspecs_dataset_prediction_shift):
    """Priduce a matrix out of the outer product of
    ``dataspecs_dataset_prediction_shift``. The matrix will be a
    pandas DataFrame, indexed similarly to ``experiments_index``."""
    all_shifts = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_prediction_shift])
    mat = np.outer(all_shifts, all_shifts)
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
