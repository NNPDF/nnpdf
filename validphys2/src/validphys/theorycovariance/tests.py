#-*- coding: utf-8 -*-
"""
tests.py
Tools for testing theory covariance matrices and their properties.
"""
from __future__ import generator_stop

import logging

from collections import namedtuple
from itertools import product
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd

from reportengine.checks import make_argcheck, check
from reportengine.figure import figure
from reportengine.table import table
from reportengine import collect

from validphys.results import results
from validphys.checks import check_two_dataspecs

from validphys.theorycovariance.construction import _check_correct_theory_combination
from validphys.theorycovariance.construction import combine_by_type, process_starting_points
from validphys.theorycovariance.construction import theory_corrmat, process_lookup
from validphys.theorycovariance.construction import commondata_experiments, results_bytheoryids
from validphys.theorycovariance.construction import experiments_results_theory, data_theory_diff
from validphys.theorycovariance.construction import theoryids_experiments_central_values
from validphys.theorycovariance.construction import each_dataset_results_bytheory
from validphys.theorycovariance.construction import each_dataset_results, dataset_names
from validphys.theorycovariance.construction import covmap, covs_pt_prescrip, theory_covmat_custom
from validphys.theorycovariance.construction import chi2_impact_custom, chi2_diag_only
from validphys.theorycovariance.construction import total_covmat_experiments, total_covmat_datasets
from validphys.theorycovariance.construction import theory_diagcovmat, theory_covmat

from validphys.theorycovariance.output import matrix_plot_labels
from validphys.theorycovariance.output import plot_covmat_heatmap, plot_corrmat_heatmap

log = logging.getLogger(__name__)

matched_dataspecs_results = collect('results', ['dataspecs'])

@make_argcheck
def _check_valid_shift_matrix_threshold_method(shift_threshold:(int, float, None) = None,
                                               method:(int, None) = None):
    """Checks that a valid method 1 or 2 is chosen where a threshold for
    removing elements of the shift correlation matrix has been specified"""
    opts = {1,2}
    if shift_threshold is not None:
        check(method is not None, "A threshold for removing elements of the "
               "shift correlation matrix has been specified but no choice of "
               "method (1 or 2) was provided")
        check(method in opts,
              "Invalid choice of method for removing shift correlation matrix "
              "elements. Please choose 1 or 2.")

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
   # expnames = np.concatenate([
   #     np.full(len(val.shifts), val.experiment_name, dtype=object)
   #     for val in matched_dataspecs_dataset_prediction_shift
   # ])
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
        names=["Dataset name", "Point"])

    return pd.DataFrame(mat, columns=index, index=index)

def shift_vector(matched_dataspecs_dataset_prediction_shift,
                 matched_dataspecs_dataset_theory):
    all_shifts = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_prediction_shift])
    all_theory = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_theory])
    norm_shifts = all_shifts/all_theory
     #build index
 #   expnames = np.concatenate([
 #       np.full(len(val.shifts), val.experiment_name, dtype=object)
 #       for val in matched_dataspecs_dataset_prediction_shift
 #   ])
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
        names=["Dataset name", "Point"])
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
#    expnames = np.concatenate([
#        np.full(len(val.shifts), val.experiment_name, dtype=object)
#        for val in matched_dataspecs_dataset_theory
#    ])
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
        names=["Dataset name", "Point"])
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
#    expnames = np.concatenate([
#        np.full(len(val.shifts),
#        val.experiment_name, dtype=object)
#        for val in matched_dataspecs_dataset_theory
#    ])
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
        names=["Dataset name", "Point"])
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
 #   expnames = np.concatenate([
 #       np.full(l, ename, dtype=object)
 #       for (l, ename) in zip(lens, enames)
 #   ])
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
        names=["Dataset name", "Point"])
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
    # Initialise array of zeros and set precision to same as FK tables
    newmat = np.zeros((len(mat),len(mat)), dtype=np.float32)
    #coarse graining for comparison
    newmat[mat>=theory_threshold]=1
    newmat[mat<=-theory_threshold]=-1
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
    # Initialise array of zeros and set precision to same as FK tables
    newthmat = np.zeros((len(thmat),len(thmat)), dtype=np.float32)
    newthmat[thmat>=theory_threshold]=1
    newthmat[thmat<=-theory_threshold]=-1
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

@_check_correct_theory_combination_theoryconfig
def evals_nonzero_basis(allthx_vector, thx_covmat, thx_vector,
                        collected_theoryids,
                        fivetheories:(str, type(None)) = None,
                        seventheories:(str, type(None)) = None,
                        eigenvalue_cutoff:(bool, type(None)) = None,
                        use_analytic:(bool, type(None)) = None):
    def shuffle_list(l, shift):
        i=0
        newlist = l.copy()
        while i <= (shift-1):
            newlist.append(newlist.pop(0))
            i = i + 1
        return newlist
    if eigenvalue_cutoff == True:
        w = None
        v = None
    else:
        orig_matrix = (thx_covmat[0]/(np.outer(thx_vector[0], thx_vector[0])))#.reorder_levels(['Dataset name',
    									#'Experiment name',
    									#'Point'])
        # constructing shift vectors
        diffs = [((thx_vector[0] - scalevarvector)/thx_vector[0])
                                        for scalevarvector in allthx_vector[0]]
        num_pts = len(diffs) + 1
        indexlist = list(diffs[0].index.values)
        procdict = {}
        for index in indexlist:
            name = index[0]
            proc = _process_lookup(name)
            if proc not in list(procdict.keys()):
                procdict[proc] = [name]
            elif name not in procdict[proc]:
                procdict[proc].append(name)
        # creating split diffs with processes separated
        splitdiffs = []
        for process, dslist in procdict.items():
            alldatasets = [y for x in list(procdict.values()) for y in x]
            otherdatasets = [x for x in alldatasets if x not in procdict[process]]
            for diff in diffs:
                splitdiff = diff.copy()
                for ds in otherdatasets:
                    splitdiff.loc[ds] = 0
                splitdiffs.append(splitdiff)
        num_procs = len(procdict)
        if (num_pts == 3) and (num_procs == 2):
            N = (1/4)
            # defining key
            pp1 = splitdiffs[0]
            mm1 = splitdiffs[1]
            pp2 = splitdiffs[2]
            mm2 = splitdiffs[3]
            ###################
            xs = [pp1 + pp2, pp1 + mm2, mm1 + pp2, mm1 + mm2]
        elif (num_pts == 5) and (num_procs == 2)  and (fivetheories == "nobar"):
            N = (1/4)
            # defining key
            pz1 = splitdiffs[0]
            mz1 = splitdiffs[1]
            zp1 = splitdiffs[2]
            zm1 = splitdiffs[3]
            pz2 = splitdiffs[4]
            mz2 = splitdiffs[5]
            zp2 = splitdiffs[6]
            zm2 = splitdiffs[7]
            ###################
            xs = [pz1 + pz2, pz1 + pz2, mz1 + mz2, mz1 + mz2,
                  zp1 + zp2, zp1 + zm2, zm1 + zp2, zm1 + zm2 ]
        elif (num_pts == 5) and (num_procs == 2) and (fivetheories == "bar"):
            N = (1/4)
            # defining key
            pp1 = splitdiffs[0]
            mm1 = splitdiffs[1]
            pm1 = splitdiffs[2]
            mp1 = splitdiffs[3]
            pp2 = splitdiffs[4]
            mm2 = splitdiffs[5]
            pm2 = splitdiffs[6]
            mp2 = splitdiffs[7]
            ###################
            xs = [pp1 + pp2, pp1 + pm2, mm1 + mp2, mm1 + mm2,
                  pm1 + pp2, pm1 + pm2, mp1 + mp2, mp1 + mm2 ]
        elif (num_pts == 7) and (num_procs == 2) and (seventheories != "original"):
            N = (1/6)
            # defining key
            pz1 = splitdiffs[0]
            mz1 = splitdiffs[1]
            zp1 = splitdiffs[2]
            zm1 = splitdiffs[3]
            pp1 = splitdiffs[4]
            mm1 = splitdiffs[5]
            pz2 = splitdiffs[6]
            mz2 = splitdiffs[7]
            zp2 = splitdiffs[8]
            zm2 = splitdiffs[9]
            pp2 = splitdiffs[10]
            mm2 = splitdiffs[11]
            ####################
            xs = [pz1 + pz2, mz1 + mz2, zp1 + zp2, zm1 + zp2, pp1 + pp2,
                  mm1 + pp2, pz1 + pz2, mz1 + mz2, zp1 + zm2, zm1 + zm2,
                  pp1 + mm2, mm1 + mm2]
        elif (num_pts == 9) and (num_procs == 2):
            N = (1/24)
            # defining key
            pz1 = splitdiffs[0]
            mz1 = splitdiffs[1]
            zp1 = splitdiffs[2]
            zm1 = splitdiffs[3]
            pp1 = splitdiffs[4]
            mm1 = splitdiffs[5]
            pm1 = splitdiffs[6]
            mp1 = splitdiffs[7]
            pz2 = splitdiffs[8]
            mz2 = splitdiffs[9]
            zp2 = splitdiffs[10]
            zm2 = splitdiffs[11]
            pp2 = splitdiffs[12]
            mm2 = splitdiffs[13]
            pm2 = splitdiffs[14]
            mp2 = splitdiffs[15]
            ####################
            xs = [ pz1 + pz2, pz1 + pz2, pz1 + pp2, pz1 + pp2, pz1 + pm2, pz1 + pm2,
    	           pp1 + pz2, pp1 + pz2, pp1 + pp2, pp1 + pp2, pp1 + pm2, pp1 + pm2,
                   pm1 + pz2, pm1 + pz2, pm1 + pp2, pm1 + pp2, pm1 + pm2, pm1 + pm2,
                   mz1 + mz2, mz1 + mz2, mz1 + mp2, mz1 + mp2, mz1 + mm2, mz1 + mm2,
                   mp1 + mz2, mp1 + mz2, mp1 + mp2, mp1 + mp2, mp1 + mm2, mp1 + mm2,
                   mm1 + mz2, mm1 + mz2, mm1 + mp2, mm1 + mp2, mm1 + mm2, mm1 + mm2,
                   zm1 + zp2, zm1 + zp2, zm1 + zp2, zm1 + zm2, zm1 + zm2, zm1 + zm2,
                   zp1 + zp2, zp1 + zp2, zp1 + zp2, zp1 + zm2, zp1 + zm2, zp1 + zm2]
        elif (num_pts == 3) and (num_procs != 2):
            if use_analytic == True:
                xs = []
                pps = splitdiffs[::(num_pts-1)]
                mms = shuffle_list(splitdiffs,1)[::(num_pts-1)]
                # the one vector with all pluses
                xs.append(sum(pps))
                # the p vectors with one minus
                for procloc, mm in enumerate(mms):
                    newvec = pps[0].copy()
                    newvec.loc[:]=0
                    subpps = pps.copy()
                    del subpps[procloc]
                    newvec = newvec + sum(subpps) + mm
                    xs.append(newvec)
            else:
                xs = splitdiffs
        elif (num_pts == 5) and (num_procs != 2) and (fivetheories == "nobar"):
            pzs = splitdiffs[::(num_pts-1)]
            mzs = shuffle_list(splitdiffs,1)[::(num_pts-1)]
            zps = shuffle_list(splitdiffs,2)[::(num_pts-1)]
            zms = shuffle_list(splitdiffs,3)[::(num_pts-1)]
            xs = []
            if use_analytic == True:
                xs.append(sum(pzs))
                xs.append(sum(mzs))
                xs.append(sum(zps))
                # the p vectors with one minus
                for procloc, zm in enumerate(zms):
                    newvec = zps[0].copy()
                    newvec.loc[:] = 0
                    subzps = zps.copy()
                    del subzps[procloc]
                    newvec = newvec + sum(subzps) + zm
                    xs.append(newvec)
            else:
                # See Richard notes pg 20, first two vectors are just all the
                #(+,0) and (-,0) elements respectively
                xs.append(sum(pzs))
                xs.append(sum(mzs))
                # Generating the other 2^p vectors
                loccombs = [p for p in product(range(2), repeat=num_procs)]
                for loccomb in loccombs:
                    newvec = pzs[0].copy()
                    newvec.loc[:] = 0
                    for index, entry in enumerate(loccomb):
                        if entry == 0:
                            newvec = newvec + zps[index]
                        elif entry == 1:
                            newvec = newvec + zms[index]
                    xs.append(newvec)
        elif (num_pts == 5) and (num_procs != 2) and (fivetheories == "bar"):
            pps = splitdiffs[::(num_pts-1)]
            mms = shuffle_list(splitdiffs,1)[::(num_pts-1)]
            pms = shuffle_list(splitdiffs,2)[::(num_pts-1)]
            mps = shuffle_list(splitdiffs,3)[::(num_pts-1)]
            xs = []
            if use_analytic == True:
                xs.append(sum(pps))
                xs.append(sum(mps))
                # the 2p vectors with one minus
                for procloc, pm in enumerate(pms):
                    newvec = pms[0].copy()
                    newvec.loc[:] = 0
                    subpps = pps.copy()
                    del subpps[procloc]
                    newvec = newvec + sum(subpps) + pm
                    xs.append(newvec)
                for procloc, mm in enumerate(mms):
                    newvec = mms[0].copy()
                    newvec.loc[:] = 0
                    submps = mps.copy()
                    del submps[procloc]
                    newvec = newvec + sum(submps) + mm
                    xs.append(newvec)
            else:
                loccombs = [p for p in product(range(2), repeat=num_procs)]
                for loccomb in loccombs:
                    newvec = pps[0].copy()
                    newvec.loc[:] = 0
                    for index, entry in enumerate(loccomb):
                        if entry == 0:
                            newvec = newvec + pps[index]
                        elif entry == 1:
                            newvec = newvec + pms[index]
                    xs.append(newvec)
                for loccomb in loccombs:
                    newvec = pps[0].copy()
                    newvec.loc[:] = 0
                    for index, entry in enumerate(loccomb):
                        if entry == 0:
                            newvec = newvec + mps[index]
                        elif entry == 1:
                            newvec = newvec + mms[index]
                    xs.append(newvec)
        elif (num_pts == 7) and (num_procs != 2) and (seventheories != "original"):
            pzs = splitdiffs[::(num_pts-1)]
            mzs = shuffle_list(splitdiffs,1)[::(num_pts-1)]
            zps = shuffle_list(splitdiffs,2)[::(num_pts-1)]
            zms = shuffle_list(splitdiffs,3)[::(num_pts-1)]
            pps = shuffle_list(splitdiffs,4)[::(num_pts-1)]
            mms = shuffle_list(splitdiffs,5)[::(num_pts-1)]
            xs = []
            if use_analytic == True:
                # 3pt-like part
                xs.append(sum(pps))
                # the p vectors with one minus
                for procloc, mm in enumerate(mms):
                    newvec = pps[0].copy()
                    newvec.loc[:]=0
                    subpps = pps.copy()
                    del subpps[procloc]
                    newvec = newvec + sum(subpps) + mm
                    xs.append(newvec)
                # 5pt-like part
                xs.append(sum(pzs))
                xs.append(sum(mzs))
                xs.append(sum(zps))
                # the p vectors with one minus
                for procloc, zm in enumerate(zms):
                    newvec = zps[0].copy()
                    newvec.loc[:] = 0
                    subzps = zps.copy()
                    del subzps[procloc]
                    newvec = newvec + sum(subzps) + zm
                    xs.append(newvec)
            else:
                # 3pt-like part
                loccombs = [p for p in product(range(2), repeat=num_procs)]
                for loccomb in loccombs:
                    newvec = pzs[0].copy()
                    newvec.loc[:] = 0
                    for index, entry in enumerate(loccomb):
                        if entry == 0:
                            newvec = newvec + pps[index]
                        elif entry == 1:
                            newvec = newvec + mms[index]
                    xs.append(newvec)
                # 5pt-like part
                xs.append(sum(pzs))
                xs.append(sum(mzs))
                # Generating the other 2^p vectors
                loccombs = [p for p in product(range(2), repeat=num_procs)]
                for loccomb in loccombs:
                    newvec = pzs[0].copy()
                    newvec.loc[:] = 0
                    for index, entry in enumerate(loccomb):
                        if entry == 0:
                            newvec = newvec + zps[index]
                        elif entry == 1:
                            newvec = newvec + zms[index]
                    xs.append(newvec)
        elif (num_pts == 9) and (num_procs != 2):
            pzs = splitdiffs[::(num_pts-1)]
            mzs = shuffle_list(splitdiffs,1)[::(num_pts-1)]
            zps = shuffle_list(splitdiffs,2)[::(num_pts-1)]
            zms = shuffle_list(splitdiffs,3)[::(num_pts-1)]
            pps = shuffle_list(splitdiffs,4)[::(num_pts-1)]
            mms = shuffle_list(splitdiffs,5)[::(num_pts-1)]
            pms = shuffle_list(splitdiffs,6)[::(num_pts-1)]
            mps = shuffle_list(splitdiffs,7)[::(num_pts-1)]
            xs = []
            if use_analytic == True:
                xs.append(sum(pps))
                xs.append(sum(mps))
                xs.append(sum(zps))
                for procloc, zm in enumerate(zms):
                    newvec = zps[0].copy()
                    newvec.loc[:] = 0
                    subzps = zps.copy()
                    del subzps[procloc]
                    newvec = newvec + sum(subzps) + zm
                    xs.append(newvec)
                for procloc, pm in enumerate(pms):
                    newvec = pps[0].copy()
                    newvec.loc[:] = 0
                    subpps = pps.copy()
                    del subpps[procloc]
                    newvec = newvec + sum(subpps) + pm
                    xs.append(newvec)
                for procloc, mm in enumerate(mms):
                    newvec = mps[0].copy()
                    newvec.loc[:] = 0
                    submps = mps.copy()
                    del submps[procloc]
                    newvec = newvec + sum(submps) + mm
                    xs.append(newvec)
                for procloc, pz in enumerate(pzs):
                    newvec = pps[0].copy()
                    newvec.loc[:] = 0
                    subpps = pps.copy()
                    del subpps[procloc]
                    newvec = newvec + sum(subpps) + pz
                    xs.append(newvec)
                for procloc, mz in enumerate(mzs):
                    newvec = mps[0].copy()
                    newvec.loc[:] = 0
                    submps = mps.copy()
                    del submps[procloc]
                    newvec = newvec + sum(submps) + mz
                    xs.append(newvec)
            else:
                # Generating first 2^p vectors
                loccombs = [p for p in product(range(2), repeat=num_procs)]
                for loccomb in loccombs:
                    newvec = pzs[0].copy()
                    newvec.loc[:] = 0
                    for index, entry in enumerate(loccomb):
                        if entry == 0:
                            newvec = newvec + zps[index]
                        elif entry == 1:
                            newvec = newvec + zms[index]
                    xs.append(newvec)
                loccombs2 = [p for p in product(range(3), repeat=num_procs)]
                for loccomb2 in loccombs2:
                    newvec = pzs[0].copy()
                    newvec.loc[:] = 0
                    for index, entry in enumerate(loccomb2):
                        if entry == 0:
                            newvec = newvec + pzs[index]
                        elif entry == 1:
                            newvec = newvec + pps[index]
                        elif entry == 2:
                            newvec = newvec + pms[index]
                    xs.append(newvec)
                for loccomb2 in loccombs2:
                    newvec = pzs[0].copy()
                    newvec.loc[:] = 0
                    for index, entry in enumerate(loccomb2):
                        if entry == 0:
                            newvec = newvec + mzs[index]
                        elif entry == 1:
                            newvec = newvec + mps[index]
                        elif entry == 2:
                            newvec = newvec + mms[index]
                    xs.append(newvec)
        A = pd.concat(xs, axis=1)
        if num_procs == 2:
            covmat = N*A.dot(A.T)
        else:
            covmat = orig_matrix
        ys = [x/np.linalg.norm(x) for x in xs]
        for i in range(1, len(xs)):
            for j in range(0,i):
                ys[i] = ys[i] - (ys[i].T.dot(ys[j]))[0][0]*ys[j]/np.linalg.norm(ys[j])
                ys[i] = ys[i]/np.linalg.norm(ys[i])
        P = pd.concat(ys, axis=1)
       # P = scipy.linalg.orth(A)
        projected_matrix = (P.T).dot(covmat.dot(P))
        w, v_projected = la.eigh(projected_matrix)
        v = P.dot(v_projected)
    return w, v

def theory_shift_test(thx_covmat, shx_vector, thx_vector, evals_nonzero_basis,
		     eigenvalue_cutoff:(bool, type(None)) = None):
    if eigenvalue_cutoff == True:
        matrix = thx_covmat[0]/(np.outer(thx_vector[0], thx_vector[0]))
        w, v = la.eigh(matrix)
    else:
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
    if eigenvalue_cutoff == True:
        mod_larg_neg_eval = np.abs(w[0])
        nonzero_locs = np.nonzero(w>10*mod_larg_neg_eval)[0]
        w_nonzero = []
        for loc in nonzero_locs:
            if loc >=0:
                w_nonzero.append(w[loc])
    # ^ taking 0th element to extract list from tuple
    else:
        nonzero_locs = range(len(w))
        w_nonzero = w[nonzero_locs]
    v_nonzero = []
    for loc in nonzero_locs:
        if loc >=0:
            v_nonzero.append(v[:,loc])
    projectors = np.sum(f*v_nonzero, axis=1)
    # Initialise array of zeros and set precision to same as FK tables
    projected_evectors = np.zeros((len(projectors), (len(f))), dtype=np.float32)
    for i in range(len(projectors)):
        projected_evectors[i] = projectors[i]*v_nonzero[i]
    fmiss = f - np.sum(projected_evectors, axis=0)
    return w_nonzero, v_nonzero, projectors, f, fmiss, w_max, w, all_projectors

def cutoff(theory_shift_test, eigenvalue_cutoff:(bool, type(None)) = None,
           use_analytic:(bool, type(None)) = None):
    if eigenvalue_cutoff == True:
        cutoff = "10 times modulus of largest 'negative' eigenvalue"
    else:
        cutoff = "Eigenvalues determined by projection onto space of non-zero eigenvalues."
    if use_analytic == True:
        cutoff = cutoff + "\n Linearly independent vectors determined analytically"
    print(f"cutoff = {cutoff}")
    return cutoff

@table
def theory_covmat_eigenvalues(theory_shift_test):
    w_nonzero, v_nonzero, projectors = theory_shift_test[:3]
    s_scrambled = np.sqrt(np.abs(w_nonzero))
    projectors_scrambled = np.ndarray.tolist(projectors)
    ratio_scrambled = projectors_scrambled/s_scrambled
    table = pd.DataFrame([s_scrambled[::-1], projectors_scrambled[::-1], ratio_scrambled[::-1]],
         		index = [r'$s_a$', r'$\delta_a$', r'$\delta_a/s_a$'],
                 columns = np.arange(1,len(s_scrambled)+1,1))
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
def projector_eigenvalue_ratio(theory_shift_test,
                               eigenvalue_cutoff:(bool, type(None)) = None,
                               use_analytic:(bool, type(None)) = None):
    surviving_evals = theory_shift_test[0][::-1]
    all_projectors = theory_shift_test[7][::-1]
    all_evals = theory_shift_test[6][::-1]
    fmiss = theory_shift_test[4]
    fmiss_mod = np.sqrt(np.sum(fmiss**2))
    ratio = np.abs(all_projectors)/np.sqrt(np.abs(all_evals))
    # Initialise array of zeros and set precision to same as FK tables
    masked_evals = np.zeros((len(all_evals)), dtype=np.float32)
    for loc, eval in enumerate(all_evals):
        if eval in surviving_evals:
            masked_evals[loc] = eval
     # Ordering according to shift vector
    mask = np.argsort(np.abs(all_projectors))[::-1]
    all_evals = np.asarray(all_evals)[mask]
    all_projectors = all_projectors[mask]
    ratio = ratio[mask]
    masked_evals = masked_evals[mask]
    xvals = np.arange(1,len(masked_evals)+1,1)
    fig, (ax1, ax2) = plt.subplots(2, figsize=(5,5))
    ax1.plot(xvals, np.abs(all_projectors), 's', label = r'|$\delta_a$|')
    ax1.plot(xvals, np.sqrt(np.abs(all_evals)), 'o', label = r'$|s_a|$')
    if use_analytic == None:
        ax1.plot(xvals, np.sqrt(np.abs(masked_evals)), 'o', label = r'surviving $|s_a|$', color='k')
    ax1.plot(0, fmiss_mod, '*', label=r'$|\delta_{miss}|$', color='b')
    ax2.plot(xvals,ratio, 'D', color="red")
    ax2.plot(0,0, '.', color="w")
    ax1.set_title(f"Number of eigenvalues = {len(surviving_evals)}", fontsize=10)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    if eigenvalue_cutoff == True:
        ax1.set_xscale('log')
        ax2.set_xscale('log')
    ax1.legend()
    labels = [item.get_text() for item in ax1.get_xticklabels()]
    ax1.set_xticklabels(labels)
    ax2.set_xticklabels(labels)
    ax2.axhline(y=3, color='k', label=r'|$\delta_a$/$s_a$| = 3')
    ax2.legend()
    ax2.set_ylabel(r"|$\delta_a$/$s_a$|")
    print(f"Subspace dimension = {len(all_evals)}")
    return fig

@figure
def shift_diag_cov_comparison(shx_vector, thx_covmat, thx_vector):
    matrix = thx_covmat[0]/(np.outer(thx_vector[0], thx_vector[0]))
    fnorm = -shx_vector[0]
    indexlist = list(matrix.index.values)
    # adding process index for plotting
    dsnames = []
    processnames= []
    ids = []
    for index in indexlist:
        name = index[0]
        i = index[1]
        dsnames.append(name)
        ids.append(i)
        proc = _process_lookup(name)
        processnames.append(proc)
    tripleindex = pd.MultiIndex.from_arrays([processnames, dsnames, ids],
                        names = ("process", "dataset", "id"))
    matrix = pd.DataFrame(matrix.values, index=tripleindex, columns=tripleindex)
    matrix.sort_index(0, inplace=True)
    matrix.sort_index(1, inplace=True)
    oldindex=matrix.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    matrix = matrix.reindex(newindex)
    matrix = (matrix.T.reindex(newindex)).T
    sqrtdiags = np.sqrt(np.diag(matrix))
    fnorm = pd.DataFrame(fnorm.values, index=tripleindex)
    fnorm.sort_index(0, inplace=True)
    fnorm = fnorm.reindex(newindex)
    fig, ax = plt.subplots(figsize=(20,10))
    ax.plot(sqrtdiags*100, '.-', label="Theory", color = "red")
    ax.plot(-sqrtdiags*100, '.-', color = "red")
    ax.plot(fnorm.values*100, '.-', label="NNLO-NLO Shift", color = "black")
    ticklocs, ticklabels, startlocs = matrix_plot_labels(matrix)
    plt.xticks(ticklocs, ticklabels, rotation=45, fontsize=20)
    ax.vlines(startlocs, -70, 70, linestyles='dashed')
    ax.margins(x=0, y=0)
    ax.set_ylabel("% of central theory", fontsize=20)
    ax.legend(fontsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    return fig

@figure
def plot_shift_scaleavg_comparison(shx_vector, thx_vector,
                                   allthx_vector, thx_covmat):
    diffs = [((thx_vector[0] - scalevarvector)/thx_vector[0])
                                        for scalevarvector in allthx_vector[0]]
    diffsconcat = pd.concat(diffs, axis=1)
    avgdiffs = diffsconcat.mean(axis=1)
    fig, ax = plt.subplots(figsize=(20,10))
    ax.plot(100*avgdiffs.values, '.-',
            label=r"Average of $\Delta$s",color = "blue")
    ax.plot(-100*shx_vector[0].values, '.-', label="NNLO-NLO shift", color = "black")
    ticklocs, ticklabels, startlocs = matrix_plot_labels(thx_covmat[0])
    plt.xticks(ticklocs, ticklabels, rotation=45, fontsize=20)
    ax.set_ylabel("% of central theory", fontsize=20)
    ax.legend(fontsize=20)
    return fig
