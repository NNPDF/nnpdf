# -*- coding: utf-8 -*-
"""
theorycovariance_tests.py
Tools for testing theory covariance matrices and their properties.
"""
from __future__ import generator_stop

import logging

from collections import namedtuple
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd

from reportengine.checks import make_argcheck
from reportengine.figure import figure
from reportengine.table import table
from reportengine import collect

from validphys.results import results
from validphys.checks import check_two_dataspecs
from validphys.theorycovariance import process_lookup, combine_by_type, process_starting_points
from validphys.theorycovariance import covmap, covs_pt_prescrip, theory_covmat_custom
from validphys.theorycovariance_output import matrix_plot_labels

log = logging.getLogger(__name__)

matched_dataspecs_results = collect('results', ['dataspecs'])

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
    xs = [(thx_vector[0] - scalevarvector)/thx_vector[0] for scalevarvector in allthx_vector[0]]
    # iteratively orthogonalising deltas
    ys = [x/np.linalg.norm(x) for x in xs]
    xdashs = [None]*len(ys)
    for n, x in enumerate(xs):
        sub_ys = ys[:n]
        subtract_terms = [None]*len(sub_ys)
        xlist = [x]*len(sub_ys)
        for i in range(len(sub_ys)):
            subtract_terms[i] = np.dot(sub_ys[i], np.dot(sub_ys[i].T, xlist[i]))
        xdashs[n] = x - np.sum(subtract_terms, axis=0)
        ys[n] = xdashs[n]/np.linalg.norm(xdashs[n])
#    xdash = xs[1] - ys[0]*np.dot(ys[0].T, xs[1])[0]
#    ys[1] = xdash/np.linalg.norm(xdash)
    P = np.column_stack(ys)
    projected_matrix = np.dot(P.T, np.dot(orig_matrix, P))
    w, v_projected = la.eigh(projected_matrix)
    v = np.dot(P, v_projected)
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
#    if num_evals is not None:
#        w_nonzero = w[-num_evals:]
#        nonzero_locs = range(len(w)-num_evals, len(w))
#    elif evalue_cutoff is not None:
#        w_nonzero = w[w>evalue_cutoff*w_max]
#        nonzero_locs = np.nonzero(w>evalue_cutoff*w_max)[0]
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
    projected_evectors = np.zeros((len(projectors), (len(f))))
    for i in range(len(projectors)):
        projected_evectors[i] = projectors[i]*v_nonzero[i]
    fmiss = f - np.sum(projected_evectors, axis=0)
    return w_nonzero, v_nonzero, projectors, f, fmiss, w_max, w, all_projectors

def cutoff(theory_shift_test, eigenvalue_cutoff:(bool, type(None)) = None):
    w_max = theory_shift_test[5]
    if eigenvalue_cutoff == True:
        cutoff = "10 times modulus of largest 'negative' eigenvalue"
    else:
        cutoff = "Eigenvalues determined by projection onto space of non-zero eigenvalues"
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
