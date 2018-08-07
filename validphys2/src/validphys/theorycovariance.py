# -*- coding: utf-8 -*-
"""
theorycovariance.py

Tools for constructing and studying theory covariance matrices.
"""
from __future__ import generator_stop

import logging

import numpy as np
import scipy.linalg as la

import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors, rcParams as rc
import pandas as pd
from collections import defaultdict

from reportengine.figure import figure
from reportengine.checks import make_argcheck, CheckError
from reportengine.table import table
from reportengine import collect

from validphys.results import results, experiment_results, experiments_central_values
from validphys.results import Chi2Data, experiments_chi2_table
from validphys.calcutils import all_chi2_theory, central_chi2_theory
from validphys.plotoptions import get_info
from validphys import plotutils

from IPython import embed

log = logging.getLogger(__name__)

theoryids_experiments_central_values = collect(experiments_central_values, ('theoryids',))

@make_argcheck
def _check_three_or_seven_theories(theoryids):
    l = len(theoryids)
    if l!=3 and l!=7:
        raise CheckError(f"Expecting exactly 3 or 7 theories, but got {l}.")

@make_argcheck
def _check_three_theories(theoryids):
    l = len(theoryids)
    if l!=3:
        raise CheckError(f"Expecting exactly 3 theories, but got {l}.")

def make_scale_var_covmat(predictions):
    """Takes N theory predictions at different scales and applies N-pt scale variations
    to produce a covariance matrix."""
    central, *others = predictions
    diffs = (other - central for other in others)
    s = sum(np.outer(d,d) for d in diffs)/len(others)  
    return s

@table
@_check_three_or_seven_theories
def theory_covmat(theoryids_experiments_central_values, experiments_index):
    """Calculates the theory covariance matrix for scale variations.
    The matrix is a dataframe indexed by experiments_index."""
    s = make_scale_var_covmat(theoryids_experiments_central_values)
    df = pd.DataFrame(s, index=experiments_index, columns=experiments_index)
    return df

results_theoryids = collect(results,('theoryids',))
each_dataset_results_theory = collect('results_theoryids', ('experiments', 'experiment'))

@_check_three_or_seven_theories
def theory_covmat_datasets(each_dataset_results_theory):
    """Produces an array of theory covariance matrices. Each matrix corresponds
    to a different dataset, which must be specified in the runcard. """
    dataset_covmats=[]
    for dataset in each_dataset_results_theory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        dataset_covmats.append(s)
    return dataset_covmats

@_check_three_or_seven_theories
def total_covmat_datasets(each_dataset_results_theory):
    """Produces an array of total covariance matrices; the sum of experimental
    and scale-varied theory covariance matrices. Each matrix corresponds
    to a different dataset, which must be specified in the runcard.
    These are needed for calculation of chi2 per dataset. """
    dataset_covmats=[]
    for dataset in each_dataset_results_theory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        sigma = dataset[0][0].covmat
        cov = s + sigma
        dataset_covmats.append(cov)
    return dataset_covmats

def theory_block_diag_covmat(theory_covmat_datasets, experiments_index):
    """Takes the theory covariance matrices for individual datasets and
    returns a data frame with a block diagonal theory covariance matrix
    by dataset"""
    s  = la.block_diag(*theory_covmat_datasets)
    df = pd.DataFrame(s, index=experiments_index, columns=experiments_index)   
    return df

experiments_results = collect(experiment_results, ('experiments',))
experiments_results_theory = collect('experiments_results', ('theoryids',))

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

def process_lookup(experiments_xq2map, each_dataset_results_theory):
    """Produces a dictionary with keys corresponding to dataset names
    and values corresponding to process types"""
    dict = {}
    for experiment, commondata, fitted, masked in experiments_xq2map:
        info = get_info(commondata)
        key  = commondata.name
        if key in dict:
            pass
        else:
            dict.setdefault(key, [])
        dict[key].append(info.process_description)
    return dict

def dataset_names(experiments_xq2map):
    """Returns a list of the names of the datasets, in the same order as 
    they are inputted in the runcard"""
    names = []
    for experiment, commondata, fitted, masked in experiments_xq2map:
        info = get_info(commondata)
        name = commondata.name
        names.append(commondata.name)
    return names

def combine_by_type(process_lookup, each_dataset_results_theory, dataset_names):
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
    for dataset, name in zip(each_dataset_results_theory, dataset_names):
        theory_centrals = [x[1].central_value for x in dataset]
        dataset_size[name] = len(theory_centrals[0])
        proc_type = process_lookup[name][0]
        current_value = theories_by_process[proc_type]
        ordered_names[proc_type].append(name)
        if current_value == []:
            new_value = theory_centrals
        else:
            new_value = np.concatenate((current_value, theory_centrals), axis=1)
        theories_by_process[proc_type] = new_value 
    return theories_by_process, ordered_names, dataset_size

@_check_three_or_seven_theories
def theory_covmat_by_type(combine_by_type, theory_block_diag_covmat, experiments_index, dataset_names):
    """
    Calculates the theory covariance matrix for scale variations 
    with variations by process type.  Scale variations are carried out,
    correlating by process type. The resulting covariance matrix is then 
    reshuffled so it is ordered by dataset as they appear
    in the runcard, rather than by process type (i.e. it is then indexed
    in the same was as the experiment covariance matrix so can be easily
    combined and compared).

    Notes: 
    ------
    - cov_by_proc is the covariance matrix ordered by process type, 
      cov_by_exp is the covariance matrix orderd by experiment.
    - to reshuffle cov_by_proc to cov_by_exp, a mapping "map" is 
      constructed which maps the index of a datapoint in cov_by_proc
      onto its corresponding index in cov_by_exp. This is done by comparing
      the ordering of the datasets in cov_by_proc (according to ordered_names)
      with the ordering of datasets in cov_by_exp (according to dataset_names), 
      and using knowledge of the size of each dataset, stored in dataset_size.
    """
    theories_by_process, ordered_names, dataset_size = combine_by_type
    covmats = defaultdict(list)
    for process, theory_centrals in theories_by_process.items():
        s = make_scale_var_covmat(theory_centrals)
        covmats[process] = s
    covmats_list = covmats.values()
    cov_by_proc = la.block_diag(*covmats_list)
    start_exp = defaultdict(list)
    map = defaultdict(list)
    running_index = 0
    for dataset in dataset_names:
        size = dataset_size[dataset]
        start_exp[dataset] = running_index
        running_index += size
    start = 0
    names_by_proc_list = [item for sublist in ordered_names.values() for item in sublist]
    for dataset in names_by_proc_list:
        for i in range(dataset_size[dataset]):
            map[start+i] = start_exp[dataset] + i 
        start += dataset_size[dataset]
    cov_by_exp = np.zeros((len(cov_by_proc), len(cov_by_proc)))
    for i in range(len(cov_by_proc)):
        for j in range(len(cov_by_proc)):
            cov_by_exp[map[i]][map[j]] = cov_by_proc[i][j]
    df = pd.DataFrame(cov_by_exp, index=experiments_index, columns=experiments_index)
    return df

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
def theory_corrmat_by_type(theory_covmat_by_type):
    """Calculates the theory correlation matrix for scale variations 
    with variations by process type"""
    mat = theory_corrmat(theory_covmat_by_type)
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
def theory_normcovmat_by_type(theory_covmat_by_type, experiments_data):
    """Calculates the theory covariance matrix for scale variations 
    normalised to data, with variations by process type."""
    df = theory_covmat_by_type
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
def experimentsplusblocktheory_covmat(experiments_covmat, theory_block_diag_covmat):
    """Calculates the experiment + theory covariance matrix for scale variations."""
    df = experiments_covmat + theory_block_diag_covmat
    return df

@table
def experimentsplustheory_covmat_by_type(experiments_covmat, theory_covmat_by_type):
    """Calculates the experiment + theory covariance matrix for scale variations correlated
    by process type."""
    df = experiments_covmat + theory_covmat_by_type
    return df

@table
def experimentsplustheory_normcovmat(experiments_covmat, theory_covmat, experiments_data):
    """Calculates the experiment + theory covariance matrix for scale
       variations normalised to data."""
    df = experiments_covmat + theory_covmat
    experiments_data_array = np.array(experiments_data)
    mat = df/np.outer(experiments_data_array, experiments_data_array)
    return mat

@table
def experimentsplusblocktheory_normcovmat(experiments_covmat, theory_block_diag_covmat, experiments_data, experimentsplustheory_normcovmat):
    """Calculates the experiment + theory covariance matrix for scale
       variations normalised to data, block diagonal by data set."""
    mat = experimentsplustheory_normcovmat(experiments_covmat, theory_block_diag_covmat, experiments_data)
    return mat

@table
def experimentsplustheory_normcovmat_by_type(experiments_covmat, theory_covmat_by_type, experiments_data, experimentsplustheory_normcovmat):
    """Calculates the experiment + theory covariance matrix for scale
       variations normalised to data, correlations by process type."""
    mat = experimentsplustheory_normcovmat(experiments_covmat, theory_covmat_by_type, experiments_data)

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
def experimentsplusblocktheory_corrmat(experiments_covmat, theory_block_diag_covmat):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices, block diagonal by dataset."""
    corrmat = experimentsplustheory_corrmat(experiments_covmat, theory_block_diag_covmat)
    return corrmat

@table
def experimentsplustheory_corrmat_by_type(experiments_covmat, theory_covmat_by_type):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices, correlations by process type."""
    corrmat = experimentsplustheory_corrmat(experiments_covmat, theory_covmat_by_type)
    return corrmat

def chi2_impact(theory_covmat, experiments_covmat, experiments_results):
    """ Returns total chi2 including theory cov mat """
    dataresults, theoryresults = zip(*experiments_results)
    dat_central_list = [x.central_value for x in dataresults]
    th_central_list = [x.central_value for x in theoryresults]
    dat_central = np.concatenate(dat_central_list)
    th_central  = np.concatenate([x for x in th_central_list])
    central_diff = dat_central - th_central
    cov = theory_covmat.values + experiments_covmat.values
    elements = central_diff.T@(la.inv(cov)@central_diff)
    chi2 = (1/len(central_diff))*np.sum(elements)
    return chi2

def data_theory_diff(experiments_results):
    """Returns (D-T) for central theory, for use in chi2 calculations"""
    dataresults, theoryresults = zip(*experiments_results)
    dat_central_list = [x.central_value for x in dataresults]
    th_central_list = [x.central_value for x in theoryresults]
    dat_central = np.concatenate(dat_central_list)
    th_central  = np.concatenate(th_central_list)
    central_diff = dat_central - th_central
    return central_diff

def chi2_block_impact(theory_block_diag_covmat, experiments_covmat, experiments_results):
    """ Returns total chi2 including theory cov mat """
    chi2 = chi2_impact(theory_block_diag_covmat, experiments_covmat, experiments_results)
    return chi2

def chi2_impact_by_type(theory_covmat_by_type, experiments_covmat, experiments_results):
    """ Returns total chi2 including theory cov mat """
    chi2 = chi2_impact(theory_covmat_by_type, experiments_covmat, experiments_results)
    return chi2

def chi2_diag_only(theory_covmat, experiments_covmat, data_theory_diff):
    """ Returns total chi2 including only diags of theory cov mat """
    s = theory_covmat.values
    s_diag = np.zeros((len(data_theory_diff),len(data_theory_diff)))
    np.fill_diagonal(s_diag, np.diag(s))
    cov = s_diag + experiments_covmat.values
    elements = np.dot(data_theory_diff.T,np.dot(la.inv(cov),data_theory_diff))
    chi2 = (1/len(data_theory_diff))*np.sum(elements)
    return chi2

each_dataset_results = collect(results, ('experiments', 'experiment'))

def abs_chi2_data_theory_dataset(each_dataset_results, total_covmat_datasets):
    """ Returns an array of tuples (member_chi², central_chi², numpoints)
    corresponding to each data set, where theory errors are included"""
    chi2data_array = []
    for results, covmat in zip(each_dataset_results, total_covmat_datasets):
        data_result, th_result = results
        chi2s = all_chi2_theory(results,covmat)
        central_result = central_chi2_theory(results, covmat)
        chi2data_array.append(Chi2Data(th_result.stats_class(chi2s[:,np.newaxis]),
                                   central_result, len(data_result)))
    return chi2data_array

def abs_chi2_data_theory_experiment(experiments_results, total_covmat_experiments):
    """ Like abs_chi2_data_theory_dataset but for experiments not datasets"""
    chi2data_array = []
    for results, covmat in zip(experiments_results, total_covmat_experiments):
        data_result, th_result = results
        chi2s = all_chi2_theory(results, covmat)
        central_result = central_chi2_theory(results, covmat)
        chi2data_array.append(Chi2Data(th_result.stats_class(chi2s[:,np.newaxis]),
                              central_result, len(data_result)))
    return chi2data_array

@table
def experiments_chi2_table_theory(experiments, pdf, abs_chi2_data_theory_experiment,
                                  abs_chi2_data_theory_dataset):
    """Same as experiments_chi2_table but including theory covariance matrix"""
    return experiments_chi2_table(experiments, pdf, abs_chi2_data_theory_experiment,
                                abs_chi2_data_theory_dataset)

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
    matrixplot = ax.matshow(100*matrix, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.01,
                            linscale=10, vmin=-100*matrix.max(), vmax=100*matrix.max()))
    cbar = fig.colorbar(matrixplot, label="% of data")
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
    cbar = fig.colorbar(matrixplot)
    ax.set_title(title)
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks(ticklocs, ticklabels, rotation=30, ha="right")
    plt.gca().xaxis.tick_bottom()
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_normexpcovmat_heatmap(experiments_normcovmat):
    """Matrix plot of the experiment covariance matrix normalised to data."""
    fig = plot_covmat_heatmap(experiments_normcovmat, "Experiment covariance matrix")
    return fig

@figure
def plot_expcorrmat_heatmap(experiments_corrmat):
    """Matrix plot of the experiment correlation matrix"""
    fig = plot_corrmat_heatmap(experiments_corrmat, "Experiment correlation matrix")
    return fig

@figure
def plot_normthcovmat_heatmap(theory_normcovmat):
    """Matrix plot of the theory covariance matrix for
     3/7-point scale variations normalised to data."""
    fig = plot_covmat_heatmap(theory_normcovmat, "Theory covariance matrix")
    return fig

@figure
def plot_normthblockcovmat_heatmap(theory_normblockcovmat):
    """Matrix plot for block diagonal theory covariance matrix"""
    fig = plot_covmat_heatmap(theory_normblockcovmat, "Block diagonal theory covariance matrix by dataset")
    return fig

@figure
def plot_normthcovmat_heatmap_by_type(theory_normcovmat_by_type):
    """Matrix plot for block diagonal theory covariance matrix by process type"""
    fig = plot_covmat_heatmap(theory_normcovmat_by_type, "Block diagonal theory covariance matrix by process type")
    return fig

@figure
def plot_thcorrmat_heatmap(theory_corrmat):
    """Matrix plot of the theory correlation matrix"""
    fig = plot_corrmat_heatmap(theory_corrmat, "Theory correlation matrix")
    return fig

@figure
def plot_thblockcorrmat_heatmap(theory_blockcorrmat):
    """Matrix plot of the theory correlation matrix"""
    fig = plot_corrmat_heatmap(theory_blockcorrmat, "Theory correlation matrix block diagonal by dataset")
    return fig

@figure
def plot_thcorrmat_heatmap_by_type(theory_corrmat_by_type):
    """Matrix plot of the theory correlation matrix, correlations by process type"""
    fig = plot_corrmat_heatmap(theory_corrmat_by_type, "Theory correlation matrix block diagonal by process type")
    return fig

@figure
def plot_normexpplusthcovmat_heatmap(experimentsplustheory_normcovmat):
    """Matrix plot of the exp + theory covariance matrix normalised to data"""
    fig = plot_covmat_heatmap(experimentsplustheory_normcovmat, "Experiment + theory covariance matrix")
    return fig

@figure
def plot_normexpplusblockthcovmat_heatmap(experimentsplusblocktheory_normcovmat):
    """Matrix plot of the exp + theory covariance matrix normalised to data"""
    fig = plot_covmat_heatmap(experimentsplusblocktheory_normcovmat, "Experiment + theory (block diagonal by dataset) covariance matrix")
    return fig

@figure
def plot_normexpplusthcovmat_heatmap_by_type(experimentsplustheory_normcovmat_by_type):
    """Matrix plot of the exp + theory covariance matrix normalised to data"""
    fig = plot_covmat_heatmap(experimentsplustheory_normcovmat_by_type, "Experiment + theory (block diagonal by process type) covariance matrix")
    return fig

@figure
def plot_expplusthcorrmat_heatmap(experimentsplustheory_corrmat):
    """Matrix plot of the exp + theory correlation matrix"""
    fig = plot_corrmat_heatmap(experimentsplustheory_corrmat, "Experiment + theory correlation matrix")
    return fig

@figure
def plot_expplusblockthcorrmat_heatmap(experimentsplusblocktheory_corrmat):
    """Matrix plot of the exp + theory correlation matrix"""
    fig = plot_corrmat_heatmap(experimentsplusblocktheory_corrmat,"Experiment + theory (block diagonal by dataset) correlation matrix")
    return fig

@figure
def plot_expplusthcorrmat_heatmap_by_type(experimentsplustheory_corrmat_by_type):
    """Matrix plot of the exp + theory correlation matrix"""
    fig = plot_corrmat_heatmap(experimentsplustheory_corrmat_by_type,"Experiment + theory (block diagonal by process type) correlation matrix")
    return fig

@figure
def plot_covdiff_heatmap(theory_covmat, experiments_covmat):
    """Matrix plot (thcov + expcov)/expcov"""
    fig = plot_covmat_heatmap((theory_covmat.values+experiments_covmat.values)/np.mean(experiments_covmat.values),"(Theory + experiment)/mean(experiment) covariance matrices'")
    return fig

@figure
def plot_blockcovdiff_heatmap(theory_block_diag_covmat, experiments_covmat):
    """Matrix plot (thcov + expcov)/expcov"""
    fig = plot_covmat_heatmap((theory_block_diag_covmat.values+experiments_covmat.values)/np.mean(experiments_covmat.values),"(Theory + experiment)/mean(experiment) covariance matrices for block diagonal theory covmat by dataset")
    return fig

@figure
def plot_covdiff_heatmap_by_type(theory_covmat_by_type, experiments_covmat):
    """Matrix plot (thcov + expcov)/expcov"""
    fig = plot_covmat_heatmap((theory_covmat_by_type.values+experiments_covmat.values)/np.mean(experiments_covmat.values),"(Theory + experiment)/mean(experiment) covariance matrices for block diagonal theory covmat by process type")
    return fig

@figure
def plot_diag_cov_comparison(theory_covmat, experiments_covmat, experiments_data):
    """Plot of sqrt(cov_ii)/|data_i| for cov = exp, theory, exp+theory"""
    data = np.abs(experiments_data)
    df_theory = theory_covmat
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
    ax.set_title("Square root of diagonal elements of covariances matrices, " 
                 + "normalised to absolute value of data")
    ax.legend()
    return fig

@figure
def plot_diag_cov_impact(theory_covmat, experiments_covmat, experiments_data):
    """Plot ((expcov)^-1_ii)^-0.5 versus ((expcov + thcov)^-1_ii)^-0.5"""
    data = experiments_data
    df_theory = theory_covmat
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
    ax.set_title("Diagonal impact of adding theory covariance matrix")
    ax.legend()
    return fig

@_check_three_theories
@figure
def plot_theory_error_test(theory_covmat, experiments_covmat, experiments_data,
                           theoryids_experiments_central_values):
    """This is a test function which works for 3 point scale variations only"""
    rc.update({'font.size': 30})
    data = experiments_data.values
    df_theory = theory_covmat
    df_experiment = experiments_covmat
    matrix_theory = df_theory.values
    matrix_experiment = df_experiment.values
    central, low, high = np.array(theoryids_experiments_central_values)
    experrors = np.sqrt(np.diag(matrix_experiment))
    theoryerrors = np.sqrt(np.diag(matrix_theory))
    fig,ax = plt.subplots(figsize=(20, 10))
    ax.plot(central/data, label="central", color="red")
    ax.plot(low/data, label="low", color="blue")
    ax.plot(high/data, label="high", color="blue")
    ax.errorbar(np.arange(len(data)), data/data, yerr=experrors/data,fmt='--o',
                label="experiment errors", color="black")
    ax.errorbar(np.arange(len(data))+0.25, data/data, yerr=theoryerrors/data,fmt='none',
                label="theory errors", color="green")
    ax.set_ylabel("Observable normalised to experiment")
    ax.set_title("Theory error comparison")
    ax.legend()
    return fig

@figure
def plot_datasets_chi2_theory(experiments, experiments_chi2, each_dataset_chi2,
                              abs_chi2_data_theory_experiment, abs_chi2_data_theory_dataset):
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
    for experiment, expres in zip(experiments, abs_chi2_data_theory_experiment):
        for dataset, dsres in zip(experiment, dstheory):
            dschi2theory.append(dsres.central_result/dsres.ndata)
    plotvalues = np.stack((dschi2theory, dschi2))
    fig,ax = plotutils.barplot(plotvalues, collabels=xticks,
                               datalabels=["experiment + theory", "experiment"])
    ax.set_title(r"$\chi^2$ distribution for datasets")
    ax.legend(fontsize=14)
    return fig
