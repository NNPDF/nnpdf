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

from reportengine.figure import figure
from reportengine.checks import make_argcheck, CheckError
from reportengine.table import table
from reportengine import collect

from validphys.results import results, experiment_results, experiments_central_values
from validphys.results import Chi2Data, experiments_chi2_table
from validphys.calcutils import all_chi2_theory, central_chi2_theory
from validphys import plotutils

log = logging.getLogger(__name__)

theoryids_experiments_central_values = collect(experiments_central_values, ('theoryids',))

@make_argcheck
def _check_three_or_seven_theories(theoryids):
    l = len(theoryids)
    if l!=3 and l!=7:
        raise CheckError(f"Expecting exactly 3 or 7 theories, but got {l}.")


def make_scale_var_covmat(theory_centrals):
    """Takes N theories at different scales and applies N-pt scale variations
    to produce a covariance matrix """
    central, *others = theory_centrals
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

<<<<<<< HEAD
def theory_block_diag_covmat(theory_covmat_datasets, experiments_index):
    """Takes the theory covariance matrices for individual datasets and
    returns a data frame with a block diagonal theory covariance matrix
    by dataset"""
    s  = la.block_diag(*theory_covmat_datasets)
    df = pd.DataFrame(s, index=experiments_index, columns=experiments_index)   
    return df

experiments_results = collect(experiment_results, ('experiments',))
experiments_results_theory = collect('experiments_results', ('theoryids',))

def theory_covmat_experiments(experiments_results_theory, make_scale_var_covmat):
    """Same as theory_covmat_datasets but per experiment rather than
    per dataset. Needed for calculation of chi2 per experiment."""
    exp_covmats = []
    for exp_result in zip(*experiments_results_theory):
        theory_centrals = [x[1].central_value for x in exp_result]
        s = make_scale_var_covmat(theory_centrals)
        sigma = exp_result[0][0].covmat
        cov = s + sigma
        exp_result_covmats.append(cov)
    return exp_result_covmats

@table
def theory_corrmat(theory_covmat):
    """Calculates the theory correlation matrix for scale variations."""
    df = theory_covmat
    covmat = df.as_matrix()
    diag_minus_half = (np.diagonal(covmat))**(-0.5)
    mat = diag_minus_half[:,np.newaxis]*df*diag_minus_half
    return mat

@table
def theory_blockcorrmat(theory_block_diag_covmat):
    """Calculates the theory correlation matrix for scale variations 
    with block diagonal entries by dataset only"""
    df = theory_block_diag_covmat
    covmat = df.as_matrix()
    diag_minus_half = (np.diagonal(covmat))**(-0.5)
    mat = diag_minus_half[:,np.newaxis]*df*diag_minus_half
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
def experimentsplustheory_normcovmat(experiments_covmat, theory_covmat, experiments_data):
    """Calculates the experiment + theory covariance matrix for scale
       variations normalised to data."""
    df = experiments_covmat + theory_covmat
    experiments_data_array = np.array(experiments_data)
    mat = df/np.outer(experiments_data_array, experiments_data_array)
    return mat

@table
def experimentsplusblocktheory_normcovmat(experiments_covmat, theory_block_diag_covmat, experiments_data):
    """Calculates the experiment + theory covariance matrix for scale
       variations normalised to data, block diagonal by data set."""
    df = experiments_covmat + theory_block_diag_covmat
    experiments_data_array = np.array(experiments_data)
    mat = df/np.outer(experiments_data_array, experiments_data_array)
    return mat

@table
def experimentsplustheory_corrmat(experiments_covmat, theory_covmat):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices."""
    total_df = experiments_covmat + theory_covmat
    total_cov = (experiments_covmat + theory_covmat).as_matrix()
    diag_minus_half = (np.diagonal(total_cov))**(-0.5)
    corrmat = diag_minus_half[:,np.newaxis]*total_df*diag_minus_half
    return corrmat

@table
def experimentsplusblocktheory_corrmat(experiments_covmat, theory_block_diag_covmat):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices, block diagonal by dataset."""
    total_df = experiments_covmat + theory_block_diag_covmat
    total_cov = (experiments_covmat + theory_block_diag_covmat).as_matrix()
    diag_minus_half = (np.diagonal(total_cov))**(-0.5)
    corrmat = diag_minus_half[:,np.newaxis]*total_df*diag_minus_half
    return corrmat

def chi2_impact(theory_covmat, experiments_covmat, experiments_results):
    """ Returns total chi2 including theory cov mat """
    dataresults = [ x[0] for x in experiments_results ]
    theoryresults = [ x[1] for x in experiments_results ]
    dat_central_list = [x.central_value for x in dataresults]
    th_central_list = [x.central_value for x in theoryresults]
    dat_central = np.concatenate([x for x in dat_central_list])
    th_central  = np.concatenate([x for x in th_central_list])
    central_diff = dat_central - th_central
    cov = theory_covmat.as_matrix() + experiments_covmat.as_matrix()
    elements = np.dot(central_diff.T,np.dot(la.inv(cov),central_diff))
    chi2 = (1/len(central_diff))*np.sum(elements)
    return chi2

def chi2_block_impact(theory_block_diag_covmat, experiments_covmat, experiments_results):
    """ Returns total chi2 including theory cov mat """
    dataresults = [ x[0] for x in experiments_results ]
    theoryresults = [ x[1] for x in experiments_results ]
    dat_central_list = [x.central_value for x in dataresults]
    th_central_list = [x.central_value for x in theoryresults]
    dat_central = np.concatenate([x for x in dat_central_list])
    th_central  = np.concatenate([x for x in th_central_list])
    central_diff = dat_central - th_central
    cov = theory_block_diag_covmat.as_matrix() + experiments_covmat.as_matrix()
    elements = np.dot(central_diff.T,np.dot(la.inv(cov),central_diff))
    chi2 = (1/len(central_diff))*np.sum(elements)
    return chi2

def chi2_diag_only(theory_covmat, experiments_covmat, experiments_results):
    """ Returns total chi2 including only diags of theory cov mat """
    dataresults = [ x[0] for x in experiments_results ]
    theoryresults = [ x[1] for x in experiments_results ]
    dat_central_list = [x.central_value for x in dataresults]
    th_central_list = [x.central_value for x in theoryresults]
    dat_central = np.concatenate(dat_central_list)
    th_central  = np.concatenate(th_central_list)
    central_diff = dat_central - th_central
    s = theory_covmat.as_matrix()
    s_diag = np.zeros((len(central_diff),len(central_diff)))
    np.fill_diagonal(s_diag, np.diag(s))
    cov = s_diag + experiments_covmat.as_matrix()
    elements = np.dot(central_diff.T,np.dot(la.inv(cov),central_diff))
    chi2 = (1/len(central_diff))*np.sum(elements)
    return chi2

each_dataset_results = collect(results, ('experiments', 'experiment'))

def abs_chi2_data_theory_dataset(each_dataset_results, theory_covmat_datasets):
    """ Returns an array of tuples (member_chi², central_chi², numpoints)
    corresponding to each data set, where theory errors are included"""
    chi2data_array = []
    for i, results in enumerate(each_dataset_results):
        data_result, th_result = results
        covmat = theory_covmat_datasets[i]
        chi2s = all_chi2_theory(results, covmat)
        central_result = central_chi2_theory(results, covmat)
        chi2data_array.append(Chi2Data(th_result.stats_class(chi2s[:,np.newaxis]),
                                   central_result, len(data_result)))
    return chi2data_array

def abs_chi2_data_theory_experiment(experiments_results, theory_covmat_experiments):
    """ Like abs_chi2_data_theory_dataset but for experiments not datasets"""
    chi2data_array = []
    for i, results in enumerate(experiments_results):
        data_result, th_result = results
        covmat = theory_covmat_experiments[i]
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
    explabels  = [list(df)[x][0] for x in range(len(list(df)))]
    points     = [list(df)[x][2] for x in range(len(list(df)))]
    unique_exp = [[0 for x in range(2)] for y in range(len(explabels))]
    unique_exp[0] = [explabels[0],points[0]]
    i=1
    for x in range(len(explabels)-1):
        if explabels[x+1] != explabels[x]:
            unique_exp[i] = [explabels[x+1],x+1]
            i=i+1
    unique_exp = [sublist for i, sublist in enumerate(unique_exp) if sublist[0] != 0]
    ticklabels = [unique_exp[x][0] for x in range(len(unique_exp))]
    startlocs = [unique_exp[x][1] for x in range(len(unique_exp))]
    startlocs += [len(explabels)]
    ticklocs = [0 for x in range(len(startlocs)-1)]
    for i in range(len(startlocs)-1):
        ticklocs[i] = 0.5*(startlocs[i+1]+startlocs[i])
    print("Experiment names:   " + str(ticklabels))
    print("Datapoint start locations:   " + str(startlocs))
    return ticklocs, ticklabels

@figure
def plot_normexpcovmat_heatmap(experiments_normcovmat):
    """Matrix plot of the experiment covariance matrix normalised to data."""
    df = experiments_normcovmat
    matrix = df.as_matrix()
    fig,ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix*100, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.01,
                            linscale=10, vmin=-100*matrix.max(), vmax=100*matrix.max()))
    cbar = fig.colorbar(matrixplot, label="% of data")
    ax.set_title('Experiment covariance matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig


@figure
def plot_expcorrmat_heatmap(experiments_corrmat):
    """Matrix plot of the experiment correlation matrix"""
    df = experiments_corrmat
    matrix = df.as_matrix()
    fig, ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    cbar = fig.colorbar(matrixplot)
    ax.set_title('Experiment correlation matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_normthcovmat_heatmap(theory_normcovmat):
    """Matrix plot of the theory covariance matrix for
     3/7-point scale variations normalised to data."""
    df = theory_normcovmat
    matrix = df.as_matrix()
    fig,ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix*100, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.1,
                            linscale=10, vmin=-100*matrix.max(), vmax=100*matrix.max()))
    cbar = fig.colorbar(matrixplot, label="% of data")
    ax.set_title('Theory covariance matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_normthblockcovmat_heatmap(theory_normblockcovmat):
    df = theory_normblockcovmat
    matrix = df.as_matrix()
    fig,ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix*100, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.1, linscale=10, vmin=-100*matrix.max(), vmax=100*matrix.max()))
    cbar = fig.colorbar(matrixplot, label="% of data")
    ax.set_title('Block diagonal theory covariance matrix by dataset')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_thcorrmat_heatmap(theory_corrmat):
    """Matrix plot of the theory correlation matrix"""
    df = theory_corrmat
    matrix = df.as_matrix()
    fig, ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    cbar = fig.colorbar(matrixplot)
    ax.set_title('Theory correlation matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_thblockcorrmat_heatmap(theory_blockcorrmat):
    """Matrix plot of the theory correlation matrix"""
    df = theory_blockcorrmat
    matrix = df.as_matrix()
    fig, ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    cbar = fig.colorbar(matrixplot)
    ax.set_title('Theory correlation matrix block diagonal by dataset')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_normexpplusthcovmat_heatmap(experimentsplustheory_normcovmat):
    """Matrix plot of the exp + theory covariance matrix normalised to data"""
    df = experimentsplustheory_normcovmat
    matrix = experimentsplustheory_normcovmat.as_matrix()
    fig, ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix*100, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.1,
                            linscale=10, vmin=-100*matrix.max(), vmax=100*matrix.max()))
    cbar = fig.colorbar(matrixplot, label="% of data")
    ax.set_title('Experiment + theory covariance matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_normexpplusblockthcovmat_heatmap(experimentsplusblocktheory_normcovmat):
    """Matrix plot of the exp + theory covariance matrix normalised to data"""
    df = experimentsplusblocktheory_normcovmat
    matrix = df.as_matrix()
    fig, ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix*100, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.1,
                            linscale=10, vmin=-100*matrix.max(), vmax=100*matrix.max()))
    cbar = fig.colorbar(matrixplot, label="% of data")
    ax.set_title('Experiment + theory (block diagonal by dataset) covariance matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_expplusthcorrmat_heatmap(experimentsplustheory_corrmat):
    """Matrix plot of the exp + theory correlation matrix"""
    df = experimentsplustheory_corrmat
    matrix = experimentsplustheory_corrmat.as_matrix()
    fig, ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    cbar = fig.colorbar(matrixplot)
    ax.set_title('Experiment + theory correlation matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_expplusblockthcorrmat_heatmap(experimentsplusblocktheory_corrmat):
    """Matrix plot of the exp + theory correlation matrix"""
    df = experimentsplusblocktheory_corrmat
    matrix = df.as_matrix()
    fig, ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    cbar = fig.colorbar(matrixplot)
    ax.set_title('Experiment + theory (block diagonal by dataset) correlation matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_covdiff_heatmap(theory_covmat, experiments_covmat):
    """Matrix plot (thcov + expcov)/expcov"""
    df_theory = theory_covmat
    df_experiment = experiments_covmat
    matrix_theory = df_theory.as_matrix()
    matrix_experiment = df_experiment.as_matrix()
    matrix = (matrix_theory+matrix_experiment)/np.mean(matrix_experiment)
    fig,ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.1,
                            linscale=10, vmin=-matrix.max(), vmax=matrix.max()))
    cbar = fig.colorbar(matrixplot)
    ax.set_title('(Theory + experiment)/mean(experiment) covariance matrices')
    ticklocs, ticklabels = matrix_plot_labels(df_experiment)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_blockcovdiff_heatmap(theory_block_diag_covmat, experiments_covmat):
    """Matrix plot (thcov + expcov)/expcov"""
    df_theory = theory_block_diag_covmat
    df_experiment = experiments_covmat
    matrix_theory = df_theory.as_matrix()
    matrix_experiment = df_experiment.as_matrix()
    matrix = (matrix_theory+matrix_experiment)/np.mean(matrix_experiment)
    fig,ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.1,
                            linscale=10, vmin=-matrix.max(), vmax=matrix.max()))
    cbar = fig.colorbar(matrixplot)
    ax.set_title('(Theory + experiment)/mean(experiment) covariance matrices for block diagonal theory covmat by dataset')
    ticklocs, ticklabels = matrix_plot_labels(df_experiment)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_diag_cov_comparison(theory_covmat, experiments_covmat, experiments_data):
    """Plot of sqrt(cov_ii)/|data_i| for cov = exp, theory, exp+theory"""
    data = np.abs(experiments_data)
    df_theory = theory_covmat
    df_experiment = experiments_covmat
    df_total = df_theory + df_experiment
    sqrtdiags1 = np.sqrt(np.diag(df_theory.as_matrix()))
    sqrtdiags2 = np.sqrt(np.diag(df_experiment.as_matrix()))
    sqrtdiags3 = np.sqrt(np.diag(df_total.as_matrix()))
    fig,ax = plt.subplots(figsize=(20,10))
    ax.plot((sqrtdiags2/data).as_matrix(), '.', label="Experiment", color="orange")
    ax.plot((sqrtdiags1/data).as_matrix(), '.', label="Theory", color = "red")
    ax.plot((sqrtdiags3/data).as_matrix(), '.', label="Total", color = "blue")
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
    matrix_theory = df_theory.as_matrix()
    matrix_experiment = df_experiment.as_matrix()
    a = (np.diag(la.inv(matrix_experiment)))**(-0.5)
    b = (np.diag(la.inv(matrix_theory+matrix_experiment)))**(-0.5)
    fig,ax = plt.subplots()
    ax.plot((a/data).as_matrix(), '.', label="Experiment", color="orange")
    ax.plot((b/data).as_matrix(), '.', label="Experiment + Theory", color="mediumseagreen")
    ticklocs, ticklabels = matrix_plot_labels(df_experiment)
    plt.xticks(ticklocs, ticklabels, rotation="vertical")
    ax.set_ylabel(r"$\frac{1}{D_i}\frac{1}{\sqrt{[cov^{-1}_]{ii}}}$")
    ax.set_title("Diagonal impact of adding theory covariance matrix")
    ax.legend()
    return fig

@figure
def plot_theory_error_test(theory_covmat, experiments_covmat, experiments_data,
                           theoryids_experiments_central_values):
    rc.update({'font.size': 30})
    data = experiments_data.as_matrix()
    df_theory = theory_covmat
    df_experiment = experiments_covmat
    matrix_theory = df_theory.as_matrix()
    matrix_experiment = df_experiment.as_matrix()
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
    for experiment, expres in zip(experiments, experiments_chi2):
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
