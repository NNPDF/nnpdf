# -*- coding: utf-8 -*-
"""
theorycovariance.py

Tools for constructing and studying theory covariance matrices.
"""
from __future__ import generator_stop

import logging

from IPython import embed

import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors, rcParams as rc
import pandas as pd

from reportengine.figure import figure
from reportengine.checks import make_argcheck, CheckError
from reportengine.table import table
from reportengine import collect

from validphys.results import results, experiment_results, experiments_central_values, Chi2Data, experiments_chi2_table
from validphys.calcutils import all_chi2, central_chi2, calc_chi2, all_chi2_theory, central_chi2_theory
from validphys import plotutils

log = logging.getLogger(__name__)



theoryids_experiments_central_values = collect(experiments_central_values, ('theoryids',))

@make_argcheck
def check_have_three_theories(theoryids):
    l = len(theoryids)
    if l!=3:
        raise CheckError(f"Expecting exactly 3 theories, but got {l}.")

def abs_chi2_data_theory_dataset(each_dataset_results, theory_covmat_datasets_3pt):
    """ Returns an array of tuples (member_chi², central_chi², numpoints)
    corresponding to each data set, where theory errors are included"""
    for i, results in enumerate(each_dataset_results):
        data_result, th_result = results
        covmat = theory_covmat_datasets_3pt[i]
        chi2s = all_chi2_theory(results, covmat)

        central_result = central_chi2_theory(results, covmat)

        if i==0:
            chi2data_array = [Chi2Data(th_result.stats_class(chi2s[:,np.newaxis]),
                                        central_result, len(data_result))]
        else:
            chi2data_array.append(
                          Chi2Data(th_result.stats_class(chi2s[:,np.newaxis]),
                                  central_result, len(data_result)))
    return chi2data_array 

def abs_chi2_data_theory_experiment(experiments_results, theory_covmat_experiments_3pt):
    """ Like abs_chi2_data_theory_dataset but for experiments not datasets"""
    for i, results in enumerate(experiments_results):
        data_result, th_result = results
        covmat = theory_covmat_experiments_3pt[i]

        chi2s = all_chi2_theory(results, covmat)

        central_result = central_chi2_theory(results, covmat)

        if i==0:
            chi2data_array = [Chi2Data(th_result.stats_class(chi2s[:,np.newaxis]),
                                        central_result, len(data_result))]
        else:
            chi2data_array.append(
                          Chi2Data(th_result.stats_class(chi2s[:,np.newaxis]),
                                  central_result, len(data_result)))                         

    embed()
    return chi2data_array 

@table
def experiments_chi2_table_theory(experiments, pdf, abs_chi2_data_theory_experiment,
                                  abs_chi2_data_theory_dataset):
    """Same as experiments_chi2_table but including theory covariance matrix"""
    return experiments_chi2_table(experiments, pdf, abs_chi2_data_theory_experiment,
                                abs_chi2_data_theory_dataset)

@table
@check_have_three_theories
def theory_covmat_3pt(theoryids_experiments_central_values, experiments, experiments_index):
    """Calculates the theory covariance matrix for 3-point scale variations."""
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
    print(np.shape(each_dataset_results_theory))
    for dataset in each_dataset_results_theory:
        data_centrals = [x[0].central_value for x in dataset]
        theory_centrals = [x[1].central_value for x in dataset]
        central, low, high = theory_centrals
        lowdiff = low - central
        highdiff = high - central
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
    experiments_results_theory = np.swapaxes(experiments_results_theory, 0, 1)
    for experiment in experiments_results_theory:
        data_centrals = [x[0].central_value for x in experiment]
        theory_centrals = [x[1].central_value for x in experiment]
        central, low, high = theory_centrals
        lowdiff = low - central
        highdiff = high - central
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

def matrix_plot_labels(df):
    explabels  = [list(df)[x][0] for x in range(len(list(df)))]
    datlabels  = [list(df)[x][1] for x in range(len(list(df)))]
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
    fig,ax = plt.subplots()
    matrixplot = ax.matshow(matrix*100, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.01, linscale=10, vmin=-100*matrix.max(), vmax=100*matrix.max()))
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
    fig, ax = plt.subplots()
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    cbar = fig.colorbar(matrixplot)
    ax.set_title('Experiment correlation matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_normthcovmat_heatmap(theory_normcovmat_3pt):
    """Matrix plot of the theory covariance matrix for 3-point scale variations normalised to data."""
    df = theory_normcovmat_3pt
    matrix = df.as_matrix()
    fig,ax = plt.subplots()
    matrixplot = ax.matshow(matrix*100, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.1, linscale=10, vmin=-100*matrix.max(), vmax=100*matrix.max()))
    cbar = fig.colorbar(matrixplot, label="% of data")
    ax.set_title('Theory covariance matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_thcorrmat_heatmap(theory_corrmat_3pt):
    """Matrix plot of the theory correlation matrix"""
    df = theory_corrmat_3pt
    matrix = df.as_matrix()
    fig, ax = plt.subplots()
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    cbar = fig.colorbar(matrixplot)
    ax.set_title('Theory correlation matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_normexpplusthcovmat_heatmap(experimentsplustheory_normcovmat_3pt):
    """Matrix plot of the exp + theory covariance matrix normalised to data"""
    df = experimentsplustheory_normcovmat_3pt
    matrix = experimentsplustheory_normcovmat_3pt.as_matrix()
    fig, ax = plt.subplots()
    matrixplot = ax.matshow(matrix*100, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.1, linscale=10, vmin=-100*matrix.max(), vmax=100*matrix.max()))
    cbar = fig.colorbar(matrixplot, label="% of data")
    ax.set_title('Experiment + theory covariance matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_expplusthcorrmat_heatmap(experimentsplustheory_corrmat_3pt):
    """Matrix plot of the exp + theory correlation matrix"""
    df = experimentsplustheory_corrmat_3pt
    matrix = experimentsplustheory_corrmat_3pt.as_matrix()
    fig, ax = plt.subplots()
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    cbar = fig.colorbar(matrixplot)
    ax.set_title('Experiment + theory correlation matrix')
    ticklocs, ticklabels = matrix_plot_labels(df)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_covdiff_heatmap(theory_covmat_3pt, experiments_covmat):
    """Matrix plot (thcov + expcov)/expcov"""
    df_theory = theory_covmat_3pt
    df_experiment = experiments_covmat
    matrix_theory = df_theory.as_matrix()
    matrix_experiment = df_experiment.as_matrix()
    matrix = (matrix_theory+matrix_experiment)/np.mean(matrix_experiment)
    fig,ax = plt.subplots()
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, norm=mcolors.SymLogNorm(linthresh=0.1, linscale=10, vmin=-matrix.max(), vmax=matrix.max()))
    cbar = fig.colorbar(matrixplot)
    ax.set_title('(Theory + experiment)/mean(experiment) covariance matrices')
    ticklocs, ticklabels = matrix_plot_labels(df_experiment)
    plt.xticks([])
    plt.yticks(ticklocs, ticklabels)
    return fig

@figure
def plot_diag_cov_comparison(theory_covmat_3pt, experiments_covmat, experiments_data):
    """Plot of sqrt(cov_ii)/data_i for cov = exp, theory, exp+theory"""
    data = experiments_data
    df_theory = theory_covmat_3pt
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
    ax.set_ylabel(r"$\frac{\sqrt{cov_{ii}}}{D_i}$")
    ax.set_title("Diagonal elements of normalised covariance matrices")
    ax.legend()
    return fig

@figure
def plot_diag_cov_impact(theory_covmat_3pt, experiments_covmat, experiments_data):
    """Plot ((expcov)^-1_ii)^-0.5 versus ((expcov + thcov)^-1_ii)^-0.5"""
    data = experiments_data
    df_theory = theory_covmat_3pt
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
def plot_theory_error_test(theory_covmat_3pt, experiments_covmat, experiments_data, theoryids_experiments_central_values):
    rc.update({'font.size': 30})
    data = experiments_data.as_matrix()
    df_theory = theory_covmat_3pt
    df_experiment = experiments_covmat
    matrix_theory = df_theory.as_matrix()
    matrix_experiment = df_experiment.as_matrix()
    central, low, high = np.array(theoryids_experiments_central_values)
    experrors = np.sqrt(np.diag(matrix_experiment))
    theoryerrors = np.sqrt(np.diag(matrix_theory))
    fig,ax = plt.subplots(figsize=(20, 10))
    ax.plot(central/data, label="central", color="red")
    ax.plot(low/data, label="low",color="blue")
    ax.plot(high/data, label="high",color="blue")
    ax.errorbar(np.arange(len(data)),data/data, yerr=experrors/data,fmt='--o', label="experiment errors",color="black")
    ax.errorbar(np.arange(len(data))+0.25,data/data, yerr=theoryerrors/data,fmt='none', label="theory errors",color="green")
    ticklocs, ticklabels = matrix_plot_labels(df_experiment)
#    plt.xticks(ticklocs, ticklabels, rotation="vertical")
    ax.set_ylabel("Observable normalised to experiment")
    ax.set_title("Theory error comparison")
    ax.legend()
 #   ax.set_ylim(-1,3)
 #   ax.set_xlim(2145,2167)
    plt.show()
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
#            xticks.append(dataset.name)
    plotvalues = np.stack((dschi2theory, dschi2))
    fig,ax = plotutils.barplot(plotvalues, collabels=xticks,
                               datalabels=["experiment + theory","experiment"])

    ax.set_title(r"$\chi^2$ distribution for datasets")
    ax.legend(fontsize=14)

    return fig
