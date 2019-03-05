# -*- coding: utf-8 -*-
"""
output.py
Basic tools for plotting theory covariance matrices and their properties.
"""
from __future__ import generator_stop

import logging

import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors

from reportengine.figure import figure
from reportengine.table import table

from validphys import plotutils
from validphys.theorycovariance.construction import abs_chi2_data_theory_experiment
from validphys.theorycovariance.construction import abs_chi2_data_theory_dataset

from validphys.results import experiments_covmat, experiments_normcovmat, experiments_corrmat
from validphys.results import experiments_chi2_table

from validphys.theorycovariance.construction import theory_block_diag_covmat, theory_normblockcovmat
from validphys.theorycovariance.construction import theory_covmat_custom, theory_normcovmat_custom
from validphys.theorycovariance.construction import theory_blockcorrmat, theory_corrmat_custom
from validphys.theorycovariance.construction import experimentsplusblocktheory_normcovmat
from validphys.theorycovariance.construction import experimentsplustheory_normcovmat_custom
from validphys.theorycovariance.construction import experimentsplusblocktheory_corrmat
from validphys.theorycovariance.construction import experimentsplustheory_corrmat_custom

log = logging.getLogger(__name__)

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
