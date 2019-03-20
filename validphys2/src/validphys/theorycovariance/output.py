# -*- coding: utf-8 -*-
"""
output.py
Basic tools for plotting theory covariance matrices and their properties.
"""
from __future__ import generator_stop

import logging

from math import inf
import pandas as pd
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors

from reportengine.figure import figure
from reportengine.table import table

from validphys import plotutils
from validphys.results import experiments_chi2_table

from validphys.theorycovariance.construction import dataset_index_byprocess

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
    """Returns the tick locations and labels, and the starting
    point values for each category,  based on a dataframe
    to be plotted. The dataframe is assumed to be multiindexed by
    (process, dataset, points) or else (dataset, points). The tick
    location is in the centre of the dataset, and labelling is by
    the outermost index of the multiindex."""
    if len(df.index[0]) == 3:
        proclabels = [x[0] for x in df.index]
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
    return ticklocs, ticklabels, startlocs

def plot_covmat_heatmap(covmat, title, dataset_index_byprocess):
    """Matrix plot of a covariance matrix.
    WARNING: Plotting is ordered by process and dataset in the custom
    order given in _procorder and _dsorder. In order for this to work
    the covmat provided MUST be indexed and ordered the same as
    experiments_covmat, that is by (exp, dataset, point) in the order
    of experiments and datasets listed in the runcard. Otherwise the
    labels will not correspond to the correct points."""
    df = pd.DataFrame(covmat.values, index=dataset_index_byprocess,
                      columns=dataset_index_byprocess)
    df.sort_index(0, inplace=True)
    df.sort_index(1, inplace=True)
    oldindex = df.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    # reindex index
    newdf = df.reindex(newindex)
    # reindex columns by transposing, reindexing, then transposing back
    newdf = (newdf.T.reindex(newindex)).T
    matrix = newdf.values
    fig,ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(100*matrix,
                            cmap=cm.Spectral_r,
                            norm=mcolors.SymLogNorm(linthresh=0.01,
                            linscale=10,
                            vmin=-100*matrix.max(),
                            vmax=100*matrix.max()))
    cbar=fig.colorbar(matrixplot)
    cbar.set_label(label="% of data", fontsize=20)
    cbar.ax.tick_params(labelsize=20)
    ax.set_title(title, fontsize=25)
    ticklocs, ticklabels, startlocs = matrix_plot_labels(newdf)
    plt.xticks(ticklocs, ticklabels, rotation=30, ha="right", fontsize=20)
    plt.gca().xaxis.tick_bottom()
    plt.yticks(ticklocs, ticklabels, fontsize=20)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x-0.5 for x in startlocs]
    ax.vlines(startlocs_lines, -0.5, len(matrix)-0.5, linestyles='dashed')
    ax.hlines(startlocs_lines, -0.5, len(matrix)-0.5, linestyles='dashed')
    ax.margins(x=0, y=0)
    return fig

_procorder = ('DIS NC', 'DIS CC', 'DY', 'JETS', 'TOP')

_dsorder = (   'BCDMSP', 'BCDMSD', 'SLACP', 'SLACD', 'NMC', 'NMCPD',
               'HERAF2CHARM', 'HERACOMBNCEP460', 'HERACOMBNCEP575',
               'HERACOMBNCEP820', 'HERACOMBNCEP920', 'HERACOMBNCEM',
               'CHORUSNU', 'CHORUSNB', 'NTVNUDMN', 'NTVNBDMN',
               'HERACOMBCCEP', 'HERACOMBCCEM', 'CDFZRAP', 'D0ZRAP',
               'D0WEASY', 'D0WMASY', 'ATLASWZRAP36PB', 'ATLASZHIGHMASS49FB',
               'ATLASLOMASSDY11EXT', 'ATLASWZRAP11', 'ATLASZPT8TEVMDIST',
               'ATLASZPT8TEVYDIST', 'CMSWEASY840PB', 'CMSWMASY47FB',
               'CMSWCHARMRAT', 'CMSDY2D11', 'CMSWMU8TEV', 'CMSWCHARMTOT',
	       'CMSZDIFF12', 'LHCBZ940PB', 'LHCBWZMU7TEV', 'LHCBWZMU8TEV',
               'LHCBZEE2FB', 'ATLAS1JET11', 'CMSJETS11', 'CDFR2KT',
               'ATLASTTBARTOT', 'ATLASTOPDIFF8TEVTRAPNORM', 'CMSTTBARTOT',
               'CMSTOPDIFF8TEVTTRAPNORM')

def _get_key(element):
    """The key used to sort covariance matrix dataframes according to
    the ordering of processes and datasets specified in _procorder and
    _dsorder."""
    x1, y1, z1 = element
    x2 = _procorder.index(x1) if x1 in _procorder else inf
    y2 = _dsorder.index(y1) if y1 in _dsorder else inf
    z2 = z1
    newelement = (x2, y2, z2)
    return newelement

@figure
def plot_corrmat_heatmap(corrmat, title, dataset_index_byprocess):
    """Matrix plot of a correlation matrix"""
    df = pd.DataFrame(corrmat.values, index=dataset_index_byprocess,
		 	columns=dataset_index_byprocess)
    df.sort_index(0, inplace=True)
    df.sort_index(1, inplace=True)
    oldindex = df.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    # reindex index
    newdf = df.reindex(newindex)
    # reindex columns by transposing, reindexing, then transposing back
    newdf = (newdf.T.reindex(newindex)).T
    matrix = newdf.values
    fig, ax = plt.subplots(figsize=(15,15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    cbar=fig.colorbar(matrixplot)
    cbar.ax.tick_params(labelsize=20)
    ax.set_title(title, fontsize=25)
    ticklocs, ticklabels, startlocs = matrix_plot_labels(newdf)
    plt.xticks(ticklocs, ticklabels, rotation=30, ha="right", fontsize=20)
    plt.gca().xaxis.tick_bottom()
    plt.yticks(ticklocs, ticklabels, fontsize=20)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x-0.5 for x in startlocs]
    ax.vlines(startlocs_lines, -0.5, len(matrix)-0.5, linestyles='dashed')
    ax.hlines(startlocs_lines, -0.5, len(matrix)-0.5, linestyles='dashed')
    ax.margins(x=0, y=0)
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
    fig = plot_covmat_heatmap(df,
                              "(Theory + experiment)/mean(experiment)" +
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
    oldindex = sqrtdiags_th.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    sqrtdiags_th = sqrtdiags_th.reindex(newindex)
    sqrtdiags_exp = np.sqrt(np.diag(experiments_covmat))/data
    sqrtdiags_exp = pd.DataFrame(sqrtdiags_exp.values, index=dataset_index_byprocess)
    sqrtdiags_exp.sort_index(0,inplace=True)
    sqrtdiags_exp = sqrtdiags_exp.reindex(newindex)
    df_total = theory_covmat_custom + experiments_covmat
    sqrtdiags_tot = np.sqrt(np.diag(df_total))/data
    sqrtdiags_tot = pd.DataFrame(sqrtdiags_tot.values, index=dataset_index_byprocess)
    sqrtdiags_tot.sort_index(0,inplace=True)
    sqrtdiags_tot = sqrtdiags_tot.reindex(newindex)
    fig,ax = plt.subplots(figsize=(20,10))
    ax.plot(sqrtdiags_exp.values, '.', label="Experiment", color="orange")
    ax.plot(sqrtdiags_th.values, '.', label="Theory", color = "red")
    ax.plot(sqrtdiags_tot.values, '.', label="Total", color = "blue")
    ticklocs, ticklabels, startlocs = matrix_plot_labels(sqrtdiags_th)
    plt.xticks(ticklocs, ticklabels, rotation=45, fontsize=20)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x-0.5 for x in startlocs]
    ax.vlines(startlocs_lines, 0, len(data), linestyles='dashed')
    ax.set_ylabel(r"$\frac{\sqrt{cov_{ii}}}{|D_i|}$", fontsize=30)
    ax.yaxis.set_tick_params(labelsize=20)
    ax.set_ylim([0,0.5])
    ax.set_title(f"Square root of diagonal elements of covariances matrices for {l} points, "
                 + "normalised to absolute value of data",
                 fontsize=20)
    ax.legend(fontsize=20)
    ax.margins(x=0)
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
    oldindex = df_inv_exp.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    df_inv_exp = df_inv_exp.reindex(newindex)
    df_inv_tot = pd.DataFrame(inv_tot, index=dataset_index_byprocess)
    df_inv_tot.sort_index(0,inplace=True)
    df_inv_tot = df_inv_tot.reindex(newindex)
    fig,ax = plt.subplots()
    ax.plot(df_inv_exp.values, '.', label="Experiment", color="orange")
    ax.plot(df_inv_tot.values, '.', label="Experiment + Theory", color="mediumseagreen")
    ticklocs, ticklabels, startlocs = matrix_plot_labels(df_inv_exp)
    plt.xticks(ticklocs, ticklabels, rotation="vertical", fontsize=20)
    ax.vlines(startlocs, 0, len(matrix_theory), linestyles='dashed')
    ax.set_ylabel(r"$\frac{1}{D_i}\frac{1}{\sqrt{[cov^{-1}_]{ii}}}$", fontsize=30)
    ax.yaxis.set_tick_params(labelsize=20)
    ax.set_title(f"Diagonal impact of adding theory covariance matrix for {l} points",
                 fontsize=20)
    ax.legend(fontsize=20)
    ax.margins(x=0)
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
