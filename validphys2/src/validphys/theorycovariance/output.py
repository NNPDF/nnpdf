"""
output.py
Basic tools for plotting theory covariance matrices and their properties.
"""

import logging
from math import inf

from matplotlib import cm
from matplotlib import colors as mcolors
import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.figure import figure
from validphys import plotutils
from validphys.theorycovariance.construction import theory_covmat_custom

log = logging.getLogger(__name__)


def matrix_plot_labels(df):
    """Returns the tick locations and labels, and the starting
    point values for each category,  based on a dataframe
    to be plotted. The dataframe is assumed to be multiindexed by
    (process, dataset, points) or else (dataset, points). The tick
    location is in the centre of the dataset, and labelling is by
    the outermost index of the multiindex."""
    if len(df.index[0]) == 3:
        proclabels = [x[0] for x in df.index]
        labels = proclabels
    elif len(df.index[0]) == 2:
        dslabels = [x[0] for x in df.index]
        labels = dslabels
    unique_ds = []
    unique_ds.append([labels[0], 0])
    for x in range(len(labels) - 1):
        if labels[x + 1] != labels[x]:
            unique_ds.append([labels[x + 1], x + 1])
    ticklabels = [unique_ds[x][0] for x in range(len(unique_ds))]
    startlocs = [unique_ds[x][1] for x in range(len(unique_ds))]
    startlocs += [len(labels)]
    ticklocs = [0 for x in range(len(startlocs) - 1)]
    for i in range(len(startlocs) - 1):
        ticklocs[i] = 0.5 * (startlocs[i + 1] + startlocs[i])
    return ticklocs, ticklabels, startlocs


def plot_covmat_heatmap(covmat, title):
    """Matrix plot of a covariance matrix."""
    df = covmat
    df.sort_index(axis=0, inplace=True)
    df.sort_index(axis=1, inplace=True)
    oldindex = df.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    # reindex index
    newdf = df.reindex(newindex)
    # reindex columns by transposing, reindexing, then transposing back
    newdf = (newdf.T.reindex(newindex)).T
    matrix = newdf.values
    fig, ax = plotutils.subplots(figsize=(15, 15))
    matrixplot = ax.matshow(
        100 * matrix,
        cmap=cm.Spectral_r,
        norm=mcolors.SymLogNorm(
            linthresh=0.00001, linscale=1, vmin=-100 * matrix.max(), vmax=100 * matrix.max()
        ),
    )
    cbar = fig.colorbar(matrixplot, fraction=0.046, pad=0.04)
    cbar.set_label(label="% of data value", fontsize=20)
    cbar.ax.tick_params(labelsize=20)
    ax.set_title(title, fontsize=25)
    ticklocs, ticklabels, startlocs = matrix_plot_labels(newdf)
    ax.set_xticks(ticklocs)
    ax.set_xticklabels(ticklabels, rotation=30, ha="right", fontsize=20)
    ax.xaxis.tick_bottom()
    ax.set_yticks(ticklocs)
    ax.set_yticklabels(ticklabels, fontsize=20)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ax.vlines(startlocs_lines, -0.5, len(matrix) - 0.5, linestyles="dashed")
    ax.hlines(startlocs_lines, -0.5, len(matrix) - 0.5, linestyles="dashed")
    ax.margins(x=0, y=0)
    return fig


_procorder = ("DIS NC", "DIS CC", "TOP", "DY NC", "DY CC", "SINGLETOP", "JETS", "PHOTON", "DIJET")

_dsorder = (
    "BCDMSP",
    "BCDMSD",
    "SLACP",
    "SLACD",
    "NMC",
    "NMCPD",
    "HERAF2CHARM",
    "HERACOMBNCEP460",
    "HERACOMBNCEP575",
    "HERACOMBNCEP820",
    "HERACOMBNCEP920",
    "HERACOMBNCEM",
    "CHORUSNU",
    "CHORUSNB",
    "NTVNUDMN",
    "NTVNBDMN",
    "HERACOMBCCEP",
    "HERACOMBCCEM",
    "CDFZRAP",
    "D0ZRAP",
    "D0WEASY",
    "D0WMASY",
    "ATLASWZRAP36PB",
    "ATLASZHIGHMASS49FB",
    "ATLASLOMASSDY11EXT",
    "ATLASWZRAP11",
    "ATLASZPT8TEVMDIST",
    "ATLASZPT8TEVYDIST",
    "CMSWEASY840PB",
    "CMSWMASY47FB",
    "CMSWCHARMRAT",
    "CMSDY2D11",
    "CMSWMU8TEV",
    "CMSWCHARMTOT",
    "CMSZDIFF12",
    "LHCBZ940PB",
    "LHCBWZMU7TEV",
    "LHCBWZMU8TEV",
    "LHCBZEE2FB",
    "ATLAS1JET11",
    "CMSJETS11",
    "CDFR2KT",
    "ATLASTTBARTOT",
    "ATLASTOPDIFF8TEVTRAPNORM",
    "CMSTTBARTOT",
    "CMSTOPDIFF8TEVTTRAPNORM",
)


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


def plot_corrmat_heatmap(corrmat, title):
    """Matrix plot of a correlation matrix"""
    df = corrmat
    df.sort_index(axis=0, inplace=True)
    df.sort_index(axis=1, inplace=True)
    oldindex = df.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    # reindex index
    newdf = df.reindex(newindex)
    # reindex columns by transposing, reindexing, then transposing back
    newdf = (newdf.T.reindex(newindex)).T
    matrix = newdf.values
    fig, ax = plotutils.subplots(figsize=(15, 15))
    matrixplot = ax.matshow(matrix, cmap=cm.Spectral_r, vmin=-1, vmax=1)
    cbar = fig.colorbar(matrixplot, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=20)
    ax.set_title(title, fontsize=25)
    ticklocs, ticklabels, startlocs = matrix_plot_labels(newdf)
    ax.set_xticks(ticklocs)
    ax.set_xticklabels(ticklabels, rotation=30, ha="right", fontsize=20)
    ax.xaxis.tick_bottom()
    ax.set_yticks(ticklocs)
    ax.set_yticklabels(ticklabels, fontsize=20)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ax.vlines(startlocs_lines, -0.5, len(matrix) - 0.5, linestyles="dashed")
    ax.hlines(startlocs_lines, -0.5, len(matrix) - 0.5, linestyles="dashed")
    ax.margins(x=0, y=0)
    return fig


@figure
def plot_normexpcovmat_heatmap(procs_normcovmat):
    """Matrix plot of the experiment covariance matrix normalised to data."""
    fig = plot_covmat_heatmap(procs_normcovmat, "Experimental Covariance Matrix")
    return fig


@figure
def plot_expcorrmat_heatmap(procs_corrmat):
    """Matrix plot of the experiment correlation matrix"""
    fig = plot_corrmat_heatmap(procs_corrmat, "Experimental Correlation Matrix")
    return fig


@figure
def plot_normthcovmat_heatmap_custom(theory_normcovmat_custom, theoryids, point_prescription):
    """Matrix plot for block diagonal theory covariance matrix by process type"""
    l = len(theoryids)
    if l == 5:
        if point_prescription == "5bar point":
            l = r"$\bar{5}$"
    fig = plot_covmat_heatmap(theory_normcovmat_custom, f"Theory Covariance matrix ({l} pt)")
    return fig


@figure
def plot_thcorrmat_heatmap_custom(theory_corrmat_custom, theoryids, point_prescription):
    """Matrix plot of the theory correlation matrix, correlations by process type"""
    l = len(theoryids)
    if l == 5:
        if point_prescription == "5bar point":
            l = r"$\bar{5}$"
    fig = plot_corrmat_heatmap(theory_corrmat_custom, f"Theory Correlation matrix ({l} pt)")
    return fig


@figure
def plot_expplusthcorrmat_heatmap_custom(
    experimentplustheory_corrmat_custom, theoryids, point_prescription
):
    """Matrix plot of the exp + theory correlation matrix"""
    l = len(theoryids)
    if l == 5:
        if point_prescription == "5bar point":
            l = r"$\bar{5}$"
    fig = plot_corrmat_heatmap(
        experimentplustheory_corrmat_custom, f"Experimental + Theory Correlation Matrix ({l} pt)"
    )
    return fig


@figure
def plot_diag_cov_comparison(
    theory_covmat_custom, procs_covmat, procs_data_values, theoryids, point_prescription
):
    """Plot of sqrt(cov_ii)/|data_i| for cov = exp, theory, exp+theory"""
    l = len(theoryids)
    if l == 5:
        if point_prescription == "5bar point":
            l = r"$\bar{5}$"
    data = np.abs(procs_data_values)
    sqrtdiags_th = np.sqrt(np.diag(theory_covmat_custom)) / data

    # sort by theory covmat processes
    newindex = sorted(sqrtdiags_th.index, key=_get_key)
    sqrtdiags_th = sqrtdiags_th.reindex(newindex)

    sqrtdiags_exp = np.sqrt(np.diag(procs_covmat)) / data
    sqrtdiags_exp = sqrtdiags_exp.reindex(newindex)

    df_total = theory_covmat_custom + procs_covmat
    sqrtdiags_tot = np.sqrt(np.diag(df_total)) / data
    sqrtdiags_tot = sqrtdiags_tot.reindex(newindex)

    fig, ax = plotutils.subplots(figsize=(20, 10))
    ax.plot(sqrtdiags_exp.values, ".", label="Experiment", color="orange")
    ax.plot(sqrtdiags_th.values, ".", label="Theory", color="red")
    ax.plot(sqrtdiags_tot.values, ".", label="Total", color="blue")
    ticklocs, ticklabels, startlocs = matrix_plot_labels(sqrtdiags_th)
    ax.set_xticks(ticklocs)
    ax.set_xticklabels(ticklabels, rotation=45, fontsize=20)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ax.vlines(startlocs_lines, 0, len(data), linestyles="dashed")
    ax.set_ylabel(r"$\frac{\sqrt{S_{ii}}}{|D_i|}$", fontsize=30)
    ax.yaxis.set_tick_params(labelsize=20)
    ax.set_ylim([0, 0.5])
    ax.set_title(
        f"Square root of diagonal elements of covariance matrices ({l} pt), "
        + "normalised to absolute value of data",
        fontsize=20,
    )
    ax.legend(fontsize=20)
    ax.margins(x=0)
    return fig


theory_covmat_custom_dataspecs = collect(theory_covmat_custom, ("dataspecs",))


@figure
def plot_diag_cov_comparison_by_process(
    theory_covmat_custom_dataspecs, procs_covmat, procs_data_values, dataspecs
):
    """Plot of sqrt(cov_ii)/|data_i| for cov = exp, theory, exp+theory, by process

    dataspecs here is used also to specify whether to plot or not the total uncertainty,
    defined as the experimental uncertainty plus the uncertainty coming from the theory
    covariance matrix for the relevant dataspec entry.
    """

    fig, ax = plotutils.subplots(figsize=(20, 10))

    # sort by theory covmat processes
    newindex = sorted(procs_data_values.index, key=_get_key)

    data = np.abs(procs_data_values)
    sqrtdiags_exp = np.sqrt(np.diag(procs_covmat)) / data
    sqrtdiags_exp = sqrtdiags_exp.reindex(newindex)
    ax.plot(sqrtdiags_exp.values, "*", markersize=4, label="Experimental uncertanties")

    # loop on th covmat
    for specs, theory_covmat_custom in zip(dataspecs, theory_covmat_custom_dataspecs):
        label = specs["speclabel"]

        sqrtdiags_th = np.sqrt(np.diag(theory_covmat_custom)) / data
        sqrtdiags_th = sqrtdiags_th.reindex(newindex)
        ax.plot(sqrtdiags_th.values, "o", markersize=4, label=label)

        if specs["plot_total"]:
            df_total = theory_covmat_custom + procs_covmat
            sqrtdiags_tot = np.sqrt(np.diag(df_total)) / data
            sqrtdiags_tot = sqrtdiags_tot.reindex(newindex)
            ax.plot(
                df_total.values, "v", markersize=4, label=f"Experimental uncertanties + {label}"
            )

    ticklocs, ticklabels, startlocs = matrix_plot_labels(sqrtdiags_th)
    ax.set_xticks(ticklocs)
    ax.set_xticklabels(ticklabels, rotation=45, fontsize=20)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ax.vlines(startlocs_lines, 0, len(data), linestyles="dashed")
    ax.set_ylabel(r"$\frac{\sqrt{S_{ii}}}{|D_i|}$", fontsize=30)
    ax.yaxis.set_tick_params(labelsize=20)
    ax.set_ylim([0, 0.5])
    ax.legend(fontsize=20)
    ax.margins(x=0)
    return fig


@figure
def plot_diag_cov_comparison_by_experiment(
    theory_covmat_custom_dataspecs,
    experiments_covmat_no_table,
    procs_data_values_experiment,
    dataspecs,
):
    """Plot of sqrt(cov_ii)/|data_i| for cov = exp, theory, exp+theory, by experiment"""

    fig, ax = plotutils.subplots(figsize=(20, 10))

    procs_data_values_experiment.sort_index(level=0, inplace=True)
    data = np.abs(procs_data_values_experiment)
    plot_index = procs_data_values_experiment.index

    # plot exp values, take diagonal first
    sqrtdiags_exp = pd.DataFrame(
        np.sqrt(np.diag(experiments_covmat_no_table)), index=experiments_covmat_no_table.index
    )
    sqrtdiags_exp.sort_index(level=0, inplace=True)
    sqrtdiags_exp = sqrtdiags_exp[0] / data.values
    ax.plot(sqrtdiags_exp.values, "*", markersize=4, label="Experimental uncertanties")

    # loop on th covmat
    for specs, theory_covmat_custom in zip(dataspecs, theory_covmat_custom_dataspecs):
        label = specs["speclabel"]

        sqrtdiags_th = pd.DataFrame(
            np.sqrt(np.diag(theory_covmat_custom)), index=theory_covmat_custom.index.droplevel(0)
        )
        # sort the diag th covmat by experiments
        temp_index_nogroup = plot_index.droplevel(0)
        sqrtdiags_th = sqrtdiags_th.reindex(temp_index_nogroup)
        sqrtdiags_th = sqrtdiags_th[0].values / data
        ax.plot(sqrtdiags_th.values, "o", markersize=4, label=label)

        if specs["plot_total"]:
            df_total = np.sqrt(sqrtdiags_th**2 + sqrtdiags_exp**2)
            ax.plot(
                df_total.values, "v", markersize=4, label=f"Experimental uncertanties + {label}"
            )

    ticklocs, ticklabels, startlocs = matrix_plot_labels(sqrtdiags_th)
    ax.set_xticks(ticklocs)
    ax.set_xticklabels(ticklabels, rotation=45, fontsize=20)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ax.vlines(startlocs_lines, 0, len(data), linestyles="dashed")
    ax.set_ylabel(r"$\frac{\sqrt{S_{ii}}}{|D_i|}$", fontsize=30)
    ax.yaxis.set_tick_params(labelsize=20)
    ax.set_ylim([0, 0.5])
    ax.legend(fontsize=20)
    ax.margins(x=0)
    return fig
