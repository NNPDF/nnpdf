# -*- coding: utf-8 -*-
"""
mc_gen.py

Tools to check the pseudo-data MC generation.
"""
# The functions in this module have been ported to not use libNNPDF
# but they should not be used as an example as they follow the libNNPDF logic
import logging
import matplotlib.patches as mpatches
from matplotlib.figure import Figure
import numpy as np
import pandas as pd
from scipy.stats import moment as mom

from reportengine.table import table
from reportengine.figure import figure

log = logging.getLogger(__name__)


def art_rep_generation(groups_data, make_replicas):
    """Generates the nreplica pseudodata replicas"""
    real_data_list = []

    for group in groups_data:
        # Load all the commondata
        real_group = group.load_commondata()
        real_data = np.concatenate([i.get_cv() for i in real_group])
        real_data_list.append(real_data)

    real_data = np.concatenate(real_data_list)

    art_replicas = np.stack(make_replicas)
    art_data = np.mean(art_replicas, axis=0)
    normart_replicas = art_replicas / real_data

    return real_data, art_replicas, normart_replicas, art_data


@figure
def art_data_residuals(art_rep_generation, color="green"):
    """
    Plot the residuals distribution of pseudodata compared to experiment.
    """
    real_data, _, _, art_data = art_rep_generation

    residuals = real_data - art_data
    normresiduals = residuals / real_data
    fig = Figure()
    ax = fig.subplots()

    ax.hist(
        normresiduals, bins=50, histtype="step", stacked=True, fill=False, color=color
    )

    ax.set_ylabel(r"Data points")
    ax.set_xlabel(r"$(D^0-<D^{(r)}>)/D^0$")
    ax.set_title(r"Residuals distribution")

    return fig


@figure
def art_data_distribution(
    art_rep_generation, title="Artificial Data Distribution", color="green"
):
    """
    Plot of the distribution of pseudodata.
    """
    real_data, _, _, art_data = art_rep_generation

    normart_data = art_data / real_data
    fig = Figure()
    ax = fig.subplots()

    ax.hist(
        normart_data, bins=50, histtype="step", stacked=True, fill=False, color=color
    )

    ax.set_ylabel(r"Data points")
    ax.set_xlabel(r"$<D^{(r)}>/D^0$")
    ax.set_title(title)

    return fig


@figure
def art_data_moments(art_rep_generation, color="green"):
    """
    Returns the moments of the distributions per data point, as a histogram.
    """
    _, _, normart_replicas, _ = art_rep_generation

    artrep_array = np.asarray(normart_replicas)
    
    fig = Figure(figsize=(10, 12))
    axes = [fig.add_subplot(3, 1, i+1) for i in range(3)]
    # Plot histogram of moments
    for momno, ax in zip(range(1, 4), axes.flatten()):
        # Calculate moments
        moms = []
        for _, datapoint in zip(range(len(artrep_array.T)), artrep_array.T):
            moment = mom(datapoint, moment=momno)
            moms.append(moment)
        ax.hist(moms, bins=50, histtype="step", stacked=True, fill=False, color=color)
        ax.set_ylabel("Data points")
        ax.set_xlabel(f"Moment {momno}")

    return fig


@figure
def art_data_comparison(art_rep_generation, nreplica: int):
    """
    Plots per datapoint of the distribution of replica values.
    """
    real_data, _, normart_replicas, art_data = art_rep_generation

    artrep_array = np.asarray(normart_replicas)
    normart_data = art_data / real_data

    nrows=len(artrep_array.T)
    fig = Figure(figsize=(4, 2 * len(artrep_array.T)))
    axes = [fig.add_subplot(nrows, 1, i+1) for i in range(nrows)]

    for i, ax, datapoint, normartdatapoint in zip(
        range(len(artrep_array.T)), axes.flatten(), artrep_array.T, normart_data
    ):
        ax.hist(datapoint, bins=10, histtype="step", stacked=True, fill=False)
        extraString = f"Datapoint number = {i}"
        handles, _ = ax.get_legend_handles_labels()
        handles.append(mpatches.Patch(color="none", label=extraString))
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylim(0, 0.5 * nreplica)
        ax.vlines(1, ax.get_ylim()[0], ax.get_ylim()[1])
        ax.vlines(
            normartdatapoint,
            ax.get_ylim()[0],
            ax.get_ylim()[1],
            linestyle="-",
            color="darkorchid",
        )
        ax.vlines(
            0, ax.get_ylim()[0], ax.get_ylim()[1], linestyle="-", color="dodgerblue"
        )
        ax.vlines(
            2, ax.get_ylim()[0], ax.get_ylim()[1], linestyle="-", color="dodgerblue"
        )
        ax.legend(handles=handles)
        ax.set_xlabel(r"$D^{(r)}/D^0$")
        ax.set_ylabel("Frequency")

    return fig


@figure
def one_art_data_residuals(groups_data, indexed_make_replicas):
    """
    Residuals plot for the first datapoint.
    """
    one_data_index = 0
    all_replicas = pd.concat(indexed_make_replicas, axis=1)
    group_level = all_replicas.index.get_level_values("group")

    group_level = indexed_make_replicas[0].index.get_level_values("group")

    all_normresidual = []
    for group in groups_data:
        real_group = group.load_commondata()
        real_data = np.concatenate([i.get_cv() for i in real_group])
        one_art_data = all_replicas[group_level == group.name].iloc[one_data_index]

        residual = one_art_data - real_data[one_data_index]
        all_normresidual.append(residual / real_data[one_data_index])

    fig = Figure()
    ax = fig.subplots()

    ax.hist(all_normresidual, bins=50, histtype="step", stacked=True, fill=False)

    ax.set_ylabel(r"replicas")
    ax.set_xlabel(r"$(D^{(r)}_{0} - D^0_{0})/D^0_{0}$")
    ax.set_title(r"Residual for Data Point 0")

    return fig


@table
def art_data_mean_table(art_rep_generation, groups_data):
    """Generate table for artdata mean values"""
    real_data, _, _, art_data = art_rep_generation

    data = []
    for group in groups_data:
        for dataset in group.datasets:
            Ndata = dataset.load_commondata().ndata
            for i in range(Ndata):
                line = [
                    dataset.name,
                    art_data[i],
                    real_data[i],
                    abs(art_data[i] - real_data[i]),
                ]
                data.append(line)

    df = pd.DataFrame(data, columns=["DataSet", "ArtData", "ExpData", "abs(residual)"])

    return df
