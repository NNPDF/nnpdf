# -*- coding: utf-8 -*-
"""
mc_gen_checks.py

Tools to check the pseudo-data MC generation.
"""
from __future__ import generator_stop

import logging
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import moment as mom

from NNPDF import Experiment, RandomGenerator
from reportengine.table import table
from reportengine.figure import figure
log = logging.getLogger(__name__)


def art_rep_generation(experiments, nreplica:int, experiments_index):
    """Generates the nreplica pseudodata replicas for a given experiment"""

    RandomGenerator.InitRNG(0,0)

    for exp in experiments:
        #Since we are going to modify the experiments, we copy them
        #(and work on the copies) to avoid all
        #sorts of weirdness with other providers. We don't want this to interact
        #with ExperimentSpec at all, because it could do funny things with the
        #cache when calling load(). We need to copy this yet again, for each
        # of the noisy replicas.
        real_exp = Experiment(exp.load())

        art_replicas = []
        normart_replicas = []
        real_data = real_exp.get_cv()

        # producing replicas
        for i in range(nreplica):
            replica_exp = Experiment(real_exp)
            replica_exp.MakeReplica()
            artrep = replica_exp.get_cv()
            normartrep = artrep/real_data
            art_replicas.append(artrep)
            normart_replicas.append(normartrep)

        art_data = np.mean(art_replicas, axis=0)

        return real_data, art_replicas, normart_replicas, art_data

def per_point_art_rep_generation(experiments, nreplica:int, experiments_index):
    """Generates the nreplica pseudodata replicas for a given experiment"""

    RandomGenerator.InitRNG(0,0)

    for exp in experiments:
        #Since we are going to modify the experiments, we copy them
        #(and work on the copies) to avoid all
        #sorts of weirdness with other providers. We don't want this to interact
        #with ExperimentSpec at all, because it could do funny things with the
        #cache when calling load(). We need to copy this yet again, for each
        # of the noisy replicas.
        real_exp = Experiment(exp.load())

        art_replicas = []
        normart_replicas = []
        real_data = real_exp.get_cv()

        # producing replicas
        for i in range(nreplica):
            replica_exp = Experiment(real_exp)
            for point in range(len(replica_exp.get_cv())):
                replica_exp.MakePerPointReplica(point)
            artrep = replica_exp.get_cv()
            normartrep = artrep/real_data
            art_replicas.append(artrep)
            normart_replicas.append(normartrep)

        art_data = np.mean(art_replicas, axis=0)

        return real_data, art_replicas, normart_replicas, art_data

@figure
def art_data_residuals(art_rep_generation, nreplica:int, color="green"):

    #pass
    """
    Plot the residuals distribution of pseudodata compared to experiment.
    """
    real_data, art_replicas, normart_replicas, art_data = art_rep_generation

    residuals=real_data-art_data
    normresiduals = residuals/real_data
    fig, ax = plt.subplots()

    ax.hist(normresiduals,bins=50,histtype='step', stacked=True, fill=False, color=color)

    ax.set_ylabel(r'Data points')
    ax.set_xlabel(r'$(D^0-<D^{(r)}>)/D^0$')
    ax.set_title(r'Residuals distribution')

    return fig

@figure
def per_point_art_data_residuals(per_point_art_rep_generation, nreplica:int):
    return art_data_residuals(per_point_art_rep_generation, nreplica, color="orange")


@figure
def art_data_distribution(art_rep_generation, nreplica:int, title='Artificial Data Distribution', color="green"):
    """
    Plot of the distribution of pseudodata.
    """
    real_data, art_replicas, normart_replicas, art_data = art_rep_generation

    normart_data = art_data/real_data
    fig, ax = plt.subplots()

    ax.hist(normart_data, bins=50, histtype='step', stacked=True, fill=False, color=color)

    ax.set_ylabel(r'Data points')
    ax.set_xlabel(r'$<D^{(r)}>/D^0$')
    ax.set_title(title)

    return fig

@figure
def per_point_art_data_distribution(per_point_art_rep_generation, nreplica:int):
    return art_data_distribution(per_point_art_rep_generation, nreplica, title='Uncorrelated Artificial Data Distribution', color="orange")


@figure
def art_data_moments(art_rep_generation, nreplica:int, color="green"):
    """
    Returns the moments of the distributions per data point, as a histogram.
    """
    real_data, art_replicas, normart_replicas, art_data = art_rep_generation

    artrep_array = np.asarray(normart_replicas)
    normart_data = art_data/real_data


    fig, axes = plt.subplots(nrows=3, figsize=(10,12))
    # Plot histogram of moments
    for momno, ax in zip(range(1,4), axes.flatten()):
        # Calculate moments
        moms = []
        for i, datapoint, normartdatapoint in zip(range(len(artrep_array.T)), artrep_array.T, normart_data):
            moment = mom(datapoint, moment=momno)
            moms.append(moment)
        ax.hist(moms, bins=50, histtype='step', stacked=True, fill=False, color=color)
        ax.set_ylabel("Data points")
        ax.set_xlabel(f"Moment {momno}")

    return fig

@figure
def per_point_art_data_moments(per_point_art_rep_generation, nreplica:int):
    return art_data_moments(per_point_art_rep_generation, nreplica, color="orange")

@figure
def art_data_comparison(art_rep_generation, nreplica:int):

    #pass
    """
    Plots per datapoint of the distribution of replica values.
    """
    real_data, art_replicas, normart_replicas, art_data = art_rep_generation

    artrep_array = np.asarray(normart_replicas)
    normart_data = art_data/real_data

    fig, axes = plt.subplots(nrows=len(artrep_array.T), figsize=(4,2*len(artrep_array.T)))

    for i, ax, datapoint, normartdatapoint in zip(range(len(artrep_array.T)), axes.flatten(), artrep_array.T, normart_data):
        ax.hist(datapoint, bins=10, histtype="step", stacked=True, fill=False)
        extraString = f"Datapoint number = {i}"
        handles, labels = ax.get_legend_handles_labels()
        handles.append(mpatches.Patch(color="none", label=extraString))
        ax.set_xlim(-0.5,2.5)
        ax.set_ylim(0,0.5*nreplica)
        ax.vlines(1, ax.get_ylim()[0], ax.get_ylim()[1])
        ax.vlines(normartdatapoint, ax.get_ylim()[0], ax.get_ylim()[1], linestyle="-", color="darkorchid")
        ax.vlines(0, ax.get_ylim()[0], ax.get_ylim()[1], linestyle="-", color="dodgerblue")
        ax.vlines(2, ax.get_ylim()[0], ax.get_ylim()[1], linestyle="-", color="dodgerblue")
        ax.legend(handles=handles)
        ax.set_xlabel(r"$D^{(r)}/D^0$")
        ax.set_ylabel("Frequency")

    return fig


@figure
def one_art_data_residuals(art_rep_generation, nreplica:int, experiments):

    #pass
    """
    Residuals plot for the first datapoint.
    """
    RandomGenerator.InitRNG(0,0)
    for exp in experiments:

        real_exp = Experiment(exp.load())
        real_data = real_exp.get_cv()
        one_art_data = np.zeros(nreplica)
        one_data_index=0

        #producing replicas
        for i in range(nreplica):
            replica_exp = Experiment(real_exp)
            replica_exp.MakeReplica()
            one_art_data[i]=replica_exp.get_cv()[one_data_index]

    fig, ax = plt.subplots()

    residual = one_art_data-real_data[one_data_index]
    normresidual = residual/real_data[one_data_index]
    ax.hist(normresidual,bins=50,histtype='step', stacked=True, fill=False)

    ax.set_ylabel(r'replicas')
    ax.set_xlabel(r'$(D^{(r)}_{0} - D^0_{0})/D^0_{0}$')
    ax.set_title(r'Residual for Data Point 0')

    return fig

@figure
def plot_deviation_from_mean(art_rep_generation, per_point_art_rep_generation, nreplica:int, experiments):
    real_data, art_replicas, normart_replicas, art_data = art_rep_generation
    ppreal_data, ppart_replicas, ppnormart_replicas, ppart_data = per_point_art_rep_generation

    fig, ax = plt.subplots()

    residuals = (real_data - art_data)/np.std(art_replicas, axis=0)
    ppresiduals = (real_data - ppart_data)/np.std(ppart_replicas, axis=0)

    ax.plot(ppresiduals, color="orange", label="Uncorrelated")
    ax.plot(residuals, color="green", label="Correlated")
    ax.legend()
    ax.set_xlabel("Datapoint index")
    ax.set_ylabel(r"$(<D> - D_0)/\sigma$")
    ax.set_title("Deviation from the mean")
    ax.hlines(0, ax.get_xlim()[0], ax.get_xlim()[1], linestyle="-", color="black")

    return fig

@table
def art_data_mean_table(art_rep_generation, nreplica:int, experiments):
    """Generate table for artdata mean values
    """
    real_data, art_replicas, normart_replicas, art_data = art_rep_generation

    #residuals=real_data-art_data
    data=[]
    for experiment in experiments:
        for dataset in experiment.datasets:
            ds = dataset.load()
            Ndata = ds.GetNData()
            for i in range(Ndata):
                line=[dataset.name,art_data[i],real_data[i],abs(art_data[i]-real_data[i])]
                data.append(line)

    df =  pd.DataFrame(data,columns=["DataSet","ArtData","ExpData","abs(residual)"])

    return df

