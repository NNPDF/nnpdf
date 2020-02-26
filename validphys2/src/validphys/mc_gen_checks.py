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
        art_data = np.zeros(real_data.shape)

        # producing replicas
        for i in range(nreplica):
            replica_exp = Experiment(real_exp)
            replica_exp.MakeReplica()
            artrep = replica_exp.get_cv()
            normartrep = artrep/real_data
            art_replicas.append(artrep) 
            normart_replicas.append(normartrep)
            
       # mean of the replicas
        for j in range(nreplica):
            art_data+=art_replicas[j]
        art_data/=nreplica
        
        return real_data, art_replicas, normart_replicas, art_data
    
@figure
def art_data_residuals(art_rep_generation, nreplica:int):

    #pass
    """
    Plot the residuals distribution of pseudodata compared to experiment.
    """
    real_data, art_replicas, normart_replicas, art_data = art_rep_generation

    residuals=real_data-art_data
    normresiduals = residuals/real_data
    fig, ax = plt.subplots()

    ax.hist(normresiduals,bins=50,histtype='step', stacked=True, fill=False)

    ax.set_ylabel(r'Data points')
    ax.set_xlabel(r'$(D^0-<D^{(r)}>)/D^0$')
    ax.set_title(r'Residuals distribution')
    
    return fig


@figure
def art_data_distribution(art_rep_generation, nreplica:int):
    """
    Plot of the distribution of pseudodata.
    """
    real_data, art_replicas, normart_replicas, art_data = art_rep_generation

    normart_data = art_data/real_data
    fig, ax = plt.subplots()

    ax.hist(normart_data, bins=50, histtype='step', stacked=True, fill=False)

    ax.set_ylabel(r'Data points')
    ax.set_xlabel(r'$<D^{(r)}>/D^0$')
    ax.set_title(r'Artificial Data Distribution')

    return fig

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
        ax.legend(handles=handles)
        ax.set_xlim(0,2)
        ax.set_ylim(0,50)
        ax.vlines(1, ax.get_ylim()[0], ax.get_ylim()[1])
        ax.vlines(normartdatapoint, ax.get_ylim()[0], ax.get_ylim()[1], linestyle="-", color="darkorchid")
        ax.set_xlabel(r"$D^{(r)}/D^0$")
        ax.set_ylabel("Frequency")

    return fig

@figure
def one_art_data_residuals(art_rep_generation, nreplica:int):

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

