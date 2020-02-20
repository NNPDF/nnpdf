# -*- coding: utf-8 -*-
"""
MCgen_checks.py

Tools to check the pseudo-data MC generation.
"""
from __future__ import generator_stop

import itertools
import logging


from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NNPDF import Experiment
from reportengine.table import table
from reportengine.figure import figure

log = logging.getLogger(__name__)

@figure
def ArtDataResiduals(experiments, nreplica:int, experiments_index):

    #pass
    """Generate artificial data for a given experiment.
    (inspired by: closure_pseudodata_replicas in results.py)
    nreplica: Number of pseudodata replicas.

    """
    from  NNPDF import RandomGenerator
    RandomGenerator.InitRNG(0,0)
    for exp in experiments:
        #Since we are going to modify the experiments, we copy them
        #(and work on the copies) to avoid all
        #sorts of weirdness with other providers. We don't want this to interact
        #with ExperimentSpec at all, because it could do funny things with the
        #cache when calling load(). We need to copy this yet again, for each
        # of the noisy replicas.
        real_exp = Experiment(exp.load())

        exp_location = experiments_index.get_loc(real_exp.GetExpName())

        index = itertools.count()

        art_replicas = []
        real_data = real_exp.get_cv()
        art_data = np.zeros(real_data.shape)

        #producing replicas
        for i in range(nreplica):
            replica_exp = Experiment(real_exp)
            replica_exp.MakeReplica()
            art_replicas.append(replica_exp.get_cv())
        
        #mean of the replicas
        for j in range(nreplica):
            art_data+=art_replicas[j]
        art_data/=nreplica

    residuals=real_data-art_data
    normresiduals = residuals/real_data
    fig, ax = plt.subplots()

    ax.hist(normresiduals,bins=50,histtype='step', stacked=True, fill=False)

    ax.set_ylabel(r'data points')
    ax.set_xlabel(r'$(D^0-<D^{(r)}>)/D^0$')
    ax.set_title(r'Residuals distribution')
    
    return fig


@figure
def ArtDataDistribution(experiments, nreplica: int, experiments_index):

    #pass
    """Generate artificial data for a given experiment.
    (inspired by: closure_pseudodata_replicas in results.py)
    nreplica: Number of pseudodata replicas.

    """
    from NNPDF import RandomGenerator
    RandomGenerator.InitRNG(0, 0)
    for exp in experiments:
        #Since we are going to modify the experiments, we copy them
        #(and work on the copies) to avoid all
        #sorts of weirdness with other providers. We don't want this to interact
        #with ExperimentSpec at all, because it could do funny things with the
        #cache when calling load(). We need to copy this yet again, for each
        # of the noisy replicas.
        real_exp = Experiment(exp.load())

        exp_location = experiments_index.get_loc(real_exp.GetExpName())

        index = itertools.count()

        art_replicas = []
        real_data = real_exp.get_cv()
        art_data = np.zeros(real_data.shape)

        #producing replicas
        for i in range(nreplica):
            replica_exp = Experiment(real_exp)
            replica_exp.MakeReplica()
            art_replicas.append(replica_exp.get_cv())

        #mean of the replicas
        for j in range(nreplica):
            art_data += art_replicas[j]
        art_data /= nreplica

    
    normart_data = art_data/real_data
    fig, ax = plt.subplots()

    ax.hist(normart_data, bins=50, histtype='step', stacked=True, fill=False)

    ax.set_ylabel(r'data points')
    ax.set_xlabel(r'$<D^{(r)}>/D^0$')
    ax.set_title(r'Artificial Data Distribution')

    return fig

@figure
def ArtDataComparison(experiments, nreplica: int, experiments_index):

    #pass
    """Generate artificial data for a given experiment.
    (inspired by: closure_pseudodata_replicas in results.py)
    nreplica: Number of pseudodata replicas.

    """
    from NNPDF import RandomGenerator
    RandomGenerator.InitRNG(0, 0)
    for exp in experiments:
        #Since we are going to modify the experiments, we copy them
        #(and work on the copies) to avoid all
        #sorts of weirdness with other providers. We don't want this to interact
        #with ExperimentSpec at all, because it could do funny things with the
        #cache when calling load(). We need to copy this yet again, for each
        # of the noisy replicas.
        real_exp = Experiment(exp.load())

        exp_location = experiments_index.get_loc(real_exp.GetExpName())

        index = itertools.count()

        art_replicas = []
        real_data = real_exp.get_cv()
        normart_data = np.zeros(real_data.shape)

        #producing replicas
        for i in range(nreplica):
            replica_exp = Experiment(real_exp)
            replica_exp.MakeReplica()
            artrep = replica_exp.get_cv()
            normartrep = artrep/real_data
            art_replicas.append(normartrep)  
                #mean of the replicas
        for j in range(nreplica):
            normart_data += art_replicas[j]
        normart_data /= nreplica
    
    artrep_array = np.asarray(art_replicas)
    from IPython import embed
  #  embed()
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
        
   # binnedreps = np.apply_along_axis(lambda a: np.histogram(a, bins=4)[0], 0, artrep_array)
   # fig = plt.figure()
   # ax = Axes3D(fig)
    
  #  embed()
   # x = np.outer(range(len(binnedreps)), np.ones(len(binnedreps.T)))
   # y  = np.outer(np.ones(len(binnedreps)), range(len(binnedreps.T)))
   # ax.plot_wireframe(x, y, binnedreps)   # ax.set_xlabel(r'Data point index')
   # ax.set_title(r'Experimental replicas')

    return fig

@figure
def OneArtDataResiduals(experiments, nreplica:int, experiments_index):

    #pass
    """Generate artificial data for a given experiment.
    (inspired by: closure_pseudodata_replicas in results.py)
    nreplica: Number of pseudodata replicas.

    """
    from  NNPDF import RandomGenerator
    RandomGenerator.InitRNG(0,0)
    for exp in experiments:
        #Since we are going to modify the experiments, we copy them
        #(and work on the copies) to avoid all
        #sorts of weirdness with other providers. We don't want this to interact
        #with ExperimentSpec at all, because it could do funny things with the
        #cache when calling load(). We need to copy this yet again, for each
        # of the noisy replicas.
        real_exp = Experiment(exp.load())

        exp_location = experiments_index.get_loc(real_exp.GetExpName())

        index = itertools.count()

        art_replicas = []
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
def ArtDataMeanTable(experiments, nreplica: int, experiments_index):
    """Generate table for artdata mean values
    """
    from NNPDF import RandomGenerator
    RandomGenerator.InitRNG(0, 0)
    art_data = np.zeros(1)
    real_data = np.zeros(1)

    for exp in experiments:
        #Since we are going to modify the experiments, we copy them
        #(and work on the copies) to avoid all
        #sorts of weirdness with other providers. We don't want this to interact
        #with ExperimentSpec at all, because it could do funny things with the
        #cache when calling load(). We need to copy this yet again, for each
        # of the noisy replicas.
        real_exp = Experiment(exp.load())

        art_replicas = []
        real_data = real_exp.get_cv()
        art_data = np.zeros(real_data.shape)

        #producing replicas
        for i in range(nreplica):
            replica_exp = Experiment(real_exp)
            replica_exp.MakeReplica()
            art_replicas.append(replica_exp.get_cv())
            

        #mean of the replicas
        for j in range(nreplica):
            art_data += art_replicas[j]
        art_data /= nreplica

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

