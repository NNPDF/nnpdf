# -*- coding: utf-8 -*-
"""
MCgen_checks.py

Tools to check the pseudo-data MC generation.
"""
from __future__ import generator_stop

from collections import OrderedDict, namedtuple, Sequence
import itertools
import logging

import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la
import pandas as pd

from NNPDF import CommonData, Experiment
from reportengine.checks import require_one, remove_outer, check_not_empty
from reportengine.table import table
from reportengine.figure import figure
from reportengine import collect

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
    fig, ax = plt.subplots()

    ax.hist(residuals,bins=50,histtype='step', stacked=True, fill=False)

    ax.set_ylabel(r'data points')
    ax.set_xlabel(r'$\mu_{art} - \mu_{exp}$')
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


    fig, ax = plt.subplots()

    ax.hist(art_data, bins=50, histtype='step', stacked=True, fill=False)

    ax.set_ylabel(r'data points')
    ax.set_xlabel(r'$\mu_{art} - \mu_{exp}$')
    ax.set_title(r'Residuals distribution')

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
    ax.hist(residual,bins=50,histtype='step', stacked=True, fill=False)

    ax.set_ylabel(r'replicas')
    ax.set_xlabel(r'$dist(art) - \mu_{exp}$')
    ax.set_title(r'one point residual')
    
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

