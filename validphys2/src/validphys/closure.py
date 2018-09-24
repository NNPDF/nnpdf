# -*- coding: utf-8 -*-
"""
Functions and Plots relating to Closure Test 
Statistical Estimators.
"""
from __future__ import generator_stop

import logging

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors, ticker as mticker

from reportengine.figure import figure
from reportengine.checks import make_argcheck, CheckError
from reportengine import collect

from validphys.results import total_experiments_chi2data, experiment_results
from validphys import plotutils
from validphys.calcutils import calc_chi2

log = logging.getLogger(__name__)


@figure
def plot_biases(experiments, bias_experiments, closures_speclabel):
    """Plot the biases of all experiments with bars."""
    biases = np.array(bias_experiments).T
    print(biases)
    labels = closures_speclabel
    xticks = [experiment.name for experiment in experiments]
    fig, ax = plotutils.barplot(biases, collabels=xticks, datalabels=labels)
    ax.set_title("biases for experiments")
    return fig

def bias_experiment(experiment_result_closure,
                    experiment_result_fakepdf):
    bias_out = np.zeros(len(experiment_result_closure))
    for i, (ct, fake) in enumerate(zip(experiment_result_closure,
                                       experiment_result_fakepdf)):
        ((dt_ct, th_ct), (dt_fake, th_fake)) = ct, fake
        central_diff = th_ct.central_value - th_fake.central_value
        bias_out[i] = calc_chi2(dt_ct.sqrtcovmat, central_diff)
    return bias_out

closures_speclabel = collect('speclabel', ('closures',), element_default=None)

experiment_result_closure = collect(experiment_results, ('closures',))
experiment_result_fakepdf = collect(experiment_results, ('closures', 'fakepdf',))

bias_experiments = collect(bias_experiment, ('experiments',))
