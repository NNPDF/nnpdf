# -*- coding: utf-8 -*-
"""
Functions and Plots relating to Closure Test 
Statistical Estimators.
"""

import logging

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors, ticker as mticker

from reportengine.figure import figure
from reportengine import collect

from validphys.results import experiment_results
from validphys import plotutils
from validphys.calcutils import calc_chi2

log = logging.getLogger(__name__)


@figure
def plot_biases(experiments, bias_experiments, closures_speclabel):
    """Plot the biases of all experiments with bars."""
    biases = np.array(bias_experiments).T
    labels = closures_speclabel
    xticks = [experiment.name for experiment in experiments]
    fig, ax = plotutils.barplot(biases, collabels=xticks, datalabels=labels)
    ax.set_title("biases for experiments")
    ax.legend()
    return fig

def bias_experiment(exp_result_closure,
                    exp_result_t0):
    """Calculates the bias for all closure fit specified in runcard for
    one experiment. The bias is the chi2 between the level zero closure
    replica and the level zero of the PDF used to generate the data.
    The underlying law is taken to be the same as the PDF used to generate
    the t0 covariance matrix
    """
    bias_out = np.zeros(len(exp_result_closure))
    for i, (ct, ul) in enumerate(zip(exp_result_closure,
                                       exp_result_t0)):
        ((dt_ct, th_ct), (_, th_ul)) = ct, ul
        central_diff = th_ct.central_value - th_ul.central_value
        bias_out[i] = calc_chi2(dt_ct.sqrtcovmat, central_diff)/len(dt_ct)
    return bias_out

#Closure test collect functions

closures_speclabel = collect('speclabel', ('closures',), element_default=None)

exp_result_closure = collect(experiment_results, ('closures',))
exp_result_t0 = collect(experiment_results, ('closures', 'fitunderlyinglaw',))

bias_experiments = collect(bias_experiment, ('experiments',))
