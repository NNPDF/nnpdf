"""
TODO
"""

import numpy as np
from reportengine.figure import figure
import matplotlib.pyplot as plt

import sys

sys.path.insert(0, "./")
from multiclosure_nsigma import ALPHA_RANGE


@figure
def plot_probability_inconsistent(
    probability_inconsistent, set_1_alpha, weighted_dataset, n_fits
):
    """
    The set of inconsistent fits I_alpha can be defined in different ways, two possible cases are:

    1. I_alpha_1 = (1_alpha intersect 2_alpha) union (3_alpha)

    2. I_alpha_2 = I_alpha_1 union (~1_alpha intersect 2_alpha)

    The probability of a dataset being inconsistent is defined as:
            P(inconsistent) = |I_alpha| / N
    where N is the total number of fits.

    """
    rates, rates_cons = probability_inconsistent

    fig, ax = plt.subplots()

    rates_set1 = []
    for alpha in ALPHA_RANGE:
        rates_set1.append(len(set_1_alpha[alpha]) / n_fits)

    ax.plot(ALPHA_RANGE, rates_cons, label="P(inconsistent) (conservative)")
    ax.plot(ALPHA_RANGE, rates, label="P(inconsistent)")
    ax.plot(ALPHA_RANGE, rates_set1, label="Reference fit (set 1)")
    ax.set_xlabel("alpha")
    ax.set_title(
        f"Probability of classifying {weighted_dataset} dataset as inconsistent"
    )
    ax.legend()
    return fig


@figure
def plot_probability_consistent(
    probability_inconsistent, comp_set_1_alpha, weighted_dataset, n_fits
):
    """
    TODO
    """
    rates, rates_cons = probability_inconsistent
    fig, ax = plt.subplots()

    rates_comp_set1 = []
    for alpha in ALPHA_RANGE:
        rates_comp_set1.append(len(comp_set_1_alpha[alpha]) / n_fits)

    ax.plot(ALPHA_RANGE, 1 - np.array(rates_cons), label="P(consistent) (conservative)")
    ax.plot(ALPHA_RANGE, 1 - np.array(rates), label="P(consistent)")
    ax.plot(ALPHA_RANGE, rates_comp_set1, label="Reference fit (set 1)")
    ax.set_xlabel("alpha")
    ax.set_title(f"Probability of classifying {weighted_dataset} dataset as consistent")
    ax.legend()
    return fig


@figure
def plot_set1_alpha(set_1_alpha, n_fits):
    """
    TODO
    """
    fig, ax = plt.subplots()
    set1_vals = []
    for alpha in set_1_alpha.keys():
        set1_vals.append(len(set_1_alpha[alpha]) / n_fits)

    ax.plot(set_1_alpha.keys(), set1_vals, label="Set 1 alpha")
    ax.set_title(r"$1_{\alpha}$")
    ax.set_xlabel(r"$\alpha$")
    ax.legend()
    return fig


@figure
def plot_set3_alpha(set_3_alpha, n_fits):
    """
    TODO
    """
    fig, ax = plt.subplots()
    set3_vals = []
    for alpha in set_3_alpha.keys():
        set3_vals.append(len(set_3_alpha[alpha]) / n_fits)

    ax.plot(set_3_alpha.keys(), set3_vals, label="Set 3 alpha")
    ax.set_title(r"$3_{\alpha}$")
    ax.set_xlabel(r"$\alpha$")
    ax.legend()
    return fig


@figure
def plot_set2_alpha(set_2_alpha, n_fits):
    """
    TODO
    """
    fig, ax = plt.subplots()
    set2_vals = []
    for alpha in set_2_alpha.keys():
        set2_vals.append(len(set_2_alpha[alpha]) / n_fits)

    ax.plot(set_2_alpha.keys(), set2_vals, label="Set 2 alpha")
    ax.set_title(r"$2_{\alpha}$")
    ax.set_xlabel("alpha")
    ax.legend()
    return fig


@figure
def plot_set1_vs_set3_alpha(set_1_alpha, set_3_alpha, weighted_dataset, n_fits):
    """
    TODO
    """
    fig, ax = plt.subplots()
    set1_vals = []
    set3_vals = []
    for alpha in set_1_alpha.keys():
        set1_vals.append(len(set_1_alpha[alpha]) / n_fits)
        set3_vals.append(len(set_3_alpha[alpha]) / n_fits)

    ax.plot(set_1_alpha.keys(), set1_vals, label="TPR, reference fit (set 1)")
    ax.plot(set_3_alpha.keys(), set3_vals, label="TPR, weighted fit (set 3)")
    ax.set_title(f"TPR for {weighted_dataset}")

    ax.legend()
    return fig


@figure
def plot_set2(set_2_alpha, n_fits):
    """
    TODO
    """
    fig, ax = plt.subplots()
    set2_vals = []
    for alpha in set_2_alpha.keys():
        set2_vals.append(len(set_2_alpha[alpha]) / n_fits)

    ax.plot(set_2_alpha.keys(), set2_vals, label="Set 2 alpha")

    ax.legend()
    return fig
