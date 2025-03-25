"""
Module for plotting the results of the multiclosure_nsigma.py script.

Can be used to reproduce the plots in Sec. 4 of arXiv: 2503.17447
"""

import numpy as np
from reportengine.figure import figure
import matplotlib.pyplot as plt

from validphys.closuretest.multiclosure_nsigma import Z_ALPHA_RANGE


@figure
def plot_all_sets(set_1, set_3, set_2, n_fits):
    """
    Plots S_1, S_2 and S_3.
    """
    fig, ax = plt.subplots()
    s1, s2, s3 = [], [], []
    for z_alpha in set_1.keys():
        s1.append(len(set_1[z_alpha]) / n_fits)
        s2.append(len(set_2[z_alpha]) / n_fits)
        s3.append(len(set_3[z_alpha]) / n_fits)

    ax.plot(set_1.keys(), s1, linewidth=3, label=r"$P_{\rm flag}, S_1$")
    ax.plot(set_2.keys(), s2, linestyle='--', linewidth=3, label=r"$P_{\rm flag}, S_2$")
    ax.plot(set_3.keys(), s3, linestyle=':', linewidth=3, label=r"$P_{\rm flag}, S_3$")
    ax.set_title(r"HERA I + II $\sigma_{e^+p}$, $E_p$=575 GeV")
    ax.set_xlabel(r"$Z$", fontsize='large')
    ax.set_ylabel(r"$P_{\rm flag}$", fontsize='large')
    ax.legend(fontsize='large')

    return fig


@figure
def plot_1_minus_all_sets(set_1, set_3, set_2, n_fits):
    """
    Plots complement of S_1, S_2 and S_3.
    """
    fig, ax = plt.subplots()
    s1, s2, s3 = [], [], []
    for z_alpha in set_1.keys():
        s1.append(len(set_1[z_alpha]) / n_fits)
        s2.append(len(set_2[z_alpha]) / n_fits)
        s3.append(len(set_3[z_alpha]) / n_fits)

    ax.plot(set_1.keys(), 1 - np.array(s1), linewidth=3, label=r"$P_{\rm flag}, S_1$")
    ax.plot(
        set_2.keys(), 1 - np.array(s2), linewidth=3, linestyle="--", label=r"$P_{\rm flag}, S_2$"
    )
    ax.plot(
        set_3.keys(), 1 - np.array(s3), linewidth=3, linestyle=":", label=r"$P_{\rm flag}, S_3$"
    )

    ax.set_title(r"HERA I + II $\sigma_{e^+p}$, $E_p$=575 GeV")
    ax.set_xlabel(r"$Z$", fontsize='large')
    ax.set_ylabel(r"$1 - P_{\rm flag}$", fontsize='large')
    ax.legend(fontsize='large')

    return fig


@figure
def plot_probability_inconsistent(probability_inconsistent, set_1, weighted_dataset, n_fits):
    """
    The set of inconsistent fits:

    1. C_1 = S_1

    2. C_2 = (S_1 intersect S_3) union (S_2)

    3. C_3 = S_1 union (~S_1 intersect S_3)

    The probability of a dataset being inconsistent is defined as:
            P(inconsistent) = |I_alpha| / N
    where N is the total number of fits.

    """
    c_3_rates, c_2_rates = probability_inconsistent

    fig, ax = plt.subplots()

    rates_set1 = []
    for z_alpha in Z_ALPHA_RANGE:
        rates_set1.append(len(set_1[z_alpha]) / n_fits)

    ax.plot(Z_ALPHA_RANGE, rates_set1, linewidth=3, label=r"$P_{\rm flag}, C_1$")
    ax.plot(Z_ALPHA_RANGE, c_2_rates, linewidth=3, linestyle="--", label=r"$P_{\rm flag}, C_2$")
    ax.plot(Z_ALPHA_RANGE, c_3_rates, linewidth=3, linestyle=":", label=r"$P_{\rm flag}, C_3$")

    ax.set_xlabel(r"$Z$", fontsize='large')
    ax.set_ylabel(r"$P_{\rm flag}$", fontsize='large')
    ax.set_title(r"HERA I + II $\sigma_{e^+p}$, $E_p$=575 GeV")
    ax.legend(fontsize='large')
    return fig


@figure
def plot_probability_consistent(probability_inconsistent, comp_set_1, weighted_dataset, n_fits):
    """
    Plots the probability of dataset being flagged as consistent.
    """
    c_3_rates, c_2_rates = probability_inconsistent
    fig, ax = plt.subplots()

    rates_comp_set1 = []
    for z_alpha in Z_ALPHA_RANGE:
        rates_comp_set1.append(len(comp_set_1[z_alpha]) / n_fits)

    ax.plot(Z_ALPHA_RANGE, rates_comp_set1, linewidth=3, label=r"$1 - P_{\rm flag}, C_1$")
    ax.plot(
        Z_ALPHA_RANGE,
        1 - np.array(c_2_rates),
        linewidth=3,
        linestyle="--",
        label=r"$1 - P_{\rm flag}, C_2$",
    )
    ax.plot(
        Z_ALPHA_RANGE,
        1 - np.array(c_3_rates),
        linewidth=3,
        linestyle=":",
        label=r"$1 - P_{\rm flag}, C_3$",
    )

    ax.set_xlabel(r"$Z$", fontsize='large')
    ax.set_ylabel(r"$1 - P_{\rm flag}$", fontsize='large')
    ax.set_title(r"HERA I + II $\sigma_{e^+p}$, $E_p$=575 GeV")
    ax.legend(fontsize='large')
    return fig
