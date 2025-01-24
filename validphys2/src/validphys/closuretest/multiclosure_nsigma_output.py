"""
TODO
"""

import numpy as np
from reportengine.figure import figure
import matplotlib.pyplot as plt

import sys

sys.path.insert(0, "./")
from multiclosure_nsigma import Z_ALPHA_RANGE


@figure
def plot_all_alpha_sets(set_1_alpha, set_2_alpha, set_3_alpha, n_fits):
    """
    TODO
    """
    fig, ax = plt.subplots()
    set_1, set_2, set_3 = [], [], []
    for z_alpha in set_1_alpha.keys():
        set_1.append(len(set_1_alpha[z_alpha]) / n_fits)
        set_2.append(len(set_2_alpha[z_alpha]) / n_fits)
        set_3.append(len(set_3_alpha[z_alpha]) / n_fits)

    ax.plot(set_1_alpha.keys(), set_1, label=r"$P_{\rm flag}, S_1$")
    ax.plot(set_2_alpha.keys(), set_2, label=r"$P_{\rm flag}, S_2$")
    ax.plot(set_3_alpha.keys(), set_3, label=r"$P_{\rm flag}, S_3$")
    ax.set_title(
        r"$S_1 = \{i | n_{\sigma}^i > Z_{\alpha} \}$"
        + "\n"
        + r"$S_2 = \{j \neq i | n_{\sigma}^{{\rm weighted}, j} - n_{\sigma}^{j} > Z_{\alpha} \}$"
        + "\n"
        + r"$S_3 = \{i | n_{\sigma}^{{\rm weighted}, i} > Z_{\alpha} \}$"
    )

    ax.set_xlabel(r"$Z_{\alpha}$")
    ax.set_ylabel(r"$P_{\rm flag}$")
    ax.legend()

    return fig


@figure
def plot_probability_inconsistent(
    probability_inconsistent, set_1_alpha, weighted_dataset, n_fits
):
    """
    The set of inconsistent fits I_alpha:

    1. I_alpha = 1_alpha

    2. I_alpha_1 = (1_alpha intersect 2_alpha) union (3_alpha)

    3. I_alpha_2 = I_alpha_1 union (~1_alpha intersect 2_alpha)

    The probability of a dataset being inconsistent is defined as:
            P(inconsistent) = |I_alpha| / N
    where N is the total number of fits.

    """
    rates, rates_cons = probability_inconsistent

    fig, ax = plt.subplots()

    rates_set1 = []
    for z_alpha in Z_ALPHA_RANGE:
        rates_set1.append(len(set_1_alpha[z_alpha]) / n_fits)

    ax.plot(Z_ALPHA_RANGE, rates_set1, label=r"$P_{\rm flag}, C_1$")
    ax.plot(Z_ALPHA_RANGE, rates_cons, label=r"$P_{\rm flag}, C_2$")
    ax.plot(Z_ALPHA_RANGE, rates, label=r"$P_{\rm flag}, C_3$")

    ax.set_xlabel(r"$Z_{\alpha}$")
    ax.set_ylabel(r"$P_{\rm flag}$")
    ax.set_title(
        f"Probability of flagging {weighted_dataset} dataset as inconsistent \n"
        + r"$C_1 = S_1$"
        + "\n"
        r"$C_2 = (S_1 \cap S_2) \cup S_3$"
        + "\n"
        + r"$C_3 = ((S_1 \cap S_2) \cup S_3) \cap (S_1^c \cap S_2)$"
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
    for z_alpha in Z_ALPHA_RANGE:
        rates_comp_set1.append(len(comp_set_1_alpha[z_alpha]) / n_fits)

    ax.plot(Z_ALPHA_RANGE, rates_comp_set1, label=r"$1 - P_{\rm flag}, C_1$")
    ax.plot(Z_ALPHA_RANGE, 1 - np.array(rates), label=r"$1 - P_{\rm flag}, C_2$")
    ax.plot(Z_ALPHA_RANGE, 1 - np.array(rates_cons), label=r"$1 - P_{\rm flag}, C_3$")

    ax.set_xlabel(r"$Z_{\alpha}$")
    ax.set_ylabel(r"$1 - P_{\rm flag}$")
    ax.set_title(
        f"Probability of flagging {weighted_dataset} dataset as consistent"
        + r"$((1_{\alpha} \cap 2_{\alpha}) \cup 3_{\alpha})^c$"
    )

    ax.set_title(
        f"Probability of flagging {weighted_dataset} dataset as consistent \n"
        + r"$C_1 = S_1$"
        + "\n"
        r"$C_2 = (S_1 \cap S_2) \cup S_3$"
        + "\n"
        + r"$C_3 = ((S_1 \cap S_2) \cup S_3) \cap (S_1^c \cap S_2)$"
    )
    ax.legend()
    return fig


@figure
def plot_set1_alpha(set_1_alpha, n_fits):
    """
    TODO
    """
    fig, ax = plt.subplots()
    set1_vals = []
    for z_alpha in set_1_alpha.keys():
        set1_vals.append(len(set_1_alpha[z_alpha]) / n_fits)

    ax.plot(set_1_alpha.keys(), set1_vals, label="Set 1 alpha")
    ax.set_title(r"$1_{\alpha}$")
    ax.set_xlabel(r"$Z_{\alpha}$")
    ax.legend()
    return fig


@figure
def plot_set3_alpha(set_3_alpha, n_fits):
    """
    TODO
    """
    fig, ax = plt.subplots()
    set3_vals = []
    for z_alpha in set_3_alpha.keys():
        set3_vals.append(len(set_3_alpha[z_alpha]) / n_fits)

    ax.plot(set_3_alpha.keys(), set3_vals, label="Set 3 alpha")
    ax.set_title(r"$3_{\alpha}$")
    ax.set_xlabel(r"$Z_{\alpha}$")
    ax.legend()
    return fig


@figure
def plot_set2_alpha(set_2_alpha, n_fits):
    """
    TODO
    """
    fig, ax = plt.subplots()
    set2_vals = []
    for z_alpha in set_2_alpha.keys():
        set2_vals.append(len(set_2_alpha[z_alpha]) / n_fits)

    ax.plot(set_2_alpha.keys(), set2_vals, label="Set 2 alpha")
    ax.set_title(r"$2_{\alpha}$")
    ax.set_xlabel(r"$Z_{\alpha}$")
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
    for z_alpha in set_1_alpha.keys():
        set1_vals.append(len(set_1_alpha[z_alpha]) / n_fits)
        set3_vals.append(len(set_3_alpha[z_alpha]) / n_fits)

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
    for z_alpha in set_2_alpha.keys():
        set2_vals.append(len(set_2_alpha[z_alpha]) / n_fits)

    ax.plot(set_2_alpha.keys(), set2_vals, label="Set 2 alpha")

    ax.legend()
    return fig
