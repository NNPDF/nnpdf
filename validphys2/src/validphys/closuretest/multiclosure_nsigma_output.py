"""
TODO
"""
from reportengine.figure import figure
import matplotlib.pyplot as plt


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