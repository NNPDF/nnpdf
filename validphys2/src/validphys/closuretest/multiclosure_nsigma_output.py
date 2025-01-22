"""
TODO
"""
from reportengine.figure import figure
import matplotlib.pyplot as plt


@figure
def plot_probability_inconsistent(set_2_alpha, set_1_alpha, set_3_alpha, comp_set_1_alpha, n_fits):
    """
    P(inconsistent) = P(set 2 | set 1) + P(set 3) + P(set 2 | ~ set 1)
    """

    fig, ax = plt.subplots()

    tagged_rates = []
    for alpha in set_2_alpha.keys():
        set_2_inters_1 = set_2_alpha[alpha].intersection(set_1_alpha[alpha])
        set_3 = set_3_alpha[alpha]
        set_2_inters_comp_1 = set_2_alpha[alpha].intersection(comp_set_1_alpha[alpha])

        set_tagged_fits = set_2_inters_1.union(set_3).union(set_2_inters_comp_1)
        tagged_rates.append(len(set_tagged_fits) / n_fits)

    ax.plot(set_2_alpha.keys(), tagged_rates, label="P(inconsistent)")
    ax.legend()
    return fig

@figure
def plot_probability_consistent(set_2_alpha, set_1_alpha, set_3_alpha, comp_set_1_alpha, comp_set_2_alpha, comp_set_3_alpha, n_fits):
    """
    P(consistent) = 1 - P(inconsistent) ?= P(~set 2 | ~set 1) + P(~set 2 and ~set 3| set 1)
    """
    fig, ax = plt.subplots()

    untagged_rates = []
    untagged_rates_method2 = []
    for alpha in set_2_alpha.keys():
        set_2_inters_1 = set_2_alpha[alpha].intersection(set_1_alpha[alpha])
        set_3 = set_3_alpha[alpha]
        set_2_inters_comp_1 = set_2_alpha[alpha].intersection(comp_set_1_alpha[alpha])

        comp_set_2_inters_comp_set1 = comp_set_2_alpha[alpha].intersection(comp_set_1_alpha[alpha])
        comp_set_2_inters_comp_set3 = comp_set_2_alpha[alpha].intersection(comp_set_3_alpha[alpha])
        comp_set_2_inters_comp_set3_inter_set_1 = comp_set_2_inters_comp_set3.intersection(set_1_alpha[alpha])

        set_untagged_fits = comp_set_2_inters_comp_set1.union(comp_set_2_inters_comp_set3_inter_set_1)


        set_tagged_fits = set_2_inters_1.union(set_3).union(set_2_inters_comp_1)
        untagged_rates.append(1 - len(set_tagged_fits) / n_fits)

        untagged_rates_method2.append(len(set_untagged_fits) / n_fits)


    ax.plot(set_2_alpha.keys(), untagged_rates, label="1 - P(inconsistent)")
    ax.plot(set_2_alpha.keys(), untagged_rates_method2, label="P(~set 2 | ~set 1) + P(~set 2 and ~set 3| set 1)")
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