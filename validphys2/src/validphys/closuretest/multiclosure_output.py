"""
multiclosure_output

Module containing the actions which produce some output in validphys
reports i.e figures or tables for multiclosure estimators in the space of
data.

"""

import numpy as np
import pandas as pd
import scipy.special as special
import scipy.stats

from reportengine.figure import figure, figuregen
from reportengine.table import table
from validphys import plotutils


@figure
def plot_data_central_diff_histogram(experiments_replica_central_diff):
    """Histogram of the difference between central prediction
    and underlying law prediction normalised by the corresponding replica
    standard deviation by concatenating the difference across all data. plot a
    scaled gaussian for reference. Total xi is the number of central differences
    which fall within the 1-sigma confidence interval of the scaled gaussian.

    """
    scaled_diffs = np.concatenate(
        [
            (central_diff / sigma).flatten()
            for sigma, central_diff in experiments_replica_central_diff
        ]
    )
    fig, ax = plotutils.subplots()
    ax.hist(scaled_diffs, bins=50, density=True, label="Central prediction distribution")
    xlim = (-5, 5)
    ax.set_xlim(xlim)

    x = np.linspace(*xlim, 100)
    ax.plot(x, scipy.stats.norm.pdf(x), "-k", label="Normal distribution")
    ax.legend()
    ax.set_xlabel("Difference to underlying prediction")
    return fig


@table
def experiments_bootstrap_sqrt_ratio_table(experiments_bootstrap_sqrt_ratio, experiments_data):
    """Given experiments_bootstrap_sqrt_ratio, which a bootstrap
    resampling of the sqrt(bias/variance) for each experiment and the total
    across all data, tabulate the mean and standard deviation across bootstrap
    samples.

    """
    indices = list(map(str, experiments_data)) + ["Total"]
    records = []
    for i, exp in enumerate(indices):
        ratio_boot = experiments_bootstrap_sqrt_ratio[i]
        records.append(
            dict(
                experiment=exp,
                mean_ratio=np.mean(ratio_boot, axis=0),
                std_ratio=np.std(ratio_boot, axis=0),
            )
        )
    df = pd.DataFrame.from_records(
        records, index="experiment", columns=("experiment", "mean_ratio", "std_ratio")
    )
    df.columns = ["Bootstrap mean sqrt(bias/variance)", "Bootstrap std. dev. sqrt(bias/variance)"]
    return df


@table
def groups_bootstrap_sqrt_ratio_table(groups_bootstrap_sqrt_ratio, groups_data):
    """Like :py:func:`experiments_bootstrap_sqrt_ratio_table` but for
    metadata groups.
    """
    df = experiments_bootstrap_sqrt_ratio_table(groups_bootstrap_sqrt_ratio, groups_data)
    idx = df.index.rename("group")
    return df.set_index(idx)


@table
def experiments_bootstrap_expected_xi_table(experiments_bootstrap_expected_xi, experiments_data):
    """Tabulate the mean and standard deviation across bootstrap samples of the
    expected xi calculated from the ratio of bias/variance. Returns a table with
    two columns, for the bootstrap mean and standard deviation
    and a row for each experiment plus the total across all experiments.

    """
    df = experiments_bootstrap_sqrt_ratio_table(experiments_bootstrap_expected_xi, experiments_data)
    # change the column headers
    df.columns = [
        r"Bootstrap mean expected $\xi_{1\sigma}$ from ratio",
        r"Bootstrap std. dev. expected $\xi_{1\sigma}$ from ratio",
    ]
    return df


@table
def groups_bootstrap_expected_xi_table(groups_bootstrap_expected_xi, groups_data):
    """Like :py:func:`experiments_bootstrap_expected_xi_table` but for metadata
    groups.
    """
    df = experiments_bootstrap_expected_xi_table(groups_bootstrap_expected_xi, groups_data)
    idx = df.index.rename("group")
    return df.set_index(idx)


@table
def experiments_bootstrap_xi_table(experiments_bootstrap_xi, experiments_data, total_bootstrap_xi):
    """Tabulate the mean and standard deviation of xi_1sigma across bootstrap
    samples. Note that the mean has already be taken across data points
    (or eigenvectors in the basis which diagonalises the covariance
    matrix) for each individual bootstrap sample.

    Tabulate the results for each experiment and for the total xi across all data.

    """
    # take mean across data points for each experiment
    xi_1sigma = [np.mean(exp_xi, axis=1) for exp_xi in experiments_bootstrap_xi]
    # take mean across all data
    xi_1sigma.append(np.mean(total_bootstrap_xi, axis=1))
    df = experiments_bootstrap_sqrt_ratio_table(xi_1sigma, experiments_data)
    df.columns = [r"Bootstrap mean $\xi_{1\sigma}$", r"Bootstrap std. dev. $\xi_{1\sigma}$"]
    return df


@table
def groups_bootstrap_xi_table(groups_bootstrap_xi, groups_data, total_bootstrap_xi):
    """Like :py:func:`experiments_bootstrap_xi_table` but for metadata groups."""
    df = experiments_bootstrap_xi_table(groups_bootstrap_xi, groups_data, total_bootstrap_xi)
    idx = df.index.rename("group")
    return df.set_index(idx)


@table
def experiments_bootstrap_xi_comparison(
    experiments_bootstrap_xi_table, experiments_bootstrap_expected_xi_table
):
    """Table comparing the mean and standard deviation across bootstrap samples of
    the measured xi_1sigma and the expected xi_1sigma calculated from
    bias/variance.

    """
    return pd.concat(
        (experiments_bootstrap_xi_table, experiments_bootstrap_expected_xi_table), axis=1
    )


@table
def groups_bootstrap_xi_comparison(groups_bootstrap_xi_table, groups_bootstrap_expected_xi_table):
    """Like :py:func:`experiments_bootstrap_xi_comparison` but for metadata
    groups.
    """
    return experiments_bootstrap_xi_comparison(
        groups_bootstrap_xi_table, groups_bootstrap_expected_xi_table
    )


@figuregen
def plot_experiments_sqrt_ratio_bootstrap_distribution(
    experiments_bootstrap_sqrt_ratio, experiments_data
):
    """Plots a histogram for each experiment and the total, showing the
    distribution of bootstrap samples. Takes the mean and std deviation of the
    bootstrap sample and plots the corresponding scaled normal distribution
    for comparison. The limits are set to be +/- 3 std deviations of the mean.

    """
    # experiments_bootstrap_sqrt_ratio includes total. str(exp) is only used to
    # generate title, so appending string is fine.
    for sqrt_ratio_sample, exp in zip(
        experiments_bootstrap_sqrt_ratio, experiments_data + ["Total"]
    ):
        fig, ax = plotutils.subplots()
        ax.hist(sqrt_ratio_sample, bins=20, density=True)
        mean = np.mean(sqrt_ratio_sample)
        std = np.std(sqrt_ratio_sample)

        xlim = (mean - 3 * std, mean + 3 * std)
        ax.set_xlim(xlim)

        x = np.linspace(*xlim, 100)
        ax.plot(
            x,
            scipy.stats.norm.pdf(x, mean, std),
            "-r",
            label=f"Corresponding normal distribution: mean = {mean:.2g}, std = {std:.2g}",
        )
        ax.legend()
        ax.set_title(f"Bootstrap distribution of sqrt(bias/variance) for {exp}")
        ax.set_xlabel("Sqrt(bias/variance)")
        yield fig


@figuregen
def plot_experiments_xi_bootstrap_distribution(
    experiments_bootstrap_xi, total_bootstrap_xi, experiments_data
):
    """Similar to :py:func:`plot_sqrt_ratio_bootstrap_distribution` except plots
    the bootstrap distribution of xi_1sigma, along with a corresponding
    scaled gaussian, for each experiment (and for all data).

    """
    # take mean across data points for each experiment
    xi_1sigma = [np.mean(exp_xi, axis=1) for exp_xi in experiments_bootstrap_xi]
    # take mean across all data
    xi_1sigma.append(np.mean(total_bootstrap_xi, axis=1))
    # use plotting function from above
    xi_plots = plot_experiments_sqrt_ratio_bootstrap_distribution(xi_1sigma, experiments_data)
    # Update the title and x label on each plot to reflect that we're plotting
    # \xi_1sigma, don't forget Total plot.
    for fig, exp in zip(xi_plots, experiments_data + ["Total"]):
        ax = fig.gca()
        ax.set_title(r"Bootstrap distribution of $\xi_{1\sigma}$ for " + str(exp))
        ax.set_xlabel(r"$\xi_{1\sigma}$")
        yield fig


@figuregen
def xq2_data_prcs_maps(xq2_data_map, each_dataset):
    """
    Heat map of the ratio bias variance (and xi, quantile estimator) for each datapoint
    in a dataset. The x and y axis are the x and Q2 coordinates of the datapoints.
    The color of each point is determined by the value of the ratio bias variance (and xi, quantile estimator).
    """
    keys = ["std_devs", "xi"]
    for j, elem in enumerate(xq2_data_map):

        for k in keys:
            if k == "std_devs":
                title = r"$R_{bv}$"
            if k == "xi":
                title = r"$\xi$"
            fig, ax = plotutils.subplots()
            im = ax.scatter(
                elem['x_coords'], elem['Q_coords'], c=(np.asarray(elem[k])), cmap='viridis'
            )
            fig.colorbar(im, label=title)
            ax.set_xscale('log')  # Set x-axis to log scale
            ax.set_yscale('log')  # Set y-axis to log scale
            ax.set_xlabel('x')
            ax.set_ylabel('Q2')
            ax.set_title(each_dataset[j].commondata.metadata.plotting.dataset_label)
            yield fig
