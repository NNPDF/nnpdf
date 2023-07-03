"""
multiclosure_inconsistent_output

Module containing the actions which produce some output in validphys
reports i.e figures or tables for (inconsistent) multiclosure 
estimators in the space of data

"""

from reportengine.figure import figure
from validphys import plotutils


@figure
def plot_variance_distribution_multi(multi_dataset_fits_bias_replicas_variance_samples, dataspecs):
    """
    histogram of the distribution of the variances (k) defined as the
    distance between central replica and replica (k) in units of the
    experimental covariance matrix.
    If more than one group of dataspecs (e.g. consistent and inconsistent)
    fits are defined, than plot comparison of these.
    """
    fig, ax = plotutils.subplots()
    for (_, variance_dist, _), spec in zip(
        multi_dataset_fits_bias_replicas_variance_samples, dataspecs
    ):
        label = spec['speclabel']

        ax.hist(variance_dist, bins='auto', density=True, alpha=0.5, label=label)

    ax.legend()
    return fig


@figure
def plot_bias_distribution_multi(multi_dataset_fits_bias_replicas_variance_samples, dataspecs):
    """
    histogram of the distribution of the biases (l) defined as the
    distance between central replica and underlying law in units of the
    experimental covariance matrix.
    If more than one group of dataspecs (e.g. consistent and inconsistent)
    fits are defined, than plot comparison of these.
    """
    fig, ax = plotutils.subplots()
    for (bias_dist, _, _), spec in zip(
        multi_dataset_fits_bias_replicas_variance_samples, dataspecs
    ):
        label = spec['speclabel']

        ax.hist(bias_dist, bins='auto', density=True, alpha=0.5, label=label)

    ax.legend()
    return fig


@figure
def plot_multi_bias_distribution_bootstrap(multi_fits_bootstrap_dataset_bias_variance, dataspecs):
    """
    histogram of the distribution of the biases obtained by bootstrapping with
    `fits_bootstrap_dataset_bias_variance` over the existing fits.

    If more than one group of dataspecs (e.g. consistent and inconsistent)
    fits are defined, than plot comparison of these.
    """
    fig, ax = plotutils.subplots()
    for (bias, _), spec in zip(multi_fits_bootstrap_dataset_bias_variance, dataspecs):
        label = spec['speclabel']

        ax.hist(bias, bins='auto', density=True, alpha=0.5, label=label)

    ax.legend()
    return fig


@figure
def plot_multi_ratio_bias_variance_distribution_bootstrap(
    multi_fits_bootstrap_dataset_bias_variance, dataspecs
):
    """
    histogram of the ratio bias variance distribution obtained by bootstrapping bias
    and variance with `fits_bootstrap_dataset_bias_variance`.

    If more than one group of dataspecs (e.g. consistent and inconsistent)
    fits are defined, than plot comparison of these.
    """
    fig, ax = plotutils.subplots()
    for (bias, variance), spec in zip(multi_fits_bootstrap_dataset_bias_variance, dataspecs):
        ratio = bias / variance
        label = spec['speclabel']

        ax.hist(ratio, bins='auto', density=True, alpha=0.5, label=label)

    ax.legend()
    return fig
