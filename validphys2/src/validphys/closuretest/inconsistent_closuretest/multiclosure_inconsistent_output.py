"""
multiclosure_inconsistent_output

Module containing the actions which produce some output in validphys
reports i.e figures or tables for (inconsistent) multiclosure 
estimators in the space of data

"""

from reportengine.figure import figure
from validphys import plotutils


@figure
def plot_replica_central_diff_multi(dataset_replica_minus_central_multi, dataspecs):
    """
    histogram of the distribution of the differences between central replica
    and replicas.
    Datapoints are averaged over.
    """
    fig, ax = plotutils.subplots()

    for ds_replica_minus_central, spec in zip(dataset_replica_minus_central_multi, dataspecs):
        diff, _, _ = ds_replica_minus_central

        label = spec['speclabel']

        ax.hist(
            diff,
            bins='auto',
            density=True,
            alpha=0.5,
            label=f"{label}",
        )

    ax.legend()

    return fig


@figure
def plot_variance_distribution_multi(dataset_variance_per_replica_multi, dataspecs):
    """
    histogram of the distribution of the variances (k) defined as the
    distance between central replica and replica (k) in units of the
    experimental covariance matrix.
    If more than one group of dataspecs (e.g. consistent and inconsistent)
    fits are defined, than plot comparison of these.
    """
    fig, ax = plotutils.subplots()
    for variance_dist, spec in zip(dataset_variance_per_replica_multi, dataspecs):
        label = spec['speclabel']

        ax.hist(variance_dist, bins='auto', density=True, alpha=0.5, label=label)

    ax.legend()
    return fig
