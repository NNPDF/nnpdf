"""
multiclosure_inconsistent_output

Module containing the actions which produce some output in validphys
reports i.e figures or tables for (inconsistent) multiclosure 
estimators in the space of data

"""

from reportengine.figure import figure
from validphys import plotutils


@figure
def plot_replica_central_diff_consistent_vs_inconsistent(
    consistent_dataset_replica_minus_central,
    inconsistent_dataset_replica_minus_central,
    dataset
):
    """ 
    For a given dataset plot the distributions of the replica theory predictions around
    the central replica prediction.
    Note:
        - Comparison is done between consistent and inconsistent fits
        - average is done over fits and datapoints of one dataset

    """
    diff_consistent = consistent_dataset_replica_minus_central
    diff_inconsistent = inconsistent_dataset_replica_minus_central

    fig, ax = plotutils.subplots()
    ax.hist(
        diff_consistent, density=True, alpha=0.5, label=f"Consistent, Central-Replica Distribution"
    )
    ax.hist(
        diff_inconsistent,
        density=True,
        alpha=0.5,
        label=f"Inconsistent, Central-Replica Distribution",
    )
    ax.set_title(f"Dataset: {dataset}")
    ax.legend()

    return fig
