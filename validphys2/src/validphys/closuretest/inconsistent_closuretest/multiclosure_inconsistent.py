"""
closuretest/inconsistent_closuretest/multiclosure_inconsistent.py

Module containing all of the statistical estimators which are
averaged across multiple inconsistent fits. The actions
in this module are used to produce results which are plotted in
``multiclosure_inconsistent_output.py``

"""

import numpy as np
from validphys.calcutils import calc_chi2

from reportengine import collect


""" To load several multiclosure fits. Useful for inconsistent closure test analysis """
multi_dataset_loader = collect("internal_multiclosure_dataset_loader", ("dataspecs",))


def dataset_replica_minus_central(internal_multiclosure_dataset_loader):
    """For a given dataset calculate the difference between theory prediction
    of each replica and central replica.

    Returns
    -------
    tuple
        3d tuple containing:
        - diff of shape (Nrep,): average is done over fits and datapoints
        - diff of shape (Ndat, Nrep): average is done over fits
        - sqrt covariance matrix of shape (Ndat, Ndat, )
    """
    closures_th, _, _, sqrt_cov = internal_multiclosure_dataset_loader
    replicas = np.asarray([th.error_members for th in closures_th])
    centrals = np.mean(replicas, axis=-1)

    diff = centrals[:, :, np.newaxis] - replicas
    # average over fits
    diff_ndat_nrep = np.mean(diff, axis=0)
    # average over fits and data points
    diff_nrep = np.mean(diff, axis=(0, 1))
    return diff_nrep, diff_ndat_nrep, sqrt_cov


def dataset_replica_minus_central_multi(multi_dataset_loader):
    """
    like `dataset_replica_minus_central` but for several groups of fits

    Returns
    -------
    list
        list each element of which is a tuple from `dataset_replica_minus_central`
    """
    multi_diffs = []
    for internal_mct_ds_loader in multi_dataset_loader:
        diff_nrep, diff_ndat_nrep, sqrt_cov = dataset_replica_minus_central(internal_mct_ds_loader)
        multi_diffs.append((diff_nrep, diff_ndat_nrep, sqrt_cov))
    return multi_diffs


def dataset_variance_per_replica(dataset_replica_minus_central):
    """
    Compute the variance for each replica for one or multiple fits.
    The average over fits is done in `dataset_replica_minus_central`.

    Returns
    -------
    np.ndarray
            array of shape (Nreplica,) for the distribution of the
            variance over replicas
    """
    _, diff, sqrt_cov = dataset_replica_minus_central

    variances = calc_chi2(sqrt_cov, diff)
    return variances


def dataset_variance_per_replica_multi(dataset_replica_minus_central_multi):
    """
    like `dataset_variance_per_replica` but for different groups of fits.

    Returns
    -------
    list
        list each element of which is an array obtained with `dataset_variance_per_replica`
    """
    multi_variances = []
    for ds_rep_minus_central in dataset_replica_minus_central_multi:
        variances = dataset_variance_per_replica(ds_rep_minus_central)
        multi_variances.append(variances)
    return multi_variances
