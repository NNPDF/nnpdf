"""
closuretest/inconsistent_closuretest/multiclosure_inconsistent.py

Module containing all of the statistical estimators which are
averaged across multiple inconsistent fits. The actions
in this module are used to produce results which are plotted in
``multiclosure_inconsistent_output.py``

"""

import numpy as np

from reportengine import collect


""" To load several multiclosure fits. Useful for inconsistent closure test analysis """
multi_dataset_loader = collect("internal_multiclosure_dataset_loader", ("dataspecs",))

multi_dataset_fits_bias_replicas_variance_samples = collect(
    "dataset_fits_bias_replicas_variance_samples", ("dataspecs",)
)

multi_dataset_fits_bias_replicas_variance_samples_pdf_covmat = collect(
    "dataset_fits_bias_replicas_variance_samples_pdf_covmat", ("dataspecs",)
)

multi_fits_bootstrap_dataset_bias_variance = collect(
    "fits_bootstrap_dataset_bias_variance", ("dataspecs",)
)

multi_bias_variance_resampling_dataset = collect("bias_variance_resampling_dataset", ("dataspecs",))



## PDF SPACE
multi_fits_pdf_total_bias_variance = collect(
    "fits_pdf_total_bias_variance", ("dataspecs", )
)


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
