"""
closuretest/inconsistent_closuretest/multiclosure_inconsistent.py

Module containing all of the statistical estimators which are
averaged across multiple inconsistent fits. The actions
in this module are used to produce results which are plotted in
``multiclosure_inconsistent_output.py``

"""

import numpy as np
from validphys.closuretest.multiclosure import internal_multiclosure_dataset_loader



def consistent_internal_multiclosure_dataset_loader(
    dataset,
    consistent_pdfs,
    consistent_multiclosure_underlyinglaw,
    consistent_fits,
    t0_covmat_from_systematics,
):
    """like `internal_multiclosure_dataset_loader` but explicitly for consistent fits only"""
    return internal_multiclosure_dataset_loader(
        dataset,
        consistent_pdfs,
        consistent_multiclosure_underlyinglaw,
        consistent_fits,
        t0_covmat_from_systematics,
    )


def inconsistent_internal_multiclosure_dataset_loader(
    dataset,
    inconsistent_pdfs,
    inconsistent_multiclosure_underlyinglaw,
    inconsistent_fits,
    t0_covmat_from_systematics,
):
    """like `internal_multiclosure_dataset_loader` but explicitly for inconsistent fits only"""
    return internal_multiclosure_dataset_loader(
        dataset,
        inconsistent_pdfs,
        inconsistent_multiclosure_underlyinglaw,
        inconsistent_fits,
        t0_covmat_from_systematics,
    )



def dataset_replica_minus_central(internal_multiclosure_dataset_loader):
    """For a given dataset calculate the difference between theory prediction
    of each replica and central replica.
    Average over different fits and over all datapoints.

    Returns
    -------
    numpy.ndarray
            array of dimension Nreplicas
    """
    closures_th, law_th, covmat, _ = internal_multiclosure_dataset_loader
    replicas = np.asarray([th.error_members for th in closures_th])
    centrals = np.mean(replicas, axis=-1)

    diff = centrals[:, :, np.newaxis] - replicas
    diff = np.mean(diff, axis=(0, 1))
    return diff


def consistent_dataset_replica_minus_central(consistent_internal_multiclosure_dataset_loader):
    """like `dataset_replica_minus_central` but for explicitly consistent fits"""
    return dataset_replica_minus_central(consistent_internal_multiclosure_dataset_loader)


def inconsistent_dataset_replica_minus_central(inconsistent_internal_multiclosure_dataset_loader):
    """like `dataset_replica_minus_central` but for explicitly inconsistent fits"""
    return dataset_replica_minus_central(inconsistent_internal_multiclosure_dataset_loader)

