"""
covmat.py

Module containing actions which return constructed covariance matrices for
datasets/groups of datasets.

"""
from validphys.results_providers.covmat_construction import (
    covmat_from_systematics,
    datasets_covmat_from_systematics,
)


def experimental_covmat(loaded_commondata_with_cuts):
    """Returns the experimental covariance matrix. Details of how
    the covmat is constructed can be found in :py:func:`covmat_from_systematics`.
    The experimental covariance matrix uses the experimental central values
    to calculate the absolute uncertainties from the multiplicative systematics.

    Parameters
    ----------
    loaded_commondata_with_cuts: validphys.coredata.CommonData

    Returns
    -------
    covmat: np.array

    """
    return covmat_from_systematics(loaded_commondata_with_cuts)


def t0_covmat(loaded_commondata_with_cuts, dataset_t0_predictions):
    """Like :py:func:`experimental_covmat` except uses the t0 predictions
    to calculate the absolute constributions to the covmat from multiplicative
    uncertainties. For more info on the t0 predictions see
    :py:func:`validphys.results_providers.dataset_t0_predictions`.

    Parameters
    ----------
    loaded_commondata_with_cuts: validphys.coredata.CommonData
        commondata object for which to generate the covmat.
    dataset_t0_predictions: np.array
        1-D array with t0 predictions.

    Returns
    -------
    t0_covmat: np.array
        t0 covariance matrix

    """
    return covmat_from_systematics(
        loaded_commondata_with_cuts, dataset_t0_predictions)


def dataset_inputs_experimental_covmat(dataset_inputs_loaded_cd_with_cuts):
    """Like :py:func:`experimental_covmat` except for all data

    Parameters
    ----------
    dataset_inputs_loaded_cd_with_cuts: list[validphys.coredata.CommonData]
        The CommonData for all datasets defined in ``dataset_inputs``.

    Returns
    -------
    covmat: np.array
        Covariance matrix for list of datasets.
    """
    return datasets_covmat_from_systematics(dataset_inputs_loaded_cd_with_cuts)

def dataset_inputs_t0_covmat(
    dataset_inputs_loaded_cd_with_cuts, dataset_inputs_t0_predictions):
    """Like :py:func:`t0_covmat` except for all data

    Parameters
    ----------
    dataset_inputs_loaded_cd_with_cuts: list[validphys.coredata.CommonData]
        The CommonData for all datasets defined in ``dataset_inputs``.
    dataset_inputs_t0_predictions: list[np.array]
        The t0 predictions for all datasets.

    Returns
    -------
    t0_covmat: np.array
        t0 covariance matrix matrix for list of datasets.
    """
    return datasets_covmat_from_systematics(
        dataset_inputs_loaded_cd_with_cuts, dataset_inputs_t0_predictions)
