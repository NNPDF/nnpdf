"""
This module contains some helper functions that are used for the computation of nsigma in the
context of a multi-closure test.
"""

import dataclasses
import logging

import numpy as np
import pandas as pd

from nnpdf_data.coredata import CommonData
from reportengine import collect
from validphys.calcutils import calc_chi2
from validphys.convolution import central_predictions
from validphys.core import DataSetSpec

log = logging.getLogger(__name__)


@dataclasses.dataclass
class CentralChi2Data:
    value: float
    ndata: int
    dataset: DataSetSpec

    @property
    def reduced(self):
        return self.value / self.ndata


def central_member_chi2(
    central_predictions: pd.DataFrame,
    sqrt_covmat: np.ndarray,
    dataset: DataSetSpec,
    loaded_commondata_with_cuts: CommonData,
) -> CentralChi2Data:
    """
    Computes the chi2 value for a dataset.

    Parameters
    ----------
    central_predictions:
        The central predictions for the dataset.
    sqrt_covmat: np.ndarray
        The square root of the covariance matrix.
    dataset: DataSetSpec
        The dataset.
    loaded_commondata_with_cuts: nnpdf_data.coredata.CommonData

    Returns
    -------
    CentralChi2Data
    """
    diff = (
        central_predictions.values.ravel()
        - loaded_commondata_with_cuts.central_values.values.ravel()
    )

    value = calc_chi2(sqrt_covmat, diff)
    ndata = len(central_predictions)
    return CentralChi2Data(value=value, ndata=ndata, dataset=dataset)


def chi2_nsigma_deviation(central_member_chi2: CentralChi2Data) -> float:
    """
    Computes n_sigma as: (chi2 - ndata) / sqrt(2 * ndata)

    Parameters
    ----------
    central_member_chi2: CentralChi2Data

    Returns
    -------
    float
        The deviation in units of sigma.
    """
    diff = central_member_chi2.value - central_member_chi2.ndata
    return diff / np.sqrt(2 * central_member_chi2.ndata)


datasets_chi2_nsigma_deviation = collect("chi2_nsigma_deviation", ("data_input",))
"""
Collect the n_sigma values over list of ``dataset_input``.
"""


fits_datasets_chi2_nsigma_deviation = collect(
    "datasets_chi2_nsigma_deviation", ("fits", "fitcontext")
)
"""
Collects over fits and for all datasets the n_sigma values.
"""


fits_data = collect("data", ("fits", "fitinputcontext"))
"""
Collects the data for each fit.
"""


def is_weighted(fits_data: list) -> bool:
    """
    Returns whether the considered multiclosure tests has been weighted or not.
    If the weighted datasets are not the same for all fits,
    or there is more than one weighted dataset, an error is raised.

    Parameters
    ----------
    fits_data: list
        List of data for each fit.

    Returns
    -------
    bool
    """
    # Extract the set of unique weighted dataset names from all fits
    weighted_ds_sets = [{ds.name for ds in data.datasets if ds.weight != 1} for data in fits_data]

    # Ensure all fits have the same set of weighted datasets
    if len(set(frozenset(ds_set) for ds_set in weighted_ds_sets)) > 1:
        error_msg = "Weighted datasets are not the same for all fits in the same multiclosure test (dataspec)."
        log.error(error_msg)
        raise ValueError(error_msg)

    # Extract the single weighted dataset set (all should be identical)
    weighted_ds = next(iter(weighted_ds_sets))

    # Ensure there is exactly one weighted dataset
    if len(weighted_ds) > 1:
        error_msg = "Only one dataset can be weighted in a multiclosure test."
        log.error(error_msg)
        raise ValueError(error_msg)

    return bool(weighted_ds)


def n_fits(dataspecs):
    """
    Computes the total number of fits in the multiclosure test.
    If the number of fits is not the same across dataspecs it raises an error.
    """
    n_fits = set()
    for dataspec in dataspecs:
        n_fits.add(len(dataspec["fits"]))

    if len(n_fits) > 1:
        error_msg = "The number of fits is not the same across dataspecs."
        log.error(error_msg)
        raise ValueError(error_msg)

    return next(iter(n_fits))
