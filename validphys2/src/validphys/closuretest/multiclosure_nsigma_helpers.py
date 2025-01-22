"""
This module contains some helper functions that are used for the computation of nsigma in the
context of a multi-closure test.
"""

import dataclasses
import pandas as pd
import numpy as np
from scipy.stats import norm
from typing import Any


from validphys.core import DataSetSpec, PDF
from validphys.coredata import CommonData
from validphys.calcutils import calc_chi2
from validphys import convolution


@dataclasses.dataclass
class CentralChi2Data:
    value: float
    ndata: int
    dataset: DataSetSpec

    @property
    def reduced(self):
        return self.value / self.ndata


def central_predictions(dataset: DataSetSpec, pdf: PDF) -> pd.DataFrame:
    """
    Computes the central prediction (central PDF member) for a dataset.

    Parameters
    ----------
    dataset: validphys.core.DataSetSpec
    pdf: validphys.core.PDF

    Returns
    -------
    pd.DataFrame
        index is datapoints, column is the central prediction.
    """
    return convolution.central_predictions(dataset, pdf)


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
    loaded_commondata_with_cuts: validphys.coredata.CommonData

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


def compute_nsigma_critical_value(
    samples: Any, mu_0: float = 0, sigma: float = 1, alpha: float = 0.05
) -> tuple:
    """
    Computes the critical value for a 1-sided composite hypothesis test defined as:

    H0: mu = mu_0;
    H1: mu > mu_0;

    under the assumption that the test statistic (n_sigma) follows a normal distribution
    with mean `mu` and standard deviation `sigma`.

    Parameters
    ----------
    samples: np.ndarray
        The sample values of the test statistic.
    mu_0: float, default=0
        The null hypothesis value.
    sigma: float, default=1
        The standard deviation of the test statistic.
    alpha: float, default=0.05
        The significance level.

    Returns
    -------
    tuple: (c, z_alpha)
        c: float, deviation of samples mean from H0 mean in units of sigma.
        z_alpha: float, the critical value.
    """
    mu = np.mean(samples)
    c_sigma = sigma / np.sqrt(len(samples))
    c = (mu - mu_0) / c_sigma

    z_alpha = norm.ppf(1 - alpha)
    return c, z_alpha
