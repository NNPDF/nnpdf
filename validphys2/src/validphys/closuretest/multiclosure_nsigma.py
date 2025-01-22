"""
This module contains the functions to compute the consistency / inconsistency sets.

TODO
"""

import dataclasses
import pandas as pd
import numpy as np
import logging
from matplotlib.figure import Figure
from scipy.stats import norm
from typing import Any, Union, Generator

import reportengine
from reportengine.figure import figuregen
from reportengine import collect

from validphys.core import DataSetSpec, PDF
from validphys.coredata import CommonData
from validphys.calcutils import calc_chi2
from validphys import convolution
from validphys.closuretest.closure_checks import (
    check_fits_areclosures,
    check_fits_underlying_law_match,
)
from validphys import plotutils
from validphys.api import API

import sys

sys.path.insert(0, "./")
from multiclosure_nsigma_helpers import CentralChi2Data


log = logging.getLogger(__name__)

"""
Quantile range for computing the true positive rate and true negative rate.
"""
ALPHA_RANGE = np.linspace(0, 1.0, 100)


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


"""
Collect the n_sigma values over list of ``dataset_input``.
"""
datasets_chi2_nsigma_deviation = collect("chi2_nsigma_deviation", ("data_input",))


"""
Collects over fits and for all datasets the n_sigma values.
"""
fits_datasets_chi2_nsigma_deviation = collect(
    "datasets_chi2_nsigma_deviation", ("fits", "fitcontext")
)


"""
Collects the data for each fit.
"""
fits_data = collect("data", ("fits", "fitinputcontext"))


@dataclasses.dataclass
class MulticlosureNsigma:
    """
    Dataclass containing nsigma values for all datasets and fits,
    also used to keep track on whether the multiclosure fit is weighted or not.

    Attributes
    ----------
    nsigma_table: pd.DataFrame
        A table containing n_sigma values.
    is_weighted: bool
        Whether the fit was weighted.
    """

    nsigma_table: pd.DataFrame
    is_weighted: bool


@check_fits_areclosures
@check_fits_underlying_law_match
def multiclosurefits_nsigma(
    fits: reportengine.namespaces.NSList,
    fits_data: list,
    fits_datasets_chi2_nsigma_deviation: list,
    is_weighted: bool,
) -> MulticlosureNsigma:
    """
    Returns a table (dataframe) containing n_sigma values.
    Index: dataset names, Columns: Level 1 seeds (filterseed).

    Parameters
    ----------
    fits: NSList
        List of fits.
    fits_data: list
        List of data for each fit.
    fits_datasets_chi2_nsigma_deviation: list
        List of n_sigma values for each dataset for each fit.
    is_weighted: bool
        Used to keep track of whether the fit was weighted.

    Returns
    -------
    MulticlosureNsigma
    """
    res_dict = {}
    for fit, input_data, nsigma_data in zip(fits, fits_data, fits_datasets_chi2_nsigma_deviation):
        res_dict[fit.as_input()["closuretest"]["filterseed"]] = d = {}
        for ds, nsigma in zip(input_data.datasets, nsigma_data):
            d[ds.name] = nsigma
    return MulticlosureNsigma(is_weighted=is_weighted, nsigma_table=pd.DataFrame(res_dict))


"""
Collect the multiclosurefits_nsigma over dataspecs.
"""
dataspecs_multiclosurefits_nsigma = collect("multiclosurefits_nsigma", ("dataspecs",))


@dataclasses.dataclass
class NsigmaAlpha:
    """
    Dataclass storing the set 1 alpha values (can be used both for the set 1 and its complement).

    Attributes
    ----------
    alpha_dict: dict
        A dictionary containing the set 1 alpha values.
    is_weighted: bool
        Whether the fit was weighted.
    """

    alpha_dict: dict
    is_weighted: bool


def def_of_nsigma_alpha(
    multiclosurefits_nsigma: pd.DataFrame, weighted_dataset: str, complement: bool = False
) -> NsigmaAlpha:
    """
    Defines how the set 1 alpha values are computed.
    It allows to compute both the set 1 and its complement.

    Parameters
    ----------
    multiclosurefits_nsigma: pd.DataFrame
        The nsigma table.
    weighted_dataset: str
        The name of the weighted dataset.
    complement: bool, default=False
        Whether to compute the complement set 1 alpha values.

    Returns
    -------
    NsigmaAlpha
    """
    df = multiclosurefits_nsigma.nsigma_table
    nsigma_values = df[df.index == weighted_dataset].values.flatten()
    set1_alpha = {}
    for alpha in ALPHA_RANGE:
        z_alpha = norm.ppf(1 - alpha)
        if complement:
            fit_idxs = np.where(nsigma_values < z_alpha)[0]
        else:
            fit_idxs = np.where(nsigma_values > z_alpha)[0]
        # save it as set to allow for easy intersection with other sets
        set1_alpha[alpha] = set(df.columns[fit_idxs])

    return NsigmaAlpha(alpha_dict=set1_alpha, is_weighted=multiclosurefits_nsigma.is_weighted)


def nsigma_alpha(multiclosurefits_nsigma: pd.DataFrame, weighted_dataset: str) -> NsigmaAlpha:
    """
    Computes the set 1 alpha values.
    """
    return def_of_nsigma_alpha(multiclosurefits_nsigma, weighted_dataset, complement=False)


"""
Collect set 1 alpha over dataspecs.
"""
dataspecs_nsigma_alpha = collect("nsigma_alpha", ("dataspecs",))


def comp_nsigma_alpha(multiclosurefits_nsigma: pd.DataFrame, weighted_dataset: str) -> NsigmaAlpha:
    """
    Computes the complement set 1 alpha values.
    """
    return def_of_nsigma_alpha(multiclosurefits_nsigma, weighted_dataset, complement=True)


"""
Collect complement set 1 alpha over dataspecs.
"""
dataspecs_comp_nsigma_alpha = collect("comp_nsigma_alpha", ("dataspecs",))


def set_1_alpha(dataspecs_nsigma_alpha: list) -> dict:
    """
    Returns the set 1 alpha values, these are defined as

    1_{\alpha} = {i | n_{\sigma}^{i} > Z_{\alpha}}

    where i is the index of the fit and n_{\sigma}^{i} is the n-sigma value computed
    for fit i.

    Parameters
    ----------
    dataspecs_nsigma_alpha: list
        List of NsigmaAlpha dataclasses.

    Returns
    -------
    dict
    """
    for dataspec_nsigma in dataspecs_nsigma_alpha:
        if not dataspec_nsigma.is_weighted:
            return dataspec_nsigma.alpha_dict


def set_3_alpha(dataspecs_nsigma_alpha: list) -> dict:
    """
    Same as the set 1 alpha values, but for the weighted fits.

    3_{\alpha} = {i | n_{weighted, \sigma}^{i} > Z_{\alpha}}

    where i is the index of the fit and n_{weighted, \sigma}^{i} is the n-sigma value computed
    on the weighted dataset for fit i.

    Parameters
    ----------
    dataspecs_nsigma_alpha: list
        List of NsigmaAlpha dataclasses.

    Returns
    -------
    dict
    """
    for dataspec_nsigma in dataspecs_nsigma_alpha:
        if dataspec_nsigma.is_weighted:
            return dataspec_nsigma.alpha_dict


def comp_set_1_alpha(dataspecs_comp_nsigma_alpha: list) -> dict:
    """
    Returns the complement set 1 alpha values.
    """
    for dataspec_nsigma in dataspecs_comp_nsigma_alpha:
        if not dataspec_nsigma.is_weighted:
            return dataspec_nsigma.alpha_dict


def comp_set_3_alpha(dataspecs_comp_nsigma_alpha: list) -> dict:
    """
    Returns the complement set 3 alpha values.
    """
    for dataspec_nsigma in dataspecs_comp_nsigma_alpha:
        if dataspec_nsigma.is_weighted:
            return dataspec_nsigma.alpha_dict


def def_set_2(
    dataspecs_multiclosurefits_nsigma: list, weighted_dataset: str, complement: bool = False
) -> dict:
    """
    Defines how the set 2 alpha values are computed.
    It allows to compute both the set 2 and its complement.

    Parameters
    ----------
    dataspecs_multiclosurefits_nsigma: list
        List of MulticlosureNsigma dataclasses.
    weighted_dataset: str
        The name of the weighted dataset.
    complement: bool, default=False
        Whether to compute the complement set 2 alpha values.

    Returns
    -------
    dict
    """
    # Order the dataspecs so that the weighted dataset is the first one
    dataspecs_mct = []
    for mct_nsigma in dataspecs_multiclosurefits_nsigma:
        if mct_nsigma.is_weighted:
            dataspecs_mct.insert(0, mct_nsigma)
        else:
            dataspecs_mct.append(mct_nsigma)

    df_weight = dataspecs_mct[0].nsigma_table
    df_weight = df_weight[df_weight.index != weighted_dataset]

    df_ref = dataspecs_mct[1].nsigma_table
    df_ref = df_ref[df_ref.index != weighted_dataset]

    # ensure that weighted and reference dfs have the columns in the same order
    # (needed to properly compare fits)
    df_ref = df_ref[df_weight.columns]

    set2_alpha = {}

    for alpha in ALPHA_RANGE:
        z_alpha = norm.ppf(1 - alpha)

        if complement:
            columns_bools = np.any((df_weight - df_ref).values < z_alpha, axis=0)
        else:
            columns_bools = np.any((df_weight - df_ref).values > z_alpha, axis=0)

        set2_alpha[alpha] = set(df_weight.columns[columns_bools])

    return set2_alpha


def set_2_alpha(dataspecs_multiclosurefits_nsigma: list, weighted_dataset: str) -> dict:
    """
    Computes the set 2 alpha values. The set 2 is defined as:

    2_{\alpha} = {i | n_{weighted, \sigma}^{i} - n_{ref, \sigma}^{i}>  + Z_{\alpha}}

    where the n-sigma is computed on all datasets that are not the weighted dataset.
    Moreover if for a fit i any dataset has a n-sigma value greater than Z_{\alpha}, then
    the fit i is included in the set.
    """
    return def_set_2(dataspecs_multiclosurefits_nsigma, weighted_dataset, complement=False)


def comp_set_2_alpha(dataspecs_multiclosurefits_nsigma: list, weighted_dataset: str) -> dict:
    """
    Computes the complement set 2 alpha values.
    """
    return def_set_2(dataspecs_multiclosurefits_nsigma, weighted_dataset, complement=True)
