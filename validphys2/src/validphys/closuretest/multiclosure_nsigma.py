"""
This module contains the functions to compute the consistency / inconsistency sets.

Assuming that we have two datasets A and B, and that we are investigating whether A is consistent or not
we can define the following sets:


1⍺ = {i | nσi > Z⍺} ...

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
    str or None
        Name of the weighted dataset.
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


def n_fits(dataspecs):
    """
    Computes the total number of fits in the multiclosure test.
    If the number of fits is not the same across dataspecs it raises an error.
    """
    n_fits = set()
    for dataspec in dataspecs:        
        n_fits.add(len(dataspec['fits']))
    
    if len(n_fits) > 1:
        error_msg = "The number of fits is not the same across dataspecs."
        log.error(error_msg)
        raise ValueError(error_msg)

    return next(iter(n_fits))


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
        set1_alpha[alpha] = df.columns[fit_idxs].tolist()

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
    Returns the set 1 alpha values.

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
    Same as the set 1 alpha values, but for the weighted datasets.

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
    Same as the complement set 1 alpha values, but for the weighted datasets.
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

        columns = df_weight.columns[columns_bools].to_list()

        set2_alpha[alpha] = columns

    return set2_alpha


def set_2_alpha(dataspecs_multiclosurefits_nsigma: list, weighted_dataset: str) -> dict:
    """
    Computes the set 2 alpha values.
    """
    return def_set_2(dataspecs_multiclosurefits_nsigma, weighted_dataset, complement=False)


def comp_set_2_alpha(dataspecs_multiclosurefits_nsigma: list, weighted_dataset: str) -> dict:
    """
    Computes the complement set 2 alpha values.
    """
    return def_set_2(dataspecs_multiclosurefits_nsigma, weighted_dataset, complement=True)
