r"""
This module contains the functions used in Sec. 4 of paper: arXiv: 2503.17447

`set_1`, `set_2`, and `set_3` correspond to \(\S_1\), \(\S_2\), and \(\S_3\) in Eq. 4.3, 4.6, and 4.7.
"""

import dataclasses

import numpy as np
import pandas as pd
from scipy.stats import norm

import reportengine
from reportengine import collect
from validphys.closuretest.closure_checks import (
    check_fits_areclosures,
    check_fits_underlying_law_match,
)

Z_ALPHA_RANGE = norm.ppf(1 - np.linspace(0, 0.5, 100))
"""
Quantile range for computing the true positive rate and true negative rate.
"""


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


dataspecs_multiclosurefits_nsigma = collect("multiclosurefits_nsigma", ("dataspecs",))
"""
Collect the multiclosurefits_nsigma over dataspecs.
"""


@dataclasses.dataclass
class NsigmaAlpha:
    """
    Dataclass storing the set 1 values (can be used both for the set 1 and its complement).

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
    for z_alpha in Z_ALPHA_RANGE:
        if complement:
            fit_idxs = np.where(nsigma_values < z_alpha)[0]
        else:
            fit_idxs = np.where(nsigma_values > z_alpha)[0]
        # save it as set to allow for easy intersection with other sets
        set1_alpha[z_alpha] = set(df.columns[fit_idxs])

    return NsigmaAlpha(alpha_dict=set1_alpha, is_weighted=multiclosurefits_nsigma.is_weighted)


def nsigma_alpha(multiclosurefits_nsigma: pd.DataFrame, weighted_dataset: str) -> NsigmaAlpha:
    """
    Computes the set 1 alpha values.
    """
    return def_of_nsigma_alpha(multiclosurefits_nsigma, weighted_dataset, complement=False)


dataspecs_nsigma_alpha = collect("nsigma_alpha", ("dataspecs",))
"""
Collect set 1 alpha over dataspecs.
"""


def comp_nsigma_alpha(multiclosurefits_nsigma: pd.DataFrame, weighted_dataset: str) -> NsigmaAlpha:
    """
    Computes the complement set 1 alpha values.
    """
    return def_of_nsigma_alpha(multiclosurefits_nsigma, weighted_dataset, complement=True)


dataspecs_comp_nsigma_alpha = collect("comp_nsigma_alpha", ("dataspecs",))
"""
Collect complement set 1 alpha over dataspecs.
"""


def set_1(dataspecs_nsigma_alpha: list) -> dict:
    r"""
    Returns the set 1 alpha values, these are defined as

    S_1 = {j | n_{\sigma}^{j} > Z_{\alpha}}

    where j is the index of the fit and n_{\sigma}^{j} is the n-sigma value computed
    for fit j.

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


def set_2(dataspecs_nsigma_alpha: list) -> dict:
    r"""
    Same as the set 1 alpha values, but for the weighted fits.

    S_2 = {i | n_{weighted, \sigma}^{i} > Z_{\alpha}}

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


def comp_set_1(dataspecs_comp_nsigma_alpha: list) -> dict:
    """
    Returns the complement set 1 alpha values.
    """
    for dataspec_nsigma in dataspecs_comp_nsigma_alpha:
        if not dataspec_nsigma.is_weighted:
            return dataspec_nsigma.alpha_dict


def def_set_3(
    dataspecs_multiclosurefits_nsigma: list, weighted_dataset: str, complement: bool = False
) -> dict:
    """
    Defines how the set 3 values are computed.
    It allows to compute both the set 3 and its complement.

    Parameters
    ----------
    dataspecs_multiclosurefits_nsigma: list
        List of MulticlosureNsigma dataclasses.
    weighted_dataset: str
        The name of the weighted dataset.
    complement: bool, default=False
        Whether to compute the complement set 3 alpha values.

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

    set3_alpha = {}

    for z_alpha in Z_ALPHA_RANGE:
        if complement:
            columns_bools = np.any((df_weight - df_ref).values < z_alpha, axis=0)
        else:
            columns_bools = np.any((df_weight - df_ref).values > z_alpha, axis=0)

        set3_alpha[z_alpha] = set(df_weight.columns[columns_bools])

    return set3_alpha


def set_3(dataspecs_multiclosurefits_nsigma: list, weighted_dataset: str) -> dict:
    r"""
    Computes the set 3 alpha values. The set 3 is defined as:

    S_3 = {i | n_{weighted, \sigma}^{i} - n_{ref, \sigma}^{i}>  + Z_{\alpha}}

    where the n-sigma is computed on all datasets that are not the weighted dataset.
    Moreover if for a fit i any dataset has a n-sigma value greater than Z_{\alpha}, then
    the fit i is included in the set.
    """
    return def_set_3(dataspecs_multiclosurefits_nsigma, weighted_dataset, complement=False)


def probability_inconsistent(set_1, set_2, set_3, comp_set_1, n_fits, weighted_dataset):
    """
    The set of inconsistent fits can be defined in different ways, two possible cases are:

    1. C_2: (S_1 intersect S_2) union (S_3)

    2. C_3 = S_1 union (~S_1 intersect S_3)

    The probability of a dataset being inconsistent is defined as:
            P(inconsistent) = |I_alpha| / N
    where N is the total number of fits.
    """

    C_2 = []
    C_3 = []
    for z_alpha in set_3.keys():
        set_3_inters_1 = set_3[z_alpha].intersection(set_1[z_alpha])
        s2 = set_2[z_alpha]

        set_tagged_fits = set_3_inters_1.union(s2)

        C_2.append(len(set_tagged_fits) / n_fits)

        set_tagged_fits = set_tagged_fits.union(comp_set_1[z_alpha].intersection(set_3[z_alpha]))
        C_3.append(len(set_tagged_fits) / n_fits)

    return C_3, C_2
