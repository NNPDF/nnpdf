"""
Module containing functions for computation and representation of the true positive rate (TPR) and true negative rate (TNR)
computed from the nsigma test in a multiclosure test.

The module is structured as follows:

Single Multiclosure Test
------------------------

1. Compute the nsigma deviation for all fits in a multiclosure test (`multiclosurefits_nsigma`).

2. Given the nsigma values, compute the true positive rate (TPR) and true negative rate (TNR) for a range of alpha values `nsigma_true_negatives_positives`.

   For example, the true positive rate (TPR) for an inconsistent dataset is defined as:

   .. math::

       TPR = \frac{m}{n}, \quad m = \lvert \{n_{\sigma_i} : n_{\sigma_i} > Z_\alpha\} \rvert

   where:
   - \( n_{\sigma_i} \) is the nsigma value (of the inconsistent dataset) for the \( i \)-th fit,
   - \( Z_\alpha \) is the critical value corresponding to the significance level \( \alpha \),
   - \( n \) is the total number of fits.


Weighted vs. Reference Multiclosure Test
----------------------------------------
This comparison collects results across dataspecs.
One of the dataspecs should be the reference fit and the other the weighted fit.


1. Collect over dataspecs `multiclosurefits_nsigma` (`dataspecs_multiclosurefits_nsigma`). This allows for the computation
   of the weighted fit discriminator, which is a list of indexes of inconsistent fits based on the weighted fit discriminator.

2. Collect over dataspecs `nsigma_true_negatives_positives` (`dataspecs_nsigma_true_negatives_positives`).

3. Update the TPR and TNR specs based on the weighted fit discriminator. This is done by processing the TPR and TNR specs.

4. Evaluate the TPR and TNR values `tpr_tnr_weighted_dataspecs`.
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


log = logging.getLogger(__name__)

"""
Quantile range for computing the true positive rate and true negative rate.
"""
ALPHA_RANGE = np.linspace(0, 1.0, 100)


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
class MulticlosureNSigmaData:
    """
    Attributes
    ----------
    weighted: bool
        Flag to indicate if there is at least one dataset in fits that is fitted
        with an extra weight.
    nsigma_table: pd.DataFrame
        Table with n_sigma values.
    """

    weighted: bool
    nsigma_table: pd.DataFrame


@check_fits_areclosures
@check_fits_underlying_law_match
def multiclosurefits_nsigma(
    fits: reportengine.namespaces.NSList,
    fits_data: list,
    fits_datasets_chi2_nsigma_deviation: list,
    weighted_datasets: list,
) -> MulticlosureNSigmaData:
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

    Returns
    -------
    MulticlosureNSigmaData
    """
    res_dict = {}
    for fit, input_data, nsigma_data in zip(fits, fits_data, fits_datasets_chi2_nsigma_deviation):
        res_dict[fit.as_input()["closuretest"]["filterseed"]] = d = {}
        for ds, nsigma in zip(input_data.datasets, nsigma_data):
            d[ds.name] = nsigma
    return MulticlosureNSigmaData(
        weighted=bool(weighted_datasets), nsigma_table=pd.DataFrame(res_dict)
    )


"""
Collect nsigma table over all dataspecs.
"""
dataspecs_multiclosurefits_nsigma = collect("multiclosurefits_nsigma", ("dataspecs",))


def weighted_datasets(fits_data: list) -> list:
    """
    Returns the names of the datasets that are weighted.
    If the weighted datasets are not the same for all fits, an error is raised.

    Parameters
    ----------
    fits_data: list
        List of data for each fit.

    Returns
    -------
    list
        List of dataset names that are weighted.
    """
    # extract datasets that are weighted for each fit
    # need to use frozenset since sets are not hashable
    weighted_ds_fits = {
        frozenset(ds.name for ds in data.datasets if ds.weight != 1) for data in fits_data
    }

    # log warning if weighed datasets are not the same for all fits
    if len(weighted_ds_fits) > 1:
        log.error(
            "WARNING: weighted datasets are not the same for all fits in the same multiclosure test (dataspec)"
        )
        raise ValueError(
            "Weighted datasets are not the same for all fits in the same multiclosure test (dataspec)"
        )

    return list(list(weighted_ds_fits)[0])


# TODO: @check_fits_areinconsistent
@check_fits_areclosures
def inconsistent_datasets(fits: reportengine.namespaces.NSList) -> list:
    """
    NOTE: will need to be adjusted once new inconsistency code is in use in validphys.
        once the new inconsistency code has been tested and merged add the check.

    Returns the names of the datasets that are inconsistent.
    If the inconsistent datasets are not the same for all fits, an error is raised.

    Parameters
    ----------
    fits: NSList
        List of fits.

    Returns
    -------
    list
        List of dataset names that are inconsistent.
    """
    inconsistent_ds = {
        frozenset(ds for ds in fit.as_input()["inconsistent_datasets"]) for fit in fits
    }

    if len(inconsistent_ds) > 1:
        log.error(
            "WARNING: inconsistent datasets are not the same for all fits in the same multiclosure test (dataspec)"
        )
        raise ValueError(
            "Inconsistent datasets are not the same for all fits in the same multiclosure test (dataspec)"
        )

    return list(list(inconsistent_ds)[0])


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


@dataclasses.dataclass
class TPosNegSpec:
    """

    Attributes
    ----------
    m: int
        size of ~critical region: |{n_sigma^i; n_sigma^i < Z_alpha}| .
    n: int
        total number of fits.
    alpha: float
        quantile value.
    weighted_dataset: bool
        flag to indicate if the dataset is weighted.
    dataset: str
        name of the dataset.
    """

    m: int
    n: int
    alpha: float
    weighted_dataset: bool
    dataset: str


@dataclasses.dataclass
class TPRSpec(TPosNegSpec):
    """
    Attributes
    ----------
    inconsistent_fit_idxs: list
        indexes of fits {i; n_sigma^i > Z_alpha}.
    """

    inconsistent_fit_idxs: list


@dataclasses.dataclass
class TNRSpec(TPosNegSpec):
    """
    Attributes
    ----------
    consistent_fit_idxs: list
        indexes of fits {i; n_sigma^i < Z_alpha}.
    """

    consistent_fit_idxs: list


def compute_tpr_specs(
    nsigma_values: Any,
    mu_0: float = 0,
    sigma: float = 1,
    alpha: float = 0.05,
    dataset: str = None,
    weighted_dataset: bool = False,
) -> TPRSpec:
    """
    Computes the specs needed to compute the True positive rate.

    Parameters
    ----------
    nsigma_values: np.ndarray (shape=(Nfits,))
        The sample values of the test statistic.
    mu_0: float, default=0
        The null hypothesis value.
    sigma: float, default=1
        The standard deviation of the test statistic.
    alpha: float, default=0.05
        The significance level.
    dataset: str, default=None
        The name of the dataset.
    weighted_dataset: bool, default=False
        Flag to indicate if the dataset is weighted.
    Returns
    -------
    TPRSpec
    """
    _, z_alpha = compute_nsigma_critical_value(nsigma_values, mu_0=mu_0, sigma=sigma, alpha=alpha)
    # datasets flagged as inconsistent by the hypothesis test
    flagged_datasets = np.where(nsigma_values > z_alpha)[0]

    return TPRSpec(
        m=len(flagged_datasets),
        n=len(nsigma_values),
        alpha=alpha,
        dataset=dataset,
        weighted_dataset=weighted_dataset,
        inconsistent_fit_idxs=list(flagged_datasets),
    )


def compute_tnr_specs(
    nsigma_values: Any,
    mu_0: float = 0,
    sigma: float = 1,
    alpha: float = 0.05,
    dataset: str = None,
    weighted_dataset: bool = False,
) -> TNRSpec:
    """
    Computes the specs needed to compute the True negative rate.

    Parameters
    ----------
    nsigma_values: np.ndarray, (shape=(Nfits,))
        The sample values of the test statistic.
    mu_0: float, default=0
        The null hypothesis value.
    sigma: float, default=1
        The standard deviation of the test statistic.
    alpha: float, default=0.05
        The significance level.
    dataset: str, default=None
        The name of the dataset.
    weighted_dataset: bool, default=False
        Flag to indicate if the dataset is weighted.

    Returns
    -------
    float
    """
    _, z_alpha = compute_nsigma_critical_value(nsigma_values, mu_0=mu_0, sigma=sigma, alpha=alpha)
    unflagged_datasets = np.where(nsigma_values < z_alpha)[0]

    return TNRSpec(
        m=len(unflagged_datasets),
        n=len(nsigma_values),
        alpha=alpha,
        dataset=dataset,
        weighted_dataset=weighted_dataset,
        consistent_fit_idxs=list(unflagged_datasets),
    )


def compute_tpr(tpr_specs: TPRSpec) -> float:
    """
    Computes the True positive rate.

    Parameters
    ----------
    tpr_specs: TPRSpec

    Returns
    -------
    float
    """
    return tpr_specs.m / tpr_specs.n


def compute_tnr(tnr_specs: TNRSpec) -> float:
    """
    Computes the True negative rate.

    Parameters
    ----------
    tnr_specs: TNRSpec

    Returns
    -------
    float
    """
    return tnr_specs.m / tnr_specs.n


@dataclasses.dataclass
class NsigmaTPNR:
    """
    Attributes
    ----------
    weighted: bool
    rates_table: pd.DataFrame
    """

    weighted: bool
    rates_table: pd.DataFrame


@check_fits_areclosures
@check_fits_underlying_law_match
def nsigma_true_negatives_positives(
    multiclosurefits_nsigma: MulticlosureNSigmaData,
    inconsistent_datasets: list,
    weighted_datasets: list,
) -> NsigmaTPNR:
    """
    True Negative Rate (TNR) is computed for consistent datasets.
    True Positive Rate (TPR) is computed for inconsistent datasets.

    The TNR and TPR are computed for a range of alpha values and for each dataset.

    Parameters
    ----------
    multiclosurefits_nsigma: MulticlosureNSigmaData
        The nsigma values.
    inconsistent_datasets: list
        List of inconsistent datasets.
    weighted_datasets: list
        List of weighted datasets.
    Returns
    -------
    NsigmaTPNR
    """
    true_false_positives_dict = {}

    nsigma_table = multiclosurefits_nsigma.nsigma_table

    for ds_name in nsigma_table.T:

        rate_true_negative = (
            []
        )  # is the amount of consistent datasets that are not classified as inconsistent
        rate_true_positive = (
            []
        )  # is the amount of inconsistent datasets that are classified as inconsistent

        nsigma_values = nsigma_table.loc[ds_name].values

        for alpha in ALPHA_RANGE:

            if ds_name in inconsistent_datasets:
                tpr = compute_tpr_specs(
                    nsigma_values,
                    mu_0=0,
                    sigma=1,
                    alpha=alpha,
                    dataset=ds_name,
                    weighted_dataset=ds_name in weighted_datasets,
                )
                rate_true_positive.append(tpr)

            else:
                tnr = compute_tnr_specs(
                    nsigma_values,
                    mu_0=0,
                    sigma=1,
                    alpha=alpha,
                    dataset=ds_name,
                    weighted_dataset=ds_name in weighted_datasets,
                )
                rate_true_negative.append(tnr)

        true_false_positives_dict[ds_name] = {"TPR": rate_true_positive, "TNR": rate_true_negative}

    rows = []
    for dataset, metrics in true_false_positives_dict.items():
        for metric_name, values in metrics.items():
            for alpha, value in zip(ALPHA_RANGE, values):
                rows.append(
                    {"Dataset": dataset, "Alpha": alpha, "Metric": metric_name, "Value": value}
                )
    rates_table = pd.DataFrame(rows)

    return NsigmaTPNR(weighted=bool(weighted_datasets), rates_table=rates_table)


"""
Collect the true positive and negative rates over each dataspec.
"""
dataspecs_nsigma_true_negatives_positives = collect(
    "nsigma_true_negatives_positives", ("dataspecs",)
)


def weighted_fit_discriminator(
    dataspecs_multiclosurefits_nsigma: list, nsigma_weighted_fit_threshold: Union[float, int] = 1.0
) -> list:
    """
    In the weighted fit procedure a dataset is flagged as inconsistent if the chi2
    of a dataset in the weighted fit exceeds the chi2 of the same dataset in the reference fit
    by more that a certain threshold n_sigma_th.

    The condition for a dataset to be flagged as inconsistent is:

    chi2_weighted - chi2_reference > sqrt(2 * Ndata) * n_sigma_th

    and can be rewritten as:

    n_sigma_weighted - n_sigma_reference > n_sigma_th

    Parameters
    ----------
    dataspecs_multiclosurefits_nsigma: list
        List of nsigma values for different dataspecs.
    nsigma_weighted_fit_threshold: float, default=3.0

    Returns
    -------
    list
        List of indexes of inconsistent fits.
    """
    log.info("Computing weighted fit discriminator")
    log.info(f"nsigma_weighted_fit_threshold = {nsigma_weighted_fit_threshold}")

    # order so that weighted fit is first
    dat = []
    for mct_nsigma in dataspecs_multiclosurefits_nsigma:
        if mct_nsigma.weighted:
            dat.insert(0, mct_nsigma.nsigma_table)
        else:
            dat.append(mct_nsigma.nsigma_table)

    if len(dat) != 2:
        raise ValueError(f"Expected two dataspecs, got {len(dat)} instead")
    # ensure that the fits for different dataspecs in the runcard are ordered with the same filterseeds
    diffs = dat[0][dat[1].columns].values - dat[1].values

    diffs = dat[0].values - dat[1].values
    fit_idxs_bool = np.any(diffs > nsigma_weighted_fit_threshold, axis=0)

    # Return indexes of inconsistent fits
    return list(np.where(fit_idxs_bool)[0])


def process_TNR(TNR: TNRSpec, discriminator_list: list) -> None:
    """
    Helper function to process the TNR specs so as to reduce the number of fits rightly classified as consistent
    based on the weighted fit discriminator.

    Parameters
    ----------
    TNR: TNRSpec
    discriminator_list: list
        List of indexes of inconsistent fits.

    Returns
    -------
    None
    """
    # only process if the dataset is weighted
    if TNR.weighted_dataset:
        new_idxs = set(TNR.consistent_fit_idxs).difference(discriminator_list)
        new_idxs = list(new_idxs)
        TNR.m = len(new_idxs)
        TNR.consistent_fit_idxs = new_idxs

    return TNR


def process_TPR(TPR: TPRSpec, discriminator_list: list) -> None:
    """
    Helper function to process the TPR specs so as to augment the number of fits rightly classified as inconsistent
    based on the weighted fit discriminator.

    Parameters
    ----------
    TPR: TPRSpec
    discriminator_list: list
        List of indexes of inconsistent fits.

    Returns
    -------
    None
    """
    # only process if the dataset is weighted
    if TPR.weighted_dataset:
        new_idxs = set(TPR.inconsistent_fit_idxs).union(discriminator_list)
        new_idxs = list(new_idxs)
        TPR.m = len(new_idxs)
        TPR.inconsistent_fit_idxs = new_idxs

    return TPR


def evaluate_nsigma_true_negatives_positives(
    nsigma_true_negatives_positives: NsigmaTPNR, weighted_fit_discriminator: list
) -> NsigmaTPNR:
    """
    Replace the TPR and TNR specs with the actual rates, taking into account the weighted fit discriminator
    when the weighted flag is set to True.

    Parameters
    ----------
    nsigma_true_negatives_positives: NsigmaTPNR

    Returns
    -------
    NsigmaTPNR
    """
    rates_table = nsigma_true_negatives_positives.rates_table

    if nsigma_true_negatives_positives.weighted:
        log.info("Computing TPR and TNR values using weighted fit discriminator")
        # process TPR and TNR specs so as to reduce / augment m based on the weighted fit discriminator
        rates_table.loc[rates_table["Metric"] == "TPR", "Value"] = rates_table.loc[
            rates_table["Metric"] == "TPR", "Value"
        ].apply(process_TPR, discriminator_list=weighted_fit_discriminator)
        rates_table.loc[rates_table["Metric"] == "TNR", "Value"] = rates_table.loc[
            rates_table["Metric"] == "TNR", "Value"
        ].apply(process_TNR, discriminator_list=weighted_fit_discriminator)

    rates_table.loc[rates_table["Metric"] == "TPR", "Value"] = rates_table["Value"][
        rates_table["Metric"] == "TPR"
    ].apply(compute_tpr)
    rates_table.loc[rates_table["Metric"] == "TNR", "Value"] = rates_table["Value"][
        rates_table["Metric"] == "TNR"
    ].apply(compute_tnr)

    return nsigma_true_negatives_positives


def tpr_tnr_weighted_dataspecs(
    dataspecs_nsigma_true_negatives_positives: list, weighted_fit_discriminator: list
) -> list:
    """
    Computes the TPR and TNR values for different dataspecs using the weighted fit discriminator
    when the weighted flag is set to True.

    Parameters
    ----------
    dataspecs_nsigma_true_negatives_positives: list

    weighted_fit_discriminator: list

    Returns
    -------
    list
    """
    log.info(f"weighted_fit_discriminator = {weighted_fit_discriminator}")
    # make sure that the weighted fits are the first in the dataspecs
    dataspec_nsigma = []
    for dataspec in dataspecs_nsigma_true_negatives_positives:
        if dataspec.weighted:
            dataspec_nsigma.insert(0, dataspec)
        else:
            dataspec_nsigma.append(dataspec)

    for dataspec in dataspec_nsigma:
        evaluate_nsigma_true_negatives_positives(dataspec, weighted_fit_discriminator)

    return dataspecs_nsigma_true_negatives_positives


@figuregen
def plot_false_true_positives_nsigma_weighted_fits(
    tpr_tnr_weighted_dataspecs: list, dataspecs: list
) -> Generator[Figure, None, None]:
    """
    TODO
    Compares TPR (and TNR) for different dataspecs for the same dataset.

    Parameters
    ----------
    tpr_tnr_weighted_dataspecs:

    Yields
    ------
    matplotlib.Figure
    """
    ict_datasets = inconsistent_datasets(dataspecs[0]["fits"])
    if set(ict_datasets) != set(inconsistent_datasets(dataspecs[1]["fits"])):
        raise ValueError(
            "Inconsistent datasets are not the same among different multiclosure tests"
        )

    list_dfs = [nsigmatpnr.rates_table for nsigmatpnr in tpr_tnr_weighted_dataspecs]
    # datasets should be all the same for different dataspecs
    datasets = list_dfs[0].Dataset.unique()

    for grp_w, grp_ref in zip(list_dfs[0].groupby("Dataset"), list_dfs[1].groupby("Dataset")):
        dataset = grp_w[0]

        fig, ax = plotutils.subplots()
        ax.set_ylim(0, 1.1)
        ax.set_xlim(0, 1.1)
        ax.set_xlabel(r"$\alpha$")

        if dataset in ict_datasets:
            ax.set_title(f"Inconsistent dataset: {dataset}")
            ax.plot(grp_w[1]["Alpha"], grp_w[1]["Value"], label=f"TPR, weighted")
            ax.plot(grp_ref[1]["Alpha"], grp_ref[1]["Value"], label=f"TPR, reference")
            ax.set_ylabel("True Positive Rate")
        else:
            ax.set_title(f"Consistent dataset: {dataset}")
            ax.plot(grp_w[1]["Alpha"], grp_w[1]["Value"], label=f"TNR, weighted")
            ax.plot(grp_ref[1]["Alpha"], grp_ref[1]["Value"], label=f"TNR, reference")
            ax.set_ylabel("True Negative Rate")
        ax.legend()
        yield fig
