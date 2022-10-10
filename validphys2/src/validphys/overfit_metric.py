"""
overfit_metric.py

This module contains the functions used to calculate the overfit metric and 
produce the corresponding tables and figures.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats

from reportengine import collect
from reportengine.figure import figure
from reportengine.table import table

from validphys.checks import check_at_least_two_replicas

log = logging.getLogger(__name__)

preds = collect("predictions",("dataset_inputs",))

def _create_new_val_pseudodata(pdf_data_index, fit_data_indices_list):
    """Loads all validation pseudodata replicas used during the fiting of the
    pdf replicas

    Returns
    -------
    np.ndarray
        (nrep,ndata) sized numpy array containing the validation data used to
        fit the pdfs.
    """
    vl_data_fitrep = []
    for fitreplica_info in fit_data_indices_list:
        vl_data_fitrep.append(
            fitreplica_info.pseudodata.loc[pdf_data_index.val_idx]
        )
    return np.array(vl_data_fitrep)[:, :, 0]


@check_at_least_two_replicas
def calculate_chi2s_per_replica(
    pdf, # for the check
    fit_code_version,
    recreate_pdf_pseudodata_no_table,
    preds,
    dataset_inputs,
    groups_covmat_no_table,
):
    """Calculates, for each PDF replica, the chi2 of the validation with the
    pseudodata generated for all other replicas in the fit

    Parameters
    ----------
    recreate_pdf_pseudodata_no_table : list[namedtuple]
        List of namedtuples, each of which contains a dataframe
        containing all the data points, the training indices, and
        the validation indices.
    preds : list[pd.core.frame.DataFrame]
        List of pandas dataframes, each containing the predictions of the pdf
        replicas for a dataset_input
    dataset_inputs : list[DatasetInput]
    groups_covmat_no_table : pdf.core.frame.DataFrame

    Returns
    -------
    np.ndarray
        (Npdfs, Npdfs) sized matrix containing the chi2 of a pdf replica
        calculated to a given psuedodata replica. The diagonal values correspond
        to the cases where the PDF replica has been fitted to the coresponding
        pseudodata replica
    """
    fit_name = fit_code_version.columns[0]
    nnpdf_version = fit_code_version[fit_name]['nnpdf']
    if nnpdf_version>='4.0.5':
        pp = []
        for i, dss in enumerate(dataset_inputs):
            preds_witout_cv = preds[i].drop(0, axis=1)
            df = pd.concat({dss.name: preds_witout_cv}, names=["dataset"])
            pp.append(df)

        PDF_predictions = pd.concat(pp)

        chi2s_per_replica = []
        for enum, pdf_data_index in enumerate(recreate_pdf_pseudodata_no_table):

            prediction_filter = pdf_data_index.val_idx.droplevel(level=0)
            prediction_filter.rename(["dataset", "data"], inplace=True)
            PDF_predictions_val = PDF_predictions.loc[prediction_filter]
            PDF_predictions_val = PDF_predictions_val.values[:, enum]

            new_val_pseudodata_list = _create_new_val_pseudodata(
                pdf_data_index, recreate_pdf_pseudodata_no_table
            )

            invcovmat_vl = np.linalg.inv(
                groups_covmat_no_table[pdf_data_index.val_idx].T[
                    pdf_data_index.val_idx
                ]
            )

            tmp = PDF_predictions_val - new_val_pseudodata_list

            chi2 = np.einsum("ij,jk,ik->i", tmp, invcovmat_vl, tmp) / tmp.shape[1]
            chi2s_per_replica.append(chi2)
            ret = np.array(chi2s_per_replica)
    else:
        log.warning(f"""Since {fit_name} pseudodata generation has changed,
            hence the overfit metric cannot be determined.""")
        ret = np.array(np.nan)

    return ret


def array_expected_overfitting(
    calculate_chi2s_per_replica,
    replica_data,
    number_of_resamples=1000,
    resampling_fraction=0.95,
):
    """Calculates the expected difference in chi2 between:
    1. The chi2 of a PDF replica calculated using the corresponding pseudodata
        replica used during the fit
    2. The chi2 of a PDF replica calculated using an alternative i.i.d random
        pseudododata replicas

    The expected difference along with an error estimate is obtained through a
    bootstrapping consisting of ``number_of_resamples`` resamples per pdf replica
    where each resampling contains a fraction ``resampling_fraction`` of all
    replicas.

    Parameters
    ----------
    calculate_chi2s_per_replica : np.ndarray
        validation chi2 per pdf replica
    replica_data : list(vp.fitdata.FitInfo)
    number_of_resamples : int, optional
        number of resamples per pdf replica, by default 1000
    resampling_fraction : float, optional
        fraction of replicas used in the bootstrap resampling, by default 0.95

    Returns
    -------
    np.ndarray
        (number_of_resamples*Npdfs,) sized array containing the mean delta chi2
        values per resampled list.
    """
    # calculate_chi2s_per_replica is set to NaN if the pseudodata generation 
    # has changed sinc the fit has been performed. As a result the overfitting
    # metric can no longer be determined.
    if (calculate_chi2s_per_replica != calculate_chi2s_per_replica).all():
        list_expected_overfitting = calculate_chi2s_per_replica
    else:
        fitted_val_erf = np.array([info.validation for info in replica_data])

        number_pdfs = calculate_chi2s_per_replica.shape[0]
        list_expected_overfitting = []
        for _ in range(number_pdfs * number_of_resamples):
            mask = np.random.randint(
                0, number_pdfs, size=int(resampling_fraction * number_pdfs)
            )
            res_tmp = calculate_chi2s_per_replica[mask][:, mask]

            fitted_val_erf_tmp = fitted_val_erf[mask]
            expected_val_chi2 = res_tmp.mean(axis=0)
            delta_chi2 = fitted_val_erf_tmp - expected_val_chi2
            expected_delta_chi2 = delta_chi2.mean()

            list_expected_overfitting.append(expected_delta_chi2)
    return np.array(list_expected_overfitting)


@figure
def plot_overfitting_histogram(fit, array_expected_overfitting):
    """Plots the bootrap error and central value of the overfittedness in a
    historgram"""
    mean = array_expected_overfitting.mean()
    std = array_expected_overfitting.std()

    f, ax = plt.subplots(1, 1)

    ax.hist(array_expected_overfitting, bins=50, density=True)
    ax.axvline(x=mean, color="black")
    ax.axvline(x=0, color="black", linestyle="--")
    xrange = [
        array_expected_overfitting.min(),
        array_expected_overfitting.max(),
    ]
    xgrid = np.linspace(xrange[0], xrange[1], num=100)
    ax.plot(xgrid, stats.norm.pdf(xgrid, mean, std))
    ax.set_xlabel(r"$\mathcal{R}_O$")
    ax.set_ylabel("density")
    ax.set_title(f"{fit.label}")
    plt.tight_layout()
    return f


fits_overfitting_summary = collect(
    "fit_overfitting_summary", ("fits", "fitcontext")
)


@table
def fit_overfitting_summary(fit, array_expected_overfitting):
    """Creates a table containing the overfitting information:
    - mean chi2 difference
    - bootstrap error
    - sigmas away from 0
    """
    mean = array_expected_overfitting.mean()
    std = array_expected_overfitting.std()
    return pd.DataFrame(
        [mean, std, mean / std],
        columns=[fit.label],
        index=["mean", "bootstrap error", "sigmas away from 0"],
    )


@table
def summarise_overfitting(fits_overfitting_summary):
    """Same as `fit_overfitting_summary`, but collected over all `fits` in the
    runcard and put in a single table.
    """
    return pd.concat(fits_overfitting_summary, axis=1)
