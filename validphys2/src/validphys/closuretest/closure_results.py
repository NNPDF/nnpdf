"""
closuretest/closure_results.py

Module containing actiosn to calculate sigle closure test estimators.
This is useful for quickly checking the bias of a fit without having 
to run the full multiclosure analysis.

"""

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.table import table
from validphys.closuretest.closure_checks import (
    check_fit_isclosure,
    check_fits_areclosures,
    check_fits_same_filterseed,
    check_fits_underlying_law_match,
    check_use_fitcommondata,
)


fits_experiments = collect("experiments_data", ("fits", "fitcontext"))


experiments_bootstrap_chi2_central = collect(
    "dataset_inputs_bootstrap_chi2_central", ("group_dataset_inputs_by_experiment",)
)


fits_exps_bootstrap_chi2_central = collect(
    "experiments_bootstrap_chi2_central", ("fits", "fitcontext")
)


fits_level_1_noise = collect("total_chi2_data", ("fits", "fitinputcontext", "fitunderlyinglaw"))


@check_use_fitcommondata
@check_fits_areclosures
@check_fits_same_filterseed
@check_fits_underlying_law_match
def delta_chi2_bootstrap(
    fits_level_1_noise, fits_exps_bootstrap_chi2_central, fits, use_fitcommondata
):
    """Bootstraps delta chi2 for specified fits.
    Delta chi2 measures whether the level one data is fitted better by
    the underlying law or the specified fit, it is a measure of overfitting.

    delta chi2 = (chi2(T[<f>], D_1) - chi2(T[f_in], D_1))/chi2(T[f_in], D_1)

    where T[<f>] is central theory prediction from fit, T[f_in] is theory
    prediction from t0 pdf (input) and D_1 is level 1 closure data

    Exact details on delta chi2 can be found in 1410.8849 eq (28).
    """
    closure_total_chi2_boot = np.sum(fits_exps_bootstrap_chi2_central, axis=1)
    t0_pseudodata_chi2 = np.array([chi2.central_result for chi2 in fits_level_1_noise])
    deltachi2boot = (
        closure_total_chi2_boot - t0_pseudodata_chi2[:, np.newaxis]
    ) / t0_pseudodata_chi2[:, np.newaxis]
    return deltachi2boot


# Note that these collect over the experiments as specified in fit in case of
# TEST set
fits_exps_level_1_noise = collect(
    "experiments_chi2_data", ("fits", "fitinputcontext", "fitunderlyinglaw")
)
fits_exps_chi2 = collect("experiments_chi2_data", ("fits", "fitcontext"))


@table
@check_use_fitcommondata
@check_fits_areclosures
@check_fits_same_filterseed
@check_fits_underlying_law_match
def delta_chi2_table(
    fits_exps_chi2,
    fits_exps_level_1_noise,
    fits_name_with_covmat_label,
    fits_experiments,
    fits,
    use_fitcommondata,
):
    """Calculated delta chi2 per experiment and put in table
    Here delta chi2 is just normalised by ndata and is equal to

    delta_chi2 = (chi2(T[<f>], D_1) - chi2(T[f_in], D_1))/ndata
    """
    dfs = []
    cols = ("ndata", r"$\Delta_{chi^2}$ (normalised by ndata)")
    for label, experiments, exps_chi2, exps_level_1_noise in zip(
        fits_name_with_covmat_label, fits_experiments, fits_exps_chi2, fits_exps_level_1_noise
    ):
        records = []
        for experiment, exp_chi2, level_1_noise in zip(experiments, exps_chi2, exps_level_1_noise):
            delta_chi2 = (exp_chi2.central_result - level_1_noise.central_result) / exp_chi2.ndata
            npoints = exp_chi2.ndata
            records.append(dict(experiment=str(experiment), npoints=npoints, delta_chi2=delta_chi2))
        df = pd.DataFrame.from_records(
            records, columns=("experiment", "npoints", "delta_chi2"), index=("experiment",)
        )
        df.columns = pd.MultiIndex.from_product(([label], cols))
        dfs.append(df)
    res = pd.concat(dfs, axis=1)
    return res


@check_fit_isclosure
def fit_underlying_pdfs_summary(fit, fitunderlyinglaw):
    """Returns a table with a single column for the `fit` with a row indication
    the PDF used to generate the data and the t0 pdf
    """
    t0name = fit.as_input()["datacuts"]["t0pdfset"]
    df = pd.DataFrame(
        [str(fitunderlyinglaw["pdf"]), t0name],
        columns=[fit.label],
        index=["underlying PDF", "t0 PDF"],
    )
    return df


fits_underlying_pdfs_summary = collect("fit_underlying_pdfs_summary", ("fits",))


@table
def summarise_closure_underlying_pdfs(fits_underlying_pdfs_summary):
    """Collects the underlying pdfs for all fits and concatenates them into a single table"""
    return pd.concat(fits_underlying_pdfs_summary, axis=1)
