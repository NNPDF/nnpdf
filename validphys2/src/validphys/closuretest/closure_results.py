"""
closuretest/results.py

underlying actions to calculate closure test estimators plus some table actions
"""
from collections import namedtuple

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.table import table
from validphys.calcutils import bootstrap_values, calc_chi2
from validphys.checks import check_pdf_is_montecarlo
from validphys.closuretest.closure_checks import (
    check_fit_isclosure,
    check_fits_areclosures,
    check_fits_same_filterseed,
    check_fits_underlying_law_match,
    check_use_fitcommondata,
)

BiasData = namedtuple("BiasData", ("bias", "ndata"))

underlying_results = collect("results", ("fitunderlyinglaw",))


@check_fit_isclosure
@check_use_fitcommondata
def bias_dataset(results, underlying_results, fit, use_fitcommondata):
    """Calculate the bias for a given dataset and fit. The bias is defined as
    chi2 between the prediction from the underlying PDF (which was used to
    generate the closure pseudodata), also known as level zero closure data, and
    the central prediction calculated from the fitted PDF.

    we require that use_fitcommondata is true because the generated closure data
    is used to generate the multiplicative contributions to the covariance
    matrix
    """
    dt_ct, th_ct = results
    # does collect need to collect a list even with one element?
    ((_, th_ul),) = underlying_results
    central_diff = th_ct.central_value - th_ul.central_value
    bias_out = calc_chi2(dt_ct.sqrtcovmat, central_diff)  # unnormalised
    return BiasData(bias_out, len(th_ct))


underlying_dataset_inputs_results = collect("dataset_inputs_results", ("fitunderlyinglaw",))


@check_fit_isclosure
@check_use_fitcommondata
def bias_experiment(
    dataset_inputs_results, underlying_dataset_inputs_results, fit, use_fitcommondata
):
    """Like `bias_dataset` but for a whole experiment."""
    return bias_dataset(
        dataset_inputs_results, underlying_dataset_inputs_results, fit, use_fitcommondata
    )


experiments_bias = collect("bias_experiment", ("group_dataset_inputs_by_experiment",))
fits_experiments_bias = collect("experiments_bias", ("fits", "fitcontext"))
fits_experiments = collect("experiments_data", ("fits", "fitcontext"))


@table
@check_fits_same_filterseed
@check_fits_underlying_law_match
def biases_table(fits_experiments, fits_experiments_bias, fits, show_total: bool = False):
    """Creates a table with fits as the columns and the experiments from both
    fits as the row index.
    """
    col = ["bias"]
    dfs = []
    for fit, experiments, biases in zip(fits, fits_experiments, fits_experiments_bias):
        total = 0
        total_points = 0
        records = []
        for biasres, experiment in zip(biases, experiments):
            records.append(
                dict(
                    experiment=str(experiment), bias=biasres.bias / biasres.ndata  # normalised bias
                )
            )
            if show_total:
                total += biasres.bias
                total_points += biasres.ndata
        if show_total:
            total /= total_points
            records.append(dict(experiment="Total", bias=total))

        df = pd.DataFrame.from_records(
            records, columns=("experiment", "bias"), index=("experiment")
        )
        df.columns = pd.MultiIndex.from_product(([str(fit)], col))
        dfs.append(df)
    return pd.concat(dfs, axis=1, sort=True)


@check_pdf_is_montecarlo
@check_fit_isclosure
def bootstrap_bias_experiment(
    dataset_inputs_results, underlying_dataset_inputs_results, bootstrap_samples=500
):
    """Calculates bias as per `bias_experiment` but performs bootstrap sample
    across replicas. note that bias_experiment returns a named tuple like
    (unnormalised_bias, ndata) whereas this actions simply returns an array
    `boostrap_bias` with length bootstrap_samples. Each element of
    returned array is bias/n_data (bias normalised by number of datapoints)
    """
    dt_ct, th_ct = dataset_inputs_results
    ((_, th_ul),) = underlying_dataset_inputs_results
    th_ct_boot_cv = bootstrap_values(th_ct.error_members, bootstrap_samples)
    boot_diffs = th_ct_boot_cv - th_ul.central_value[:, np.newaxis]
    boot_bias = calc_chi2(dt_ct.sqrtcovmat, boot_diffs) / len(dt_ct)
    return boot_bias


experiments_bootstrap_bias = collect(
    "bootstrap_bias_experiment", ("group_dataset_inputs_by_experiment",)
)
fits_experiments_bootstrap_bias = collect("experiments_bootstrap_bias", ("fits", "fitcontext"))


@table
@check_use_fitcommondata
@check_fits_areclosures
@check_fits_same_filterseed
@check_fits_underlying_law_match
def fits_bootstrap_bias_table(
    fits_experiments_bootstrap_bias,
    fits_name_with_covmat_label,
    fits_experiments,
    fits,
    use_fitcommondata,
):
    """Produce a table with bias for each experiment for each fit, along with
    variance calculated from doing a bootstrap sample
    """
    col = ["bias", "bias std. dev."]
    dfs = []
    for fit, experiments, biases in zip(
        fits_name_with_covmat_label, fits_experiments, fits_experiments_bootstrap_bias
    ):
        records = []
        for biasres, experiment in zip(biases, experiments):
            records.append(
                dict(experiment=str(experiment), bias=biasres.mean(), bias_std=biasres.std())
            )
        df = pd.DataFrame.from_records(
            records, columns=("experiment", "bias", "bias_std"), index=("experiment")
        )
        df.columns = pd.MultiIndex.from_product(([str(fit)], col))
        dfs.append(df)
    return pd.concat(dfs, axis=1, sort=True)


VarianceData = namedtuple("VarianceData", ("variance", "ndata"))


@check_fit_isclosure
@check_use_fitcommondata
def variance_dataset(results, fit, use_fitcommondata):
    """calculate the variance for a given dataset, which is the spread of
    replicas measured in the space of the covariance matrix. Given by:

        E_rep [ (T - E_rep[T])_i C^{-1}_ij (T - E_rep[T])_j ]

    where E_rep is the expectation value across replicas. The quantity is the
    same as squaring `phi_data`, however it is redefined here in a way which can
    be made fully independent of the closure data. This is useful when checking
    the variance of data which was not included in the fit.

    """
    dt, th = results
    diff = th.central_value[:, np.newaxis] - th.error_members
    var_unnorm = calc_chi2(dt.sqrtcovmat, diff).mean()
    return VarianceData(var_unnorm, len(th))


@check_fit_isclosure
@check_use_fitcommondata
def variance_experiment(dataset_inputs_results, fit, use_fitcommondata):
    """Like variance_dataset but for a whole experiment"""
    return variance_dataset(dataset_inputs_results, fit, use_fitcommondata)


@check_fit_isclosure
def bootstrap_variance_experiment(dataset_inputs_results, bootstrap_samples=500):
    """Calculate the variance as in `variance_experiment` but performs bootstrap
    sample of the estimator. Returns an array of variance for each resample,
    normalised to the number of data in the experiment.
    """
    dt_ct, th_ct = dataset_inputs_results
    diff = th_ct.central_value[:, np.newaxis] - th_ct.error_members
    var_unnorm_boot = bootstrap_values(
        diff, bootstrap_samples, apply_func=(lambda x, y: calc_chi2(y, x)), args=[dt_ct.sqrtcovmat]
    ).mean(
        axis=0
    )  # mean across replicas
    return var_unnorm_boot / len(th_ct)  # normalise by n_data


experiments_boostrap_variance = collect(
    "bootstrap_variance_experiment", ("group_dataset_inputs_by_experiment",)
)

fits_exps_bootstrap_var = collect("experiments_boostrap_variance", ("fits", "fitcontext"))


@table
@check_use_fitcommondata
@check_fits_areclosures
@check_fits_same_filterseed
@check_fits_underlying_law_match
def fits_bootstrap_variance_table(
    fits_exps_bootstrap_var, fits_name_with_covmat_label, fits_experiments, fits, use_fitcommondata
):
    """Produce a table with variance and its error. Variance is defined as

        var = sum_ij E_rep[(E_rep[T_i], T_i) invcov_ij (E_rep[T_j], T_j)] / N_data

    which is the expectation value across replicas of the chi2 between the central
    theory predictions and replica theory predictions. It is the same as phi^2
    and gives the variance of the theory predictions in units of the covariance
    matrix, normalised by the number of data points.

    The error is the standard deviation across bootstrap samples.

    """
    col = ["variance", "variance std. dev."]
    dfs = []
    for fit, experiments, fit_vars in zip(
        fits_name_with_covmat_label, fits_experiments, fits_exps_bootstrap_var
    ):
        records = []
        for var_res, experiment in zip(fit_vars, experiments):
            records.append(
                dict(
                    experiment=str(experiment), variance=var_res.mean(), variance_std=var_res.std()
                )
            )
        df = pd.DataFrame.from_records(
            records, columns=("experiment", "variance", "variance_std"), index=("experiment")
        )
        df.columns = pd.MultiIndex.from_product(([str(fit)], col))
        dfs.append(df)
    return pd.concat(dfs, axis=1, sort=True)


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
