"""
multiclosure_pseudodata

actions which load fit pseudodata and compute actions related to overfitting.
Estimators here can only be calculated on data used in the fit.

"""
import numpy as np
import pandas as pd
from reportengine import collect
from reportengine.table import table

from validphys.calcutils import calc_chi2
from validphys.closuretest.closure_checks import check_use_fitcommondata
from validphys.core import cut_mask
from validphys.tableloader import sane_load


# collect datasets across fits, with use_fitcommondata=True
fits_exps_data_inputs = collect("experiments_data", ("fits",))
fits_replica_indexes = collect("fitted_replica_indexes", ("fits", "fitpdf"))

PSEUDODATA_TABLE_NAME = "datacuts_theory_closuretest_fitting_pseudodata_table.csv"

@check_use_fitcommondata
#check savepseudodata
#check fitted data? This is implictly checked by the loader.
def internal_multiclosure_pseudodata_loader(
    fits_exps_data_inputs, fits, fits_replica_indexes, experiments_index):
    """Load the pseudodata and central values for each fit. The pseudodata
    is loaded from CSV for each of the fitted replicas and then concatenated
    for each fit. The central values are loaded and then concatenated
    across fits.

    Returns
    -------
    fits_cvs: pd.DataFrame
        indexed by :py:func`experiments_index` each column corresponds to
        a different ``fit``
    fits_replicas_pseudodata: list[pd.DataFrame]
        list of dataframes, each one indexed by :py:func`experiments_index`
        with each column corresponding to a different replica pseudodata. The
        replicas are ordered in the same order as the PDF replicas.

    """
    fits_pseudodata = []
    fits_central_values = []
    for fit_exps_data, fit, replica_indexes in zip(fits_exps_data_inputs, fits, fits_replica_indexes):
        reps_pseudodata = []
        # load replicas from CSV
        for rep in replica_indexes:
            tab_path = fit.path / "nnfit" / f"replica_{rep}" / PSEUDODATA_TABLE_NAME
            # contains pseudodata for all datasets in fit.
            df = sane_load(tab_path, index_col=[0, 1, 2], header=0)
            # reindex by experiments_index in case report is on subset of data.
            reps_pseudodata.append(df.reindex(experiments_index))
        # keep replicas in correct order.
        fits_pseudodata.append(pd.concat(reps_pseudodata, axis=1, sort=False))

        # load commondata using python cd loader
        dsets_cv = []
        for exp in fit_exps_data:
            for ds in exp.datasets:
                # using the official loader is really slow, open the CSV
                # and then cut the central values manually.
                # TODO: Save central values in nice table like pseudodata
                # but this should be done beyond NNPDF4.0
                cd_df = pd.read_csv(ds.commondata.datafile, sep=r'\s+', skiprows=1, header=None)
                # based on columns from python cd reader:
                # ['entry', 'process', 'kin1', 'kin2', 'kin3', 'data', 'stat']
                dsets_cv.append(cd_df.iloc[cut_mask(ds.cuts), 5].to_numpy())

        fits_central_values.append(
            pd.DataFrame(
                np.concatenate(dsets_cv),
                index=experiments_index,
                columns=[str(fit)]
            )
        )
    return pd.concat(fits_central_values, axis=1, sort=False), fits_pseudodata


# align theory predictions with experiment grouping.
experiments_internal_multiclosure_loader = collect(
    "internal_multiclosure_data_loader", ("group_dataset_inputs_by_experiment",))


@table
def experiments_closure_pseudodata_estimators_table(
    group_dataset_inputs_by_experiment,
    experiments_internal_multiclosure_loader,
    internal_multiclosure_pseudodata_loader
):
    """Tabulate for each experiment the pseudodata relates estimators for
    multiple closure fits. The estimators are shift, noise, delta replica chi2,
    delta_chi2 and delta_epsilon. The first two are the shift and noise applied
    the the underlying law to generate the pseudodata. delta replica chi2 is
    the replica chi2 - (shift + noise) and tells us if the replica predictions
    overfit the pseudodata replicas if delta replica chi2 < 0, or underfit
    the data if delta replica chi2 > 0. Then the last
    two estimators break down delta replica chi2 into the component which is
    how much the average replica overfits the level 1 data and how much each
    replicas overfits the level 2 noise, which in theory can be independent.

    The expectation value is taken over fits and replicas for each of the
    estimators.

    """
    fits_cv, fits_reps_pseudo = internal_multiclosure_pseudodata_loader
    records = []
    total_ndata = 0
    total_erep_delta_chi2 = 0
    total_delta_chi2 = 0
    total_shift = 0
    total_noise = 0
    for (exp, exp_internal_loader_tuple) in zip(
        group_dataset_inputs_by_experiment, experiments_internal_multiclosure_loader
    ):
        exp_name = exp["group_name"]
        closures_th, law_th, _, sqrt_covmat = exp_internal_loader_tuple
        law_central = law_th.central_value

        fits_shift = []
        fits_noise = []
        fits_delta_chi2 = []
        fits_erep_delta_chi2 = []
        for i_fit, fit_th in enumerate(closures_th):
            # some of these could be done outside of loop, but easier to do here.
            th_replicas = fit_th.error_members
            th_central = np.mean(th_replicas, axis=-1)
            dt_replicas = fits_reps_pseudo[i_fit].xs(exp_name, axis=0, level=0).to_numpy()
            dt_central = fits_cv.xs(exp_name, axis=0, level=0).iloc[:, i_fit].to_numpy()

            shift = calc_chi2(sqrt_covmat, law_central - dt_central)
            e_rep_noise = np.mean(calc_chi2(sqrt_covmat, dt_central[:, np.newaxis] - dt_replicas))
            chi2_cent = calc_chi2(sqrt_covmat, th_central - dt_central)
            e_rep_chi2_rep = np.mean(calc_chi2(sqrt_covmat, th_replicas - dt_replicas))
            fits_erep_delta_chi2.append(e_rep_chi2_rep - (e_rep_noise + shift))
            fits_delta_chi2.append(chi2_cent - shift)
            fits_shift.append(shift)
            fits_noise.append(e_rep_noise)

        total_ndata += len(law_central)
        total_erep_delta_chi2 += np.mean(fits_erep_delta_chi2)
        total_delta_chi2 += np.mean(fits_delta_chi2)
        total_shift += np.mean(fits_shift)
        total_noise += np.mean(fits_noise)


        records.append(
            dict(
                experiment=exp_name,
                ndata=len(law_central),
                shift=np.mean(fits_shift),
                noise=np.mean(fits_noise),
                e_rep_delta_chi2=np.mean(fits_erep_delta_chi2),
                delta_chi2=np.mean(fits_delta_chi2),
                delta_eps=np.mean(fits_erep_delta_chi2) - np.mean(fits_delta_chi2)
            )
        )
    records.append(
        dict(
            experiment="Total",
            ndata=total_ndata,
            shift=total_shift,
            noise=total_noise,
            e_rep_delta_chi2=total_erep_delta_chi2,
            delta_chi2=total_delta_chi2,
            delta_eps=total_erep_delta_chi2 - total_delta_chi2
        )
    )
    df = pd.DataFrame.from_records(records, index="experiment")
    # normalise by ndata
    df.iloc[:, 1:] = df.iloc[:, 1:] / df.iloc[:, [0]].to_numpy()
    df.columns = [
        "ndata",
        "shift",
        "noise",
        r"$\Delta_{{\chi^2}^{(k)}}$",
        r"$\Delta_{\chi^2}$",
        r"$\Delta_{\epsilon}$"
    ]
    return df

@table
def compare_delta_chi2_bias_variance_table(
    experiments_bias_variance_table,
    experiments_closure_pseudodata_estimators_table,
):
    """Convenience function, which joins relevant columns from
    ``experiments_bias_variance_table`` and
    ``experiments_closure_pseudodata_estimators_table`` to generate a table
    which has bias, delta_chi2, variance and delta_epsilon.

    Bias and variance
    contextualise delta_chi2 and delta_epsilon. An alternative form of
    delta_chi2 is

        delta_chi2 = bias - shift cross term,

    where the shift cross term is the covariance between the level one shift
    and the difference between central prediction and underlying law (in units
    of covariance). By comparing delta_chi2 and bias you can get an idea of
    how correlated the shift and the difference between the central prediction
    and underlying law are.

    """
    bias_var_tab = experiments_bias_variance_table.drop("ndata", axis=1)
    total_df = pd.concat(
        [
            bias_var_tab,
            experiments_closure_pseudodata_estimators_table,
        ],
        axis=1
    )
    return total_df[
        ["ndata", "bias", r"$\Delta_{\chi^2}$", "variance", r"$\Delta_{\epsilon}$"]
    ]
