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


# NOTE: for some reason the fit doesn't get properly resolved if you try to
# collect data over fits
fits_dataset = collect("dataset", ("fits",))

@check_use_fitcommondata
def fits_dataset_cvs(fits_dataset):
    """Internal function for loading the level one data for all fits
    for a single dataset. This function avoids using the c++ loading of
    commondata which is very slow and also avoids the stringent metadata
    checks of the newer python commondata parser.

    """
    fits_cv = []
    for ds in fits_dataset:
        # using the official loader is really slow, open the CSV
        # and then cut the central values manually.
        # TODO: Save central values in nice table like pseudodata
        # but this should be done beyond NNPDF4.0
        cd_df = pd.read_csv(ds.commondata.datafile, sep=r'\s+', skiprows=1, header=None)
        # based on columns from python cd reader:
        # ['entry', 'process', 'kin1', 'kin2', 'kin3', 'data', 'stat']
        fits_cv.append(cd_df.iloc[cut_mask(ds.cuts), 5].to_numpy())
    return fits_cv

data_fits_cv = collect(fits_dataset_cvs, ("data",))

def expected_data_delta_chi2(
    data_fits_cv,
    internal_multiclosure_data_loader
):
    """For ``data``, calculate the mean of delta chi2 across all fits, returns
    a tuple of number of data points and unnormalised delta chi2.
    """
    closures_th, law_th, _, sqrt_covmat = internal_multiclosure_data_loader
    law_central = law_th.central_value
    fits_delta_chi2 = []
    for i_fit, fit_th in enumerate(closures_th):
        # transpose the datasets fits cvs into the cvs for all datasets for single fit
        dt_central = np.concatenate([fits_cvs[i_fit] for fits_cvs in data_fits_cv])
        th_replicas = fit_th._rawdata
        th_central = np.mean(th_replicas, axis=-1)
        shift = calc_chi2(sqrt_covmat, law_central - dt_central)
        chi2_cent = calc_chi2(sqrt_covmat, th_central - dt_central)
        fits_delta_chi2.append(chi2_cent - shift)
    ndata = len(law_central)
    return ndata, np.mean(fits_delta_chi2)


exps_expected_delta_chi2 = collect(
    "expected_data_delta_chi2", ("group_dataset_inputs_by_experiment",))


def total_expected_data_delta_chi2(exps_expected_delta_chi2):
    """Takes :py:func:`expected_data_delta_chi2` evaluated for each experiment
    and then sums across experiments. Returns the total number of datapoints
    and unnormalised delta chi2.
    """
    ndata, delta_chi2 = np.sum(exps_expected_delta_chi2, axis=0)
    return ndata, delta_chi2


groups_expected_delta_chi2 = collect(
    "expected_data_delta_chi2", ("group_dataset_inputs_by_metadata",))

@table
def expected_delta_chi2_table(
    groups_expected_delta_chi2,
    group_dataset_inputs_by_metadata,
    total_expected_data_delta_chi2,
):
    """Tabulate the expectation value of delta chi2 across fits for groups
    with an additional row with the total across all data at the bottom.
    """
    records = []
    for group, delta_chi2_res in zip(
        group_dataset_inputs_by_metadata,
        groups_expected_delta_chi2,
    ):
        name = group["group_name"]
        ndata, delta_chi2 = delta_chi2_res

        records.append(
            dict(
                group=name,
                ndata=ndata,
                delta_chi2=delta_chi2 / ndata,
            )
        )
    ndata, delta_chi2 = total_expected_data_delta_chi2
    records.append(
        dict(
            group="Total",
            ndata=ndata,
            delta_chi2=delta_chi2 / ndata,
        )
    )
    df = pd.DataFrame.from_records(records, index="group")
    df.columns = [
        "ndata",
        r"$\Delta_{\chi^2}$",
    ]
    return df
