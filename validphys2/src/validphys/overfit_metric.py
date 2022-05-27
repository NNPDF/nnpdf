"""
overfit_metric.py

This calculates the overfit metric
"""

import logging


from reportengine import collect
from reportengine.figure import figure
from reportengine.table import table

import scipy.stats as stats
import pandas as pd
import numpy as np


log = logging.getLogger(__name__)

preds = collect("predictions", ("pdfs","dataset_inputs",))

def foo(get_val_erf, calculate_chi2s_per_replica):
    import ipdb; ipdb.set_trace()

def array_expected_delta_chi2(calculate_chi2s_per_replica, replica_data):

    res_over = calculate_chi2s_per_replica

    fitted_val_erf = np.array([info.validation for info in replica_data])

    number_pdfs = res_over.shape[0]
    list_expected_delta_chi2 = []
    for i in range(number_pdfs*1000):
        mask = np.random.randint(0, number_pdfs, size=int(0.95*number_pdfs))
        res_tmp = res_over[mask][:,mask]

        fitted_val_erf_tmp = fitted_val_erf[mask]
        expected_val_chi2 = res_tmp.mean(axis=0)
        delta_chi2 = fitted_val_erf_tmp - expected_val_chi2
        expected_delta_chi2 = delta_chi2.mean()

        list_expected_delta_chi2.append(expected_delta_chi2)
    return np.array(list_expected_delta_chi2)


def _create_new_val_pseudodata(pdf_data_index, fit_data_indices_list):
    vl_data_fitrep = []
    for fitreplica_info in fit_data_indices_list:
        vl_data_fitrep.append(fitreplica_info.pseudodata.loc[pdf_data_index.val_idx])
    return np.array(vl_data_fitrep).squeeze()


def calculate_chi2s_per_replica(
    recreate_pdf_pseudodata_no_table,
    preds,
    dataset_inputs,
    groups_covmat_no_table,
    ):

    groups_covmat = groups_covmat_no_table

    pp = []
    for i, dss in enumerate(dataset_inputs):
        preds_witout_cv = preds[i].drop(0, axis=1)
        df = pd.concat({dss.name: preds_witout_cv}, names=['dataset'])
        pp.append(df)

    PDF_predictions = pd.concat(pp)

    chi2s_per_replica=[]
    for enum, pdf_data_index in enumerate(recreate_pdf_pseudodata_no_table):

        prediction_filter=pdf_data_index.val_idx.droplevel(level=0)
        prediction_filter.rename(["dataset","data"], inplace=True)
        PDF_predictions_val = PDF_predictions.loc[prediction_filter]
        PDF_predictions_val = PDF_predictions_val.values[:,enum]

        new_val_pseudodata_list = _create_new_val_pseudodata(pdf_data_index, recreate_pdf_pseudodata_no_table)

        invcovmat_vl = np.linalg.inv(groups_covmat[pdf_data_index.val_idx].T[pdf_data_index.val_idx])
        
        tmp = ( PDF_predictions_val - new_val_pseudodata_list)

        chi2 = np.einsum("ij,jk,ik->i", tmp, invcovmat_vl, tmp) / tmp.shape[1]
        chi2s_per_replica.append(chi2)

    # array of chi2 per replica
    return np.array(chi2s_per_replica)

import matplotlib.pyplot as plt

@figure
def plot_deltachi2_histogram(array_expected_delta_chi2):
    mean = array_expected_delta_chi2.mean()
    std = array_expected_delta_chi2.std()

    f, ax = plt.subplots(1,1)

    ax.hist(array_expected_delta_chi2, bins=50, density=True);
    ax.axvline(x=mean,color='black')
    xrange=[array_expected_delta_chi2.min(),array_expected_delta_chi2.max()]
    xgrid = np.linspace(xrange[0], xrange[1], num=100)
    ax.set_xlim(xrange[0],xrange[1])
    ax.plot(xgrid, stats.norm.pdf(xgrid, mean, std))
    ax.set_xlabel(r"$\Delta \chi^2_{\mathrm{overfit}}$")
    ax.set_ylabel("density")
    ax.set_title(f"fitname")
    plt.tight_layout()
    return f

fits_deltachi2_summary = collect('fit_deltachi2_summary', ('fits','fitcontext'))

def fit_deltachi2_summary(fit,array_expected_delta_chi2):
    mean = array_expected_delta_chi2.mean()
    std = array_expected_delta_chi2.std()
    return pd.DataFrame(
        [mean, std, abs(mean/std)],
        columns=[fit.label],
        index=["mean", "bootstrap error", "abs(mean/bootsrap error)"]
        )

@table
def summarise_deltachi2(fits_deltachi2_summary):
    return pd.concat(fits_deltachi2_summary,axis=1)