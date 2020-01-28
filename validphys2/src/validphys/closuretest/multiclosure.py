"""
closuretest/multiclosure.py

containing all of the statistical estimators which are
averaged across multiple fits or a single replica proxy fit
"""

import numpy as np
import scipy.linalg as la
from scipy.special import erf
import pandas as pd
import matplotlib.pyplot as plt

from reportengine import collect
from reportengine.table import table
from reportengine.figure import figuregen, figure
from reportengine.checks import make_argcheck

from validphys.pdfgrids import xplotting_grid, check_basis
from validphys.pdfbases import Basis
from validphys.core import PDF
from validphys.checks import check_scale

fits_bias_experiment = collect("bias_experiment", ("fits", "fitcontext"))

def expected_bias(fits_bias_experiment):
    bias_centrals = [
        biasdata.bias*biasdata.ndata for biasdata in fits_bias_experiment]
    #ndata = fits_bias_experiment[0].ndata # ndata should be same for all fits
    bias_mean = np.mean(bias_centrals)
    return bias_mean

exps_expected_bias = collect("expected_bias", ("experiments",))

fits_phi_experiment = collect("phi_data_experiment", ("fits", "fitcontext"))

fits_variance_experiment = collect("variance_experiment", ("fits", "fitcontext"))

def expected_variance(fits_variance_experiment):
    var_mean = np.mean(fits_variance_experiment)
    return var_mean

exps_expected_var = collect("expected_variance", ("experiments",))

@table
def experiments_bias_variance_ratio(
    exps_expected_bias, exps_expected_var, experiments,):
    records = []
    bias_tot = 0
    var_tot = 0
    for exp, bias, var in zip(experiments, exps_expected_bias, exps_expected_var):
        records.append(dict(
            experiment=str(exp),
            ratio=bias/var
        ))
        bias_tot += bias
        var_tot += var
    records.append(dict(
            experiment="total",
            ratio=bias_tot/var_tot
    ))
    df = pd.DataFrame.from_records(
        records,
        index="experiment",
        columns=("experiment", "ratio")
    )
    df.columns = ["bias/variance"]
    return df

@table
def sqrt_experiments_bias_variance_ratio(experiments_bias_variance_ratio):
    df_in = experiments_bias_variance_ratio
    return pd.DataFrame(np.sqrt(df_in.values), index=df_in.index, columns=["sqrt(bias/variance)"])

@table
def expected_xi_from_bias_variance(sqrt_experiments_bias_variance_ratio):
    df_in = sqrt_experiments_bias_variance_ratio
    n_sigma_in_variance = 1/df_in.values
    res = erf(n_sigma_in_variance/np.sqrt(2))
    return pd.DataFrame(res, index=df_in.index, columns=[r"$\xi$ from ratio"])

@make_argcheck(check_basis)
def histogram_pdfgrid(
    pdf:PDF, Q:(float,int), basis:(str, Basis)='flavour',
    flavours:(list, tuple, type(None))=None):
    xgrid = np.array([0.05, 0.1, 0.2]) # internal grid, defined in 3.0 paper
    return xplotting_grid(
        pdf, Q, xgrid=xgrid, basis=basis, flavours=flavours)

@make_argcheck(check_basis)
def xi_pdfgrid(
    pdf:PDF, Q:(float,int), basis:(str, Basis)='flavour',
    flavours:(list, tuple, type(None))=None):
    xgrid_log = np.logspace(-5, -1, num=10, endpoint=False)
    xgrid_lin = np.linspace(0.1, 1, num=10, endpoint=False)
    xgrid = np.concatenate((xgrid_log, xgrid_lin), axis=0)
    return xplotting_grid(
        pdf, Q, xgrid=xgrid, basis=basis, flavours=flavours)

fits_xi_pdfgrid = collect("pdf_bias_grid", ("fits", "fitpdf"))
underlying_xi_pdfgrid = collect("pdf_bias_grid", ("fitunderlyinglaw",))

@table
def xi_1sigma(fits_xi_pdfgrid, underlying_xi_pdfgrid):
    fit_replica_grids = np.array([xpg.grid_values for xpg in fits_xi_pdfgrid])
    fit_central_grids = fit_replica_grids.mean(axis=1, keepdims=True)
    law_grid = underlying_xi_pdfgrid[0].grid_values[0][np.newaxis, np.newaxis, ...]
    #print(law_grid.shape, fit_replica_grids.shape, fit_central_grids.shape)
    central_diff = law_grid - fit_central_grids
    sigma_diff = fit_replica_grids - fit_central_grids
    variance = ((sigma_diff)**2).mean(axis=1, keepdims=True)#.mean(axis=0)
    sigma = np.sqrt(variance)
    indicator = np.array(abs(central_diff) < sigma, dtype=float) #fit, 1, flav, x
    flavs = pd.Index([
        underlying_xi_pdfgrid[0].basis.elementlabel(fl)
        for fl in underlying_xi_pdfgrid[0].flavours],
        name="flavour"
    )
    xgridindex = pd.Index(underlying_xi_pdfgrid[0].xgrid, name="x")
    df = pd.DataFrame(
        indicator.mean(axis=0)[0, ...],
        index=flavs,
        columns=xgridindex,
    )
    return df

@figuregen
@check_scale('xscale', allow_none=True)
def plot_xi_1sigma(xi_1sigma, Q, xscale):
    x = xi_1sigma.columns.values
    xi_data = xi_1sigma.values
    for i, fl in enumerate(xi_1sigma.index.values):
        fig, ax = plt.subplots()
        ax.plot(x, xi_data[i, :], '*', label=r"$\xi_{1\sigma}$ = "+f"{xi_data[i, :].mean()}")
        ax.axhline(0.68, linestyle=":", color="k", label="expected value")
        ax.set_ylim([0, 1])
        ax.set_title(r"$\xi_{1\sigma}$"+f" for Q={Q}, ${fl}$ PDF.")
        ax.set_xlabel("x")
        ax.set_ylabel(r"$\xi_{1\sigma}$")
        ax.set_xscale(xscale)
        ax.legend()
        yield fig

def xi_pdfgrid_proxy(
    proxy_pdf:PDF, Q:(float,int), basis:(str, Basis)='flavour',
    flavours:(list, tuple, type(None))=None):
    return xi_pdfgrid(proxy_pdf, Q, basis, flavours)

def xi_pdfgrid_replicas(
    replica_pdf:PDF, Q:(float,int), basis:(str, Basis)='flavour',
    flavours:(list, tuple, type(None))=None):
    return xi_pdfgrid(replica_pdf, Q, basis, flavours)

@table
def xi_1sigma_proxy(xi_pdfgrid_replicas, underlying_xi_pdfgrid, xi_pdfgrid_proxy):
    fit_replica_grids = xi_pdfgrid_replicas.grid_values
    fit_central_grids = fit_replica_grids.mean(axis=0, keepdims=True)

    law_grid = underlying_xi_pdfgrid[0].grid_values[0][np.newaxis, ...]
    proxy_grid = xi_pdfgrid_proxy.grid_values

    central_diff = law_grid - proxy_grid
    sigma_diff = fit_replica_grids - fit_central_grids

    variance = ((sigma_diff)**2).mean(axis=0, keepdims=True)
    sigma = np.sqrt(variance)

    indicator = np.array(abs(central_diff) < sigma, dtype=float) #fit, flav, x
    flavs = pd.Index([
        underlying_xi_pdfgrid[0].basis.elementlabel(fl)
        for fl in underlying_xi_pdfgrid[0].flavours],
        name="flavour"
    )
    xgridindex = pd.Index(underlying_xi_pdfgrid[0].xgrid, name="x")
    df = pd.DataFrame(
        indicator.mean(axis=0),
        index=flavs,
        columns=xgridindex,
    )
    return df

@figuregen
@check_scale('xscale', allow_none=True)
def plot_xi_1sigma_proxy(xi_1sigma_proxy, Q, xscale):
    yield from plot_xi_1sigma(xi_1sigma_proxy, Q, xscale)

@figuregen
@check_scale('xscale', allow_none=True)
def plot_xi_1sigma_comparison(xi_1sigma_proxy, xi_1sigma, Q, xscale):
    x = xi_1sigma.columns.values
    xi_data = xi_1sigma.values
    xi_proxy = xi_1sigma_proxy.values
    for i, fl in enumerate(xi_1sigma.index.values):
        fig, ax = plt.subplots()
        ax.plot(x, xi_data[i, :], '*', label=r"multiple fits $\xi_{1\sigma}$ = "+f"{xi_data[i, :].mean()}")
        ax.plot(x, xi_proxy[i, :], '*', label=r"single rep proxy $\xi_{1\sigma}$ = "+f"{xi_proxy[i, :].mean()}")
        ax.axhline(0.68, linestyle=":", color="k", label="expected value")
        ax.set_ylim([0, 1])
        ax.set_title(r"$\xi_{1\sigma}$"+f" for Q={Q}, ${fl}$ PDF.")
        ax.set_xlabel("x")
        ax.set_ylabel(r"$\xi_{1\sigma}$")
        ax.set_xscale(xscale)
        ax.legend()
        yield fig

@make_argcheck(check_basis)
def pdf_bias_grid(
    pdf:PDF, Q:(float,int), basis:(str, Basis)='flavour',
    flavours:(list, tuple, type(None))=None,
    logspace_lims=(-5, -1),
    logspace_num=3,
    linspace_lims=(0.1, 1),
    linspace_num=3):
    xgrid_log = np.logspace(*logspace_lims, num=logspace_num, endpoint=False)
    xgrid_lin = np.linspace(*linspace_lims, num=linspace_num, endpoint=False)
    xgrid = np.concatenate((xgrid_log, xgrid_lin), axis=0)
    return xplotting_grid(
        pdf, Q, xgrid=xgrid, basis=basis, flavours=flavours)

fits_pdf_bias_grid = collect("pdf_bias_grid", ("fits", "fitpdf"))
underlying_pdf_bias_grid = collect("pdf_bias_grid", ("fitunderlyinglaw",))

def multifits_flavours_covmat(fits_pdf_bias_grid):
    fits_replica_grids = np.array([xpg.grid_values for xpg in fits_pdf_bias_grid])
    fits_central_grids = fits_replica_grids.mean(axis=0, keepdims=True)
    diffs = np.concatenate((fits_central_grids - fits_replica_grids), axis=0) # superreps, flav, x
    covs = []
    for fli in range(diffs.shape[1]): # loop over fl index
        covs.append(np.cov(diffs[:, fli, :], rowvar=False))
    return covs


def multifits_flavours_bias_variance(
    fits_pdf_bias_grid, underlying_pdf_bias_grid, multifits_flavours_covmat):
    fits_replica_grids = np.array([xpg.grid_values for xpg in fits_pdf_bias_grid])
    fits_central_grids = fits_replica_grids.mean(axis=1, keepdims=True)
    underlying_grid = underlying_pdf_bias_grid[0].grid_values[0][np.newaxis, np.newaxis, ...]
    flavours_invs = np.array([la.inv(cov) for cov in multifits_flavours_covmat])

    bias_diff = (fits_central_grids - underlying_grid)[:, 0, ...] #fit, flav, x
    var_diff = fits_central_grids - fits_replica_grids # fit, rep, flav, x

    bias_xM = (bias_diff[..., np.newaxis] * flavours_invs[np.newaxis, ...]) #fit, flav, x, x
    bias = (bias_xM * bias_diff[:, :, np.newaxis, :]).sum(axis=(-2, -1)) #fit, flav
    expected_bias = np.mean(bias, axis=0)

    var_xM = (var_diff[..., np.newaxis] * flavours_invs[np.newaxis, np.newaxis, ...]) #fit, rep, flav, x, x
    var = (var_xM * var_diff[:, :, :, np.newaxis, :]).sum(axis=(-2, -1)).mean(axis=1) # fit, flav
    expected_var = np.mean(var, axis=0)

    return expected_bias/expected_var

def multifits_total_covmat(fits_pdf_bias_grid):
    fits_replica_grids = np.array([xpg.grid_values for xpg in fits_pdf_bias_grid])
    fits_central_grids = fits_replica_grids.mean(axis=0, keepdims=True)
    diffs = np.concatenate((fits_central_grids - fits_replica_grids), axis=0) # superreps, flav, x
    diffs_total = diffs.reshape(-1, diffs.shape[1]*diffs.shape[2]) # superreps, (flav * x)
    return np.cov(diffs_total, rowvar=False)

def multifits_total_pdf_bias_variance(
    multifits_total_covmat, fits_pdf_bias_grid, underlying_pdf_bias_grid):
    underlying_grid = underlying_pdf_bias_grid[0].grid_values[0][np.newaxis, np.newaxis, ...]
    flavx = underlying_grid.shape[2] * underlying_grid.shape[3]
    underlying_grid = underlying_grid.reshape(1, 1, flavx)
    fits_replica_grids = np.array([xpg.grid_values.reshape(-1, flavx) for xpg in fits_pdf_bias_grid])
    fits_central_grids = fits_replica_grids.mean(axis=1, keepdims=True)

    bias_diff = (fits_central_grids - underlying_grid)[:, 0, :] #fit, (flav*x)
    var_diff = fits_central_grids - fits_replica_grids # fit, rep, (flav*x)
    total_inv = la.inv(multifits_total_covmat) # flav*x, flav*x

    bias_xM = (bias_diff[..., np.newaxis] * total_inv[np.newaxis, ...]) #fit, flav*x, flav*x
    bias = (bias_xM * bias_diff[:, np.newaxis, :]).sum(axis=(-2, -1))/flavx #fit
    expected_bias = np.mean(bias, axis=0)

    var_xM = (var_diff[..., np.newaxis] * total_inv[np.newaxis, np.newaxis, ...]) #fit, rep, flav*x
    var = (var_xM * var_diff[:, :, np.newaxis, :]).sum(axis=(-2, -1)).mean(axis=1)/flavx # fit
    expected_var = np.mean(var, axis=0)

    return expected_bias/expected_var

@table
def multifits_pdf_bias_variance_table(
    multifits_total_pdf_bias_variance,
    multifits_flavours_bias_variance,
    underlying_pdf_bias_grid
):
    flavids = underlying_pdf_bias_grid[0].flavours
    basis = underlying_pdf_bias_grid[0].basis
    records = []
    for flid, ratio in zip(flavids, multifits_flavours_bias_variance):
        records.append(dict(
            flavour=f"${basis.elementlabel(flid)}$",
            ratio=ratio
        ))
    records.append(dict(
        flavour="total",
        ratio=multifits_total_pdf_bias_variance
    ))
    df = pd.DataFrame.from_records(records, index="flavour", columns=("flavour", "ratio"))
    df.columns = ("bias/variance",)
    return df


fits_central_diff_logistic_experiment = collect("central_diff_logistic_experiment", ("fits", "fitcontext"))

def data_xi_experiment(fits_central_diff_logistic_experiment):
    fits_logistic = np.array(fits_central_diff_logistic_experiment)
    return fits_logistic.mean(axis=0) # average across fits

@figure
def plot_data_xi_experiment(data_xi_experiment, experiment):
    fig, ax = plt.subplots()
    ax.plot(
        data_xi_experiment,
        "*",
        label=r"$\xi_{1\sigma}$ = "+f"{data_xi_experiment.mean()}, from multifits"
    )
    ax.axhline(0.68, linestyle=":", color="k", label=r"$1\sigma$ "+"expected value")
    ax.axhline(0.95, linestyle=":", color="k", label=r"$2\sigma$ "+"expected value")
    ax.set_ylim((0, 1))
    ax.set_xlabel("datapoint index")
    ax.set_title(r"$\xi_{1\sigma}$ for "+str(experiment))
    ax.legend()
    return fig