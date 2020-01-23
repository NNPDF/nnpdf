"""
closuretest/multiclosure.py

containing all of the statistical estimators which are
averaged across multiple fits or a single replica proxy fit
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from reportengine import collect
from reportengine.table import table
from reportengine.figure import figuregen
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

def expected_variance(fits_phi_experiment):
    var_centrals = [
        phi[0]**2 * phi[1] for phi in fits_phi_experiment]
    var_mean = np.mean(var_centrals)
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

fits_xi_pdfgrid = collect("xi_pdfgrid", ("fits", "fitpdf"))
underlying_xi_pdfgrid = collect("xi_pdfgrid", ("fitunderlyinglaw",))

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
