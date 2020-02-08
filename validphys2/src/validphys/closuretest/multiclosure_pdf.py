"""
multiclosure_pdf.py

Module containing all of the actions related to statistical estimators across
multiple closure fits or proxy fits defined in PDF space.

"""
import numpy as np
import scipy.linalg as la
import scipy.special
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from reportengine import collect
from reportengine.table import table
from reportengine.figure import figuregen
from reportengine.checks import make_argcheck

from validphys.pdfgrids import xplotting_grid, check_basis
from validphys.pdfbases import Basis
from validphys.core import PDF
from validphys.checks import check_scale
from validphys.calcutils import calc_chi2

# exclude charm
XI_FLAVOURS = (r'\Sigma', 'gluon', 'V', 'V3', 'V8', 'T3', 'T8',)

N = 24
SINGLET_GLUON_XGRID = np.concatenate(
    (np.logspace(-3, -1, int(N/2), endpoint=False), np.linspace(0.1, 0.5, int(N/2))),
    axis=0
)
NONSINGLET_XGRID = np.linspace(0.1, 0.5, N)


def xi_pdfgrids(pdf:PDF, Q:(float,int)):
    """Generate PDF grids which are required for calculating xi in PDF space
    in the NN31IC basis, exclusing the charm we want to specify different xgrids
    for different flavours to avoid sampling PDFs in deep extrapolation, the
    limits are:

        gluon and singlet: 10^-3 < x < 0.5
        other non-singlets: 0.1 < x < 0.5

    TODO: check above limits.

    Returns
    -------
    tuple of xplotting_grids, one for gluon and singlet and one for other
    non-singlets
    """
    #TODO: do we want to hardcode Q here as well?
    singlet_gluon_grid = xplotting_grid(
        pdf,
        Q,
        xgrid=SINGLET_GLUON_XGRID,
        basis="NN31IC",
        flavours=XI_FLAVOURS[:2]
    )

    nonsinglet_grid = xplotting_grid(
        pdf,
        Q,
        xgrid=NONSINGLET_XGRID,
        basis="NN31IC",
        flavours=XI_FLAVOURS[2:]
    )
    return singlet_gluon_grid, nonsinglet_grid

def xi_grid_values(xi_pdfgrids):
    """grid values from the xi_pdfgrids concatenated as single numpy array"""
    glu_sin_grid, nonsin_grid = xi_pdfgrids
    # grid values have shape: replica, flavour, x
    # concatenate along flavour
    return np.concatenate(
        (glu_sin_grid.grid_values, nonsin_grid.grid_values),
        axis=1
    )

def underlying_xi_grid_values(multiclosure_underlyinglaw: PDF, Q:(float,int)):
    """Like xi_pdfgrids but setting the PDF as the underlying law, extracted
    from a set of fits.
    """
    underlying_grid = xi_pdfgrids(multiclosure_underlyinglaw, Q)
    return xi_grid_values(underlying_grid)

def pdf_central_difference(
    xi_grid_values, underlying_xi_grid_values, multiclosure_underlyinglaw):
    """Calculate the difference between underlying law and central pdf for,
    specifically:

        underlying_grid - mean(grid_vals)

    where mean is across replicas

    Returns
    -------
    diffs: np.array
        array of diffs with shape (flavour, x)

    """
    underlying_central = multiclosure_underlyinglaw.stats_class(
        underlying_xi_grid_values).central_value()
    return underlying_central - np.mean(xi_grid_values, axis=0)

def pdf_replica_difference(xi_grid_values):
    """Calculate the difference between the central pdf and the replica pdfs,
    specifically:

        mean(grid_vals) - grid_vals

    where the mean is across replicas

    Returns
    -------
    diffs: np.array
        array of diffs with shape (replicas, flavour, x)
    """
    return xi_grid_values.mean(axis=0, keepdims=True) - xi_grid_values

fits_replica_difference = collect("pdf_replica_difference", ("fits", "fitpdf"))
fits_central_difference = collect("pdf_central_difference", ("fits", "fitpdf"))


def fits_covariance_matrix_by_flavour(fits_replica_difference):
    """Given a set of pdf grids from multiple closure tests, obtain an estimate
    of the covariance matrix for each flavour seperately, return as a list of
    covmats
    """
    # diffs want to be calculated on the per fit level
    super_diffs = np.concatenate(fits_replica_difference, axis=0)
    covmats = []
    for i in range(len(XI_FLAVOURS)):
        covmats.append(np.cov(super_diffs[:, i, :], rowvar=False))
    return covmats

def xi_flavour_x(
    fits_replica_difference,
    fits_central_difference,
    fits_covariance_matrix_by_flavour
):
    """for a set of fits calculate the indicator function

        I_{[-sigma, sigma]}(delta)

    where sigma is the RMS difference between central and replicas PDF
    and delta is the difference between central PDF and underlying law.

    differences are all rotated to basis which diagonalises the covariance
    matrix that was estimated from the super set of all fit replicas.

    finally take the mean across fits to get xi in flavour and x

    """
    rep_diff = np.asarray(fits_replica_difference)
    central_diff = np.asarray(fits_central_difference)

    xis = []
    for i in range(len(XI_FLAVOURS)):
        _, e_vec = la.eigh(fits_covariance_matrix_by_flavour[i])
        # uncomment to diagonalise
        # put x on first axis
        diag_central_diff = e_vec.T @ central_diff[:, i, :].T
        # put x on second to last axis
        diag_rep_diff = e_vec.T @ rep_diff[:, :, i, :].transpose(1, 2, 0)
        # uncoment for xbasis
        #diag_central_diff = central_diff[:, i, :].T
        #diag_rep_diff = rep_diff[:, :, i, :].transpose(1, 2, 0)
        var_diff = (diag_rep_diff)**2
        sigma = np.sqrt(var_diff.mean(axis=0)) # mean across reps
        # indicator and mean across fits
        xi = np.asarray(abs(diag_central_diff) < sigma, dtype=int).mean(axis=1)
        xis.append(xi)
    return np.asarray(xis)

def fits_covariance_matrix_totalpdf(fits_replica_difference):
    """Given a set of pdf grids from multiple closure tests, obtain an estimate
    of the covariance matrix allowing for correlations across flavours
    """
    # diffs want to be calculated on the per fit level
    super_diffs = np.concatenate(fits_replica_difference, axis=0).reshape(
        -1, N*len(XI_FLAVOURS)) # reshape to correlate flavours
    return np.cov(super_diffs, rowvar=False)

fits_indicator_function_totalpdf = collect(
    "pdf_indicator_function_totalpdf", ("fits", "fitpdf"))

def xi_totalpdf(
    fits_replica_difference,
    fits_central_difference,
    fits_covariance_matrix_totalpdf
):
    """Like xi_flavour_x except calculate the total xi across flavours and x
    accounting for correlations
    """
    # keep fits and reps then reshape flavour x to one dim
    rep_diff = np.asarray(fits_replica_difference).reshape(
        len(fits_replica_difference), -1, N*len(XI_FLAVOURS)
    )
    central_diff = np.asarray(fits_central_difference).reshape(
        -1, N*len(XI_FLAVOURS))
    _, e_vec = la.eigh(fits_covariance_matrix_totalpdf)
    # put flavourx on first axis
    diag_central_diff = e_vec.T @ central_diff.T
    # need reps on second axis
    diag_rep_diff = e_vec.T @ rep_diff.transpose(1, 2, 0)
    # uncomment for x basis
    #diag_central_diff = central_diff.T
    #diag_rep_diff = rep_diff.transpose(1, 2, 0)
    var_diff = (diag_rep_diff)**2
    sigma = np.sqrt(var_diff.mean(axis=0)) # mean across reps
    # indicator and mean across all
    return np.asarray(abs(diag_central_diff) < sigma, dtype=int).mean()

@table
def xi_flavour_table(xi_flavour_x, xi_totalpdf):
    """for each flavour take the mean of xi_flavour_x across x to get single
    number of proportion points on the central PDF which are within 1 sigma,
    calculated from the replicas, of the underlying PDF.

    Returns
    xi_flavour: pd.DataFrame
        table of xi by flavour
    """
    data = np.concatenate(
        (xi_flavour_x.mean(axis=-1), [xi_totalpdf]), axis=0)[:, np.newaxis]
    return pd.DataFrame(
        data,
        columns=[r"$\xi_{1\sigma}$"],
        index=[f"${XI_FLAVOURS[0]}$", *XI_FLAVOURS[1:], "Total"]
    )

@figuregen
def plot_xi_flavour_x(xi_flavour_x, Q):
    for i, fl in enumerate(XI_FLAVOURS):
        if i == 0:
            fl = f"${fl}$"
        fig, ax = plt.subplots()
        ax.plot(
            xi_flavour_x[i, :],
            '*',
            label=r"$\xi_{1\sigma}$ = "+f"{xi_flavour_x[i, :].mean():.2f}"
        )
        ax.axhline(0.68, linestyle=":", color="k", label="expected value")
        ax.set_ylim([0, 1])
        ax.set_title(r"$\xi_{1\sigma}$"+f" for Q={Q}, {fl} PDF.")
        ax.set_xlabel("estimated covariance eigenvectors (ascending in eigenvalue)")
        ax.set_ylabel(r"$\xi_{1\sigma}$")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend()
        yield fig

def fits_sqrt_covmat_by_flavour(fits_covariance_matrix_by_flavour):
    """For each flavour covariance matrix calculate the sqrt covmat
    (cholesky lower triangular)
    """
    return [
        la.cholesky(cov, lower=True)
        for cov in fits_covariance_matrix_by_flavour
    ]

def fits_pdf_flavour_ratio(
    fits_sqrt_covmat_by_flavour,
    fits_central_difference,
    fits_replica_difference,
):
    """calculate the bias (chi2 between central PDF and underlying PDF)
    for each flavour and the variance (mean chi2 between replica and central PDF)
    then return a numpy array with shape (flavours, 2) with second axis being
    bias, variance

    """
    central_diff = np.asarray(fits_central_difference)
    rep_diff = np.asarray(fits_replica_difference)
    ratios = []
    for i in range(len(XI_FLAVOURS)):
        bias = calc_chi2(
                fits_sqrt_covmat_by_flavour[i],
                central_diff[:, i, :].T
            )
        variance = np.mean(calc_chi2(
                fits_sqrt_covmat_by_flavour[i],
                rep_diff[:, :, i, :].transpose(2, 1, 0) # need x on first axis
            ), axis=0)
        ratios.append(np.mean(bias)/np.mean(variance))
    return ratios

def fits_pdf_total_ratio(
    fits_central_difference,
    fits_replica_difference,
    fits_covariance_matrix_totalpdf
):
    """Calculate the total bias and variance for all flavours and x allowing for
    correlations across flavour,

    Returns
    -------
    ratio_data: tuple
        required data for calculating mean(bias) over mean(variance) across fits
        in form of tuple (bias, variance)
    """
    central_diff = np.asarray(fits_central_difference).reshape(-1, N*len(XI_FLAVOURS))
    rep_diff = np.asarray(fits_replica_difference).reshape(
        len(fits_replica_difference), -1, N*len(XI_FLAVOURS))

    sqrtcov = la.cholesky(fits_covariance_matrix_totalpdf, lower=True)

    bias = calc_chi2(sqrtcov, central_diff.T)
    # need flav x on first axis
    variance = np.mean(calc_chi2(sqrtcov, rep_diff.transpose(2, 1, 0)), axis=0)
    return np.mean(bias) / np.mean(variance)

@table
def fits_pdf_bias_variance_ratio(
    fits_pdf_flavour_ratio, fits_pdf_total_ratio):
    """Returns a table with the values of mean bias / mean variance with mean
    referring to mean across fits, by flavour. Includes total across all
    flavours allowing for correlations

    """
    records = []
    for i, fl in enumerate(XI_FLAVOURS):
        if i == 0:
            fl = f"${fl}$"
        records.append(dict(
            flavour=fl,
            ratio=fits_pdf_flavour_ratio[i]
        ))
    records.append(dict(
        flavour="Total",
        ratio=fits_pdf_total_ratio
    ))
    df = pd.DataFrame.from_records(records,
        index="flavour",
        columns=["flavour", "ratio"]
    )
    df.columns = ["bias/variance"]
    return df

@table
def fits_pdf_sqrt_ratio(
    fits_pdf_bias_variance_ratio
):
    """Like fits_pdf_bias_variance_ratio except taking the sqrt, we see how
    faithful our uncertainty is in units of the standard deviation
    """
    df_in = fits_pdf_bias_variance_ratio
    data = np.sqrt(df_in.values)
    return pd.DataFrame(
        data,
        index=df_in.index,
        columns=["sqrt bias/variance"]
    )

@table
def fits_pdf_expected_xi_from_ratio(fits_pdf_sqrt_ratio):
    """Like expected_xi_from_bias_variance but in PDF space. Estimate the
    integral across central difference distribution, with domain defined by the
    replica distribution for more details see expected_xi_from_bias_variance
    """
    df_in = fits_pdf_sqrt_ratio
    n_sigma_in_variance = 1/df_in.values
    # pylint can't find erf here, disable error in this function
    #pylint: disable=no-member
    estimated_integral = scipy.special.erf(n_sigma_in_variance/np.sqrt(2))
    return pd.DataFrame(
        np.sqrt(estimated_integral),
        index=df_in.index,
        columns=[r"expected $\xi_{1\sigma}$"]
    )

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

    indicator = np.asarray(abs(central_diff) < sigma, dtype=float) #fit, flav, x
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
        ax.plot(x, xi_data[i, :], '*', label=r"multiple fits $\xi_{1\sigma}$ = "+f"{xi_data[i, :].mean():.2f}")
        ax.plot(x, xi_proxy[i, :], '*', label=r"single rep proxy $\xi_{1\sigma}$ = "+f"{xi_proxy[i, :].mean():.2f}")
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
    fits_replica_grids = np.asarray([xpg.grid_values for xpg in fits_pdf_bias_grid])
    fits_central_grids = fits_replica_grids.mean(axis=0, keepdims=True)
    diffs = np.concatenate((fits_central_grids - fits_replica_grids), axis=0) # superreps, flav, x
    covs = []
    for fli in range(diffs.shape[1]): # loop over fl index
        covs.append(np.cov(diffs[:, fli, :], rowvar=False))
    return covs


def multifits_flavours_bias_variance(
    fits_pdf_bias_grid, underlying_pdf_bias_grid, multifits_flavours_covmat):
    fits_replica_grids = np.asarray([xpg.grid_values for xpg in fits_pdf_bias_grid])
    fits_central_grids = fits_replica_grids.mean(axis=1, keepdims=True)
    underlying_grid = underlying_pdf_bias_grid[0].grid_values[0][np.newaxis, np.newaxis, ...]
    flavours_invs = np.asarray([la.inv(cov) for cov in multifits_flavours_covmat])

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
    fits_replica_grids = np.asarray([xpg.grid_values for xpg in fits_pdf_bias_grid])
    fits_central_grids = fits_replica_grids.mean(axis=0, keepdims=True)
    diffs = np.concatenate((fits_central_grids - fits_replica_grids), axis=0) # superreps, flav, x
    diffs_total = diffs.reshape(-1, diffs.shape[1]*diffs.shape[2]) # superreps, (flav * x)
    return np.cov(diffs_total, rowvar=False)

def multifits_total_pdf_bias_variance(
    multifits_total_covmat, fits_pdf_bias_grid, underlying_pdf_bias_grid):
    underlying_grid = underlying_pdf_bias_grid[0].grid_values[0][np.newaxis, np.newaxis, ...]
    flavx = underlying_grid.shape[2] * underlying_grid.shape[3]
    underlying_grid = underlying_grid.reshape(1, 1, flavx)
    fits_replica_grids = np.asarray([xpg.grid_values.reshape(-1, flavx) for xpg in fits_pdf_bias_grid])
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
