"""
multiclosure_pdf_output.py

Module containing all of the plots and tables for multiclosure estimators in
PDF space.

"""
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.special

from reportengine.figure import figure, figuregen
from reportengine.table import table
from validphys import plotutils
from validphys.closuretest.multiclosure import DEFAULT_SEED
from validphys.closuretest.multiclosure_pdf import (
    XI_FLAVOURS,
    bootstrap_pdf_differences,
    fits_covariance_matrix_by_flavour,
    fits_covariance_matrix_totalpdf,
    replica_and_central_diff_totalpdf,
    xi_flavour_x,
    xi_totalpdf,
)


@table
def xi_flavour_table(xi_flavour_x, xi_totalpdf):
    """For each flavour take the mean of xi_flavour_x across x to get a single
    number, which is the proportion of points on the central PDF which are
    within 1 sigma. This is calculated from the replicas of the underlying PDF.

    Returns
    -------

    xi_flavour: pd.DataFrame
        table of xi by flavour

    """
    data = np.concatenate((xi_flavour_x.mean(axis=-1), [xi_totalpdf]), axis=0)[:, np.newaxis]
    index = pd.Index([f"${XI_FLAVOURS[0]}$", *XI_FLAVOURS[1:], "Total"], name="flavour")
    return pd.DataFrame(data, columns=[r"measured $\xi_{1\sigma}$"], index=index)


@figuregen
def plot_xi_flavour_x(
    xi_flavour_x,
    Q,
    internal_singlet_gluon_xgrid,
    internal_nonsinglet_xgrid,
    multiclosure_nx=4,
    use_x_basis=False,
):
    """For each flavour plot xi for each x-point. By default xi is calculated
    and plotted in the basis which diagonalises the covmat, which is estimated
    from the union of all the replicas. However, if ``use_x_basis`` is ``True``
    then xi will be calculated and plotted in the x-basis.

    """
    # treat singlet and gluon separately
    if use_x_basis:
        x_for_plot = 2 * [internal_singlet_gluon_xgrid] + 5 * [internal_nonsinglet_xgrid]
        x_label = "x"
    else:
        x_for_plot = 7 * [np.arange(multiclosure_nx)]
        x_label = "estimated covariance eigenvectors (ascending in eigenvalue)"

    for i, fl in enumerate(XI_FLAVOURS):
        if i == 0:
            fl = f"${fl}$"
        fig, ax = plotutils.subplots()
        ax.plot(
            x_for_plot[i],
            xi_flavour_x[i, :],
            "*",
            label=r"$\xi_{1\sigma}$ = " + f"{xi_flavour_x[i, :].mean():.2f}",
            clip_on=False,
        )
        ax.axhline(0.68, linestyle=":", color="k", label="expected value")
        ax.set_ylim([0, 1])
        ax.set_title(r"$\xi_{1\sigma}$" + f" for Q={Q} GeV, {fl} PDF")
        ax.set_xlabel(x_label)
        ax.set_ylabel(r"$\xi_{1\sigma}$")
        ax.legend()
        yield fig


@figure
def plot_pdf_central_diff_histogram(replica_and_central_diff_totalpdf):
    """Histogram of the difference between central PDF
    and underlying law normalised by the corresponding replica
    standard deviation for all points in x and flavour alongside a scaled
    Gaussian. Total xi is proportion of the histogram which falls within the
    central 1-sigma confidence interval.

    """
    sigma, delta = replica_and_central_diff_totalpdf
    scaled_diffs = (delta / sigma).flatten()
    fig, ax = plotutils.subplots()
    ax.hist(scaled_diffs, bins=50, density=True, label="Central PDF distribution")
    xlim = (-5, 5)
    ax.set_xlim(xlim)

    x = np.linspace(*xlim, 100)
    ax.plot(
        x,
        scipy.stats.norm.pdf(x),
        "-k",
        label="Normal distribution",
    )
    ax.legend()
    ax.set_xlabel("Difference to input PDF")
    return fig


@table
def fits_pdf_bias_variance_ratio(fits_pdf_flavour_bias_variance, fits_pdf_total_bias_variance):
    """Returns a table with the values of mean bias / mean variance with mean
    referring to mean across fits, by flavour. Includes total across all
    flavours allowing for correlations.

    """
    records = []
    bias_fl, variance_fl = fits_pdf_flavour_bias_variance
    bias_tot, variance_tot = fits_pdf_total_bias_variance
    for i, fl in enumerate(XI_FLAVOURS):
        if i == 0:
            fl = f"${fl}$"
        records.append(
            dict(
                flavour=fl,
                bias=bias_fl[i],
                variance=variance_fl[i],
                ratio=bias_fl[i] / variance_fl[i],
                ratio_sqrt=np.sqrt(bias_fl[i] / variance_fl[i]),
            )
        )
    records.append(
        dict(
            flavour="Total",
            bias=bias_tot,
            variance=variance_tot,
            ratio=bias_tot / variance_tot,
            ratio_sqrt=np.sqrt(bias_tot / variance_tot),
        )
    )
    df = pd.DataFrame.from_records(
        records, index="flavour", columns=["flavour", "bias", "variance", "ratio", "ratio_sqrt"]
    )
    df.columns = ["bias", "variance", "bias/variance", "sqrt(bias/variance)"]
    return df


@table
def fits_pdf_expected_xi_from_ratio(fits_pdf_sqrt_ratio):
    """Like
    :py:func:`validphys.closuretest.multiclosure_output.expected_xi_from_bias_variance`
    but in PDF space. An estimate is made of the integral across the central
    difference distribution, with domain defined by the replica distribution.
    For more details see
    :py:func:`validphys.closuretest.multiclosure_output.expected_xi_from_bias_variance`.

    """
    df_in = fits_pdf_sqrt_ratio
    n_sigma_in_variance = 1 / df_in.values
    # pylint can't find erf here, disable error in this function
    # pylint: disable=no-member
    estimated_integral = scipy.special.erf(n_sigma_in_variance / np.sqrt(2))
    return pd.DataFrame(
        estimated_integral, index=df_in.index, columns=[r"expected $\xi_{1\sigma}$"]
    )


@table
def fits_pdf_compare_xi_to_expected(fits_pdf_expected_xi_from_ratio, xi_flavour_table):
    """Two-column table comparing the measured value of xi for each flavour to the
    value calculated from the bias/variance.

    """
    return pd.concat((xi_flavour_table, fits_pdf_expected_xi_from_ratio), axis=1)


@table
def fits_bootstrap_pdf_xi_table(
    fits_xi_grid_values,
    underlying_xi_grid_values,
    multiclosure_underlyinglaw,
    multiclosure_nx=4,
    n_boot=100,
    boot_seed=DEFAULT_SEED,
    use_x_basis=False,
):
    """Perform a bootstrap sampling across fits and replicas of xi, by flavour
    and total and then tabulate the mean and error.

    """
    rng = np.random.RandomState(seed=boot_seed)
    xi_boot = []
    for _ in range(n_boot):
        # perform single bootstrap
        boot_central_diff, boot_rep_diff = bootstrap_pdf_differences(
            fits_xi_grid_values,
            underlying_xi_grid_values,
            multiclosure_underlyinglaw,
            rng,
        )

        flav_cov = fits_covariance_matrix_by_flavour(boot_rep_diff)
        total_cov = fits_covariance_matrix_totalpdf(boot_rep_diff, multiclosure_nx)
        rep_cent_diff = replica_and_central_diff_totalpdf(
            boot_rep_diff, boot_central_diff, total_cov, multiclosure_nx, use_x_basis
        )
        xi_flav = xi_flavour_x(boot_rep_diff, boot_central_diff, flav_cov, use_x_basis)
        xi_total = xi_totalpdf(rep_cent_diff)
        xi_data = np.concatenate((xi_flav.mean(axis=-1), [xi_total]), axis=0)
        xi_boot.append(xi_data)
    # construct table in this action, since bootstrap rawdata isn't required elsewhere
    index = pd.Index([f"${XI_FLAVOURS[0]}$", *XI_FLAVOURS[1:], "Total"], name="flavour")
    res = np.concatenate(
        (
            np.mean(xi_boot, axis=0)[:, np.newaxis],
            np.std(xi_boot, axis=0)[:, np.newaxis],
        ),
        axis=1,
    )
    return pd.DataFrame(
        res,
        columns=[r"bootstrap mean $\xi_{1\sigma}$", r"bootstrap std. dev. $\xi_{1\sigma}$"],
        index=index,
    )


@table
def fits_bootstrap_pdf_sqrt_ratio_table(fits_bootstrap_pdf_sqrt_ratio):
    """Tabulate the mean and standard deviation across bootstrap samples of the
    sqrt ratio of bias/variance in PDF space, with a row for each flavour and
    the total. For more information on the bootstrap sampling see
    :py:func:`fits_bootstrap_pdf_ratio`.

    """
    index = pd.Index([f"${XI_FLAVOURS[0]}$", *XI_FLAVOURS[1:], "Total"], name="flavour")
    # add extra dimension here so that 1-D arrays appear as columns for table
    res = np.concatenate(
        (
            np.mean(fits_bootstrap_pdf_sqrt_ratio, axis=0)[:, np.newaxis],
            np.std(fits_bootstrap_pdf_sqrt_ratio, axis=0)[:, np.newaxis],
        ),
        axis=1,
    )
    return pd.DataFrame(
        res,
        columns=[r"bootstrap mean sqrt ratio", r"bootstrap std. dev. sqrt ratio"],
        index=index,
    )


@table
def fits_bootstrap_pdf_expected_xi_table(fits_bootstrap_pdf_expected_xi):
    """Tabulate the mean and standard deviation across bootstrap samples of
    :py:func:`fits_bootstrap_pdf_expected_xi` with a row for each flavour and
    the total expected xi.

    """
    df = fits_bootstrap_pdf_sqrt_ratio_table(fits_bootstrap_pdf_expected_xi)
    return pd.DataFrame(
        df.values,
        index=df.index,
        columns=[
            r"bootstrap mean expected $\xi_{1\sigma}$ from ratio",
            r"bootstrap std. dev. expected $\xi_{1\sigma}$ from ratio",
        ],
    )


@table
def fits_bootstrap_pdf_compare_xi_to_expected(
    fits_bootstrap_pdf_expected_xi_table, fits_bootstrap_pdf_xi_table
):
    """Table comparing the mean and standard deviation across bootstrap samples
    of the measured value of xi to the value calculated from bias/variance
    in PDF space. This is done for each flavour and for the total across all
    flavours accounting for correlations.

    """
    return fits_pdf_compare_xi_to_expected(
        fits_bootstrap_pdf_expected_xi_table, fits_bootstrap_pdf_xi_table
    )


def plot_pdf_matrix(matrix, n_x, **kwargs):
    """Utility function which, given a covmat/corrmat for all flavours and
    x, plots it with appropriate labels. Input matrix is expected to be
    size (n_flavours*n_x) * (n_flavours*n_x).

    Parameters
    ----------

    matrix: np.array
        square matrix which must be (n_flavours*n_x) * (n_flavours*n_x) with
        elements ordered like:
        (flavour0_x0, flavour0_x1, ..., flavourN_x0, ..., flavourN_xN)
        i.e. the points along x for flavour 0, then points along x for flavour 1
        etc.
    **kwargs:
        keyword arguments for the matplotlib.axes.Axes.imshow function

    Notes
    -----

    See matplotlib.axes.Axes.imshow for more details on the plotting function.

    """
    labels = [f"${XI_FLAVOURS[0]}$", *XI_FLAVOURS[1:]]
    # we want to centre the labels on each of the xgrids
    # ticks appear shifted right by 0.5, so we account for this here
    ticks = (n_x) / 2 + n_x * np.arange(len(XI_FLAVOURS)) - 0.5
    fig, ax = plotutils.subplots()
    im = ax.imshow(matrix, **kwargs)
    fig.colorbar(im, ax=ax)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    # for some reason the axis can be resized, so fix that here
    ax.set_ylim([n_x * len(XI_FLAVOURS) - 0.5, -0.5])
    ax.set_xlim([-0.5, n_x * len(XI_FLAVOURS) - 0.5])
    # looks more natural with xticks at top
    ax.xaxis.tick_top()
    return fig, ax


@figure
def plot_multiclosure_covariance_matrix(fits_covariance_matrix_totalpdf, multiclosure_nx=4):
    """Plot the covariance matrix for all flavours. The covariance matrix has
    shape n_flavours * n_x, where each block is the covariance of the replica
    PDFs on the x-grid defined in :py:func:`xi_pdfgrids`.

    """
    fig, ax = plot_pdf_matrix(fits_covariance_matrix_totalpdf, multiclosure_nx)
    ax.set_title("Covariance matrix estimated from multiclosure replicas")
    return fig


@figure
def plot_multiclosure_correlation_matrix(fits_correlation_matrix_totalpdf, multiclosure_nx=4):
    """Like plot_multiclosure_covariance_matrix but plots the total correlation
    matrix.

    """
    fig, ax = plot_pdf_matrix(fits_correlation_matrix_totalpdf, multiclosure_nx, vmin=0, vmax=1)
    ax.set_title("Correlation matrix estimated from multiclosure replicas")
    return fig


@figure
def plot_multiclosure_correlation_eigenvalues(fits_correlation_matrix_totalpdf):
    """Plot scatter points for each of the eigenvalues from the estimated
    correlation matrix from the multiclosure PDFs in flavour and x.

    In the legend add the ratio of the largest eigenvalue over the smallest
    eigenvalue, aka the l2 condition number of the correlation matrix.

    """
    # e_val is in ascending order
    e_val, _ = la.eigh(fits_correlation_matrix_totalpdf)
    l2_condition = e_val[-1] / e_val[0]
    fig, ax = plotutils.subplots()
    ax.plot(e_val, "*", label=f"l2-condition: {l2_condition:,.0f}")
    ax.set_ylabel(r"$\lambda_{\rm corr}$")
    ax.set_xlabel("eigenvalue index (ascending)")
    ax.set_title("Eigenvalues of PDF correlation matrix, $N_{x}=$" + f"{len(e_val)}")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend()
    return fig
