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
from reportengine.figure import figuregen, figure

from validphys.pdfgrids import xplotting_grid
from validphys.core import PDF
from validphys.calcutils import calc_chi2

# Define the NN31IC basis with the charm PDF excluded. It is excluded because
# the exercises carried out with this module are intended to be done in the
# data region and at the fitting scale, where the charm is noisy. Results
# obtained with it are therefore likely to be dominated by fluctuations
XI_FLAVOURS = (r"\Sigma", "gluon", "V", "V3", "V8", "T3", "T8")


def internal_Nx(multiclosure_nx=4):
    """Set a default value for number of points Nx to be used in multiclosure
    PDF actions, allows for studying dependence on this number but by default
    sets Nx = 4.

    This was found to be a reasonable value for 30 fits each with
    40 replicas. Increasing the number of points beyond this leads to large
    correlations across points, which would require significantly higher
    number of fits/replicas to properly estimate.

    """
    return multiclosure_nx


def internal_singlet_gluon_xgrid(internal_Nx):
    """Given the number of x points, set up the singlet and gluon xgrids,
    which are defined as half the points being logarithmically spaced
    between 10^-3 and 0.1 and the other half of the points being linearly
    spaced between 0.1 and 0.5
    """
    return np.concatenate(
        (
            np.logspace(-3, -1, int(internal_Nx / 2), endpoint=False),
            np.linspace(0.1, 0.5, int(internal_Nx / 2)),
        ),
        axis=0,
    )


def internal_nonsinglet_xgrid(internal_Nx):
    """Given the number of x points, set up the xgrid for flavours which
    are not singlet or gluon, defined as being linearly spaced points between
    0.1 and 0.5
    """
    return np.linspace(0.1, 0.5, internal_Nx)


def xi_pdfgrids(
    pdf: PDF, Q: (float, int), internal_singlet_gluon_xgrid, internal_nonsinglet_xgrid
):
    """Generate PDF grids which are required for calculating xi in PDF space
    in the NN31IC basis, excluding the charm. We want to specify different xgrids
    for different flavours to avoid sampling PDFs in deep extrapolation regions.
    The limits are chosen to achieve this and specifically they are chosen to be:

        gluon and singlet: 10^-3 < x < 0.5
        other non-singlets: 0.1 < x < 0.5

    Returns
    -------

    tuple of xplotting_grids, one for gluon and singlet and one for other
    non-singlets

    """
    # NOTE: Could we hardcode Q to the initial scale/infer from fits?
    singlet_gluon_grid = xplotting_grid(
        pdf,
        Q,
        xgrid=internal_singlet_gluon_xgrid,
        basis="NN31IC",
        flavours=XI_FLAVOURS[:2],
    )

    nonsinglet_grid = xplotting_grid(
        pdf,
        Q,
        xgrid=internal_nonsinglet_xgrid,
        basis="NN31IC",
        flavours=XI_FLAVOURS[2:],
    )
    return singlet_gluon_grid, nonsinglet_grid


def xi_grid_values(xi_pdfgrids):
    """Grid values from the xi_pdfgrids concatenated as single numpy array"""
    glu_sin_grid, nonsin_grid = xi_pdfgrids
    # grid values have shape: replica, flavour, x
    # concatenate along flavour
    return np.concatenate((glu_sin_grid.grid_values, nonsin_grid.grid_values), axis=1)


def underlying_xi_grid_values(
    multiclosure_underlyinglaw: PDF,
    Q: (float, int),
    internal_singlet_gluon_xgrid,
    internal_nonsinglet_xgrid,
):
    """Like xi_pdfgrids but setting the PDF as the underlying law, extracted
    from a set of fits
    """
    underlying_grid = xi_pdfgrids(
        multiclosure_underlyinglaw,
        Q,
        internal_singlet_gluon_xgrid,
        internal_nonsinglet_xgrid,
    )
    return xi_grid_values(underlying_grid)


def pdf_central_difference(
    xi_grid_values, underlying_xi_grid_values, multiclosure_underlyinglaw
):
    """Calculate the difference between underlying law and central PDF for,
    specifically:

        underlying_grid - mean(grid_vals)

    where mean is across replicas

    Returns:

    diffs: np.array
        array of diffs with shape (flavour, x)

    """
    underlying_central = multiclosure_underlyinglaw.stats_class(
        underlying_xi_grid_values
    ).central_value()
    return underlying_central - np.mean(xi_grid_values, axis=0)


def pdf_replica_difference(xi_grid_values):
    """Calculate the difference between the central PDF and the replica PDFs,
    specifically:

        mean(grid_vals) - grid_vals

    where the mean is across replicas.

    Returns:

    diffs: np.array
        array of diffs with shape (replicas, flavour, x)

    """
    return xi_grid_values.mean(axis=0, keepdims=True) - xi_grid_values


fits_replica_difference = collect("pdf_replica_difference", ("fits", "fitpdf"))
fits_central_difference = collect("pdf_central_difference", ("fits", "fitpdf"))


def fits_covariance_matrix_by_flavour(fits_replica_difference):
    """Given a set of PDF grids from multiple closure tests, obtain an estimate
    of the covariance matrix for each flavour separately, return as a list of
    covmats
    """
    # diffs should be calculated on the per fit level
    super_diffs = np.concatenate(fits_replica_difference, axis=0)
    covmats = []
    for i in range(len(XI_FLAVOURS)):
        covmats.append(np.cov(super_diffs[:, i, :], rowvar=False))
    return covmats


def xi_flavour_x(
    fits_replica_difference,
    fits_central_difference,
    fits_covariance_matrix_by_flavour,
    use_x_basis=False,
):
    """For a set of fits calculate the indicator function

        I_{[-sigma, sigma]}(delta)

    where sigma is the RMS difference between central and replicas PDF
    and delta is the difference between central PDF and underlying law.

    The differences are all rotated to basis which diagonalises the covariance
    matrix that was estimated from the super set of all fit replicas.

    Finally take the mean across fits to get xi in flavour and x.

    """
    rep_diff = np.asarray(fits_replica_difference)
    central_diff = np.asarray(fits_central_difference)

    xis = []
    for i in range(len(XI_FLAVOURS)):
        if use_x_basis:
            # put x on first axis
            diag_central_diff = central_diff[:, i, :].T
            # put x on second to last axis
            diag_rep_diff = rep_diff[:, :, i, :].transpose(1, 2, 0)
        else:
            _, e_vec = la.eigh(fits_covariance_matrix_by_flavour[i])
            # put x on first axis
            diag_central_diff = e_vec.T @ central_diff[:, i, :].T
            # put x on second to last axis
            diag_rep_diff = e_vec.T @ rep_diff[:, :, i, :].transpose(1, 2, 0)
        var_diff = (diag_rep_diff) ** 2
        sigma = np.sqrt(var_diff.mean(axis=0))  # mean across reps
        # indicator and mean across fits
        xi = np.asarray(abs(diag_central_diff) < sigma, dtype=int).mean(axis=1)
        xis.append(xi)
    return np.asarray(xis)


def fits_covariance_matrix_totalpdf(fits_replica_difference, internal_Nx):
    """Given a set of PDF grids from multiple closure tests, obtain an estimate
    of the covariance matrix allowing for correlations across flavours
    """
    # diffs should be calculated on the per fit level
    super_diffs = np.concatenate(fits_replica_difference, axis=0).reshape(
        -1, internal_Nx * len(XI_FLAVOURS)
    )  # reshape to correlate flavours
    return np.cov(super_diffs, rowvar=False)


fits_indicator_function_totalpdf = collect(
    "pdf_indicator_function_totalpdf", ("fits", "fitpdf")
)


def xi_totalpdf(
    fits_replica_difference,
    fits_central_difference,
    fits_covariance_matrix_totalpdf,
    internal_Nx,
    use_x_basis=False,
):
    """Like :py:func:`xi_flavour_x` except calculate the total xi across flavours and x
    accounting for correlations
    """
    # keep fits and reps then reshape flavour x to one dim
    rep_diff = np.asarray(fits_replica_difference).reshape(
        len(fits_replica_difference), -1, internal_Nx * len(XI_FLAVOURS)
    )
    central_diff = np.asarray(fits_central_difference).reshape(
        -1, internal_Nx * len(XI_FLAVOURS)
    )
    if use_x_basis:
        diag_central_diff = central_diff.T
        diag_rep_diff = rep_diff.transpose(1, 2, 0)
    else:
        _, e_vec = la.eigh(fits_covariance_matrix_totalpdf)
        # put flavourx on first axis
        diag_central_diff = e_vec.T @ central_diff.T
        # need reps on second axis
        diag_rep_diff = e_vec.T @ rep_diff.transpose(1, 2, 0)

    var_diff = (diag_rep_diff) ** 2
    sigma = np.sqrt(var_diff.mean(axis=0))  # mean across reps
    # indicator and mean across all
    return np.asarray(abs(diag_central_diff) < sigma, dtype=int).mean()


@table
def xi_flavour_table(xi_flavour_x, xi_totalpdf):
    """For each flavour take the mean of xi_flavour_x across x to get single
    number of proportion points on the central PDF which are within 1 sigma,
    calculated from the replicas, of the underlying PDF.

    Returns
    xi_flavour: pd.DataFrame
        table of xi by flavour
    """
    data = np.concatenate((xi_flavour_x.mean(axis=-1), [xi_totalpdf]), axis=0)[
        :, np.newaxis
    ]
    index = pd.Index([f"${XI_FLAVOURS[0]}$", *XI_FLAVOURS[1:], "Total"], name="flavour")
    return pd.DataFrame(data, columns=[r"Measured $\xi_{1\sigma}$"], index=index)


@figuregen
def plot_xi_flavour_x(
    xi_flavour_x,
    Q,
    internal_singlet_gluon_xgrid,
    internal_nonsinglet_xgrid,
    internal_Nx,
    use_x_basis=False,
):
    """For each flavour plot xi for each x-point. By default xi is calculated and
    plotted in the basis which diagonalises the covmat estimated from the union
    of all the replicas. However, if ``use_x_basis`` is ``True`` then xi will be
    calculated and plotted in x basis.

    """
    # treat singlet and gluon separately
    if use_x_basis:
        x_for_plot = 2 * [internal_singlet_gluon_xgrid] + 5 * [
            internal_nonsinglet_xgrid
        ]
        x_label = "x"
    else:
        x_for_plot = 7 * [np.arange(internal_Nx)]
        x_label = "estimated covariance eigenvectors (ascending in eigenvalue)"

    for i, fl in enumerate(XI_FLAVOURS):
        if i == 0:
            fl = f"${fl}$"
        fig, ax = plt.subplots()
        ax.plot(
            x_for_plot[i],
            xi_flavour_x[i, :],
            "*",
            label=r"$\xi_{1\sigma}$ = " + f"{xi_flavour_x[i, :].mean():.2f}",
            clip_on=False,
        )
        ax.axhline(0.68, linestyle=":", color="k", label="expected value")
        ax.set_ylim([0, 1])
        ax.set_title(r"$\xi_{1\sigma}$" + f" for Q={Q}, {fl} PDF.")
        ax.set_xlabel(x_label)
        ax.set_ylabel(r"$\xi_{1\sigma}$")
        # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend()
        yield fig


def fits_sqrt_covmat_by_flavour(fits_covariance_matrix_by_flavour):
    """For each flavour covariance matrix calculate the sqrt covmat
    (cholesky lower triangular)
    """
    return [la.cholesky(cov, lower=True) for cov in fits_covariance_matrix_by_flavour]


def fits_pdf_flavour_ratio(
    fits_sqrt_covmat_by_flavour, fits_central_difference, fits_replica_difference
):
    """Calculate the bias (chi2 between central PDF and underlying PDF)
    for each flavour and the variance (mean chi2 between replica and central PDF),
    then return a numpy array with shape (flavours, 2) with second axis being
    bias, variance

    """
    central_diff = np.asarray(fits_central_difference)
    rep_diff = np.asarray(fits_replica_difference)
    ratios = []
    for i in range(len(XI_FLAVOURS)):
        bias = calc_chi2(fits_sqrt_covmat_by_flavour[i], central_diff[:, i, :].T)
        variance = np.mean(
            calc_chi2(
                fits_sqrt_covmat_by_flavour[i],
                rep_diff[:, :, i, :].transpose(2, 1, 0),  # need x on first axis
            ),
            axis=0,
        )
        ratios.append(np.mean(bias) / np.mean(variance))
    return ratios


def fits_pdf_total_ratio(
    fits_central_difference,
    fits_replica_difference,
    fits_covariance_matrix_totalpdf,
    internal_Nx,
):
    """Calculate the total bias and variance for all flavours and x allowing for
    correlations across flavour.

    Returns:

    ratio_data: tuple
        required data for calculating mean(bias) over mean(variance) across fits
        in form of tuple (bias, variance)
    """
    central_diff = np.asarray(fits_central_difference).reshape(
        -1, internal_Nx * len(XI_FLAVOURS)
    )
    rep_diff = np.asarray(fits_replica_difference).reshape(
        len(fits_replica_difference), -1, internal_Nx * len(XI_FLAVOURS)
    )

    sqrtcov = la.cholesky(fits_covariance_matrix_totalpdf, lower=True)

    bias = calc_chi2(sqrtcov, central_diff.T)
    # need flav x on first axis
    variance = np.mean(calc_chi2(sqrtcov, rep_diff.transpose(2, 1, 0)), axis=0)
    return np.mean(bias) / np.mean(variance)


@table
def fits_pdf_bias_variance_ratio(fits_pdf_flavour_ratio, fits_pdf_total_ratio):
    """Returns a table with the values of mean bias / mean variance with mean
    referring to mean across fits, by flavour. Includes total across all
    flavours allowing for correlations.

    """
    records = []
    for i, fl in enumerate(XI_FLAVOURS):
        if i == 0:
            fl = f"${fl}$"
        records.append(dict(flavour=fl, ratio=fits_pdf_flavour_ratio[i]))
    records.append(dict(flavour="Total", ratio=fits_pdf_total_ratio))
    df = pd.DataFrame.from_records(
        records, index="flavour", columns=["flavour", "ratio"]
    )
    df.columns = ["bias/variance"]
    return df


@table
def fits_pdf_sqrt_ratio(fits_pdf_bias_variance_ratio):
    """Like :py:func:`fits_pdf_bias_variance_ratio` except taking the sqrt. To see how
    faithful our uncertainty is in units of the standard deviation.

    """
    df_in = fits_pdf_bias_variance_ratio
    data = np.sqrt(df_in.values)
    return pd.DataFrame(data, index=df_in.index, columns=["sqrt bias/variance"])


@table
def fits_pdf_expected_xi_from_ratio(fits_pdf_sqrt_ratio):
    """Like expected_xi_from_bias_variance but in PDF space. Estimate the
    integral across central difference distribution, with domain defined by the
    replica distribution. For more details see expected_xi_from_bias_variance.

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
    """Two column table comparing the measured value of xi for each flavour to the
    one calculated from bias/variance

    """
    return pd.concat((xi_flavour_table, fits_pdf_expected_xi_from_ratio), axis=1)


fits_xi_grid_values = collect("xi_grid_values", ("fits", "fitpdf"))


def bootstrap_pdf_differences(
    fits_xi_grid_values, underlying_xi_grid_values, multiclosure_underlyinglaw, rng
):
    """Generate a single bootstrap sample of ``pdf_central_difference`` and
    ``pdf_replica_difference`` given the multiclosure fits grid values
    (``fits_xi_grid_values``); the underlying law grid values and the underlying
    law; and a numpy random state which is used to generate random indices
    for bootstrap sample. The bootstrap does include repeats and has the same
    number of fits and replicas as the original ``fits_xi_grid_values`` which is
    being resampled.

    Returns
    -------

        pdf_difference: tuple
            a tuple of 2 lists: the central differences and the replica
            differences. Each list is n_fits long and each element is a resampled
            differences array for a randomly selected fit, randomly selected
            replicas.
    """
    fit_boot_index = rng.choice(len(fits_xi_grid_values), size=len(fits_xi_grid_values))
    boot_central_diff = []
    boot_rep_diff = []
    for i_fit in fit_boot_index:
        fit_xi_grid = fits_xi_grid_values[i_fit]
        rep_boot_index = rng.choice(fit_xi_grid.shape[0], size=fit_xi_grid.shape[0])
        xi_gv = fit_xi_grid[rep_boot_index, ...]
        boot_central_diff.append(
            pdf_central_difference(
                xi_gv, underlying_xi_grid_values, multiclosure_underlyinglaw
            )
        )
        boot_rep_diff.append(pdf_replica_difference(xi_gv))
    return boot_central_diff, boot_rep_diff


@table
def fits_bootstrap_pdf_xi_table(
    fits_xi_grid_values,
    underlying_xi_grid_values,
    multiclosure_underlyinglaw,
    internal_Nx,
    n_boot=100,
    boot_seed=None,
    use_x_basis=False,
):
    """Perform a bootstrap sampling across fits and replicas of xi, by flavour
    and total and then tabulate the mean and error

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
        total_cov = fits_covariance_matrix_totalpdf(boot_rep_diff, internal_Nx)
        xi_flav = xi_flavour_x(boot_rep_diff, boot_central_diff, flav_cov, use_x_basis)
        xi_total = xi_totalpdf(
            boot_rep_diff, boot_central_diff, total_cov, internal_Nx, use_x_basis
        )
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
        columns=[r"bootstrap mean $\xi_{1\sigma}$", r"bootstrap std. $\xi_{1\sigma}$"],
        index=index,
    )


def fits_bootstrap_pdf_ratio(
    fits_xi_grid_values,
    underlying_xi_grid_values,
    multiclosure_underlyinglaw,
    internal_Nx,
    n_boot=100,
    boot_seed=None,
):
    """Perform a bootstrap sampling across fits and replicas of the sqrt ratio,
    by flavour and total and then tabulate the mean and error

    """
    rng = np.random.RandomState(seed=boot_seed)
    ratio_boot = []
    for _ in range(n_boot):
        # perform single bootstrap
        boot_central_diff, boot_rep_diff = bootstrap_pdf_differences(
            fits_xi_grid_values,
            underlying_xi_grid_values,
            multiclosure_underlyinglaw,
            rng,
        )
        # need various dependencies for ratio actions
        flav_cov = fits_covariance_matrix_by_flavour(boot_rep_diff)
        flav_sqrt_cov = fits_sqrt_covmat_by_flavour(flav_cov)
        total_cov = fits_covariance_matrix_totalpdf(boot_rep_diff, internal_Nx)
        ratio_flav = fits_pdf_flavour_ratio(
            flav_sqrt_cov, boot_central_diff, boot_rep_diff
        )
        ratio_tot = fits_pdf_total_ratio(
            boot_central_diff, boot_rep_diff, total_cov, internal_Nx
        )
        ratio_data = np.concatenate((ratio_flav, [ratio_tot]), axis=0)
        ratio_boot.append(ratio_data)
    return ratio_boot


def fits_bootstrap_pdf_sqrt_ratio(fits_bootstrap_pdf_ratio):
    """Take the square root of fits_bootstrap_pdf_ratio"""
    return np.sqrt(fits_bootstrap_pdf_ratio)


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
        columns=[r"bootstrap mean sqrt ratio", r"bootstrap std. sqrt ratio"],
        index=index,
    )


def fits_bootstrap_pdf_expected_xi(fits_bootstrap_pdf_sqrt_ratio):
    """Using fits_bootstrap_pdf_sqrt_ratio calculate a bootstrap of the expected
    xi using the same procedure as in expected_xi_from_bias_variance

    """
    n_sigma_in_variance = 1 / fits_bootstrap_pdf_sqrt_ratio
    # pylint can't find erf here, disable error in this function
    # pylint: disable=no-member
    estimated_integral = scipy.special.erf(n_sigma_in_variance / np.sqrt(2))
    return estimated_integral


@table
def fits_bootstrap_pdf_expected_xi_table(fits_bootstrap_pdf_expected_xi):
    """Tabulate the mean and standard deviation across bootstrap samples
    of :py:func:`fits_bootstrap_pdf_expected_xi` with a row for each flavour and the
    total expected xi

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
    of the measured value of xi to the one calculated from bias/variance
    in PDF space, for each flavour and total across all flavours accounting
    for correlations

    """
    return fits_pdf_compare_xi_to_expected(
        fits_bootstrap_pdf_expected_xi_table, fits_bootstrap_pdf_xi_table
    )


def plot_pdf_matrix(matrix, n_x, **kwargs):
    """Utility function which, given a covmat/corrmat for all flavours and
    x, plots it with appropriate labels. Input matrix is expected to be
    size (n_flavours*n_x) * (n_flavours*n_x)

    Parameters:

    matrix: np.array
        square matrix which must be (n_flavours*n_x) * (n_flavours*n_x) with
        elements ordered like:
        (flavour0_x0, flavour0_x1, ..., flavourN_x0, ..., flavourN_xN)
        i.e the points along x for flavour 0, then points along x for flavour 1
        etc.
    **kwargs:
        keyword arguments for the matplotlib.axes.Axes.imshow function

    Notes:

    See matplotlib.axes.Axes.imshow for more details on plotting function

    """
    labels = [f"${XI_FLAVOURS[0]}$", *XI_FLAVOURS[1:]]
    # we want to centre the labels on each of the xgrids
    # ticks appear shifted right by 0.5, so we account for this here
    ticks = (n_x) / 2 + n_x * np.arange(len(XI_FLAVOURS)) - 0.5
    fig, ax = plt.subplots()
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
def plot_multiclosure_covariance_matrix(fits_covariance_matrix_totalpdf, internal_Nx):
    """Plot the covariance matrix for all flavours, covariance matrix
    has shape n_flavours * n_x each block is the covariance of the replica
    PDFs on the x-grid defined in xi_pdfgrids

    """
    fig, ax = plot_pdf_matrix(fits_covariance_matrix_totalpdf, internal_Nx)
    ax.set_title("Covariance matrix estimated from multiclosure replicas")
    return fig


def fits_correlation_matrix_totalpdf(fits_covariance_matrix_totalpdf):
    """Given the fits_covariance_matrix_totalpdf, returns the corresponding
    correlation matrix

    """
    d = np.sqrt(np.diag(fits_covariance_matrix_totalpdf))
    return (fits_covariance_matrix_totalpdf / d) / d[:, np.newaxis]


@figure
def plot_multiclosure_correlation_matrix(fits_correlation_matrix_totalpdf, internal_Nx):
    """Like plot_multiclosure_covariance_matrix but plots the total
    correlation matrix

    """
    fig, ax = plot_pdf_matrix(
        fits_correlation_matrix_totalpdf, internal_Nx, vmin=0, vmax=1
    )
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
    fig, ax = plt.subplots()
    ax.plot(
        e_val,
        "*",
        label=f"Eigenvalues of correlation matrix, l2-condition: {l2_condition:.2f}",
    )
    ax.set_ylabel(r"$\lambda_{\rm corr}$")
    ax.set_xlabel("eigenvalue index (ascending)")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend()
    return fig
