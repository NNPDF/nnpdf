"""
multiclosure_preprocessing.py

Module containing all of the actions related to preprocessing exponents. In
particular, comparing the next preprocessing exponents across the multiple
closure fits with the previous effective exponents, to see if there is a big
dependence on the level 1 shift.

"""
from matplotlib.figure import Figure
import pandas as pd

from reportengine import collect
from reportengine.figure import figuregen
from reportengine.table import table

from validphys.closuretest.closure_checks import (
    check_fits_areclosures,
    check_fits_have_same_basis,
    check_fits_underlying_law_match,
)
from validphys.plotutils import plot_horizontal_errorbars


def _next_multiclosure_preprocessing_table(fits_pdf, fits_preprocessing_table):
    """Internal function which tabulates the next preprocessing range of all
    fits with the pdf name as the index. The columns have a ``MultiIndex`` of
    flavour and next exponent range minimum/maxmimum.

    Parameters
    ----------
    fits_pdf: list[validphys.core.FitSpec]
        Multiclosure fit pdfs, used to build the table index
    fits_preprocessing_table: list[pd.DataFrame]
        DataFrames containing the next exponent ranges (minimum and maximum)
        for each flavour for either alpha or beta exponents.

    Returns
    -------
    multiclosure_preprocessing_table: pd.DataFrame
        Containing the next preprocessing ranges for each fit and flavour.

    """
    dfs = []
    for fit_pdf, pp_table in zip(fits_pdf, fits_preprocessing_table):
        next_pp_tab = pp_table.loc[:, f"next ({fit_pdf.label})"]
        # flatten table by creating multiindex
        dfs.append(next_pp_tab.stack().to_frame(name=fit_pdf.label))
    # transpose concatenated table to put fits on index.
    return pd.concat(dfs, axis=1).T


fits_fitbasis_alpha_lines = collect("get_alpha_lines", ("fits", "fitpdfandbasis"))
fits_basis = collect("basis", ("fits", "basisfromfit"))


@check_fits_areclosures
@check_fits_have_same_basis
@check_fits_underlying_law_match
@table
def next_multiclosure_alpha_preprocessing_table(
    fits, fits_basis, fits_pdf, fits_fitbasis_alpha_lines
):
    """Returns a table with the next alpha preprocessing exponent for each fit
    with a multiindex column of flavour and next preprocessing range limits.

    For more information see :py:func:`_next_multiclosure_preprocessing_table`
    """
    return _next_multiclosure_preprocessing_table(fits_pdf, fits_fitbasis_alpha_lines)


fits_fitbasis_beta_lines = collect("get_beta_lines", ("fits", "fitpdfandbasis"))


@check_fits_areclosures
@check_fits_have_same_basis
@check_fits_underlying_law_match
@table
def next_multiclosure_beta_preprocessing_table(
    fits, fits_basis, fits_pdf, fits_fitbasis_beta_lines
):
    """Returns a table with the next beta preprocessing exponent for each fit
    with a multiindex column of flavour and next preprocessing range limits.

    For more information see :py:func:`_next_multiclosure_preprocessing_table`
    """
    return _next_multiclosure_preprocessing_table(fits_pdf, fits_fitbasis_beta_lines)


@figuregen
def plot_next_multiclosure_alpha_preprocessing(
    fits_fitbasis_alpha_lines, fits_pdf, next_multiclosure_alpha_preprocessing_table,
):
    """Using the table produced by
    :py:func:`next_multiclosure_alpha_preprocessing_table`, plot the next
    alpha preprocessing exponent ranges. The ranges are represented by
    horizontal error bars, with vertical lines indicating the previous range
    limits of the first fit.

    """
    first_prev_ranges = fits_fitbasis_alpha_lines[0].loc[
        :, f"prev ({fits_pdf[0].label})"
    ]
    flavours = next_multiclosure_alpha_preprocessing_table.columns.get_level_values(
        0
    ).unique()
    for flavour in flavours:
        next_flavour_range = next_multiclosure_alpha_preprocessing_table.loc[:, flavour]
        next_flavour_range_vals = next_flavour_range.to_numpy()
        half_diff = (next_flavour_range_vals[:, 1] - next_flavour_range_vals[:, 0]) / 2
        cvs = (next_flavour_range_vals[:, 1] + next_flavour_range_vals[:, 0]) / 2

        fig, ax = plot_horizontal_errorbars(
            [cvs],
            [half_diff],
            next_flavour_range.index.to_numpy(),
            ["Next preprocessing ranges"],
        )
        xlims = ax.get_xlim()
        prev_flavour_range = first_prev_ranges.loc[flavour, :].to_numpy().squeeze()
        ax.vlines(
            prev_flavour_range,
            *ax.get_ylim(),
            linestyle=":",
            color="k",
            label="Previous limits",
            clip_on=False,
        )
        # use convenience method to set xlims in case vlines are off plot.
        ax.autoscale(True, "x")
        # remove markers since the CV is meaningless.
        handles, labels = ax.get_legend_handles_labels()
        _, (errorbar_points, _, _) = handles
        errorbar_points.set_marker(None)
        ax.legend(handles, labels)
        ax.set_title(f"Multiclosure fits {flavour} alpha preprocessing exponents.")
        yield fig


@figuregen
def plot_next_multiclosure_beta_preprocessing(
    fits_fitbasis_beta_lines, fits_pdf, next_multiclosure_beta_preprocessing_table,
):
    """Using the table produced by
    :py:func:`next_multiclosure_beta_preprocessing_table`, plot the next
    beta preprocessing exponent ranges. The ranges are represented by
    horizontal error bars, with vertical lines indicating the previous range
    limits of the first fit.

    """
    for fig in plot_next_multiclosure_alpha_preprocessing(
        fits_fitbasis_beta_lines, fits_pdf, next_multiclosure_beta_preprocessing_table,
    ):
        # fixup title.
        ax = fig.gca()
        new_title = ax.get_title().replace("alpha", "beta")
        ax.set_title(new_title)
        yield fig


@figuregen
def plot_next_multiclosure_alpha_preprocessing_range_width(
    fits_fitbasis_alpha_lines, fits_pdf, next_multiclosure_alpha_preprocessing_table,
):
    """Using the table produced by
    :py:func:`next_multiclosure_alpha_preprocessing_table`, plot the next
    alpha preprocessing exponent ranges width, aka max alpha - min alpha
    as a histogram over fits for each flavour. Add a vertical line of the
    previous range width of
    the first fit for reference

    """
    first_prev_ranges = fits_fitbasis_alpha_lines[0].loc[
        :, f"prev ({fits_pdf[0].label})"
    ]
    flavours = next_multiclosure_alpha_preprocessing_table.columns.get_level_values(
        0
    ).unique()
    for flavour in flavours:
        next_flavour_range = next_multiclosure_alpha_preprocessing_table.loc[:, flavour]
        next_flavour_range_vals = next_flavour_range.to_numpy()
        diffs = next_flavour_range_vals[:, 1] - next_flavour_range_vals[:, 0]
        fig = Figure()
        ax = fig.subplots()
        ax.hist(diffs, label="Next range width")
        prev_lims = first_prev_ranges.loc[flavour].to_numpy().squeeze()
        ax.axvline(
            prev_lims[1] - prev_lims[0],
            linestyle=":",
            color="k",
            label="Previous range width.",
        )
        ax.set_title(
            f"Multiclosure fits {flavour} alpha preprocessing exponents range width."
        )
        ax.legend()
        yield fig


@figuregen
def plot_next_multiclosure_beta_preprocessing_range_width(
    fits_fitbasis_beta_lines, fits_pdf, next_multiclosure_beta_preprocessing_table,
):
    """Using the table produced by
    :py:func:`next_multiclosure_beta_preprocessing_table`, plot the next
    beta preprocessing exponent ranges width, aka max beta - min beta
    as a histogram over fits for each flavour. Add a vertical line of the
    previous range width of
    the first fit for reference

    """
    for fig in plot_next_multiclosure_alpha_preprocessing_range_width(
        fits_fitbasis_beta_lines, fits_pdf, next_multiclosure_beta_preprocessing_table,
    ):
        # fixup title.
        ax = fig.gca()
        new_title = ax.get_title().replace("alpha", "beta")
        ax.set_title(new_title)
        yield fig
