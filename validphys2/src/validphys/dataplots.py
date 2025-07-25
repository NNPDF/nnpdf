"""
Plots of relations between data PDFs and fits.
"""

from __future__ import generator_stop

from collections import defaultdict
from collections.abc import Sequence
import itertools
import logging

import matplotlib as mpl
from matplotlib import cm
from matplotlib import colors as mcolors
from matplotlib import ticker as mticker
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
import scipy.stats as stats

from nnpdf_data.coredata import KIN_NAMES
from reportengine import collect
from reportengine.checks import CheckError, check, make_argcheck, make_check
from reportengine.figure import figure, figuregen
from reportengine.floatformatting import format_number
from validphys import plotutils
from validphys.checks import check_not_using_pdferr
from validphys.core import CutsPolicy, MCStats, cut_mask
from validphys.plotoptions.core import get_info, kitable, transform_result
from validphys.results import chi2_stat_labels, chi2_stats
from validphys.sumrules import POL_LIMS, partial_polarized_sum_rules
from validphys.utils import sane_groupby_iter, scale_from_grid, split_ranges

log = logging.getLogger(__name__)


@figure
def plot_chi2dist_experiments(total_chi2_data, experiments_chi2_stats, pdf):
    """Plot the distribution of chi²s of the members of the pdfset."""
    fig, ax = _chi2_distribution_plots(total_chi2_data, experiments_chi2_stats, pdf, "hist")
    ax.set_title(r"Experiments $\chi^2$ distribution")
    return fig


@figure
def kde_chi2dist_experiments(total_chi2_data, experiments_chi2_stats, pdf):
    """KDE plot for experiments chi2."""
    fig, ax = _chi2_distribution_plots(total_chi2_data, experiments_chi2_stats, pdf, "kde")
    ax.set_ylabel(r"Density")
    ax.set_title(r"Experiments $\chi^2 KDE plot$")
    return fig


@figure
def plot_chi2dist(dataset, abs_chi2_data, chi2_stats, pdf):
    """Plot the distribution of chi²s of the members of the pdfset."""
    setlabel = dataset.name
    fig, ax = _chi2_distribution_plots(abs_chi2_data, chi2_stats, pdf, "hist")
    ax.set_title(r"$\chi^2$ distribution for %s" % setlabel)
    return fig


@figure
def plot_chi2dist_sv(dataset, abs_chi2_data_thcovmat, pdf):
    """Same as ``plot_chi2dist`` considering also the theory covmat in the calculation"""
    chi2_stats_thcovmat = chi2_stats(abs_chi2_data_thcovmat)
    return plot_chi2dist(dataset, abs_chi2_data_thcovmat, chi2_stats_thcovmat, pdf)


def _chi2_distribution_plots(chi2_data, stats, pdf, plot_type):
    fig, ax = plotutils.subplots()
    label = pdf.name
    alldata, central, npoints = chi2_data
    if not isinstance(alldata, MCStats):
        ax.set_facecolor("#ffcccc")
        log.warning("Chi² distribution plots have a " "different meaning for non MC sets.")
        label += " (%s!)" % pdf.error_type
    label += '\n' + '\n'.join(str(chi2_stat_labels[k]) + (' %.2f' % v) for (k, v) in stats.items())
    ax.set_xlabel(r"Replica $\chi^2$")

    if plot_type == "hist":
        ax.hist(alldata.error_members(), label=label, zorder=100)
    elif plot_type == "kde":
        # We need the squeeze here to change shape from (x, 1) to (x,)
        ax = plotutils.kde_plot(alldata.data.squeeze(), label=label)
    else:
        raise ValueError(f"plot_type must either be hist or kde, not {plot_type}")

    l = ax.legend()
    l.set_zorder(1000)
    return fig, ax


@figure
def plot_phi(groups_data, groups_data_phi, processed_metadata_group):
    """plots phi for each group of data as a bar for a single
    PDF input

    See `phi_data` for information on how phi is calculated
    """
    phi = [exp_phi for (exp_phi, npoints) in groups_data_phi]
    xticks = [group.name for group in groups_data]
    fig, ax = plotutils.barplot(phi, collabels=xticks, datalabels=[r'$\phi$'])
    ax.set_title(rf"$\phi$ by {processed_metadata_group}")
    return fig


@figure
def plot_fits_groups_data_phi(fits_groups_phi_table, processed_metadata_group):
    """Plots a set of bars for each fit, each bar represents the value of phi for the corresponding
    group of datasets, which is defined according to the keys in the PLOTTING info file"""
    fig, ax = _plot_chis_df(fits_groups_phi_table)
    ax.set_title(rf"$\phi$ by {processed_metadata_group}")
    return fig


@figure
def plot_dataset_inputs_phi_dist(data, dataset_inputs_bootstrap_phi_data):
    """Generates a bootstrap distribution of phi and then plots a histogram
    of the individual bootstrap samples for `dataset_inputs`. By default
    the number of bootstrap samples is set to a sensible number (500) however
    this number can be changed by specifying `bootstrap_samples` in the runcard
    """
    phi = dataset_inputs_bootstrap_phi_data
    label = '\n'.join(
        [
            fr'$\phi$ mean = {format_number(phi.mean())}',
            fr'$\phi$ std dev = {format_number(phi.std())}',
        ]
    )
    fig, ax = plotutils.subplots()
    ax.hist(phi, label=label)
    ax.set_title(r"$\phi$ distribution for " + data.name)
    ax.legend()
    return fig


@make_argcheck
def _check_same_group_data_name(dataspecs_groups):
    lst = dataspecs_groups
    if not lst:
        return
    for j, x in enumerate(lst[1:]):
        if len(x) != len(lst[0]):
            raise CheckError("All dataspecs should have the same number of groups of data")
        for i, group in enumerate(x):
            if group.name != lst[0][i].name:
                raise CheckError(
                    "\n".join(
                        [
                            "All groups of data must have the same name",
                            fr"dataspec {j+1}, group {i+1}: {group.name}",
                            fr"dataspec 1, group {i+1}: {lst[0][i].name}",
                        ]
                    )
                )


@_check_same_group_data_name
@figure
def plot_phi_scatter_dataspecs(
    dataspecs_groups, dataspecs_speclabel, dataspecs_groups_bootstrap_phi
):
    """For each of the dataspecs, a bootstrap distribution of phi is generated
    for all specified groups of datasets. The distribution is then represented as a
    scatter point which is the median of the bootstrap distribution and an
    errorbar which spans the 68% confidence interval. By default the number
    of bootstrap samples is set to a sensible value, however it can be
    controlled by specifying `bootstrap_samples` in the runcard.
    """
    labels = dataspecs_speclabel
    phis = dataspecs_groups_bootstrap_phi
    exps = dataspecs_groups
    xticks = [group.name for group in exps[0]]
    x = range(1, len(xticks) + 1)
    fig, ax = plotutils.subplots()
    phi_stats = np.percentile(phis, [16, 50, 84], axis=2)
    for i, label in enumerate(labels):
        phi_errs = np.vstack(
            (phi_stats[2, i, :] - phi_stats[1, i, :], phi_stats[1, i, :] - phi_stats[0, i, :])
        )
        ax.errorbar(x, phi_stats[1, i, :], yerr=phi_errs, fmt='.', label=label)
    ax.set_xticks(x, minor=False)
    ax.set_xticklabels(xticks, minor=False, rotation=45)
    ax.legend()
    return fig


# TODO: This should be simplified if at all possible. For now some more examples
# are needed for a spec to emerge.
@make_check
def check_normalize_to(ns, **kwargs):
    """Transforn normalize_to into an index."""

    msg = (
        "normalize_to should be either 'data', a pdf id or an index of the "
        "result (0 for the data, and i for the ith pdf)"
    )

    val = ns.get('normalize_to', None)
    if val is None:
        return

    if 'pdf' in ns:
        names = ['data', ns['pdf'].name]
    else:
        names = ['data', *(pdf.name for pdf in ns['pdfs'])]
    if isinstance(val, int):
        if not val < len(names):
            raise CheckError(msg)
        return

    if isinstance(val, str):
        try:
            val = names.index(val)
        except ValueError:
            raise CheckError(msg, val, alternatives=names)
        ns['normalize_to'] = val
        return

    raise RuntimeError("Should not be here")


# TODO: This interface is horrible.
# We need to think how to adapt it to make this use case easier
def _plot_fancy_impl(
    results, commondata, cutlist, normalize_to: (int, type(None)) = None, labellist=None
):
    """Implementation of the data-theory comparison plots. Providers are
    supposed to call (yield from) this.
    Parameters
    -----------
    results : list
        A list of results, where the first one is a data result and the
        subsequent ones are theory predictions.
    commondata : ``CommonDataSpec``
        The specification corresponfing to the commondata to be plotted.
    cutlist : list
        The list of ``CutSpecs`` or ``None`` corresponding to the cuts for each
        result.
    normalize_to : int or None
        The index of the result to which ratios will be computed. If ``None``,
        plot absolute values.
    labellist : list or None
        The labesl that will appear in the plot. They sill be deduced
        (from the PDF names) if None is given.
    Returns
    -------
    A generator over figures.
    """

    info = get_info(commondata, normalize=(normalize_to is not None))

    table = kitable(commondata, info)
    nkinlabels = len(table.columns)
    ndata = len(table)

    # This is easier than cheking every time
    if labellist is None:
        labellist = [None] * len(results)

    if normalize_to is not None:
        norm_result = results[normalize_to]
        mask = cut_mask(cutlist[normalize_to])
        cv = np.full(ndata, np.nan)
        cv[mask] = norm_result.central_value

        err = np.full(ndata, np.nan)
        err[mask] = norm_result.std_error
        # We modify the table, so we pass only the label columns
        norm_cv, _ = transform_result(cv, err, table.iloc[:, :nkinlabels], info)

    cvcols = []
    for i, (result, cuts) in enumerate(zip(results, cutlist)):
        # We modify the table, so we pass only the label columns
        mask = cut_mask(cuts)
        cv = np.full(ndata, np.nan)
        cv[mask] = result.central_value
        err = np.full(ndata, np.nan)
        err[mask] = result.std_error

        cv, err = transform_result(cv, err, table.iloc[:, :nkinlabels], info)

        # By doing tuple keys we avoid all possible name collisions
        cvcol = ('cv', i)
        if normalize_to is None:
            table[cvcol] = cv
            table[('err', i)] = err
        else:
            table[cvcol] = cv / norm_cv
            table[('err', i)] = np.abs(err / norm_cv)
        cvcols.append(cvcol)

    figby = sane_groupby_iter(table, info.figure_by)

    for samefig_vals, fig_data in figby:
        # Nothing to plot if all data is cut away
        if np.all(np.isnan(fig_data[cvcols])):
            continue
        # For some reason matplotlib doesn't set the axis right
        min_vals = []
        max_vals = []
        fig, ax = plotutils.subplots()
        ax.set_title(
            "{} {}".format(info.dataset_label, info.group_label(samefig_vals, info.figure_by))
        )

        lineby = sane_groupby_iter(fig_data, info.line_by)

        first = True

        for sameline_vals, line_data in lineby:
            ax.set_prop_cycle(None)
            labels = first
            first = False

            offset_iter = plotutils.offset_xcentered(len(results), ax)

            x = info.get_xcol(line_data)

            try:
                x = np.asanyarray(x, float)
            except ValueError:
                xticklabels = x
                npoints = len(x)
                x = np.arange(npoints)
                ax.set_xticks(x)
                ax.set_xticklabels(xticklabels)
                # TODO: Remove this when mpl stops doing the wrong thing
                # (in v2?)
                ax.set_xlim(-npoints / 20, npoints - 1 + npoints / 20)

            # Use black for the first iteration (data),
            # and follow the cycle for
            # the rest.
            next_color = itertools.chain(['#262626'], plotutils.color_iter())

            for i, (res, lb, color) in enumerate(zip(results, labellist, next_color)):
                if labels:
                    if lb:
                        label = lb
                    else:
                        label = res.label
                else:
                    label = None

                cv = line_data[('cv', i)].values
                err = line_data[('err', i)].values
                ax.errorbar(
                    x,
                    cv,
                    yerr=err,
                    lw=0.25,
                    label=label,
                    # elinewidth = 2,
                    capsize=2,
                    marker='s',
                    markeredgewidth=0.25,
                    c=color,
                    zorder=1000,
                    transform=next(offset_iter),
                )

                # We 'plot' the empty lines to get the labels. But
                # if everything is empty we skip the plot.
                if np.any(np.isfinite(cv)):
                    max_vals.append(np.nanmax(cv + err))
                    min_vals.append(np.nanmin(cv - err))

            glabel = info.group_label(sameline_vals, info.line_by)

            # Use some anchor that is not in y=1 for ratio plots
            if normalize_to is not None:
                next_after_normalize = (normalize_to + 1) % len(results)
                annotate_point = x[-1], line_data[('cv', next_after_normalize)].values[-1]
            else:
                annotate_point = x[-1], line_data[('cv', 0)].values[-1]
            # This is a workaround for https://github.com/matplotlib/matplotlib/issues/12648
            if np.isfinite(annotate_point).all():
                ax.annotate(
                    glabel,
                    annotate_point,
                    xytext=(15, -10),
                    size='xx-small',
                    textcoords='offset points',
                    zorder=10000,
                )

        if info.x_scale:
            ax.set_xscale(info.x_scale)

        if normalize_to is None:
            if info.y_label:
                ax.set_ylabel(info.y_label)
            if info.y_scale:
                ax.set_yscale(info.y_scale)
        else:
            # Rebuild the limits of the plot looking at the max and min values
            margin_fraction = 0.05
            ymax = max(max_vals)
            ymin = min(min_vals)
            margin = (ymax - ymin) * margin_fraction
            ax.set_ylim((ymin - margin, ymax + margin))

            # Normalized plots should always be linear in the y axis
            ax.set_yscale("linear")

            lb = labellist[normalize_to]
            ax.set_ylabel(f"Ratio to {lb if lb else norm_result.label}")

        ax.legend().set_zorder(100000)
        ax.set_xlabel(info.xlabel)
        fig.tight_layout()
        yield fig


@check_not_using_pdferr
@check_normalize_to
@figuregen
def plot_fancy(
    one_or_more_results,
    commondata,
    cuts,
    normalize_to: (int, str, type(None)) = None,
    use_pdferr: bool = False,  # pylint: disable=unused-argument # for checks
):
    """
    Read the PLOTTING configuration for the dataset and generate the
    corrspondig data theory plot.

    The input results are assumed to be such that the first one is the data,
    and the subsequent ones are the predictions for the PDFfs. See
    ``one_or_more_results``. The labelling of the predictions can be
    influenced by setting ``label`` attribute of theories and pdfs.

    normalize_to: should be either 'data', a pdf id or an index of the
    result (0 for the data, and i for the ith pdf). None means plotting
    absolute values.


    See docs/plotting_format.md for details on the format of the PLOTTING
    files.
    """
    yield from _plot_fancy_impl(
        results=one_or_more_results,
        commondata=commondata,
        cutlist=[cuts] * len(one_or_more_results),
        normalize_to=normalize_to,
    )


@make_argcheck
def _check_same_dataset_name(dataspecs_commondata):
    lst = dataspecs_commondata
    if not lst:
        return
    ele = lst[0].name
    for x in lst[1:]:
        if x.name != ele:
            raise CheckError("All datasets must have the same name")


@make_argcheck
def _check_dataspec_normalize_to(normalize_to, dataspecs):
    if normalize_to in (0, None) or (
        isinstance(normalize_to, int) and normalize_to <= len(dataspecs)
    ):
        return
    if normalize_to == 'data':
        return {'normalize_to': 0}

    raise CheckError(
        "Unrecignized format for normalize_to. Must be either "
        "'data', 0 or the 1-indexed index of the dataspec "
        f"(<{len(dataspecs)}), not {normalize_to}"
    )


@check_not_using_pdferr
@_check_same_dataset_name
@_check_dataspec_normalize_to
@figuregen
def plot_fancy_dataspecs(
    dataspecs_results,
    dataspecs_commondata,
    dataspecs_cuts,
    dataspecs_speclabel,
    normalize_to: (str, int, type(None)) = None,
    use_pdferr: bool = False,  # pylint: disable=unused-argument # for checks
):
    """
    General interface for data-theory comparison plots.

    The user should define an arbitrary list of mappings called "dataspecs".
    In each of these, ``dataset`` must resolve to a dataset with the same name
    (but could be e.g. different theories). The production rule
    ``matched_datasets_from_datasepcs`` may be used for this purpose.

    The result will be a plot combining all the predictions from the dataspecs
    mapping (whch could vary in theory, pdf, cuts, etc).

    The user can define a "speclabel" key in each datasspec (or only on some).
    By default, the PDF label will be used in the legend (like in
    ``plot_fancy``).

    ``normalize_to must`` be either:

        - The string 'data' or the integer 0 to plot the ratio to data,

        - or the 1-based index of the dataspec to normalize to the
          corresponding prediction,


        - or None (default) to plot absolute values.

    A limitation at the moment is that the data cuts and errors will be taken
    from the first specifiaction.
    """
    # We have at least one element
    if not dataspecs_results:
        return

    # For now, simply take the first data result. We'll need to improve this.
    results = [dataspecs_results[0][0], *[r[1] for r in dataspecs_results]]
    cutlist = [dataspecs_cuts[0], *dataspecs_cuts]
    commondata = dataspecs_commondata[0]
    labellist = [None, *dataspecs_speclabel]

    yield from _plot_fancy_impl(
        results=results,
        commondata=commondata,
        cutlist=cutlist,
        labellist=labellist,
        normalize_to=normalize_to,
    )


@_check_same_dataset_name
@_check_dataspec_normalize_to
@figuregen
def plot_fancy_sv_dataspecs(
    dataspecs_results_with_scale_variations,
    dataspecs_commondata,
    dataspecs_cuts,
    dataspecs_speclabel,
    normalize_to: (str, int, type(None)) = None,
):
    """
    Exactly the same as ``plot_fancy_dataspecs`` but the theoretical results passed down
    are modified so that the 1-sigma error bands correspond to a combination of the
    PDF error and the scale variations collected over theoryids

    See: :py:func:`validphys.results.results_with_scale_variations`
    """
    return plot_fancy_dataspecs(
        dataspecs_results_with_scale_variations,
        dataspecs_commondata,
        dataspecs_cuts,
        dataspecs_speclabel,
        normalize_to=normalize_to,
    )


def _scatter_marked(ax, x, y, marked_dict, *args, **kwargs):
    kwargs['s'] = kwargs.get('s', 30) + 10
    x = np.array(x, copy=False)
    y = np.array(y, copy=False)
    for label, indexes in marked_dict.items():
        ax.scatter(
            x[indexes],
            y[indexes],
            *args,
            **kwargs,
            label=label,
            facecolors='none',
            linewidth=0.5,
            edgecolor='red',
        )
        kwargs['s'] += 10


@figure
def plot_dataspecs_groups_chi2_spider(dataspecs_groups_chi2_table):
    fig, ax = _plot_chi2s_spider_df(dataspecs_groups_chi2_table)
    return fig


@figure
def plot_fits_chi2_spider(fits, fits_groups_chi2, fits_groups_data, processed_metadata_group):
    """Plots the chi²s of all groups of datasets
    on a spider/radar diagram."""

    fig, ax = plotutils.add_subplot(figsize=(12, 12), projection='polar')

    for fit, fitchi2, fitgroup in zip(fits, fits_groups_chi2, fits_groups_data):
        exchi2 = [group_res.central_result / group_res.ndata for group_res in fitchi2]
        xticks = [group.name for group in fitgroup]

        ax = plotutils.spiderplot(xticks, exchi2, fit, ax)

    ax.set_title(rf"$\chi^2$ by {processed_metadata_group}")

    return fig


@figure
def plot_fits_phi_spider(fits, fits_groups_data, fits_groups_data_phi, processed_metadata_group):
    """Like plot_fits_chi2_spider but for phi."""

    fig, ax = plotutils.add_subplot(figsize=(12, 12), projection='polar')

    for fit, fitphi, fitgroup in zip(fits, fits_groups_data_phi, fits_groups_data):
        phi = [exp_phi for (exp_phi, _npoints) in fitphi]
        xticks = [group.name for group in fitgroup]

        ax = plotutils.spiderplot(xticks, phi, fit, ax)

    ax.set_title(rf"$\phi$ by {processed_metadata_group}")

    return fig


@figure
def plot_groups_data_chi2_spider(groups_data, groups_chi2, processed_metadata_group, pdf):
    """Plot the chi² of all groups of datasets as a spider plot."""

    exchi2 = [group_res.central_result / group_res.ndata for group_res in groups_chi2]
    xticks = [group.name for group in groups_data]

    fig, ax = plotutils.add_subplot(projection='polar')

    ax = plotutils.spiderplot(xticks, exchi2, pdf, ax)
    ax.set_title(rf"$\chi^2$ by {processed_metadata_group}")
    return fig


@figure
def plot_groups_data_phi_spider(groups_data, groups_data_phi, processed_metadata_group, pdf):
    """Plot the phi of all groups of datasets as a spider plot."""
    phi = [exp_phi for (exp_phi, _npoints) in groups_data_phi]
    xticks = [group.name for group in groups_data]
    fig, ax = plotutils.add_subplot(projection='polar')

    ax = plotutils.spiderplot(xticks, phi, pdf, ax)
    ax.set_title(rf"$\phi$ by {processed_metadata_group}")
    return fig


@figure
def plot_groups_data_chi2(groups_data, groups_chi2, processed_metadata_group):
    """Plot the chi² of all groups of datasets with bars."""
    exchi2 = [group_res.central_result / group_res.ndata for group_res in groups_chi2]
    xticks = [group.name for group in groups_data]

    fig, ax = plotutils.barplot(exchi2, collabels=xticks, datalabels=[r'$\chi^2$'])
    ax.set_title(rf"$\chi^2$ by {processed_metadata_group}")
    return fig


plot_experiments_chi2 = collect("plot_groups_data_chi2", ("group_dataset_inputs_by_experiment",))


@figure
def plot_datasets_chi2(groups_data, groups_chi2):
    """Plot the chi² of all datasets with bars."""
    dschi2 = []
    xticks = []
    for group, group_res in zip(groups_data, groups_chi2):
        xticks = [dataset.name for dataset in group]
        dschi2 = [dsres.central_result / dsres.ndata for dsres in group_res]
    fig, ax = plotutils.barplot(dschi2, collabels=xticks, datalabels=[r'$\chi^2$'])

    ax.set_title(r"$\chi^2$ distribution for datasets")

    return fig


each_dataset_chi2_pdfs = collect("each_dataset_chi2", ("pdfs",))


@figure
def plot_datasets_pdfs_chi2(data, each_dataset_chi2_pdfs, pdfs):
    """
    Plot the chi² of all datasets with bars, and for different
    pdfs.
    """

    chi2_pdfs = list(each_dataset_chi2_pdfs)

    pdf_dict = {
        dataset.name: [chi2_pdfs[i][j] for i in range(len(chi2_pdfs))]
        for j, dataset in enumerate(data)
    }

    vals = []
    collabels = []

    for ds, val in pdf_dict.items():
        vals.append([chi2.central_result / chi2.ndata for chi2 in val])
        collabels.append(ds)

    fig, ax = plotutils.barplot(
        np.array(vals).T, collabels, datalabels=[rf"$\chi^2$, {str(pdf)}" for pdf in pdfs]
    )
    ax.set_title(r"$\chi^2$ distribution for datasets")
    ax.legend()
    return fig


each_dataset_chi2_sv = collect("abs_chi2_data_thcovmat", ("data",))
each_dataset_chi2_pdfs_sv = collect("each_dataset_chi2_sv", ("pdfs",))


@figure
def plot_datasets_pdfs_chi2_sv(data, each_dataset_chi2_pdfs_sv, pdfs):
    """Same as ``plot_datasets_pdfs_chi2_sv`` with the chi²s computed including scale variations"""
    return plot_datasets_pdfs_chi2(data, each_dataset_chi2_pdfs_sv, pdfs)


@figure
def plot_datasets_chi2_spider(groups_data, groups_chi2):
    """Plot the chi² of all datasets with bars."""
    dschi2 = []
    xticks = []
    for group, group_res in zip(groups_data, groups_chi2):
        xticks = [dataset.name for dataset in group]
        dschi2 = [dsres.central_result / dsres.ndata for dsres in group_res]

    fig, ax = plotutils.add_subplot(figsize=(4, 4), projection='polar')
    ax = plotutils.spiderplot(xticks, dschi2, label=[r'$\chi^2$'], ax=ax)

    ax.set_title(r"$\chi^2$ distribution for datasets")

    return fig


def _plot_chis_df(df):
    """Takes a dataframe that is a reduced version of ``fits_dataset_chi2s_table``
    and returns a bar plot. See ``plot_fits_datasets_chi2`` for use"""
    chilabel = df.columns.get_level_values(1)[1]
    data = df.iloc[:, df.columns.get_level_values(1) == chilabel].T.values
    fitnames = df.columns.get_level_values(0).unique()
    expnames = list(df.index.get_level_values(0))
    fig, ax = plotutils.barplot(data, expnames, fitnames)
    ax.grid(False)
    ax.legend()
    return fig, ax


def _plot_chi2s_spider_df(df, size=6):
    """Like _plot_chis_df but for spider plot."""
    chilabel = df.columns.get_level_values(1)[1]
    data = df.iloc[:, df.columns.get_level_values(1) == chilabel].T.values
    fitnames = df.columns.get_level_values(0).unique()
    expnames = list(df.index.get_level_values(0))
    fig, ax = plotutils.add_subplot(figsize=(size, size), projection='polar')
    for dat, fitname in zip(data, fitnames):
        ax = plotutils.spiderplot(expnames, dat, fitname, ax)
    ax.legend(bbox_to_anchor=(0.3, -0.2), fontsize=15)
    return fig, ax


@figure
def plot_fits_datasets_chi2(fits_datasets_chi2_table):
    """Generate a plot equivalent to ``plot_datasets_chi2`` using all the
    fitted datasets as input."""
    ind = fits_datasets_chi2_table.index.droplevel(0)
    df = fits_datasets_chi2_table.set_index(ind)
    cols = fits_datasets_chi2_table.columns.get_level_values(0).unique()
    dfs = []
    for col in cols:
        dfs.append(df[col].dropna())
    df_out = pd.concat(dfs, axis=1, keys=cols, sort=False)
    fig, ax = _plot_chis_df(df_out)
    ax.set_title(r"$\chi^2$ for datasets")
    return fig


@figure
def plot_fits_datasets_chi2_spider(fits_datasets_chi2_table):
    """Generate a plot equivalent to ``plot_datasets_chi2_spider`` using all the
    fitted datasets as input."""
    ind = fits_datasets_chi2_table.index.droplevel(0)
    df = fits_datasets_chi2_table.set_index(ind)
    cols = fits_datasets_chi2_table.columns.get_level_values(0).unique()
    dfs = []
    for col in cols:
        dfs.append(df[col].dropna())
    df_out = pd.concat(dfs, axis=1, keys=cols, sort=False)
    fig, ax = _plot_chi2s_spider_df(df_out, size=14)
    ax.set_title(r"$\chi^2$ for datasets")
    return fig


@figuregen
def plot_fits_datasets_chi2_spider_bygroup(fits_datasets_chi2_table):
    """Same as plot_fits_datasets_chi2_spider but one plot for each group."""
    tab = fits_datasets_chi2_table
    groups = tab.index.unique(level=0)
    # dfs = [tab.T[group].T for group in groups]
    for group in groups:
        df = tab.T[group].T
        fig, ax = _plot_chi2s_spider_df(df)
        ax.set_title(rf"$\chi^2$ for {group}")
        yield fig


@figure
def plot_dataspecs_datasets_chi2(dataspecs_datasets_chi2_table):
    """Same as plot_fits_datasets_chi2 but for arbitrary dataspecs"""
    return plot_fits_datasets_chi2(dataspecs_datasets_chi2_table)


@figure
def plot_dataspecs_datasets_chi2_spider(dataspecs_datasets_chi2_table):
    """Same as plot_fits_datasets_chi2_spider but for arbitrary dataspecs"""
    return plot_fits_datasets_chi2_spider(dataspecs_datasets_chi2_table)


@figure
def plot_fits_groups_data_chi2(fits_groups_chi2_table, processed_metadata_group):
    """Generate a plot equivalent to ``plot_groups_data_chi2`` using all the
    fitted group of data as input."""
    fig, ax = _plot_chis_df(fits_groups_chi2_table)
    ax.set_title(rf"$\chi^2$ by {processed_metadata_group}")
    return fig


@figure
def plot_dataspecs_groups_chi2(dataspecs_groups_chi2_table, processed_metadata_group):
    """Same as plot_fits_groups_data_chi2 but for arbitrary dataspecs"""
    return plot_fits_groups_data_chi2(dataspecs_groups_chi2_table, processed_metadata_group)


@figure
def plot_training_length(replica_data, fit):
    """Generate an histogram for the distribution
    of training lengths in a given fit. Each bin is normalised by the total
    number of replicas.

    """
    fig, ax = plotutils.subplots()
    x = [x.nite for x in replica_data]
    hist, bin_edges = np.histogram(x)
    # don't plot pdf, instead proportion of replicas in each bin.
    hist = hist / np.sum(hist)
    width = np.diff(bin_edges)
    ax.bar(bin_edges[:-1], hist, width=width, align="edge", label=str(fit))
    ax.set_title("Distribution of training lengths")
    ax.legend()
    return fig


@figure
def plot_training_validation(fit, replica_data, replica_filters=None):
    """Scatter plot with the training and validation chi² for each replica
    in the fit. The mean is also displayed as well as a line y=x to easily
    identify whether training or validation chi² is larger.

    """
    training, valid = zip(*((dt.training, dt.validation) for dt in replica_data))
    fig, ax = plotutils.subplots(
        figsize=(max(mpl.rcParams.get("figure.figsize")), max(mpl.rcParams.get("figure.figsize")))
    )
    ax.plot(training, valid, marker="o", linestyle="none", markersize=5, zorder=100)
    if replica_filters:
        _scatter_marked(ax, training, valid, replica_filters, zorder=90)
        ax.legend().set_zorder(10000)

    ax.set_title(fit.label)

    ax.set_xlabel(r"$\chi^2/N_{dat}$ training")
    ax.set_ylabel(r"$\chi^2/N_{dat}$ validation")

    min_max_lims = [min([*ax.get_xlim(), *ax.get_ylim()]), max([*ax.get_xlim(), *ax.get_ylim()])]
    ax.plot(min_max_lims, min_max_lims, ":k")

    ax.plot(np.mean(training), np.mean(valid), marker="s", color="red", markersize=7, zorder=1000)

    ax.set_aspect("equal")
    return fig


@figure
def plot_trainvaliddist(fit, replica_data):
    """KDEs for the trainning and validation distributions for
    each replica in the fit."""
    training, valid = zip(*((dt.training, dt.validation) for dt in replica_data))
    fig, ax = plotutils.subplots()

    kde_train = stats.gaussian_kde(training, bw_method='silverman')
    kde_valid = stats.gaussian_kde(valid, bw_method='silverman')
    mean = (np.array(training) + np.array(valid)) * 0.5
    kde_mean = stats.gaussian_kde(mean, bw_method='silverman')

    x = np.linspace(np.min([training, valid]), np.max([training, valid]), 150)
    ax.plot(x, kde_train(x), label="Training")
    ax.plot(x, kde_valid(x), label="Validation")
    ax.plot(x, kde_mean(x), label="Mean")

    ax.set_xlabel(r"$\chi^2/N_{dat}$")
    ax.set_title(f"KDE of the fit distributions for {fit.label}")

    ax.set_ylim(0, None)
    ax.legend()
    return fig


@figure
def plot_chi2_eigs(pdf, dataset, chi2_per_eig):
    fig, ax = plotutils.subplots()
    x = np.arange(1, len(chi2_per_eig) + 1)
    ax.plot(x, chi2_per_eig, 'o', markersize=10)
    ax.yaxis.grid(False)
    ax.set_title(fr"$\chi^2/N_{{dat}}$  {dataset}")
    ax.set_xlabel("# Eigenvalue")
    return fig


@figure
def plot_replica_sum_rules(pdf, sum_rules, Q):
    """Plot the value of each sum rule as a function of the replica index"""
    fig, axes = plotutils.subplots(nrows=len(sum_rules), sharex=True)
    # TODO: Get rid of this nonsense
    ncomputed = len(sum_rules[0])
    if pdf.error_type == 'replicas':
        x = np.arange(1, ncomputed + 1)
    else:
        x = np.arange(ncomputed)
    for label, rls, ax in zip(sum_rules._fields, sum_rules, axes):
        ax.scatter(x, rls)
        ax.set_ylabel(label)
    fig.suptitle(f'Sum rules for {pdf} at Q={Q} GeV')
    return fig


def _compute_hists(preds, nbins=10):
    """Given an 1D array of predictions, compute a binned histogram. The result
    is a dictionary containing the bin edges, the bin weights, and the bin centers.

    Parameters
    ----------
    preds : np.ndarray
        A 1D array of predictions
    nbins : int, optional
        an integer specifying the number of bins to be used in the histogram.

    Returns
    -------
    dict:
        a dictionary containing the bin edges, the bin weights, and the bin centers.
    """
    binning = np.linspace(preds.min(), preds.max(), num=nbins)
    frequencies, bins = np.histogram(preds, bins=binning, density=True)
    return {"x": bins[:-1], "bins": bins, "weights": frequencies}


def _compute_ellipse(preds_x, preds_y, nstd=3):
    """Given two arrays of predictions, each with 1D of the same length, compute
    the contour ellipse. The results is a matplotlib Ellipse instance.

    Parameters
    ----------
    preds_x : np.ndarray
        A 1D array of predictions
    preds_y : np.ndarray
        A 1D array of predictions
    nstd : int, optional
        number of standard deviation to draw the ellipse, by default 3

    Returns
    -------
    matplotlib.patches.Ellipse
        a matplotlib Ellipse instance
    """
    covmat = np.cov(preds_x, preds_y)
    eigval, eigvec = np.linalg.eig(covmat)
    sqrt_eigval = np.sqrt(eigval)

    kwargs = {
        "xy": (preds_x.mean(), preds_y.mean()),
        "width": 2.0 * nstd * sqrt_eigval[0],
        "height": 2.0 * nstd * sqrt_eigval[1],
        "angle": np.degrees(np.arctan2(*eigvec[:, 0][::-1])),
        "fill": False,
        "linewidth": 2.0,
        "color": "C0",
    }
    return Ellipse(**kwargs)


@figure
def plot_polarized_momentum(pdf, Q, partial_polarized_sum_rules, angular_momentum=False):
    """
    Plot the correlated uncertainties for the truncated integrals of the polarized
    gluon and singlet distributions.
    """
    ((xmina, xmaxa), (xminb, xmaxb)) = POL_LIMS
    predictions = partial_polarized_sum_rules[-1]  # large-x

    if not angular_momentum:
        xpreds = np.array(predictions["g"])
        ypreds = np.array(predictions["singlet"]) / 2.0
    else:
        preds_low_x = partial_polarized_sum_rules[0]  # small-x
        xpreds = np.array(predictions["g"]) + np.array(predictions["singlet"]) / 2.0
        xpreds = 1 / 2 - xpreds  # substract the proton spin
        ypreds = -(np.array(preds_low_x["g"]) + np.array(preds_low_x["singlet"]) / 2.0)

    params = {
        "width_ratios": [6, 1],
        "height_ratios": [1, 6],
        "nrows": 2,
        "ncols": 2,
        "sharex": "col",
        "sharey": "row",
    }  # Parameters that define the subplots
    fig, ((ax_histx, ax_empty), (ax_scatr, ax_histy)) = plotutils.subplots(**params)
    fig.subplots_adjust(wspace=1e-4, hspace=1e-1)

    # Add the scatter plots
    ax_scatr.scatter(xpreds, ypreds, label=pdf.label)
    ax_scatr.scatter(xpreds.mean(), ypreds.mean(), marker="s", c="red")

    # Add the Ellipse plot
    ax_scatr.add_artist(_compute_ellipse(xpreds, ypreds))

    # Add the histogram plots
    ax_histx.hist(**_compute_hists(xpreds))
    ax_histy.hist(**_compute_hists(ypreds), orientation="horizontal")

    # Define the Labels according to the type of plots
    xlabel = r"$\int_{xmin}^{xmax} \Delta g(x) dx$"
    ylabel = r"$\frac{1}{2}\int_{xmin}^{xmax}\Delta\Sigma (x)dx$"

    if angular_momentum:
        # Draw contributions from Angular Momentum
        ax_scatr.axline((0, 0), slope=-1, c="#fc6464", lw=0.65)
        ax_scatr.axline((1, 1), slope=-1, c="gray", lw=0.65)
        ax_scatr.axline((-1, -1), slope=-1, c="gray", lw=0.65)

        xaxis = r"$\frac{1}{2} -$" + xlabel + r"$-$" + ylabel
        yaxis = xlabel + r"$+$" + ylabel
        ax_scatr.set_xlabel(xaxis.replace("xmin", f"{xminb}").replace("xmax", f"{xmaxb}"))
        ax_scatr.set_ylabel(yaxis.replace("xmin", f"{xminb}").replace("xmax", f"{xmina}"))
    else:
        ax_scatr.set_xlabel(xlabel.replace("xmin", f"{xminb}").replace("xmax", f"{xmaxb}"))
        ax_scatr.set_ylabel(ylabel.replace("xmin", f"{xminb}").replace("xmax", f"{xmaxb}"))
    ax_scatr.legend(title=f"Computed at Q={round(Q, 1)} GeV", title_fontsize=10)

    # Turn off visibility of some axes
    ax_empty.set_axis_off()
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
    ax_histx.spines["left"].set_visible(False)
    ax_histx.get_yaxis().set_ticks([])
    ax_histy.axes.spines["bottom"].set_visible(False)
    ax_histy.get_xaxis().set_ticks([])

    return fig


@figure
def plot_orbital_momentum(pdf, Q, partial_polarized_sum_rules):
    """In addition to plotting the correlated spin moments as in `plot_polarized_momentum`,
    it also plots the contributions from the Orbital Angular Momentum.
    """
    return plot_polarized_momentum(pdf, Q, partial_polarized_sum_rules, angular_momentum=True)


@figuregen
def plot_smpdf(pdf, dataset, obs_pdf_correlations, mark_threshold: float = 0.9):
    """
    Plot the correlations between the change in the observable and the change
    in the PDF in (x,fl) space.

    mark_threshold is the proportion of the maximum absolute correlation
    that will be used to mark the corresponding area in x in the
    background of the plot. The maximum absolute values are used for
    the comparison.

    Examples
    --------
    >>> from validphys.api import API
    >>> data_input = {
    >>>    "dataset_input" : {"dataset": "HERACOMBNCEP920"},
    >>>    "theoryid": 200,
    >>>     "use_cuts": "internal",
    >>>     "pdf": "NNPDF40_nnlo_as_01180",
    >>>     "Q": 1.6,
    >>>     "mark_threshold": 0.2
    >>> }
    >>> smpdf_gen = API.plot_smpdf(**data_input)
    >>> fig = next(smpdf_gen)
    >>> fig.show()
    """
    info = get_info(dataset)

    table = kitable(dataset, info)

    basis = obs_pdf_correlations.basis

    fullgrid = obs_pdf_correlations.grid_values.data

    fls = obs_pdf_correlations.flavours
    x = obs_pdf_correlations.xgrid
    nf = len(fls)

    plotting_var = info.get_xcol(table)

    categorical = not np.issubdtype(plotting_var.dtype, np.number)
    if categorical:
        # Plot lines using a categorical color map (for a reasonable number of
        # categories), and set up the categorical labels (used below).
        categorical_keys, values = np.unique(plotting_var, return_inverse=True)
        plotting_var = values
        num_categories = len(categorical_keys)
        if num_categories <= len(cm.Set2.colors):
            cmap = mcolors.ListedColormap(cm.Set2.colors[:num_categories])
        else:
            cmap = cm.viridis.resample(num_categories)
        bins = np.linspace(0, num_categories, num_categories + 1)
        norm = mcolors.BoundaryNorm(bins, num_categories)

    else:
        cmap = cm.viridis
        # TODO: vmin vmax should be global or by figure?
        vmin, vmax = min(plotting_var), max(plotting_var)
        if info.x_scale == 'log':
            norm = mcolors.LogNorm(vmin, vmax)
        else:
            norm = mcolors.Normalize(vmin, vmax)

    table["__plotting_var"] = plotting_var
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)

    figby = sane_groupby_iter(table, info.figure_by)

    for same_vals, fb in figby:
        grid = fullgrid[np.asarray(fb.index), ...]

        # Use the maximum absolute correlation for plotting purposes
        absgrid = np.max(np.abs(grid), axis=0)
        mark_mask = absgrid > np.max(absgrid) * mark_threshold

        label = info.group_label(same_vals, info.figure_by)
        # TODO: PY36ScalarMappable
        # TODO Improve title?
        title = f"{info.dataset_label} {label if label else ''}\n[{pdf.label}]"

        # Start plotting
        w, h = mpl.rcParams["figure.figsize"]
        h *= 2.5
        fig, axes = plotutils.subplots(nrows=nf, sharex=True, figsize=(w, h), sharey=True)
        fig.suptitle(title)
        colors = sm.to_rgba(fb["__plotting_var"])
        for flindex, (ax, fl) in enumerate(zip(axes, fls)):
            for i, color in enumerate(colors):
                ax.plot(x, grid[i, flindex, :].T, color=color)

            flmask = mark_mask[flindex, :]
            ranges = split_ranges(x, flmask, filter_falses=True)
            for r in ranges:
                ax.axvspan(r[0], r[-1], color='#eeeeff')

            ax.set_ylabel("$%s$" % basis.elementlabel(fl))
            ax.set_xscale(scale_from_grid(obs_pdf_correlations))
            ax.set_ylim(-1, 1)
            ax.set_xlim(x[0], x[-1])
        ax.set_xlabel('$x$')

        cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), label=info.xlabel, aspect=100)
        if categorical:
            cbar.set_ticks(np.linspace(0.5, num_categories - 0.5, num_categories))
            cbar.ax.set_yticklabels(categorical_keys)

        # TODO: Fix title for this
        # fig.tight_layout()
        yield fig


@figure
def plot_obscorrs(corrpair_datasets, obs_obs_correlations, pdf):
    """NOTE: EXPERIMENTAL. Plot the correlation matrix between a pair of datasets."""
    fig, ax = plotutils.subplots()

    ds1, ds2 = corrpair_datasets
    # in1,in2 = get_info(ds1), get_info(ds2)

    im = ax.imshow(obs_obs_correlations, cmap=cm.Spectral_r, vmin=-1, vmax=1)

    ax.set_ylabel(str(ds1))
    ax.set_xlabel(str(ds2))
    fig.colorbar(im, ax=ax)
    return fig


@figure
def plot_positivity(pdfs, positivity_predictions_for_pdfs, posdataset, pos_use_kin=False):
    """Plot an errorbar spanning the central 68% CI of a positivity
    observable as well as a point indicating the central value (according
    to the ``pdf.stats_class.central_value()``).

    Errorbars and points are plotted on a symlog scale as a
    function of the data point index (if pos_use_kin==False) or the first
    kinematic variable (if pos_use_kin==True).

    """
    fig, ax = plotutils.subplots()
    ax.axhline(0, color='red')

    posset = posdataset.load_commondata()
    ndata = posset.ndata
    xvals = []

    if pos_use_kin:
        kin_name = KIN_NAMES[0]
        ax.set_xlabel(kin_name)
        xvals = posset.kinematics[kin_name].values
    else:
        ax.set_xlabel('idat')
        xvals = np.arange(ndata)

    offsets = plotutils.offset_xcentered(len(pdfs), ax)
    minscale = np.inf
    for pdf, pred in zip(pdfs, positivity_predictions_for_pdfs):
        cv = pred.central_value
        lower, upper = pred.stats.errorbar68()

        try:
            ax.errorbar(
                xvals,
                cv,
                yerr=[cv - lower, upper - cv],
                linestyle='--',
                marker='s',
                label=str(pdf),
                lw=0.5,
                transform=next(offsets),
            )
        except ValueError as e:
            if any(lower > cv) or any(upper < cv):
                raise ValueError(
                    f"The central value of {pdf} for {posdataset} is outside of the error bands. This is not supported"
                ) from e
            raise e
        minscale = min(minscale, np.abs(np.min(cv)))
    ax.legend()
    ax.set_title(str(posdataset))

    ax.set_ylabel('Observable Value')
    ax.set_yscale('symlog', linthresh=minscale)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    return fig


@make_argcheck
def _check_same_posdataset_name(dataspecs_posdataset):
    """Check that the ``posdataset`` key matches for ``dataspecs``"""
    _check_same_dataset_name.__wrapped__([ds.commondata for ds in dataspecs_posdataset])


@figure
@_check_same_posdataset_name
def plot_dataspecs_positivity(
    dataspecs_speclabel, dataspecs_positivity_predictions, dataspecs_posdataset, pos_use_kin=False
):
    """Like :py:meth:`plot_positivity` except plots positivity for each
    element of dataspecs, allowing positivity predictions to be generated with
    different ``theory_id`` s as well as ``pdf`` s
    """
    # we checked the positivity set matches between dataspecs so this is fine
    posset = dataspecs_posdataset[0]
    return plot_positivity(
        dataspecs_speclabel, dataspecs_positivity_predictions, posset, pos_use_kin
    )


@make_argcheck
def _check_display_cuts_requires_use_cuts(display_cuts, use_cuts):
    check(
        not (display_cuts and use_cuts is CutsPolicy.NOCUTS),
        "The display_cuts option requires setting some cuts",
    )


@make_argcheck
def _check_marker_by(marker_by):
    markers = ('process type', 'experiment', 'dataset', 'group')
    if marker_by not in markers:
        raise CheckError("Unknown marker_by value", marker_by, markers)


# TODO: Right now this is hackish Could we turn it into a permanent interface?
@make_argcheck
def _check_highlights(data_input, highlight_datasets):
    if highlight_datasets:
        values = frozenset(highlight_datasets)
        names_set = {ds.name for ds in data_input}
        diff = values - names_set
        if diff:
            formatted_diff = '\n'.join(diff)
            raise CheckError(
                f"The following highlight elements are not dataset names:\n{formatted_diff}"
            )
        return {'highlight_datasets': values}


@make_argcheck
def _check_aspect(aspect):
    aspects = ('landscape', 'portrait', 'square')
    if aspect not in aspects:
        raise CheckError(f"Unknown aspect {aspect}", aspect, aspects)


@figure
@_check_display_cuts_requires_use_cuts
@_check_marker_by
@_check_highlights
@_check_aspect
def plot_xq2(
    dataset_inputs_by_groups_xq2map,
    use_cuts,
    data_input,
    display_cuts: bool = True,
    marker_by: str = 'process type',
    highlight_label: str = 'highlight',
    highlight_datasets: (Sequence, type(None)) = None,
    aspect: str = 'landscape',
):
    """Plot the (x,Q²) coverage based of the data based on some LO
    approximations. These are governed by the relevant kintransform.

    The representation of the filtered data depends on the `display_cuts` and
    `use_cuts` options:

     - If cuts are disabled (`use_cuts` is CutsPolicy.NOCUTS), all the data
    will be plotted (and setting `display_cuts` to True is an error).

     - If cuts are enabled (`use_cuts` is either CutsPolicy.FROMFIT or
    CutsPolicy.INTERNAL) and `display_cuts` is False, the masked points will
    be ignored.

     - If cuts are enabled and `display_cuts` is True, the filtered points
    will be displaed and marked.

    The points are grouped according to the `marker_by` option. The possible
    values are: "process type", "experiment", "group" or "dataset".

    Some datasets can be made to appear highlighted in the figure: Define a key
    called ``highlight_datasets`` containing the names of the datasets to be
    highlighted and a key `highlight_label` with a string containing the label
    of the highlight, which will appear in the legend.

    Example
    -------

    Obtain a plot with some reasonable defaults::

    >>> from validphys.api import API
    >>> inp = {
    ...     'dataset_inputs': [
    ...         {'dataset': 'NMC_NC_NOTFIXED_EM-F2', 'variant': 'legacy_dw'},
    ...         {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy'},
    ...         {'dataset': 'SLAC_NC_NOTFIXED_P_EM-F2', 'variant': 'legacy_dw'},
    ...         {'dataset': 'SLAC_NC_NOTFIXED_D_EM-F2', 'variant': 'legacy_dw'},
    ...         {'dataset': 'BCDMS_NC_NOTFIXED_P_EM-F2', 'variant': 'legacy_dw'},
    ...         {'dataset': 'BCDMS_NC_NOTFIXED_D_EM-F2', 'variant': 'legacy_dw'},
    ...         {'dataset': 'CHORUS_CC_NOTFIXED_PB_NU-SIGMARED', 'variant': 'legacy_dw'},
    ...         {'dataset': 'CHORUS_CC_NOTFIXED_PB_NB-SIGMARED', 'variant': 'legacy_dw'},
    ...         {'dataset': 'NUTEV_CC_NOTFIXED_FE_NU-SIGMARED', 'cfac': ['MAS'], 'variant': 'legacy_dw'},
    ...         {'dataset': 'NUTEV_CC_NOTFIXED_FE_NB-SIGMARED', 'cfac': ['MAS'], 'variant': 'legacy_dw'},
    ...         {'dataset': 'HERA_NC_318GEV_EM-SIGMARED', 'variant': 'legacy'},
    ...         {'dataset': 'HERA_NC_225GEV_EP-SIGMARED', 'variant': 'legacy'},
    ...         {'dataset': 'HERA_NC_251GEV_EP-SIGMARED', 'variant': 'legacy'},
    ...         {'dataset': 'HERA_NC_300GEV_EP-SIGMARED', 'variant': 'legacy'},
    ...         {'dataset': 'HERA_NC_318GEV_EP-SIGMARED', 'variant': 'legacy'},
    ...         {'dataset': 'HERA_CC_318GEV_EM-SIGMARED', 'variant': 'legacy'},
    ...         {'dataset': 'HERA_CC_318GEV_EP-SIGMARED', 'variant': 'legacy'},
    ...         {'dataset': 'HERA_NC_318GEV_EAVG_CHARM-SIGMARED', 'variant': 'legacy'},
    ...         {'dataset': 'HERA_NC_318GEV_EAVG_BOTTOM-SIGMARED', 'variant': 'legacy'},
    ...         {'dataset': 'DYE866_Z0_800GEV_DW_RATIO_PDXSECRATIO', 'variant': 'legacy'},
    ...         {'dataset': 'DYE866_Z0_800GEV_PXSEC', 'variant': 'legacy'},
    ...         {'dataset': 'DYE605_Z0_38P8GEV_DW_PXSEC', 'variant': 'legacy'},
    ...         {'dataset': 'DYE906_Z0_120GEV_DW_PDXSECRATIO', 'cfac': ['ACC'], 'variant': 'legacy'},
    ...         {'dataset': 'CDF_Z0_1P96TEV_ZRAP', 'variant': 'legacy'},
    ...         {'dataset': 'D0_Z0_1P96TEV_ZRAP', 'variant': 'legacy'},
    ...         {'dataset': 'D0_WPWM_1P96TEV_ASY', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_WPWM_7TEV_36PB_ETA', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_Z0_7TEV_36PB_ETA', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_Z0_7TEV_49FB_HIMASS', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_Z0_7TEV_LOMASS_M', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_WPWM_7TEV_46FB_CC-ETA', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_Z0_7TEV_46FB_CC-Y', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_Z0_7TEV_46FB_CF-Y', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_Z0_8TEV_HIMASS_M-Y', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_Z0_8TEV_LOWMASS_M-Y', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_Z0_13TEV_TOT', 'cfac': ['NRM'], 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_WPWM_13TEV_TOT', 'cfac': ['NRM'], 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_WJ_8TEV_WP-PT', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_WJ_8TEV_WM-PT', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_Z0J_8TEV_PT-M', 'variant': 'legacy_10'},
    ...         {'dataset': 'ATLAS_Z0J_8TEV_PT-Y', 'variant': 'legacy_10'},
    ...         {'dataset': 'ATLAS_TTBAR_7TEV_TOT_X-SEC', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_TTBAR_8TEV_TOT_X-SEC', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_TTBAR_13TEV_TOT_X-SEC', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_TTBAR_8TEV_LJ_DIF_YT-NORM', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_TTBAR_8TEV_LJ_DIF_YTTBAR-NORM', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_TTBAR_8TEV_2L_DIF_YTTBAR-NORM', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_1JET_8TEV_R06_PTY', 'variant': 'legacy_decorrelated'},
    ...         {'dataset': 'ATLAS_2JET_7TEV_R06_M12Y', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_PH_13TEV_XSEC', 'cfac': ['EWK'], 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_SINGLETOP_7TEV_TCHANNEL-XSEC', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_SINGLETOP_13TEV_TCHANNEL-XSEC', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_SINGLETOP_7TEV_T-Y-NORM', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_SINGLETOP_7TEV_TBAR-Y-NORM', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_SINGLETOP_8TEV_T-RAP-NORM', 'variant': 'legacy'},
    ...         {'dataset': 'ATLAS_SINGLETOP_8TEV_TBAR-RAP-NORM', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_WPWM_7TEV_ELECTRON_ASY'},
    ...         {'dataset': 'CMS_WPWM_7TEV_MUON_ASY', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_Z0_7TEV_DIMUON_2D'},
    ...         {'dataset': 'CMS_WPWM_8TEV_MUON_Y', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_Z0J_8TEV_PT-Y', 'cfac': ['NRM'], 'variant': 'legacy_10'},
    ...         {'dataset': 'CMS_2JET_7TEV_M12Y'},
    ...         {'dataset': 'CMS_1JET_8TEV_PTY', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_TTBAR_7TEV_TOT_X-SEC', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_TTBAR_8TEV_TOT_X-SEC', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_TTBAR_13TEV_TOT_X-SEC', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_TTBAR_8TEV_LJ_DIF_YTTBAR-NORM', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_TTBAR_5TEV_TOT_X-SEC', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_TTBAR_8TEV_2L_DIF_MTTBAR-YT-NORM', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_TTBAR_13TEV_2L_DIF_YT', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_TTBAR_13TEV_LJ_2016_DIF_YT', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_SINGLETOP_7TEV_TCHANNEL-XSEC', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_SINGLETOP_8TEV_TCHANNEL-XSEC', 'variant': 'legacy'},
    ...         {'dataset': 'CMS_SINGLETOP_13TEV_TCHANNEL-XSEC', 'variant': 'legacy'},
    ...         {'dataset': 'LHCB_Z0_7TEV_DIELECTRON_Y'},
    ...         {'dataset': 'LHCB_Z0_8TEV_DIELECTRON_Y'},
    ...         {'dataset': 'LHCB_WPWM_7TEV_MUON_Y', 'cfac': ['NRM']},
    ...         {'dataset': 'LHCB_Z0_7TEV_MUON_Y', 'cfac': ['NRM']},
    ...         {'dataset': 'LHCB_WPWM_8TEV_MUON_Y', 'cfac': ['NRM']},
    ...         {'dataset': 'LHCB_Z0_8TEV_MUON_Y', 'cfac': ['NRM']},
    ...         {'dataset': 'LHCB_Z0_13TEV_DIMUON-Y'},
    ...         {'dataset': 'LHCB_Z0_13TEV_DIELECTRON-Y'}
    ...     ],
    ...     'use_cuts': 'internal',
    ...     'display_cuts': False,
    ...     'theoryid': 40_000_000,
    ...     'highlight_label': 'Tevatron',
    ...     'highlight_datasets': ['CDF_Z0_1P96TEV_ZRAP', 'D0_Z0_1P96TEV_ZRAP', 'D0_WPWM_1P96TEV_ASY']
    ... }
    >>> API.plot_xq2(**inp)
    <Figure size 1024x768 with 1 Axes>

    """

    w, h = mpl.rcParams["figure.figsize"]
    rescaling_factor = 1.6
    w *= rescaling_factor
    h *= rescaling_factor
    if aspect == 'landscape':
        figsize = w, h
    elif aspect == 'portrait':
        figsize = h, w
    elif aspect == 'square':
        figsize = h, h
    else:
        raise ValueError(f"Unknown aspect {aspect}")
    fig, ax = plotutils.subplots(figsize=figsize)

    filteredx = []
    filteredq2 = []

    x = defaultdict(list)
    q2 = defaultdict(list)

    xh = defaultdict(list)
    q2h = defaultdict(list)

    if not highlight_datasets:
        highlight_datasets = set()

    def next_options():
        # Get the colors
        prop_settings = mpl.rcParams['axes.prop_cycle']
        # Apparently calling the object gives us an infinite cycler
        settings_cycler = prop_settings()
        # So far, I don't understand how this is done with mpl "cycler"
        # objects, or wether  I like it. So far this is godd enough
        for markeropts, settings in zip(plotutils.marker_iter_plot(), settings_cycler):
            # Override last with first
            options = {'linestyle': 'none', **markeropts, **settings}
            yield options

    next_opts = next_options()
    key_options = {}

    for experiment, commondata, fitted, masked, group in dataset_inputs_by_groups_xq2map:
        info = get_info(commondata)
        if marker_by == 'process type':
            key = info.process_description
        elif marker_by == 'experiment':
            key = str(experiment)
        elif marker_by == 'dataset':
            key = info.dataset_label
        elif marker_by == "group":
            # if group is None then make sure that shows on legend.
            key = str(group)
        else:
            raise ValueError('Unknown marker_by value')

        # TODO: This is an ugly check. Is there a way to do it with .setdefault
        # or defaultdict?
        if key not in key_options:
            key_options[key] = next(next_opts)

        if commondata.name in highlight_datasets:
            xdict = xh
            q2dict = q2h
        else:
            xdict = x
            q2dict = q2

        xdict[key].append(fitted[0])
        q2dict[key].append(fitted[1])
        if display_cuts:
            xdict[key].append(masked[0])
            q2dict[key].append(masked[1])
            filteredx.append(masked[0])
            filteredq2.append(masked[1])

    for key in key_options:
        if key in x:
            coords = np.concatenate(x[key]), np.concatenate(q2[key])
        else:
            # This is to get the label key
            coords = [], []
        ax.plot(*coords, label=key, markeredgewidth=1, markeredgecolor=None, **key_options[key])

    # Iterate again so highlights are printed on top.
    for key in xh:
        ax.plot(
            np.concatenate(xh[key]),
            np.concatenate(q2h[key]),
            markeredgewidth=0.6,
            markeredgecolor="black",
            **key_options[key],
        )
    if xh:
        # Get legend key
        ax.plot(
            [],
            [],
            marker='s',
            markeredgewidth=0.6,
            color='none',
            markersize=5,
            markeredgecolor="black",
            label=f'Black edge: {highlight_label}',
        )

    if display_cuts:
        ax.scatter(
            np.concatenate(filteredx),
            np.concatenate(filteredq2),
            marker='o',
            facecolors='none',
            edgecolor='red',
            s=40,
            lw=0.8,
            label="Cut",
        )

    ax.set_title("Kinematic coverage")
    ax.legend()
    ax.set_xlabel('$x$')
    ax.set_ylabel(r'$Q^2$ (GeV$^2$)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return fig
