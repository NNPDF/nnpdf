# -*- coding: utf-8 -*-
"""
Plots of relations between data PDFs and fits.
"""
from __future__ import generator_stop

import logging
import itertools
from collections import defaultdict
from collections.abc import Sequence

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors, ticker as mticker
import scipy.stats as stats
import pandas as pd

from reportengine.figure import figure, figuregen
from reportengine.checks import make_check, CheckError, make_argcheck, check
from reportengine.floatformatting import format_number

from validphys.core import MCStats, cut_mask, CutsPolicy
from validphys.results import chi2_stat_labels
from validphys.plotoptions import get_info, kitable, transform_result
from validphys import plotutils
from validphys.utils import sane_groupby_iter, split_ranges, scale_from_grid

log = logging.getLogger(__name__)

@figure
def plot_chi2dist_experiments(total_experiments_chi2data, experiments_chi2_stats, pdf):
    """Plot the distribution of chi²s of the members of the pdfset."""
    fig, ax = _chi2_distribution_plots(total_experiments_chi2data, experiments_chi2_stats, pdf, "hist")
    ax.set_title(r"Experiments $\chi^2$ distribution")
    return fig


@figure
def kde_chi2dist_experiments(total_experiments_chi2data, experiments_chi2_stats, pdf):
    """KDE plot for experiments chi2."""
    fig, ax = _chi2_distribution_plots(total_experiments_chi2data, experiments_chi2_stats, pdf, "kde")
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


def _chi2_distribution_plots(chi2_data, stats, pdf, plot_type):
    fig, ax = plt.subplots()
    label = pdf.name
    alldata, central, npoints = chi2_data
    if not isinstance(alldata, MCStats):
        ax.set_facecolor("#ffcccc")
        log.warning("Chi² distribution plots have a "
                "different meaning for non MC sets.")
        label += " (%s!)" % pdf.ErrorType
    label += '\n'+ '\n'.join(str(chi2_stat_labels[k])+(' %.2f' % v) for (k,v) in stats.items())
    ax.set_xlabel(r"Replica $\chi^2$")

    if plot_type == "hist":
        ax.hist(alldata.data, label=label, zorder=100)
    elif plot_type == "kde":
        # We need the squeeze here to change shape from (x, 1) to (x,)
        ax = plotutils.kde_plot(alldata.data.squeeze(), label=label)
    else:
        raise ValueError(f"plot_type must either be hist or kde, not {plot_type}")

    l = ax.legend()
    l.set_zorder(1000)
    return fig, ax


@figure
def plot_phi(experiments, experiments_phi):
    """plots phi for each experiment as a bar for a single
    PDF input

    See `phi_data` for information on how phi is calculated
    """
    phi = [exp_phi for (exp_phi, npoints) in experiments_phi]
    xticks = [experiment.name for experiment in experiments]
    fig, ax = plotutils.barplot(phi, collabels=xticks, datalabels=[r'$\phi$'])
    ax.set_title(r"$\phi$ for each experiment")
    return fig

@figure
def plot_fits_experiments_phi(fits_experiments_phi_table):
    """Plots a set of bars for each fit, each bar represents the value of phi for the corresponding
    experiment, where the experiment is a group of datasets according to the `experiment` key in
    the PLOTTING info file"""
    fig, ax = _plot_chis_df(fits_experiments_phi_table)
    ax.set_title(r"$\phi$ for each experiment")
    return fig

@figure
def plot_phi_experiment_dist(experiment, bootstrap_phi_data_experiment):
    """Generates a bootstrap distribution of phi and then plots a histogram
    of the individual bootstrap samples for a single experiment. By default
    the number of bootstrap samples is set to a sensible number (500) however
    this number can be changed by specifying `bootstrap_samples` in the runcard
    """
    phi = bootstrap_phi_data_experiment
    label = '\n'.join([fr'$\phi$ mean = {format_number(phi.mean())}',
                       fr'$\phi$ std dev = {format_number(phi.std())}'])
    fig, ax = plt.subplots()
    ax.hist(phi, label=label)
    ax.set_title(r"$\phi$ distribution for " + experiment.name)
    ax.legend()
    return fig

@make_argcheck
def _check_same_experiment_name(dataspecs_experiments):
    lst = dataspecs_experiments
    if not lst:
        return
    for j, x in enumerate(lst[1:]):
        if len(x) != len(lst[0]):
            raise CheckError("All dataspecs should have the same number "
                             "of experiments")
        for i, exp in enumerate(x):
            if exp.name != lst[0][i].name:
                raise CheckError("\n".join(["All experiments must have the "
                                            "same name",
                                            fr"dataspec {j+1}, "
                                            fr"experiment {i+1}: {exp.name}",
                                            fr"dataspec 1, experiment {i+1}: "
                                            fr"{lst[0][i].name}"]))

@_check_same_experiment_name
@figure
def plot_phi_scatter_dataspecs(dataspecs_experiments,
        dataspecs_speclabel, dataspecs_experiments_bootstrap_phi):
    """For each of the dataspecs, a bootstrap distribution of phi is generated
    for all specified experiments. The distribution is then represented as a
    scatter point which is the median of the bootstrap distribution and an
    errorbar which spans the 68% confidence interval. By default the number
    of bootstrap samples is set to a sensible value, however it can be
    controlled by specifying `bootstrap_samples` in the runcard.
    """
    labels = dataspecs_speclabel
    phis = dataspecs_experiments_bootstrap_phi
    exps = dataspecs_experiments
    xticks = [experiment.name for experiment in exps[0]]
    x = range(1, len(xticks)+1)
    fig, ax = plt.subplots()
    phi_stats = np.percentile(phis, [16, 50, 84], axis=2)
    for i, label in enumerate(labels):
        phi_errs = np.vstack((phi_stats[2, i, :] - phi_stats[1, i, :],
                              phi_stats[1, i, :] - phi_stats[0, i, :]))
        ax.errorbar(x, phi_stats[1, i, :], yerr=phi_errs, fmt='.',
                    label=label)
    ax.set_xticks(x, minor=False)
    ax.set_xticklabels(xticks, minor=False, rotation=45)
    ax.legend()
    return fig

#TODO: This should be simplified if at all possible. For now some more examples
#are needed for a spec to emerge.
@make_check
def check_normalize_to(ns, **kwargs):
    """Transforn normalize_to into an index."""

    msg = ("normalize_to should be either 'data', a pdf id or an index of the "
           "result (0 for the data, and i for the ith pdf)")

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

#TODO: This interface is horrible. We need to think how to adapt libnnpdf
#to make this use case easier
def _plot_fancy_impl(results, commondata, cutlist,
               normalize_to:(int,type(None)) = None, labellist=None):

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

    #This is easier than cheking every time
    if labellist is None:
        labellist = [None]*len(results)

    if normalize_to is not None:
        norm_result = results[normalize_to]
        mask = cut_mask(cutlist[normalize_to])
        cv = np.full(ndata, np.nan)
        cv[mask] = norm_result.central_value

        err = np.full(ndata, np.nan)
        err[mask] =  norm_result.std_error
        #We modify the table, so we pass only the label columns
        norm_cv, _ = transform_result(cv,
                                   err,
                                   table.iloc[:,:nkinlabels], info)


    cvcols = []
    for i,(result, cuts) in enumerate(zip(results, cutlist)):
        #We modify the table, so we pass only the label columns
        mask = cut_mask(cuts)
        cv = np.full(ndata, np.nan)
        cv[mask] = result.central_value
        err = np.full(ndata, np.nan)
        err[mask] = result.std_error

        cv, err = transform_result(cv, err,
                                   table.iloc[:,:nkinlabels], info)
        #By doing tuple keys we avoid all possible name collisions
        cvcol = ('cv', i)
        if normalize_to is None:
            table[cvcol] = cv
            table[('err', i)] = err
        else:
            table[cvcol] = cv/norm_cv
            table[('err', i)] = err/norm_cv
        cvcols.append(cvcol)

    figby = sane_groupby_iter(table, info.figure_by)


    for samefig_vals, fig_data in figby:
        #Nothing to plot if all data is cut away
        if np.all(np.isnan(fig_data[cvcols])):
            continue
        #For some reason matplotlib doesn't set the axis right
        min_vals = []
        max_vals = []
        fig, ax = plt.subplots()
        ax.set_title("%s %s"%(info.dataset_label,
                     info.group_label(samefig_vals, info.figure_by)))

        lineby = sane_groupby_iter(fig_data, info.line_by)

        first = True


        for (sameline_vals, line_data) in lineby:
            ax.set_prop_cycle(None)
            labels = first
            first = False

            offset_iter = plotutils.offset_xcentered(len(results), ax)

            x = info.get_xcol(line_data)

            try:
                x = np.asanyarray(x, np.float)
            except ValueError:
                xticklabels = x
                npoints = len(x)
                x = np.arange(npoints)
                ax.set_xticks(x)
                ax.set_xticklabels(xticklabels)
                #TODO: Remove this when mpl stops doing the wrong thing
                #(in v2?)
                ax.set_xlim(-npoints/20, npoints - 1+ npoints/20)


            #Use black for the first iteration (data),
            #and follow the cycle for
            #the rest.
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
                ax.errorbar(x, cv, yerr=err,
                     linestyle='--',
                     lw=0.25,
                     label= label,
                     #elinewidth = 2,
                     capsize=2,
                     marker = 's',
                     markeredgewidth=0.25,
                     c=color,
                     zorder=1000,
                     transform=next(offset_iter))


                #We 'plot' the empty lines to get the labels. But
                #if everything is rmpty we skip the plot.
                if np.any(np.isfinite(cv)):
                    max_vals.append(np.nanmax(cv+err))
                    min_vals.append(np.nanmin(cv-err))

            glabel = info.group_label(sameline_vals, info.line_by)

            #Use some anchor that is not in y=1 for ratio plots
            if normalize_to is not None:
                next_after_normalize = (normalize_to + 1) % len(results)
                annotate_point = x[-1], line_data[('cv', next_after_normalize)].values[-1]
            else:
                annotate_point = x[-1], line_data[('cv', 0)].values[-1]
            #This is a workaround for https://github.com/matplotlib/matplotlib/issues/12648
            if np.isfinite(annotate_point).all():
                ax.annotate(glabel, annotate_point, xytext=(15 ,-10),
                             size='xx-small',
                             textcoords='offset points', zorder=10000)

        if info.x_scale:
            ax.set_xscale(info.x_scale)

        if info.y_scale:
            ax.set_yscale(info.y_scale)

        if normalize_to is None:
            if info.y_label:
                ax.set_ylabel(info.y_label)
        else:
            lb = labellist[normalize_to]
            ax.set_ylabel(f"Ratio to {lb if lb else norm_result.label}")


        ax.legend().set_zorder(100000)
        ax.set_xlabel(info.xlabel)
        fig.tight_layout()
        yield fig


@check_normalize_to
@figuregen
def plot_fancy(one_or_more_results, commondata, cuts,
               normalize_to:(int,str,type(None)) = None):
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
    yield from _plot_fancy_impl(results=one_or_more_results,
                                commondata=commondata,
                                cutlist=[cuts]*len(one_or_more_results),
                                normalize_to=normalize_to)

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
    if (normalize_to in (0, None) or
            (isinstance(normalize_to, int) and normalize_to <= len(dataspecs))):
        return
    if normalize_to == 'data':
        return {'normalize_to': 0}

    raise CheckError("Unrecignized format for normalize_to. Must be either "
                     "'data', 0 or the 1-indexed index of the dataspec "
                     f"(<{len(dataspecs)}), not {normalize_to}")



@_check_same_dataset_name
@_check_dataspec_normalize_to
@figuregen
def plot_fancy_dataspecs(dataspecs_results, dataspecs_commondata,
                         dataspecs_cuts, dataspecs_speclabel,
                         normalize_to:(str, int, type(None))=None):
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
    #We have at least one element
    if not dataspecs_results:
        return

    #For now, simply take the first data result. We'll need to improve this.
    results = [dataspecs_results[0][0], *[r[1] for r in dataspecs_results]]
    cutlist = [dataspecs_cuts[0], *dataspecs_cuts]
    commondata = dataspecs_commondata[0]
    labellist = [None, *dataspecs_speclabel]
    yield from _plot_fancy_impl(results = results, commondata=commondata,
                                cutlist=cutlist, labellist=labellist,
                                normalize_to=normalize_to)




def _scatter_marked(ax, x, y, marked_dict, *args, **kwargs):
    kwargs['s'] = kwargs.get('s', 30) + 10
    x = np.array(x, copy=False)
    y = np.array(y, copy=False)
    for label, indexes in marked_dict.items():
        ax.scatter(x[indexes],y[indexes], *args, **kwargs, label=label,
                   facecolors='none', linewidth=0.5, edgecolor='red')
        kwargs['s'] += 10


#I need to use the informations contained in experiments_chi2_table



@figure
def plot_experiments_chi2(experiments, experiments_chi2):
    """Plot the chi² of all experiments with bars."""
    exchi2 = []
    xticks = []
    for experiment, expres in zip(experiments, experiments_chi2):
        exchi2.append(expres.central_result/expres.ndata)
        xticks.append(experiment.name)
    fig, ax = plotutils.barplot(exchi2, collabels=xticks, datalabels=[r'$\chi^2$'])
    ax.set_title(r"$\chi^2$ distribution for experiments")
    return fig

@figure
def plot_datasets_chi2(experiments, experiments_chi2,each_dataset_chi2):
    """Plot the chi² of all datasets with bars."""
    ds = iter(each_dataset_chi2)
    dschi2 = []
    xticks = []
    for experiment, expres in zip(experiments, experiments_chi2):
        for dataset, dsres in zip(experiment, ds):
            dschi2.append(dsres.central_result/dsres.ndata)
            xticks.append(dataset.name)
    fig,ax = plotutils.barplot(dschi2, collabels=xticks,
                               datalabels=[r'$\chi^2$'])

    ax.set_title(r"$\chi^2$ distribution for datasets")

    return fig

def _plot_chis_df(df):
    chilabel = df.columns.get_level_values(1)[1]
    data = df.iloc[:, df.columns.get_level_values(1)==chilabel].T.values
    fitnames = df.columns.get_level_values(0).unique()
    expnames = list(df.index.get_level_values(0))
    fig, ax = plotutils.barplot(data, expnames, fitnames)
    ax.grid(False)
    ax.legend()
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
def plot_dataspecs_datasets_chi2(dataspecs_datasets_chi2_table):
    """Same as plot_fits_datasets_chi2 but for arbitrary dataspecs"""
    return plot_fits_datasets_chi2(dataspecs_datasets_chi2_table)

@figure
def plot_fits_experiments_chi2(fits_experiments_chi2_table):
    """Generate a plot equivalent to ``plot_experiments_chi2`` using all the
    fitted experiments as input."""
    fig, ax = _plot_chis_df(fits_experiments_chi2_table)
    ax.set_title(r"$\chi^2$ for experiments")
    return fig

@figure
def plot_dataspecs_experiments_chi2(dataspecs_experiments_chi2_table):
    """Same as plot_fits_experiments_chi2 but for arbitrary dataspecs"""
    return plot_fits_experiments_chi2(dataspecs_experiments_chi2_table)

@figure
def plot_training_length(replica_data, fit):
    """Generate an histogram for the distribution
    of training lengths in a given fit."""
    fig, ax = plt.subplots()
    x = [x.nite for x in replica_data]
    ax.hist(x, density=True, label=str(fit))
    ax.set_title("Distribution of training lengths")
    ax.legend()
    return fig


@figure
def plot_training_validation(fit, replica_data, replica_filters=None):
    """Scatter plot with the training and validation chi² for each replica
    in the fit. The mean is also displayed"""
    training, valid = zip(*((dt.training, dt.validation) for dt in replica_data))
    fig, ax = plt.subplots()
    ax.plot(training, valid, marker='o', linestyle='none', markersize=5, zorder=100)
    if replica_filters:
        _scatter_marked(ax, training,valid, replica_filters, zorder=90)
        ax.legend().set_zorder(10000)

    ax.set_title(fit.label)

    ax.set_xlabel(r'$\chi^2/N_{dat}$ train')
    ax.set_ylabel(r'$\chi^2/N_{dat}$ valid')

    ax.plot(np.mean(training), np.mean(valid),
         marker='s', color='red', markersize=7, zorder=1000)

    return fig


@figure
def plot_trainvaliddist(fit, replica_data):
    """KDEs for the trainning and validation distributions for
    each replica in the fit."""
    training, valid = zip(*((dt.training, dt.validation) for dt in replica_data))
    fig, ax = plt.subplots()

    kde_train = stats.gaussian_kde(training, bw_method='silverman')
    kde_valid = stats.gaussian_kde(valid, bw_method='silverman')
    mean = (np.array(training) + np.array(valid))*0.5
    kde_mean = stats.gaussian_kde(mean, bw_method='silverman')

    x = np.linspace(np.min([training,valid]),
                    np.max([training, valid]), 150)
    ax.plot(x, kde_train(x), label="Training")
    ax.plot(x, kde_valid(x), label="Validation")
    ax.plot(x, kde_mean(x), label="Mean")

    ax.set_xlabel(r"$\chi^2/N_{dat}$")
    ax.set_title(f"KDE of the fit distributions for {fit.label}")

    ax.set_ylim(0, None)
    ax.legend()
    return fig

@figure
def plot_covmat_eigs(experiment):
    """Plot the eigenvalues of the covariance matrix for a given experiment."""
    eigs = la.eigvalsh(experiment.load().get_covmat())
    fig,ax = plt.subplots()
    x = np.arange(1,len(eigs) + 1)
    ax.plot(x, eigs, 'o', markersize=10)
    ax.set_yscale('log')
    ax.yaxis.grid(False)
    plt.title("Covmat eigenvalues for %s" % experiment.name)
    plt.xlabel("# Eigenvector")
    return fig

@figure
def plot_corrmat_eigs(experiment):
    """Plot the eigenvalues of the correlation matrix for a given experiment."""
    covmat = experiment.load().get_covmat()
    stds = np.sqrt(np.diag(covmat))
    corrmat = covmat/np.outer(stds,stds)
    eigs = la.eigvalsh(corrmat)
    fig,ax = plt.subplots()
    ax.plot(eigs, 'o')
    ax.set_yscale('log')
    return fig


@figure
def plot_chi2_eigs(pdf,dataset,chi2_per_eig):
    fig,ax = plt.subplots()
    x = np.arange(1,len(chi2_per_eig) + 1)
    ax.plot(x, chi2_per_eig, 'o', markersize=10)
    ax.yaxis.grid(False)
    plt.title(fr"$\chi^2/N_{{dat}}$  {dataset}")
    plt.xlabel("# Eigenvalue")
    return fig

@figure
def plot_replica_sum_rules(pdf, sum_rules, Q):
    """Plot the value of each sum rule as a function of the replica index"""
    fig, axes = plt.subplots(nrows=len(sum_rules), sharex=True)
    #TODO: Get rid of this nonsense
    ncomputed = len(sum_rules[0])
    if pdf.ErrorType == 'replicas':
        x = np.arange(1, ncomputed + 1)
    else:
        x = np.arange(ncomputed)
    for label, rls, ax in zip(sum_rules._fields, sum_rules, axes):
        ax.scatter(x, rls)
        ax.set_ylabel(label)
    fig.suptitle(f'Sum rules for {pdf} at Q={Q} GeV')
    return fig

@figuregen
def plot_smpdf(pdf, dataset, obs_pdf_correlations, mark_threshold:float=0.9):
    """
    Plot the correlations between the change in the observable and the change
    in the PDF in (x,fl) space.

    mark_threshold is the proportion of the maximum absolute correlation
    that will be used to mark the corresponding area in x in the
    background of the plot. The maximum absolute values are used for
    the comparison."""
    info = get_info(dataset)

    table = kitable(dataset, info)
    figby = sane_groupby_iter(table, info.figure_by)

    basis = obs_pdf_correlations.basis

    fullgrid = obs_pdf_correlations.grid_values

    fls = obs_pdf_correlations.flavours
    x = obs_pdf_correlations.xgrid
    nf = len(fls)

    plotting_var = info.get_xcol(table)

    #TODO: vmin vmax should be global or by figure?
    vmin,vmax = min(plotting_var), max(plotting_var)
    if info.x_scale == 'log':
        norm = mcolors.LogNorm(vmin, vmax)
    else:
        norm = mcolors.Normalize(vmin, vmax)
    #http://stackoverflow.com/a/11558629/1007990
    sm = cm.ScalarMappable(cmap=cm.viridis, norm=norm)
    sm._A = []

    for same_vals, fb in figby:
        grid = fullgrid[ np.asarray(fb.index),...]


        #Use the maximum absolute correlation for plotting purposes
        absgrid = np.max(np.abs(grid), axis=0)
        mark_mask = absgrid > np.max(absgrid)*mark_threshold

        label = info.group_label(same_vals, info.figure_by)
        #TODO: PY36ScalarMappable
        #TODO Improve title?
        title = "%s %s\n[%s]" % (info.dataset_label, '(%s)'%label if label else '' ,pdf.label)

        #Start plotting
        w,h = plt.rcParams["figure.figsize"]
        h*=2.5
        fig,axes = plt.subplots(nrows=nf ,sharex=True, figsize=(w,h), sharey=True)
        fig.suptitle(title)
        colors = sm.to_rgba(info.get_xcol(fb))
        for flindex, (ax, fl) in enumerate(zip(axes, fls)):
            for i,color in enumerate(colors):
                ax.plot(x, grid[i,flindex,:].T, color=color)


            flmask = mark_mask[flindex,:]
            ranges = split_ranges(x, flmask, filter_falses=True)
            for r in ranges:
                ax.axvspan(r[0], r[-1], color='#eeeeff')

            ax.set_ylabel("$%s$"%basis.elementlabel(fl))
            ax.set_xscale(scale_from_grid(obs_pdf_correlations))
            ax.set_ylim(-1,1)
            ax.set_xlim(x[0], x[-1])
        ax.set_xlabel('$x$')
        #fig.subplots_adjust(hspace=0)

        fig.colorbar(sm, ax=axes.ravel().tolist(), label=info.xlabel,
                     aspect=100)
        #TODO: Fix title for this
        #fig.tight_layout()
        yield fig

@figure
def plot_obscorrs(corrpair_datasets, obs_obs_correlations, pdf):
    """NOTE: EXPERIMENTAL. Plot the correlation matrix between a pair of datasets."""
    fig, ax = plt.subplots()

    ds1, ds2 = corrpair_datasets
    #in1,in2 = get_info(ds1), get_info(ds2)

    im = ax.imshow(obs_obs_correlations, cmap=cm.Spectral_r, vmin=-1, vmax=1)

    ax.set_ylabel(str(ds1))
    ax.set_xlabel(str(ds2))
    fig.colorbar(im, [ax])
    return fig

@figure
def plot_positivity(pdfs, positivity_predictions_for_pdfs, posdataset, pos_use_kin=False):
    """Plot the value of a positivity observable on a symlog scale as a
    function of the data point index (if pos_use_kin==False) or the first
    kinematic variable (if pos_use_kin==True)."""
    fig, ax = plt.subplots()
    ax.axhline(0, color='red')

    posset = posdataset.load()
    ndata  = posset.GetNData()
    xvals = []

    if pos_use_kin:
        ax.set_xlabel('kin1')
        xvals = [posset.GetKinematics(i, 0) for i in range(0, ndata)]
    else:
        ax.set_xlabel('idat')
        xvals = np.arange(ndata)

    offsets = plotutils.offset_xcentered(len(pdfs), ax)
    minscale = np.inf
    for i, (pdf, pred) in enumerate(zip(pdfs, positivity_predictions_for_pdfs)):
        cv = pred.central_value
        ax.errorbar(xvals, cv, yerr=pred.std_error,
                    linestyle='--',
                    marker='s',
                    label=pdf.label, lw=0.5, transform=next(offsets))
        minscale = min(minscale, np.abs(np.min(cv)))
    ax.legend()
    ax.set_title(str(posdataset))

    ax.set_ylabel('Observable Value')
    ax.set_yscale('symlog', linthreshy=minscale)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    return fig

@make_argcheck
def _check_display_cuts_requires_use_cuts(display_cuts, use_cuts):
    check(
        not (display_cuts
             and use_cuts not in (CutsPolicy.FROMFIT, CutsPolicy.INTERNAL)),
        "The display_cuts option requires setting use_cuts to True")

@make_argcheck
def _check_marker_by(marker_by):
    markers = ('process type', 'experiment', 'dataset')
    if marker_by not in markers:
        raise CheckError("Unknown marker_by value", marker_by, markers)

#TODO: Right now this is hackish Could we turn it into a permanent interface?
@make_argcheck
def _check_highlights(experiments, highlight_datasets):
    if highlight_datasets:
        values = frozenset(highlight_datasets)
        names_set = {ds.name for experiment in experiments for ds in experiment }
        diff = values - names_set
        if diff:
            raise CheckError(f"The following highlight elements are "
                             "not dataset names: {diff}")
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
def plot_xq2(experiments_xq2map, use_cuts ,display_cuts:bool=True,
                 marker_by:str='process type', highlight_label:str='highlight',
                 highlight_datasets:(Sequence,type(None))=None,
                 aspect:str='landscape'):
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
    values are: "process type", "experiment" or "dataset".
    """

    w,h = plt.rcParams["figure.figsize"]
    rescaling_factor = 1.6
    w *= rescaling_factor
    h *= rescaling_factor
    if aspect=='landscape':
        figsize = w, h
    elif aspect=='portrait':
        figsize = h, w
    elif aspect=='square':
        figsize = h, h
    else:
        raise ValueError(f"Unknown aspect {aspect}")
    fig, ax = plt.subplots(figsize=figsize)

    filteredx = []
    filteredq2 = []

    x = defaultdict(list)
    q2 = defaultdict(list)

    xh = defaultdict(list)
    q2h = defaultdict(list)

    if not highlight_datasets:
        highlight_datasets = set()

    def next_options():
        #Get the colors
        prop_settings = plt.rcParams['axes.prop_cycle']
        #Apparently calling the object gives us an infinite cycler
        settings_cycler = prop_settings()
        #So far, I don't understand how this is done with mpl "cycler"
        #objects, or wether  I like it. So far this is godd enough
        for  markeropts, settings in zip(plotutils.marker_iter_plot(), settings_cycler):
            #Override last with first
            options = {
                'linestyle': 'none',
                **markeropts,
                **settings,
            }
            yield options

    next_opts = next_options()
    key_options = {}

    for experiment, commondata, fitted, masked in experiments_xq2map:
        info = get_info(commondata)
        if marker_by == 'process type':
            key = info.process_description
        elif marker_by == 'experiment':
            key = str(experiment)
        elif marker_by == 'dataset':
            key = info.dataset_label
        else:
            raise ValueError('Unknown marker_by value')

        #TODO: This is an ugly check. Is there a way to do it with .setdefault
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
            #This is to get the label key
            coords = [], []
        ax.plot(*coords,
            label=key,
            markeredgewidth=1,
            markeredgecolor=None,
            **key_options[key],
        )

    #Iterate again so highlights are printed on top.
    for key in xh:
        ax.plot(np.concatenate(xh[key]), np.concatenate(q2h[key]),
            markeredgewidth=0.6,
            markeredgecolor="black",
            **key_options[key],
        )
    if xh:
        #Get legend key
        ax.plot([], [], marker='s', markeredgewidth=0.6, color='none',
            markersize=5,
            markeredgecolor="black", label= f'Black edge: {highlight_label}',
        )

    if display_cuts:
        ax.scatter(np.concatenate(filteredx), np.concatenate(filteredq2),
            marker='o',
            facecolors='none', edgecolor='red', s=40, lw=0.8, label="Cut"
        )

    ax.set_title("Kinematic coverage")
    ax.legend()
    ax.set_xlabel('$x$')
    ax.set_ylabel(r'$Q^2$ (GeV$^2$)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return fig
