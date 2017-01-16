# -*- coding: utf-8 -*-
"""
Figures for visualizing results
"""
from __future__ import generator_stop

import logging
import functools
import warnings

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import scipy.stats as stats

from reportengine.figure import figure, figuregen
from reportengine.utils import saturate
from reportengine.checks import make_check, CheckError, make_argcheck

from validphys.core import MCStats
from validphys.results import chi2_stat_labels
from validphys.pdfgrids import PDG_PARTONS
from validphys.plotoptions import get_infos, kitable, transform_result
from validphys.checks import check_scale
from validphys import plotutils

log = logging.getLogger(__name__)

@figure
def plot_chi2dist(results, dataset, abs_chi2_data, chi2_stats, pdf):
    """Plot the distribution of chi²s of the members of the pdfset."""
    setlabel = dataset.name
    fig, ax = plt.subplots()
    label = pdf.name
    alldata, central, npoints = abs_chi2_data
    if not isinstance(alldata, MCStats):
        ax.set_axis_bgcolor("#ffcccc")
        log.warn("Chi² distribution plots have a different meaning for non MC sets.")
        label += " (%s!)" % pdf.ErrorType
    label += '\n'+ '\n'.join(str(chi2_stat_labels[k])+(' %.2f' % v) for (k,v) in chi2_stats.items())
    ax.set_title("$\chi^2$ distribution for %s" % setlabel)

    ax.hist(alldata.data, label=label, zorder=100)
    l = ax.legend()
    l.set_zorder(1000)
    return fig

#TODO: This should be simplified if at all possible. For now some more examples
#are needed for a spec to emerge.
@make_check
def _check_normalize_to(ns, **kwargs):
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


@_check_normalize_to
@figuregen
def plot_fancy(one_or_more_results, dataset,
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


    results = one_or_more_results

    infos = get_infos(dataset)
    for info in infos:
        table = kitable(dataset, info)
        nkinlabels = len(table.columns)

        if normalize_to is not None:
            norm_result = results[normalize_to]
            #We modify the table, so we pass only the label columns
            norm_cv, _ = transform_result(norm_result.central_value,
                                       norm_result.std_error,
                                       table.iloc[:,:nkinlabels], info)

        for i,result in enumerate(results):
            #We modify the table, so we pass only the label columns
            cv, err = transform_result(result.central_value, result.std_error,
                                       table.iloc[:,:nkinlabels], info)
            #By doing tuple keys we avoid all possible name collisions
            if normalize_to is None:
                table[('cv', i)] = cv
                table[('err', i)] = err
            else:
                table[('cv', i)] = cv/norm_cv
                table[('err', i)] = err/norm_cv
        if info.figure_by:
            figby = table.groupby(info.figure_by)
        else:
            figby = [('', table)]

        for samefig_vals, fig_data in figby:
            #For some reason matplotlib doesn't set the axis right
            min_vals = []
            max_vals = []
            #Have consistent output for one or more groupby columns
            if not isinstance(samefig_vals, tuple):
                samefig_vals = (samefig_vals, )
            fig, ax = plt.subplots()
            plotutils.setup_ax(ax)
            ax.set_title("%s %s"%(dataset.name,
                         info.group_label(samefig_vals, info.figure_by)))
            if info.line_by:
                lineby = fig_data.groupby(info.line_by)
            else:
                lineby = [('', fig_data)]

            first = True



            #http://matplotlib.org/users/transforms_tutorial.html
            for (sameline_vals, line_data) in lineby:
                ax.set_prop_cycle(None)
                if first:
                    labels = True
                else:
                    labels = False
                first = False
                if not isinstance(sameline_vals, tuple):
                    sameline_vals = (sameline_vals, )

                nres = len(results)
                first_offset = +(nres//2)

                if info.x == 'idat':
                    x = np.array(line_data.index)
                else:
                    x = line_data[info.x].as_matrix()

                try:
                    x = np.asanyarray(x, np.float)
                except ValueError:
                    xticklabels = x
                    npoints = len(x)
                    x = np.arange(npoints)
                    ax.set_xticks(x)
                    ax.set_xticklabels(xticklabels)
                    #TODO: Remove this when mpl stops doing the idiotic thing
                    #(in v2?)
                    ax.set_xlim(-npoints/20, npoints - 1+ npoints/20)


                #Use black for the first iteration (data),
                #and follow the cycle for
                #the rest.
                color = '#262626'
                for i, res in enumerate(results):

                    if labels:
                        label = res.label
                    else:
                        label = None
                    cv = line_data[('cv', i)].as_matrix()
                    err = line_data[('err', i)].as_matrix()


                    dx, dy = 0.05*(i-first_offset), 0.
                    offset = transforms.ScaledTranslation(dx, dy,
                                                          fig.dpi_scale_trans)
                    offset_transform = ax.transData + offset

                    max_vals.append(np.nanmax(cv+err))
                    min_vals.append(np.nanmin(cv-err))


                    ax.errorbar(x, cv, yerr=err,
                         linestyle='--',
                         lw=0.25,
                         label= label,
                         #elinewidth = 2,
                         #capsize=10,
                         marker = 's',
                         #markeredgewidth=0.25,
                         c=color,
                         zorder=1000,
                         transform=offset_transform)

                    color = None

            total_extremes = min(min_vals), max(max_vals)
            small_lim, big_lim = plotutils.expand_margin(*total_extremes, 1.2)
            ax.set_ylim(small_lim, big_lim)


            glabel = info.group_label(sameline_vals, info.line_by)
            annotate_point = x[-1], line_data[('cv', 0)].as_matrix()[-1]
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
                ax.set_ylabel("Ratio to %s" % norm_result.label)

            ax.legend().set_zorder(100000)
            ax.set_xlabel(info.xlabel)




            fig.tight_layout()


            yield fig


@_check_normalize_to
@figuregen
def plot_bands(one_or_more_results, dataset,
               normalize_to:(int,str,type(None)) = None):


    results = one_or_more_results

    nnpdf_dt = dataset.load()
    if not dataset.commondata.plotfiles:
        infos = [get_info(nnpdf_dt)]
    else:
        infos = []
        for p in dataset.commondata.plotfiles:
            with p.open() as f:
                infos.append(get_info(nnpdf_dt, f, cuts=dataset.cuts))
    for info in infos:
        table = kitable(nnpdf_dt, info)
        nkinlabels = len(table.columns)

        if normalize_to is not None:
            norm_result = results[normalize_to]
            #We modify the table, so we pass only the label columns
            norm_cv, _ = transform_result(norm_result.central_value,
                                       norm_result.std_error,
                                       table.iloc[:,:nkinlabels], info)

        for i,result in enumerate(results):
            #We modify the table, so we pass only the label columns
            cv, err = transform_result(result.central_value, result.std_error,
                                       table.iloc[:,:nkinlabels], info)
            #By doing tuple keys we avoid all possible name collisions
            if normalize_to is None:
                table[('cv', i)] = cv
                table[('err', i)] = err
            else:
                table[('cv', i)] = cv/norm_cv
                table[('err', i)] = err/norm_cv
        if info.figure_by:
            figby = table.groupby(info.figure_by)
        else:
            figby = [('', table)]
        for samefig_vals, fig_data in figby:
            #Have consistent output for one or more groupby columns
            if not isinstance(samefig_vals, tuple):
                samefig_vals = (samefig_vals, )
            fig, ax = plt.subplots()
            plotutils.setup_ax(ax)
            ax.set_title("%s %s"%(dataset.name,
                         info.group_label(samefig_vals, info.figure_by)))
            if info.line_by:
                lineby = fig_data.groupby(info.line_by)
            else:
                lineby = [('', fig_data)]

            first = True


            #http://matplotlib.org/users/transforms_tutorial.html
            for (sameline_vals, line_data) in lineby:
                ax.set_prop_cycle(None)
                if first:
                    labels = True
                else:
                    labels = False
                first = False
                if not isinstance(sameline_vals, tuple):
                    sameline_vals = (sameline_vals, )

                nres = len(results)
                first_offset = -nres//2

                if info.x == 'idat':
                    x = np.array(line_data.index)
                else:
                    x = line_data[info.x].as_matrix()

                #Use black for the first iteration (data),
                #and follow the cycle for
                #the rest.
                color = '#262626'
                elw = None
                mark = 3
                markw = 0.5
                ls = '--'
                al = 1

                for i, res in enumerate(results):

                    if labels:
                        label = res.label
                    else:
                        label = None
                    cv = line_data[('cv', i)].as_matrix()
                    err = line_data[('err', i)].as_matrix()

                    dx, dy = 0.05*(i-first_offset), 0.
                    offset = transforms.ScaledTranslation(dx, dy,
                                                          fig.dpi_scale_trans)
                    offset_transform = ax.transData + offset
                    ax.errorbar(x, cv, yerr=err,
                         linestyle=ls,
                         lw=.25,
                         label= label,
                         elinewidth = elw,
                         capsize=2,
                         markersize = mark,
                         marker = 's',
                         markeredgewidth=markw,
                         c=color,
                         alpha = al,
                         zorder=1000,
                         transform=offset_transform)

                    markw = None
                    color = None
                    elw = 20
                    mark = 0.
                    ls = ''
                    al = 0.6

                glabel = info.group_label(sameline_vals, info.line_by)
                annotate_point = x[-1], line_data[('cv', 0)].as_matrix()[-1]
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
                ax.set_ylabel("Ratio to %s" % norm_result.label)

            ax.legend().set_zorder(100000)
            ax.set_xlabel(info.xlabel)

            fig.tight_layout()


            yield fig


def _scatter_marked(ax, x, y, marked_dict, *args, **kwargs):
    kwargs['s'] = kwargs.get('s', 30) + 10
    x = np.array(x, copy=False)
    y = np.array(y, copy=False)
    for label, indexes in marked_dict.items():
        ax.scatter(x[indexes],y[indexes], *args, **kwargs, label=label,
                   facecolors='none', linewidth=0.5, edgecolor='red')
        kwargs['s'] += 10


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

    ax.set_title(getattr(fit, 'label', fit.name))

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
    ax.set_title("KDE of the fit distributions for %s" % getattr(fit,
                                                                 'label', fit.name))

    plt.legend()
    return fig

@figure
def plot_covmat_eigs(experiment):
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
    covmat = experiment.load().get_covmat()
    stds = np.sqrt(np.diag(covmat))
    corrmat = covmat/np.outer(stds,stds)
    eigs = la.eigvalsh(corrmat)
    fig,ax = plt.subplots()
    ax.plot(eigs, 'o')
    ax.set_yscale('log')
    return fig


#The indexing to one instead of zero is so that we can be consistent with
#how plot_fancy works, so normalize_to: 1 would normalize to the first pdf
#for both.
@make_argcheck
def _check_pdf_normalize_to(pdfs, normalize_to):
    """Transforn normalize_to into an index."""

    msg = ("normalize_to should be, a pdf id or an index of the "
           "pdf (starting from one)")

    if normalize_to is None:
        return

    names = [pdf.name for pdf in pdfs]
    if isinstance(normalize_to, int):
        normalize_to -= 1
        if not normalize_to < len(names) or normalize_to<0:
            raise CheckError(msg)
        return {'normalize_to': normalize_to}

    if isinstance(normalize_to, str):
        try:
            normalize_to = names.index(normalize_to)
        except ValueError:
            raise CheckError(msg, normalize_to, alternatives=names)
        return {'normalize_to': normalize_to}


    raise RuntimeError("Should not be here")

def _plot_pdf_factory(draw_function, setup_function=None, legend_function=None):
    """f does the actual plotting, and returns data used to compute the axis
    limits. It is called like f(ax, pdf, flindex ,grid)"""
    @figuregen
    @check_scale('xscale', allow_none=True)
    @_check_pdf_normalize_to
    def f_(pdfs, xplotting_grids, xscale:(str,type(None))=None,
           normalize_to:(int,str,type(None))=None):
        if not xplotting_grids:
            return

        if normalize_to is not None:
            normalize_pdf = pdfs[normalize_to]
            normalize_grid = xplotting_grids[normalize_to]
            normvals = normalize_pdf.stats_class(
                            normalize_grid.grid_values).central_value()
            def fp_error(tp, flag):
                log.warn("Invalid values found computing normalization to %s: "
                 "Floating point error (%s).", normalize_pdf, tp)
                #Show warning only once
                np.seterr(all='ignore')
            newgrids = []
            with np.errstate(all='call'):
                np.seterrcall(fp_error)
                for grid in xplotting_grids:
                    newvalues = grid.grid_values/normvals
                    #newgrid is like the old grid but with updated values
                    newgrid = type(grid)(**{**grid._asdict(),
                                             'grid_values':newvalues})
                    newgrids.append(newgrid)
            xplotting_grids = newgrids
            ylabel = "Ratio to {}".format(normalize_pdf.label)
        else:
            ylabel = None

        firstgrid = xplotting_grids[0]
        if xscale is None:
            xscale = 'linear' if firstgrid.scale == 'linear' else 'log'
        Q = firstgrid.Q

        for flindex, fl in enumerate(firstgrid.flavours):
            fig, ax = plt.subplots()
            if setup_function:
                setupres =  saturate(setup_function, locals())
            ax.set_title("$%s$ at %.1f GeV" % (PDG_PARTONS[fl], Q))

            all_vals = []
            for pdf, grid in zip(pdfs, xplotting_grids):
                all_vals.append(saturate(draw_function, locals()))

            #Note these two lines do not conmute!
            ax.set_xscale(xscale)
            plotutils.frame_center(ax, firstgrid.xgrid, np.concatenate(all_vals))

            ax.set_xlabel('x')
            parton_name = PDG_PARTONS[fl]
            if ylabel:
                ax.set_ylabel(ylabel)
            else:
                ax.set_ylabel('$x{}(x)$'.format(parton_name))

            ax.set_axisbelow(True)

            if legend_function:
                saturate(legend_function, locals())
            else:
                ax.legend()
            yield fig, parton_name

    #Keep the signature of f_ instead of that of f, and also keep annotations
    #for type checking
    functools.update_wrapper(f_, draw_function,
            assigned=[a for a in
                      functools.WRAPPER_ASSIGNMENTS if a!='__annotations__'])
    del f_.__wrapped__

    return f_

@make_argcheck
def _warn_any_pdf_not_montecarlo(pdfs):
    for pdf in pdfs:
        et = pdf.ErrorType
        if et != 'replicas':
            log.warn("Plotting members of a non-Monte Carlo PDF set:"
            " %s with error type '%s'.", pdf.name, et)


@_warn_any_pdf_not_montecarlo
@_plot_pdf_factory
def plot_pdfreplicas(ax, pdf, flindex ,grid):
    """Plot the replicas of the specifid PDFs.

    - xscale sets the scale of the plot. E.g. 'linear' or 'log'. Default is
    deduced from the xplotting_grid, which in turn is 'log' by default.

    - normalize_to should be, a pdf id or an index of the pdf (starting from one).
    """
    next_prop = next(ax._get_lines.prop_cycler)
    color = next_prop['color']
    gv = grid.grid_values[:,flindex,:]


    ax.plot(grid.xgrid, gv.T, alpha=0.2, linewidth=0.5,
            color=color, zorder=1)
    stats = pdf.stats_class(gv)
    ax.plot(grid.xgrid, stats.central_value(), color=color,
            linewidth=2,
            label=pdf.label)
    return gv

#Because of how pickle works, this has to have this name, and then be redefined
#Otherwise will complain about not veing able to pickle the inner f_ in the
#decorator.
def plot_pdfs(ax, pdf, flindex, grid, setupres):
    hatchit, labels, handles = (setupres['hatchit'], setupres['labels'],
                                setupres['handles'])
    stats = pdf.stats_class(grid.grid_values[:,flindex,:])
    pcycler = ax._get_lines.prop_cycler
    #This is ugly but can't think of anything better

    next_prop = next(pcycler)
    cv = stats.central_value()
    xgrid = grid.xgrid
    #Ignore spurious normalization warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        err68down, err68up = stats.errorbar68()

    color = next_prop['color']
    ax.plot(xgrid, cv, color=color)
    alpha = 0.5
    ax.fill_between(xgrid, err68up, err68down, color=color, alpha=alpha,
                    zorder=1)
    #http://stackoverflow.com/questions/5195466/matplotlib-does-not-display-hatching-when-rendering-to-pdf
    hatch = next(hatchit)
    ax.fill_between(xgrid, err68up, err68down, facecolor='None', alpha=alpha,
                    edgecolor=color, hatch=hatch, zorder=1)
    if isinstance(stats, MCStats):
        errorstdup, errorstddown = stats.errorbarstd()
        ax.plot(xgrid, errorstdup, linestyle='--', color=color)
        ax.plot(xgrid, errorstddown, linestyle='--', color=color)
        label  = "%s ($68\%%$ c.l.+$1\sigma$)" % pdf.label
        outer = True
    else:
        outer = False
        label = "%s ($68\%%$ c.l.)" % pdf.label
    handle = plotutils.HandlerSpec(color=color, alpha=alpha,
                                           hatch=hatch, outer=outer)
    handles.append(handle)
    labels.append(label)

    return [err68down, err68up]

def _plot_pdfs_setup():
    return dict(handles=[], labels=[],
                           hatchit=plotutils.hatch_iter())

def _plot_pdfs_legend(ax, setupres):
    labels, handles = (setupres['labels'], setupres['handles'])
    return ax.legend(handles, labels, handler_map={plotutils.HandlerSpec:
                                             plotutils.ComposedHandler()})

plot_pdfs = _plot_pdf_factory(plot_pdfs, _plot_pdfs_setup, _plot_pdfs_legend)
