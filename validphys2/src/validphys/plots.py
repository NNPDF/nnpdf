# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:58:14 2016

@author: Zahari Kassabov
"""
from __future__ import generator_stop

import logging

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import scipy.stats as stats

from reportengine.figure import figure, figuregen
from reportengine.checks import make_check, CheckError

from validphys.core import MCStats
from validphys.results import chi2_stat_labels
from validphys.plotoptions import get_info, kitable, transform_result
from validphys import plotutils

log = logging.getLogger(__name__)

@figure
def plot_chi2dist(results, dataset, abs_chi2_data, chi2_stats, pdf):
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
def check_normalize_to(callspec, ns, graph, **kwargs):
    """Transforn normalize_to into an index."""

    msg = ("normalize_to should be either 'data', a pdf id or an index of the "
           "result (0 for the data, and i for the ith pdf)")

    val = ns.get('normalize_to', None)
    if val is None:
        return

    if ns.get('pdf', False):
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


@check_normalize_to
@figuregen
def plot_fancy(one_or_more_results, dataset,
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

@figure
def plot_training_validation(fit, replica_data):
    training, valid = zip(*((dt.training, dt.validation) for dt in replica_data))
    fig, ax = plt.subplots()
    s=ax.plot(training,valid,linestyle='', marker='o', markersize=5, zorder=100)

    ax.set_title(getattr(fit, 'label', fit.name))

    ax.set_xlabel(r'$\chi^2/N_{dat}$ train')
    ax.set_ylabel(r'$\chi^2/N_{dat}$ valid')

    ax.plot(np.mean(training), np.mean(valid),
         marker='s', color='red', markersize=7, zorder=1000)

    return fig


@figure
def plot_trainvaliddist(fit, replica_data):
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



