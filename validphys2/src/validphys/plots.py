# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:58:14 2016

@author: Zahari Kassabov
"""
import logging

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

from reportengine.figure import figure, figuregen
from reportengine.checks import make_check, CheckError

from validphys.core import MCStats
from validphys.plotoptions import get_info, kitable, transform_result
from validphys import plotutils

log = logging.getLogger(__name__)

@figure
def plot_chi2dist(results, dataset, chi2_data, chi2_stats, pdf):
    setlabel = dataset.name
    fig, ax = plt.subplots()
    label = pdf.name
    if not isinstance(chi2_data, MCStats):
        ax.set_axis_bgcolor("#ffcccc")
        log.warn("Chi² distribution plots have a different meaning for non MC sets.")
        label += " (%s!)" % pdf.ErrorType
    label += '\n'+ ', '.join(str(k)+(' %.2f' % v) for (k,v) in chi2_stats.items())
    ax.set_title("$\chi^2$ distribution for %s" % setlabel)
    ax.hist(chi2_data.data, label=label, zorder=10000)
    ax.legend()
    return fig

@figure
def plot_results(results, dataset, normalize_to = None):

    setlabel = dataset.name

    cvs = [r.central_value for r in results]
    errors = [r.std_error for r in results]

    y_label = "Observable value"

    if normalize_to is not None:
        y_label = "Ratio to %s" % normalize_to.label
        norm_cv = normalize_to.central_value
        cvs = [cv/norm_cv for cv in cvs]
        errors = [e/norm_cv for e in errors]

    l = len(results)
    if l < 5:
        delta = iter(np.linspace(-0.05*l, 0.05*l, l))
    else:
        delta = iter(np.linspace(-0.25, 0.25, l))

    fig, ax = plt.subplots()

    ax.set_title(setlabel)
    ax.set_ylabel(y_label)

    for result, cv, error in zip(results, cvs, errors):
        x = np.arange(1, len(result) + 1) + next(delta)
        ax.errorbar(x, cv, yerr=error,
                         linestyle='none',
                         label=result.label, elinewidth = 2,
                         capsize=10)

    ax.set_xticks(range(1,len(result) + 1), len(result) // 15 + 1)

    ax.grid(axis='x')

    #TODO: Abtract this out
    dfs = ax.get_yticks() - 1
    l = len(dfs) // 2  + 1 - ((len(dfs) // 2) % 2)
    mdiff = np.max(dfs)
    ax.set_yticks(np.linspace(-mdiff, mdiff, l) + 1)



    ax.legend()

    return fig

#TODO: This should be simplified if at all possible. For now some more examples
#are needed for a spec to emerge.
@make_check
def check_normalize_to(callspec, ns, graph):
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
def plot_fancy(one_or_more_results, dataset, normalize_to = None):


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
