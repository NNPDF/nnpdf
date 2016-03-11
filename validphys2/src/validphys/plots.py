# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:58:14 2016

@author: Zahari Kassabov
"""
import logging

import numpy as np
import matplotlib.pyplot as plt

from reportengine.figure import figure

from validphys.core import MCStats

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
