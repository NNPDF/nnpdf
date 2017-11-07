"""
arclength.py

Module for the computation of arclengths
As of yet just works in it's own dedicated basis, needs
to be extended to use the Basis class.
"""
from collections import namedtuple
import numbers
import math

import numpy as np
import pandas as pd

from reportengine.figure import figure, figuregen
from reportengine.table import table
from reportengine.checks import check_positive

from validphys.core import PDF

import matplotlib.pyplot as plt

FLAVOUR_ALIAS = {
    'gluon': 21,
    'cbar'    : -4,
    'dbar'    : -1,
    'sbar'   : -3,
    'ubar'   : -2,
    'd': 1,
    'u': 2,
    's': 3,
    'c': 4      }

ArcLengthGrid = namedtuple('ArcLengthGrid', FLAVOUR_ALIAS)
@check_positive('Q')
def arc_lengths(pdf:PDF, Q:numbers.Real):
    """Compute arc lengths at scale Q"""
    lpdf = pdf.load()
    nmembers = lpdf.GetMembers()
    res = np.zeros((len(FLAVOUR_ALIAS), nmembers))
    flavours = FLAVOUR_ALIAS.values()
    def integral(fl, a, b, irep):
        # Using the same questionable integrator as vp1
        eps = (b-a)/1000.
        bdiff = lambda x: lpdf.xfxQ(x, Q=Q, n=irep, fl=fl) - lpdf.xfxQ(x-eps, Q=Q, n=irep, fl=fl)
        xvals = np.arange(a,b,eps)
        gvals = [math.sqrt(eps*eps+math.pow(x*bdiff(x),2)) for x in xvals]
        return sum(gvals)
    for irep in range(nmembers):
        for ifl, fl in enumerate(flavours):
            res[ifl,irep] = ( integral(fl, 1e-7, 1e-5, irep)
                              +integral(fl, 1e-5, 1e-3, irep)
                              +integral(fl, 1e-3, 1, irep))
    return ArcLengthGrid(*res)

@table
def arc_length_table(arc_lengths):
    """Return a table with the descriptive statistics of the arc lengths
    over members of the PDF."""
    #We don't  really want the count, which is going to be the same for all.
    #Hence the .iloc[1:,:].
    return pd.DataFrame(arc_lengths._asdict()).describe().iloc[1:,:]

@figure
def plot_arc_lengths(arc_lengths):
    """Plot the arc lengths of a PDF set"""
    fig, ax = plt.subplots()

    aldict = arc_lengths._asdict()
    xvalues = range(len(aldict))
    yvalues = []
    yerr    = []
    labels  = []
    for fl, al in aldict.items():
        labels.append(fl)
        yvalues.append(np.mean(al))
        yerr.append(np.std(al))
    ax.errorbar(xvalues, yvalues, yerr = yerr, fmt='.')
    ax.set_xticklabels(labels)
    #l = ax.legend()
    return fig

