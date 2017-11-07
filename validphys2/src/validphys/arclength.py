"""
arclength.py

Module for the computation and presentation of arclengths
As of yet just works in it's own dedicated basis, needs
to be extended to use the Basis class.
"""
from collections import namedtuple, Sequence
import numbers
import math

import numpy as np
import pandas as pd

from reportengine.figure import figure
from reportengine.table import table
from reportengine.checks import check_positive

from validphys.plots import check_normalize_to
from validphys.core  import PDF

import matplotlib.pyplot as plt

FLAVOUR_ALIAS = {
    'gluon' : 21,
    'u'     :  2,
    'd'     :  1,
    's'     :  3,
    'ubar'  : -2,
    'dbar'  : -1,
    'sbar'  : -3}

ArcLengthGrid = namedtuple('ArcLengthGrid', FLAVOUR_ALIAS)
@check_positive('Q')
#TODO convert to use Basis (probably through xplotting_grid to concentrate xfx calls in one loc)
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
@check_normalize_to
def plot_arc_lengths(pdfs:Sequence, Q:numbers.Real, normalize_to:(type(None),int)=None):
    """Plot the arc lengths of a PDF set"""
    fig, ax = plt.subplots()

    # I'm assuming that the dictionary.items() returns always in the same order here...
    norm_cv = None
    if normalize_to is not None:
        norm_al = arc_lengths(pdfs[normalize_to], Q)
        norm_cv = [ np.mean(al) for fl, al in norm_al._asdict().items() ]

    for ipdf, pdf in enumerate(pdfs):
        arclengths = arc_lengths(pdf, Q)
        aldict = arclengths._asdict()
        xvalues = np.array(range(len(aldict)))
        xlabels  = [ fl for fl in aldict.keys()]
        yvalues = [ np.mean(al) for fl, al in aldict.items()]
        yerrors = [ np.std(al)  for fl, al in aldict.items()]
        if norm_cv is not None:
            yvalues = np.divide(yvalues, norm_cv)
            yerrors = np.divide(yerrors, norm_cv)
        #TODO should do a better job of this shift
        ax.errorbar(xvalues + ipdf/5.0, yvalues, yerr = yerrors, fmt='.', label=pdf.label)
        ax.set_xticks(xvalues)
        ax.set_xticklabels(xlabels)
        ax.legend()
    return fig
