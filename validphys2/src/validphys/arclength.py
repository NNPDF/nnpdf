"""
arclength.py

Module for the computation and presentation of arclengths.
"""
from collections import namedtuple, Sequence
import numbers
import math

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.figure import figure
from reportengine.table  import table
from reportengine.checks import check_positive, make_argcheck

from validphys.pdfbases import Basis, check_basis, PDG_PARTONS
from validphys.pdfgrids import (xgrid, xplotting_grid)
from validphys.plots    import check_pdf_normalize_to
from validphys.core     import PDF
from validphys.checks   import check_pdf_is_montecarlo

import matplotlib.pyplot as plt

ArcLengthGrid = namedtuple('ArcLengthGrid', ('pdf', 'basis','scale', 'flavours', 'values'))
@check_positive('Q')
@make_argcheck(check_basis)
@check_pdf_is_montecarlo
def arc_lengths(pdf:PDF, Q:numbers.Real,
                basis:(str, Basis)='flavour',
                flavours:(list, tuple, type(None))=None):
    """Compute arc lengths at scale Q"""
    lpdf = pdf.load()
    checked = check_basis(basis, flavours)
    basis, flavours = checked['basis'], checked['flavours']
    # Using the same questionable integration as validphys
    # x-grid points and limits in three segments
    npoints = 199 # 200 intervals
    seg_min = [1e-6, 1e-4, 1e-2]
    seg_max = [1e-4, 1e-2, 1.0 ]
    res = np.zeros((len(flavours), lpdf.GetMembers()))
    # Integrate the separate segments
    for iseg, seg in enumerate(seg_min):
        a, b = seg, seg_max[iseg]                                # Integration limits
        eps = (b-a)/npoints                                      # Integration step-size
        x1grid = xgrid(a, b, 'linear', npoints)                  # Integration grid x1
        f1grid = xplotting_grid(pdf, Q, x1grid, basis, flavours) # PDFs evaluated at x1
        dfgrid = f1grid.grid_values                              # Backwards-difference
        np.swapaxes(dfgrid, 1,2)
        for irep in range(lpdf.GetMembers()):
            for ifl,fl in enumerate(flavours):
                dslice = np.diff(dfgrid[irep][ifl] * x1grid[1])
                asqr = np.square(dslice)
                res[ifl,irep] += np.sum(np.sqrt(eps*eps+asqr))
    return ArcLengthGrid(pdf, basis, Q, flavours, res)

# Collect arc_lengths over PDF list
pdfs_arc_lengths = collect(arc_lengths,['pdfs'])

@table
def arc_length_table(arc_lengths):
    """Return a table with the descriptive statistics of the arc lengths
    over members of the PDF."""
    arc_length_transpose = np.transpose(arc_lengths.values)
    arc_length_columns = ['$'+arc_lengths.basis.elementlabel(fl)+'$' for fl in arc_lengths.flavours]
    return pd.DataFrame(arc_length_transpose, columns=arc_length_columns).describe().iloc[1:,:]

@figure
@check_pdf_normalize_to
def plot_arc_lengths(pdfs_arc_lengths:Sequence, Q:numbers.Real, normalize_to:(type(None),int)=None):
    """Plot the arc lengths of provided pdfs"""
    fig, ax = plt.subplots()
    if normalize_to is not None:
        ax.set_ylabel("Arc length $Q="+str(Q)+"$ GeV (normalised)")
    else:
        ax.set_ylabel("Arc length $Q="+str(Q)+"$ GeV")

    norm_cv = None
    if normalize_to is not None:
        norm_al = pdfs_arc_lengths[normalize_to].values
        norm_cv = [ np.mean(fl) for fl in norm_al ]

    for ipdf, arclengths in enumerate(pdfs_arc_lengths):
        alvalues = arclengths.values
        xvalues = np.array(range(len(arclengths.flavours)))
        xlabels  = [ '$'+arclengths.basis.elementlabel(fl)+'$' for fl in arclengths.flavours]
        yvalues = [ np.mean(fl) for fl in alvalues ]
        # Computation of 68% C.I - this feels like it should be somewhere else
        nmembers = np.shape(alvalues)[1]
        nrep_16 = math.floor(0.16*nmembers)
        nrep_84 = math.ceil(0.84*nmembers)
        srvalues = [np.sort(fl) for fl in alvalues] # Replicas sorted
        yupper = [fl[nrep_84] - yvalues[ifl] for ifl, fl in enumerate(srvalues)]
        ylower = [yvalues[ifl] - fl[nrep_16] for ifl, fl in enumerate(srvalues)]
        if norm_cv is not None:
            yvalues = np.divide(yvalues, norm_cv)
            yupper  = np.divide(yupper,  norm_cv)
            ylower  = np.divide(ylower,  norm_cv)
        ##TODO should do a better job of this shift
        ax.errorbar(xvalues + ipdf/5.0, yvalues, yerr = (ylower,yupper), fmt='.', label=arclengths.pdf.label)
        ax.set_xticks(xvalues)
        ax.set_xticklabels(xlabels)
        ax.legend()
    return fig
