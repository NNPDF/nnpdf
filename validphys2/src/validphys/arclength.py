"""
arclength.py

Module for the computation and presentation of arclengths.
"""
from collections import namedtuple
from collections.abc import Sequence
import numbers

import numpy as np
import pandas as pd
import scipy.integrate as integrate

from reportengine import collect
from reportengine.figure import figure
from reportengine.table import table
from reportengine.checks import check_positive, make_argcheck

from validphys.pdfbases import Basis, check_basis
from validphys.pdfgrids import xgrid, xplotting_grid
from validphys.core import PDF
from validphys.checks import check_pdf_normalize_to

import matplotlib.pyplot as plt

ArcLengthGrid = namedtuple("ArcLengthGrid", ("pdf", "basis", "flavours", "stats"))


@check_positive("Q")
@make_argcheck(check_basis)
def arc_lengths(
    pdf: PDF,
    Q: numbers.Real,
    basis: (str, Basis) = "flavour",
    flavours: (list, tuple, type(None)) = None,
    xscale:(str,type(None)) = "linear"
):
    """Compute arc lengths at scale Q"""
    checked = check_basis(basis, flavours)
    basis, flavours = checked["basis"], checked["flavours"]
    # x-grid points and limits in three segments
    npoints = 199  # 200 intervals
    seg_min = [1e-6, 1e-4, 1e-2]
    seg_max = [1e-4, 1e-2, 1.0]
    res = np.zeros((pdf.get_members(), len(flavours)))
    # Integrate the separate segments
    for a, b in zip(seg_min, seg_max):
        # Finite diff. step-size, x-grid
        eps = (b - a) / npoints
        ixgrid = xgrid(a, b, xscale, npoints)
        # PDFs evaluated on grid, use the entire thing, the Stats class will chose later
        xfgrid = xplotting_grid(pdf, Q, ixgrid, basis, flavours).grid_values.data * ixgrid[1]
        fdiff = np.diff(xfgrid) / eps  # Compute forward differences
        res += integrate.simps(1 + np.square(fdiff), ixgrid[1][1:])
    stats = pdf.stats_class(res)
    return ArcLengthGrid(pdf, basis, flavours, stats)


# Collect arc_lengths over PDF list
pdfs_arc_lengths = collect(arc_lengths, ["pdfs"])


@table
def arc_length_table(arc_lengths):
    """Return a table with the descriptive statistics of the arc lengths
    over members of the PDF."""
    arc_length_data = arc_lengths.stats.error_members()
    arc_length_columns = [
        f"${arc_lengths.basis.elementlabel(fl)}$" for fl in arc_lengths.flavours
    ]
    return (
        pd.DataFrame(arc_length_data, columns=arc_length_columns).describe().iloc[1:, :]
    )


@figure
@check_pdf_normalize_to
def plot_arc_lengths(
    pdfs_arc_lengths: Sequence, Q: numbers.Real, normalize_to: (type(None), int) = None
):
    """Plot the arc lengths of provided pdfs"""
    fig, ax = plt.subplots()
    if normalize_to is not None:
        ax.set_ylabel(f"Arc length $Q={Q}$ GeV (normalised)")
    else:
        ax.set_ylabel(f"Arc length $Q={Q}$ GeV")

    for ipdf, arclengths in enumerate(pdfs_arc_lengths):
        xvalues = np.array(range(len(arclengths.flavours)))
        xlabels = [
            "$" + arclengths.basis.elementlabel(fl) + "$" for fl in arclengths.flavours
        ]
        yvalues = arclengths.stats.central_value()
        ylower, yupper = arclengths.stats.errorbar68()
        ylower = yvalues - ylower
        yupper = yupper - yvalues
        if normalize_to is not None:
            norm_cv = pdfs_arc_lengths[normalize_to].stats.central_value()
            yvalues = np.divide(yvalues, norm_cv)
            yupper = np.divide(yupper, norm_cv)
            ylower = np.divide(ylower, norm_cv)
        shift = (ipdf - (len(pdfs_arc_lengths) - 1) / 2.0) / 5.0
        ax.errorbar(
            xvalues + shift,
            yvalues,
            yerr=(ylower, yupper),
            fmt=".",
            label=arclengths.pdf.label,
        )
        ax.set_xticks(xvalues)
        ax.set_xticklabels(xlabels)
    ax.legend()
    return fig


# TODO: this should probably go somewhere else
def integrability_number(
    pdf: PDF,
    Q: numbers.Real,
    basis: (str, Basis) = "evolution",
    flavours: (list, tuple, type(None)) = None,
):
    """Return \sum_i |x_i*f(x_i)|, x_i = {1e-9, 1e-8, 1e-7} 
    for selected flavours 
    """
    checked = check_basis(basis, flavours)
    basis, flavours = checked["basis"], checked["flavours"]
    ixgrid = xgrid(1e-9, 1e-6, "log", 3)
    xfgrid = xplotting_grid(pdf, Q, ixgrid, basis, flavours).grid_values.data
    res = np.sum(np.abs(xfgrid), axis=2)
    return res.squeeze()
