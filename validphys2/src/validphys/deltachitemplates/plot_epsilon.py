"""Code I added to validphys/pdfplots.py to plot the Gaussian 
deviation estimator for all the pdfs in a Monte Carlo set.
see https://arxiv.org/abs/1505.06736 eq. (11)
"""

import logging
import warnings

import numpy as np
import matplotlib.pyplot as plt

from reportengine.figure import figuregen

from validphys import plotutils
from validphys.checks import check_scale
from validphys.pdfplots import FlavourState

log = logging.getLogger(__name__)


@figuregen
@check_scale("xscale", allow_none=True)
def plot_epsilon(
    pdfs, xplotting_grids, xscale: (str, type(None)) = None, ymin=None, ymax=None, eps=None
):
    """Plot the discrepancy (epsilon) of the 1-sigma and 68% bands at each grid value
    for all pdfs for a given Q.

    Of course this function is meaningful only if the pdfs are sets of MC replicas.
    (Still there is no check in the implementation)

    xscale is read from pdf plotting_grid scale, which is 'log' by default.

    eps defines the value at which plot a simple hline
    """

    if not xplotting_grids:
        return

    # xplotting_grids is a collection of xplotting_grid pdfs (collect defined in
    # reportengine.resourcebuilder). xplotting_grid is a tuple subclass (XPlottingGrid)
    # built with the collections.namedtuple function (see docs). xplotting_grid fields
    # are accessible by attribute lookup as in a usual class.

    # pick XPlottingGrid of the first pdf set just to extract useful information
    firstgrid = xplotting_grids[0]
    basis = firstgrid.basis  # this is an object pdfbases.Basis defined in validphys.pdfbases
    for flindex, fl in enumerate(firstgrid.flavours):
        fig, ax = plt.subplots()
        # return the label of the flavour number key passed
        parton_name = basis.elementlabel(fl)
        # create object with useful attributes
        flstate = FlavourState(flindex=flindex, fl=fl, fig=fig, ax=ax, parton_name=parton_name)
        _setup_flavour_label(flstate)
        ax.set_title("$%s$ at %.1f GeV" % (parton_name, firstgrid.Q))

        all_vals = []
        for pdf, grid in zip(pdfs, xplotting_grids):
            limits = _draw(pdf, grid, flstate)
            if limits is not None:
                all_vals.append(np.atleast_2d(limits))

        if xscale is None:
            ax.set_xscale(firstgrid.scale)
        else:
            ax.set_xscale(xscale)

        if eps is not None:
            xmin = np.min(xplotting_grids[0].xgrid)
            line2d = ax.axhline(
                eps, xmin, 1, alpha=0.5, lw=0.6, color="#666666", label="$\epsilon$=%.2f" % eps
            )
            flstate.labels.append(line2d.get_label())

        plotutils.frame_center(ax, firstgrid.xgrid, np.concatenate(all_vals))
        if ymin is not None:
            ax.set_ylim(bottom=ymin)
        if ymax is not None:
            ax.set_ylim(top=ymax)

        ax.set_xlabel("$x$")
        ax.set_xlim(firstgrid.xgrid[0])

        ax.set_ylabel("$\epsilon(x)$")

        ax.set_axisbelow(True)

        ax.legend(flstate.labels, handler_map={plotutils.HandlerSpec: plotutils.ComposedHandler()})

        plt.grid(linestyle="-", lw=0.5, alpha=0.4)

        yield fig, parton_name  # figure and string used by figuregen to save figures


def _setup_flavour_label(flstate):
    flstate.labels = []


def _draw(pdf, grid, flstate):

    ax = flstate.ax
    flindex = flstate.flindex
    labels = flstate.labels

    # pick all replicas grid_values for the flavour flindex. stats_class is a method
    # of the PDF class, which returns the stats calculator (object) for the pdf error type.
    # Basically stats is an object which says what is the type class of the replicas:
    # MCStats, HessianStats, SymmHessianStats. In this way validphys use the right methods
    # to compute statistical values.
    stats = pdf.stats_class(grid.grid_values[:, flindex, :])

    # Ignore spurious normalization warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        error68down, error68up = stats.errorbar68()  # 68% error bands

    errorstd = stats.std_error()
    # color cycle iterable
    pcycler = ax._get_lines.prop_cycler
    next_prop = next(pcycler)
    color = next_prop["color"]
    xgrid = grid.xgrid
    # the division by 2 is equivalent to considering the complete 1-sigma band (2 * error_std)
    error68 = (error68up - error68down) / 2.0
    epsilon = abs(1 - errorstd / error68)
    ax.plot(xgrid, epsilon, linestyle="-", color=color)
    labels.append(rf"{pdf.label}")

    return [5 * epsilon]
