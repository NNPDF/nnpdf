"""
deltachi2.py

Plots and data processing that can be used in a delta chi2 analysis
"""
import logging
import warnings
from collections import namedtuple

import matplotlib.pyplot as plt
import numpy as np

from reportengine.checks import CheckError, make_argcheck
from reportengine.figure import figure, figuregen
from reportengine import collect

from validphys import plotutils
from validphys.checks import check_scale, check_pdf_normalize_to, check_pdfs_noband
from validphys.core import PDF
from validphys.pdfplots import PDFPlotter, BandPDFPlotter, FlavourState


log = logging.getLogger(__name__)


@make_argcheck
def check_pdf_is_symmhessian(pdf, **kwargs):
    """Check ``pdf`` has error type of ``symmhessian``"""
    etype = pdf.ErrorType
    if etype != "symmhessian":
        raise CheckError(
            "Error: type of PDF %s must be 'symmhessian' and not %s" % (pdf, etype)
        )


@check_pdf_is_symmhessian
def delta_chi2_hessian(pdf, total_chi2_data):
    """
    Return delta_chi2 (computed as in plot_delta_chi2_hessian) relative to
    each eigenvector of the Hessian set.
    """
    delta_chi2 = (
        np.ravel(total_chi2_data.replica_result.error_members())
        - total_chi2_data.central_result
    )
    return delta_chi2


@figure
def plot_delta_chi2_hessian_eigenv(delta_chi2_hessian, pdf):
    """
    Plot of the chi2 difference between chi2 of each eigenvector of a symmHessian set
    and the central value for all experiments in a fit.
    As a function of every eigenvector in a first plot, and as a distribution in a second plot.
    """
    delta_chi2 = delta_chi2_hessian

    x = np.arange(1, len(delta_chi2) + 1)

    fig, ax = plt.subplots()

    ax.bar(x, delta_chi2, label=pdf.label)
    ax.set_xlabel("# Hessian PDF")
    ax.set_ylabel("$\Delta\chi^2$")
    ax.set_title("$\Delta\chi^2$ each eigenvector")

    ax.legend()

    return fig


@figure
def plot_delta_chi2_hessian_distribution(delta_chi2_hessian, pdf, total_chi2_data):
    """
    Plot of the chi2 difference between chi2 of each eigenvector of a symmHessian set
    and the central value for all experiments in a fit.
    As a function of every eigenvector in a first plot, and as a distribution in a second plot.
    """
    delta_chi2 = delta_chi2_hessian

    fig, ax = plt.subplots()

    bins = np.arange(np.floor(min(delta_chi2)), np.ceil(max(delta_chi2))+1)

    ax.hist(
        delta_chi2,
        bins=bins,
        label=f"{pdf.label} - $\chi^2_{0}$={total_chi2_data.central_result:.0f}",
    )
    ax.set_xlabel("$\Delta\chi^2$")
    ax.set_title("$\Delta\chi^2$ distribution")

    return fig


XPlottingGrid = namedtuple(
    "XPlottingGrid", ("Q", "basis", "flavours", "xgrid", "grid_values", "scale")
)


def pos_neg_xplotting_grids(delta_chi2_hessian, xplotting_grid):
    """
    Generates xplotting_grids correspodning to positive and negative delta chi2s.
    """
    positive_eigenvalue_mask = delta_chi2_hessian >= 0

    # include replica 0 in both new grids
    pos_mask = np.append(True, positive_eigenvalue_mask)
    neg_mask = np.append(True, ~positive_eigenvalue_mask)

    pos_grid = xplotting_grid.grid_values[pos_mask]
    neg_grid = xplotting_grid.grid_values[neg_mask]

    pos_xplotting_grid = xplotting_grid._replace(grid_values=pos_grid)
    neg_xplotting_grid = xplotting_grid._replace(grid_values=neg_grid)

    return [xplotting_grid, pos_xplotting_grid, neg_xplotting_grid]


@figuregen
@check_pdf_normalize_to
@check_pdfs_noband
@check_scale("xscale", allow_none=True)
def plot_pos_neg_pdfs(
    pdf,
    pos_neg_xplotting_grids,
    xscale: (str, type(None)) = None,
    normalize_to: (int, str, type(None)) = None,
    ymin=None,
    ymax=None,
    pdfs_noband: (list, type(None)) = None,
):
    """
    Plot the the uncertainty of the original hessian pdfs, as well as that of the positive and
    negative subset.
    """

    original_pdf = pdf.name
    # create fake PDF objects so we can reuse BandPDFPlotter
    pos_pdf = PDF(original_pdf)
    pos_pdf.label = f"{original_pdf}_pos"
    neg_pdf = PDF(original_pdf)
    neg_pdf.label = f"{original_pdf}_neg"

    pdfs = [pdf, pos_pdf, neg_pdf]

    yield from BandPDFPlotter(
        pdfs,
        pos_neg_xplotting_grids,
        xscale,
        normalize_to,
        ymin,
        ymax,
        pdfs_noband=pdfs_noband,
    )


class PDFEpsilonPlotter(PDFPlotter):
    """Subclassing PDFPlotter in order to plot epsilon (measure of gaussanity)
    for multiple PDFs, yielding a separate figure for each flavour
    """

    def setup_flavour(self, flstate):
        flstate.labels = []
        flstate.handles = []
    
    def get_ylabel(self, parton_name):
        return '$\epsilon(x)$'

    def draw(self, pdf, grid, flstate):
        """Obtains the gridvalues of epsilon (measure of Gaussianity)"""
        ax = flstate.ax
        flindex = flstate.flindex
        labels = flstate.labels
        handles = flstate.handles

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

        handle, = ax.plot(xgrid, epsilon, linestyle="-", color=color)

        handles.append(handle)
        labels.append(pdf.label)

        return [5 * epsilon]

    def legend(self, flstate):
        return flstate.ax.legend(flstate.handles, flstate.labels,
                                 handler_map={plotutils.HandlerSpec:
                                             plotutils.ComposedHandler()
                                             }
                                 )


@make_argcheck
def check_pdfs_are_montecarlo(pdfs, **kwargs):
    """Checks that the action is applied only to a pdf consisiting of MC replicas."""
    for pdf in pdfs:
        etype = pdf.ErrorType
        if etype != "replicas":
            raise CheckError(
                "Error: type of PDF %s must be 'replicas' and not '%s'" % (pdf, etype)
            )


@figuregen
@check_pdfs_are_montecarlo
@check_scale("xscale", allow_none=True)
def plot_epsilon(
    pdfs,
    xplotting_grids,
    xscale: (str, type(None)) = None,
    ymin=None,
    ymax=None,
    eps=None,
):
    """Plot the discrepancy (epsilon) of the 1-sigma and 68% bands at each grid value
    for all pdfs for a given Q. See https://arxiv.org/abs/1505.06736 eq. (11)

    xscale is read from pdf plotting_grid scale, which is 'log' by default.

    eps defines the value at which plot a simple hline
    """
    yield from PDFEpsilonPlotter(
        pdfs, xplotting_grids, xscale, normalize_to=None, ymin=ymin, ymax=ymax
    )
