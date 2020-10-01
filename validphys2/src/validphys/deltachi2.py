"""
deltachi2.py

Plots and data processing that can be used in a delta chi2 analysis
"""
import logging
import warnings
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np

from reportengine.checks import CheckError, make_argcheck
from reportengine.figure import figure, figuregen

from validphys import plotutils
from validphys.checks import check_scale


log = logging.getLogger(__name__)


class FlavourState(SimpleNamespace):
    """This is the namespace for the pats specific for each flavour"""
    pass


@make_argcheck
def check_pdf_is_replicas(pdfs, **kwargs):
    """Checks that the action is applied only to a pdf consisiting of MC replicas.
    """
    for pdf in pdfs:
        etype = pdf.ErrorType
        if etype != "replicas":
            raise CheckError("Error: type of PDF %s must be 'replicas' and not '%s'" % (pdf, etype))


@figuregen
@check_pdf_is_replicas
@check_scale("xscale", allow_none=True)
def plot_epsilon(
    pdfs, xplotting_grids, xscale: (str, type(None)) = None, ymin=None, ymax=None, eps=None
):
    """Plot the discrepancy (epsilon) of the 1-sigma and 68% bands at each grid value
    for all pdfs for a given Q. See https://arxiv.org/abs/1505.06736 eq. (11)

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

        ax.set_ylabel("$\epsilon$")

        ax.set_axisbelow(True)

        ax.legend(flstate.labels, handler_map={plotutils.HandlerSpec: plotutils.ComposedHandler()})

        plt.grid(linestyle="-", lw=0.5, alpha=0.4)

        yield fig, parton_name  # figure and string used by figuregen to save figures


def _setup_flavour_label(flstate):
    flstate.labels = []


def _draw(pdf, grid, flstate):
    """ Obtains the gridvalues of epsilon (measure of Gaussianity)
    """
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


@make_argcheck
def check_pdf_is_symmhessian(pdf, **kwargs):
    """Check plot_delta_chi2_hessian is applied only to Hessian
    set (symmhessian) converted from MC set.
    """
    etype = pdf.ErrorType
    if etype != "symmhessian":
        raise CheckError("Error: type of PDF %s must be 'symmhessian' and not %s" % (pdf, etype))


@check_pdf_is_symmhessian
def delta_chi2_hessian(pdf, experiments, groups_chi2):
    """
    Return delta_chi2 (computed as in plot_delta_chi2_hessian) relative to
    each eigenvector of the Hessian set.
    """
    experiments_chi2 = groups_chi2
    # store for each exp the chi2 from central value and chi2 from each error member
    total_chis_exps = np.zeros(
        (len(experiments), 1 + len(experiments_chi2[0].replica_result.error_members()))
    )

    for i, ch in enumerate(experiments_chi2):
        th, central, _ = ch  # chi2: error member | central value | n data
        total_chis_exps[i] = [central, *th.error_members()]

    # sum over all exps to get chi2 total for cv and each error member
    total_chis = np.sum(total_chis_exps, axis=0)

    delta_chi2 = total_chis[1:] - total_chis[0]

    return delta_chi2


@figure
@check_pdf_is_symmhessian
def plot_delta_chi2_hessian_eigenv(pdf, experiments, groups_chi2):
    """
    Plot of the chi2 difference between chi2 of each eigenvector of a symmHessian set
    and the central value for all experiments in a fit.
    As a function of every eigenvector in a first plot, and as a distribution in a second plot.
    """

    experiments_chi2 = groups_chi2

    # store for each exp the chi2 from central value and chi2 from each error memeber
    total_chis_exps = np.zeros(
        (len(experiments), 1 + len(experiments_chi2[0].replica_result.error_members()))
    )

    for i, ch in enumerate(experiments_chi2):
        th, central, _ = ch  # chi2: error member | central value | n data
        total_chis_exps[i] = [central, *th.error_members()]  # cv + chi2 each eigenvector

    # sum over all exps to get chi2 total for cv and each error member
    total_chis = np.sum(total_chis_exps, axis=0)

    delta_chi2 = total_chis[1:] - total_chis[0]

    x = np.arange(1, len(delta_chi2) + 1)

    fig, ax = plt.subplots()

    ax.bar(x, delta_chi2, label='%s' % pdf.label)
    ax.set_xlabel('# Hessian PDF')
    ax.set_ylabel('$\Delta\chi^2$')
    ax.set_title('$\Delta\chi^2$ each eigenvector')
    ax.grid(False)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center")
    plt.subplots_adjust(top=0.85)

    return fig


@figure
@check_pdf_is_symmhessian
def plot_delta_chi2_hessian_distribution(pdf, experiments, groups_chi2):
    """
    Plot of the chi2 difference between chi2 of each eigenvector of a symmHessian set
    and the central value for all experiments in a fit.
    As a function of every eigenvector in a first plot, and as a distribution in a second plot.
    """

    experiments_chi2 = groups_chi2

    # store for each exp the chi2 from central value and chi2 from each error member
    total_chis_exps = np.zeros(
        (len(experiments), 1 + len(experiments_chi2[0].replica_result.error_members()))
    )

    for i, ch in enumerate(experiments_chi2):
        th, central, _ = ch  # chi2: error member | central value | n data
        total_chis_exps[i] = [central, *th.error_members()]  # cv + chi2 each eigenvector

    # sum over all exps to get chi2 total for cv and each error member
    total_chis = np.sum(total_chis_exps, axis=0)

    delta_chi2 = total_chis[1:] - total_chis[0]

    x = np.arange(1, len(delta_chi2) + 1)

    fig, ax = plt.subplots()

    range_min = np.int32(np.floor(np.min(delta_chi2)))
    range_max = np.int32(np.ceil(np.max(delta_chi2)))
    bins = np.asarray([i for i in range(range_min - 1, range_max + 1)])
    ax.hist(delta_chi2, bins=bins, label="%s - $\chi^2_{0}$=%d" % (pdf.label, total_chis[0]))
    ax.set_xlabel("$\Delta\chi^2$")
    ax.set_title("$\Delta\chi^2$ distribution")
    ax.grid(False)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center")
    plt.subplots_adjust(top=0.85)

    return fig
