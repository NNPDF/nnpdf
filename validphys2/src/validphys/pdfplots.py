"""
pdfplots.py

Plots of quantities that are mostly functions of the PDFs only.
"""
import abc
import logging
import functools
import warnings
import numbers
import copy
from types import SimpleNamespace


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors

from reportengine.figure import figure, figuregen
from reportengine.checks import make_argcheck
from reportengine.floatformatting import format_number

from validphys import plotutils
from validphys.core import MCStats
from validphys.gridvalues import LUMI_CHANNELS
from validphys.utils import scale_from_grid
from validphys.checks import check_pdf_normalize_to, check_scale, check_have_two_pdfs
from validphys.checks import check_pdfs_noband, check_pdfs_plot_replicas

log = logging.getLogger(__name__)

class FlavourState(SimpleNamespace):
    """This is the namespace for the pats specific for each flavour"""
    pass


class PDFPlotter(metaclass=abc.ABCMeta):
    """Stateful object breaks plotting grids by favour, as a function of x and
    for fixed Q.

    This class has a lot of state, but it should all be defined at
    initialization time. Things that change e.g. per flavour should be passed
    explicitly as arguments.
    """

    def __init__(self, pdfs, xplotting_grids, xscale, normalize_to, ymin, ymax):
        self.pdfs = pdfs
        self._xplotting_grids = xplotting_grids
        self._xscale = xscale
        self.normalize_to = normalize_to
        self.xplotting_grids = self.normalize()
        self.ymin = ymin
        self.ymax = ymax

    def setup_flavour(self, flstate):
        pass

    def normalize(self):
        normalize_to = self.normalize_to
        if normalize_to is not None:
            normalize_pdf = self.normalize_pdf
            normalize_grid = self._xplotting_grids[normalize_to]
            normvals = normalize_grid.grid_values.central_value()

            #Handle division by zero more quietly
            def fp_error(tp, flag):
                log.warning("Invalid values found computing "
                    f"normalization to {normalize_pdf}: "
                    f"Floating point error ({tp}).")
                #Show warning only once
                np.seterr(all='ignore')

            newgrids = []
            with np.errstate(all='call'):
                np.seterrcall(fp_error)
                for pdf, grid in zip(self.pdfs, self._xplotting_grids):
                    newvalues = pdf.stats_class(grid.grid_values.data/normvals)
                    newgrids.append(grid.copy_grid(grid_values=newvalues))

            return newgrids
        return self._xplotting_grids

    @property
    def normalize_pdf(self):
        if self.normalize_to is None:
            raise AttributeError("Need to set a normalize_to index.")
        return self.pdfs[self.normalize_to]

    def get_ylabel(self, parton_name):
        if self.normalize_to is not None:
            return f"Ratio to {self.normalize_pdf.label}"

        base_str = f"x{parton_name}(x)"
        # Ask the xplotting grid if it has something to add:
        final_str = self.firstgrid.process_label(base_str)

        # Wrap it in latex
        return f"${final_str}$"

    def get_title(self, parton_name):
        return f"${parton_name}$ at {self.Q} GeV"

    @property
    def xscale(self):
        if self._xscale is None:
            return scale_from_grid(self.firstgrid)
        return self._xscale

    @property
    def Q(self):
        return self.firstgrid.Q

    @property
    def firstgrid(self):
        if self.xplotting_grids:
            return self.xplotting_grids[0]
        raise AttributeError("Need at least one xgrid")


    @abc.abstractmethod
    def draw(self, pdf, grid, flstate):
        """Plot the desired function of the grid and return the array to be
        used for autoscaling"""
        pass

    def legend(self, flstate):
        return flstate.ax.legend()

    def __iter__(self):
        yield from self()


    def __call__(self,):
        if not self.xplotting_grids:
            return

        basis = self.firstgrid.basis

        for flindex, fl in enumerate(self.firstgrid.flavours):
            fig, ax = plt.subplots()
            parton_name = basis.elementlabel(fl)
            flstate = FlavourState(flindex=flindex, fl=fl, fig=fig, ax=ax,
                                    parton_name=parton_name)
            self.setup_flavour(flstate)
            ax.set_title(self.get_title(parton_name))

            all_vals = []
            for pdf, grid in zip(self.pdfs, self.xplotting_grids):
                limits = self.draw(pdf, grid, flstate)
                if limits is not None:
                    all_vals.append(np.atleast_2d(limits))

            #Note these two lines do not conmute!
            ax.set_xscale(self.xscale)
            plotutils.frame_center(ax, self.firstgrid.xgrid, np.concatenate(all_vals))
            if (self.ymin is not None):
                ax.set_ylim(ymin=self.ymin)
            if (self.ymax is not None):
                ax.set_ylim(ymax=self.ymax)

            ax.set_xlabel('$x$')
            ax.set_xlim(self.firstgrid.xgrid[0])


            ax.set_ylabel(self.get_ylabel(parton_name))

            ax.set_axisbelow(True)

            self.legend(flstate)
            yield fig, parton_name



@functools.lru_cache()
def _warn_pdf_not_montecarlo(pdf):
    et = pdf.error_type
    if et != 'replicas':
        log.warning("Plotting members of a non-Monte Carlo PDF set:"
        f" {pdf.name} with error type '{et}'.")

#Cant't add the lru_cache here because pdfs is not hashable at the moment
@make_argcheck
def _warn_any_pdf_not_montecarlo(pdfs):
    for pdf in pdfs:
        _warn_pdf_not_montecarlo(pdf)


class ReplicaPDFPlotter(PDFPlotter):
    def draw(self, pdf, grid, flstate):
        ax = flstate.ax
        next_prop = next(ax._get_lines.prop_cycler)
        color = next_prop['color']
        flavour_grid = grid.select_flavour(flstate.flindex)
        stats = flavour_grid.grid_values
        gv = stats.data
        ax.plot(grid.xgrid, gv.T, alpha=0.2, linewidth=0.5,
                color=color, zorder=1)
        ax.plot(grid.xgrid, stats.central_value(), color=color,
                linewidth=2,
                label=pdf.label)
        return gv

@figuregen
@check_pdf_normalize_to
@check_scale('xscale', allow_none=True)
@_warn_any_pdf_not_montecarlo
def plot_pdfreplicas(pdfs, xplotting_grids, xscale:(str,type(None))=None,
                      normalize_to:(int,str,type(None))=None, ymin = None, ymax = None):
    """Plot the replicas of the specified PDFs. Otherise it works the same as
    plot_pdfs.

    - xscale sets the scale of the plot. E.g. 'linear' or 'log'. Default is
    deduced from the xplotting_grid, which in turn is 'log' by default.

    - normalize_to should be, a pdf id or an index of the pdf (starting from one).
    """
    yield from ReplicaPDFPlotter(pdfs=pdfs, xplotting_grids=xplotting_grids,
                                 xscale=xscale, normalize_to=normalize_to, ymin=ymin, ymax=ymax)

class UncertaintyPDFPlotter(PDFPlotter):

    def get_ylabel(self, parton_name):
        if self.normalize_to is not None:
            return r"$\sigma($%s$)$" % super().get_ylabel(parton_name)
        return r"$\sigma$"

    def draw(self, pdf, grid, flstate):
        ax = flstate.ax
        flindex = flstate.flindex
        stats = grid.select_flavour(flindex).grid_values

        res = stats.std_error()

        ax.plot(grid.xgrid, res, label=pdf.label)

        return res

    def __call__(self):
        # Fixup y limit to not be below zero
        for fig, parton_name in super().__call__():
            ax = fig.get_axes()[0]
            ymin, _ = ax.get_ylim()
            ax.set_ylim(max(0, ymin), None)
            yield fig, parton_name


@figuregen
@check_pdf_normalize_to
@check_scale('xscale', allow_none=True)
def plot_pdf_uncertainties(pdfs, xplotting_grids, xscale:(str,type(None))=None,
                      normalize_to:(int,str,type(None))=None, ymin=None, ymax=None):
    """Plot the PDF standard deviations as a function of x.
    If normalize_to is set, the ratio to that
    PDF's central value is plotted. Otherwise it is the absolute values."""
    yield from UncertaintyPDFPlotter(pdfs, xplotting_grids, xscale, normalize_to, ymin, ymax)


class AllFlavoursPlotter(PDFPlotter):
    """Auxiliary class which groups multiple PDF flavours in one plot."""

    def setup_flavour(self, flstate):
        flstate.handles= self.handles
        flstate.labels= self.doesnothing
        flstate.hatchit= self.hatchit

    def __call__(self):
        if not self.xplotting_grids:
            return

        self.handles = []
        self.doesnothing = []
        self.labels = []
        self.hatchit = plotutils.hatch_iter()

        basis = self.firstgrid.basis
        fig, ax = plt.subplots()
        ax.set_xlabel('$x$')
        ax.set_ylabel(self.get_ylabel(None))
        ax.set_xscale(self.xscale)
        ax.set_title(self.get_title(None))

        all_vals = []
        for flindex, fl in enumerate(self.firstgrid.flavours):

            parton_name = basis.elementlabel(fl)
            self.labels.append(f'${parton_name}$')
            flstate = FlavourState(flindex=flindex, fl=fl, fig=fig, ax=ax,
                                   parton_name=parton_name)
            self.setup_flavour(flstate)

            for pdf, grid in zip(self.pdfs, self.xplotting_grids):
                limits = self.draw(pdf, grid, flstate)
                if limits is not None:
                    all_vals.append(np.atleast_2d(limits))

        #It can happen that we don't get anything to concatenate
        #e.g. because we are comparing to the base PDF several times.
        if all_vals:
            plotutils.frame_center(ax, self.firstgrid.xgrid,
                                   np.concatenate(all_vals))
        if (self.ymin is not None):
            ax.set_ylim(ymin=self.ymin)
        if (self.ymax is not None):
            ax.set_ylim(ymax=self.ymax)

        ax.set_axisbelow(True)
        ax.set_xlim(self.firstgrid.xgrid[0])
        flstate.labels = self.labels
        self.legend(flstate)
        return fig

    def get_ylabel(self, parton_name):
        return ''


class DistancePDFPlotter(PDFPlotter):
    """Auxiliary class which draws the distance plots."""

    def normalize(self):
        return self._xplotting_grids

    def get_ylabel(self, parton_name):
        return "Distance from {}".format(self.normalize_pdf.label)

    def draw(self, pdf, grid, flstate):

        if pdf == self.normalize_pdf:
            return None

        ax = flstate.ax
        flindex = flstate.flindex
        pcycler = ax._get_lines.prop_cycler
        next_prop = next(pcycler)
        color = next_prop['color']

        # The grid for the distance is (1, flavours, points)
        # take only the flavour we are interested in
        gv = grid.select_flavour(flindex).grid_values.data.squeeze()

        ax.plot(grid.xgrid, gv, color=color, label='$%s$' % flstate.parton_name)

        return gv

    def get_title(self, parton_name):
        return f'{self.pdfs[(1+self.normalize_to)%2]} Q={self.Q : .1f} GeV'


class VarDistancePDFPlotter(DistancePDFPlotter):
    """Auxiliary class which draws the variance distance plots"""

    def get_ylabel(self, parton_name):
        return "Variance distance from {}".format(self.normalize_pdf.label)

    def get_title(self, parton_name):
        return f'{self.pdfs[(1+self.normalize_to)%2]} Q={self.Q : .1f} GeV'


class FlavoursDistancePlotter(DistancePDFPlotter, AllFlavoursPlotter): pass


class FlavoursVarDistancePlotter(VarDistancePDFPlotter, AllFlavoursPlotter): pass


@figure
@check_pdf_normalize_to
@check_have_two_pdfs
@check_scale('xscale', allow_none=True)
def plot_pdfdistances(pdfs, distance_grids, *,
                      xscale:(str,type(None))=None,
                      normalize_to:(int,str),ymin=None,ymax=None):
    """Plots the distances between different PDF sets and a reference PDF set
    for all flavours. Distances are normalized such that a value of order 10
    is unlikely to be explained by purely statistical fluctuations
    """
    return FlavoursDistancePlotter(pdfs, distance_grids, xscale, normalize_to, ymin, ymax)()


@figure
@check_pdf_normalize_to
@check_have_two_pdfs
@check_scale('xscale', allow_none=True)
def plot_pdfvardistances(pdfs, variance_distance_grids, *,
                      xscale:(str,type(None))=None,
                      normalize_to:(int,str),ymin=None,ymax=None):
    """Plots the distances between different PDF sets and a reference PDF set
    for all flavours. Distances are normalized such that a value of order 10
    is unlikely to be explained by purely statistical fluctuations
    """
    return FlavoursVarDistancePlotter(pdfs, variance_distance_grids, xscale, normalize_to, ymin, ymax)()


class BandPDFPlotter(PDFPlotter):
    def __init__(
        self,
        *args,
        pdfs_noband=None,
        show_mc_errors=True,
        legend_stat_labels=True,
        **kwargs
    ):
        if pdfs_noband is None:
            pdfs_noband = []
        self.pdfs_noband = pdfs_noband
        self.show_mc_errors = show_mc_errors
        self.legend_stat_labels = legend_stat_labels
        super().__init__(*args, **kwargs)

    def setup_flavour(self, flstate):
        flstate.handles=[]
        flstate.labels=[]
        flstate.hatchit=plotutils.hatch_iter()

    def draw(self, pdf, grid, flstate):
        ax = flstate.ax
        hatchit = flstate.hatchit
        labels = flstate.labels
        handles = flstate.handles
        # Take only the flavours we are interested in
        stats = grid.select_flavour(flstate.flindex).grid_values

        next_prop = next(ax._get_lines.prop_cycler)
        cv = stats.central_value()
        xgrid = grid.xgrid
        #Ignore spurious normalization warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            err68down, err68up = stats.errorbar68()

        #http://stackoverflow.com/questions/5195466/matplotlib-does-not-display-hatching-when-rendering-to-pdf
        hatch = next(hatchit)
        color = next_prop['color']
        cvline, = ax.plot(xgrid, cv, color=color)
        if pdf in self.pdfs_noband:
            labels.append(pdf.label)
            handles.append(cvline)
            return [cv, cv]
        alpha = 0.5
        ax.fill_between(xgrid, err68up, err68down, color=color, alpha=alpha,
                        zorder=1)

        ax.fill_between(xgrid, err68up, err68down, facecolor='None', alpha=alpha,
                        edgecolor=color,
                        hatch=hatch,
                        zorder=1)
        if isinstance(stats, MCStats) and self.show_mc_errors:
            errorstdup, errorstddown = stats.errorbarstd()
            ax.plot(xgrid, errorstdup, linestyle='--', color=color)
            ax.plot(xgrid, errorstddown, linestyle='--', color=color)
            label = (
                rf"{pdf.label} ($68\%$ c.l.+$1\sigma$)"
                if self.legend_stat_labels
                else pdf.label
            )
            outer = True
        else:
            outer = False
            label = (
                rf"{pdf.label} ($68\%$ c.l.)"
                if self.legend_stat_labels
                else pdf.label
            )
        handle = plotutils.HandlerSpec(color=color, alpha=alpha,
                                               hatch=hatch,
                                               outer=outer)
        handles.append(handle)
        labels.append(label)

        return [err68down, err68up]

    def legend(self, flstate):
        return flstate.ax.legend(flstate.handles, flstate.labels,
                                 handler_map={plotutils.HandlerSpec:
                                             plotutils.ComposedHandler()
                                             }
                                 )


@figuregen
@check_pdf_normalize_to
@check_pdfs_noband
@check_scale("xscale", allow_none=True)
def plot_pdfs(
    pdfs,
    xplotting_grids,
    xscale: (str, type(None)) = None,
    normalize_to: (int, str, type(None)) = None,
    ymin=None,
    ymax=None,
    pdfs_noband: (list, type(None)) = None,
    show_mc_errors: bool = True,
    legend_stat_labels: bool = True,
):
    """Plot the central value and the uncertainty of a list of pdfs as a
    function of x for a given value of Q. If normalize_to is given, plot the
    ratios to the corresponding PDF. Otherwise, plot absolute values.
    See the help for ``xplotting_grid`` for information on how to set basis,
    flavours and x ranges. Yields one figure per PDF flavour.

    normalize_to:  Either the name of one of the PDFs or its corresponding
    index in the list, starting from one, or None to plot absolute values.

    xscale: One of the matplotlib allowed scales. If undefined, it will be
    set based on the scale in xgrid, which should be used instead.

    pdfs_noband: A list of PDFs to plot without error bands, i.e. only the
    central values of these PDFs will be plotted. The list can be formed of
    strings, corresponding to PDF IDs, integers (starting from one),
    corresponding to the index of the PDF in the list of PDFs, or a mixture
    of both.

    show_mc_errors (bool): Plot 1σ bands in addition to 68% errors for Monte Carlo
    PDF.

    legend_stat_labels (bool): Show detailed information on what kind of confidence
    interval is being plotted in the legend labels.
    """
    yield from BandPDFPlotter(
        pdfs,
        xplotting_grids,
        xscale,
        normalize_to,
        ymin,
        ymax,
        pdfs_noband=pdfs_noband,
        show_mc_errors=show_mc_errors,
        legend_stat_labels=legend_stat_labels,
    )


@figuregen
@check_pdf_normalize_to
@check_pdfs_plot_replicas
@check_pdfs_noband
@check_scale("xscale", allow_none=True)
def plot_pdfs_mixed(
    pdfs,
    xplotting_grids,
    xscale: (str, type(None)) = None,
    normalize_to: (int, str, type(None)) = None,
    ymin=None,
    ymax=None,
    pdfs_noband: (list, type(None)) = None,
    show_mc_errors: bool = True,
    legend_stat_labels: bool = True,
    pdfs_plot_replicas: (list, type(None)) = None,
):
    """This function is similar to plot_pdfs, except instead of only plotting
    the central value and the uncertainty of the PDFs, those PDFs indicated by
    pdfs_plot_replicas will be plotted as replicas without the central value.

    Inputs are the same as plot_pdfs, with the exeption of pdfs_plot_replicas,
    which only exists here.

    pdfs_plot_replicas: A list of PDFs to plot as replicas, i.e. the
    central values and replicas of these PDFs will be plotted. The list can be
    formed of strings, corresponding to PDF IDs, integers (starting from one),
    corresponding to the index of the PDF in the list of PDFs, or a mixture
    of both.
    """
    yield from MixMatchPDFPlotter(
        pdfs,
        xplotting_grids,
        xscale,
        normalize_to,
        ymin,
        ymax,
        pdfs_noband=pdfs_noband,
        show_mc_errors=show_mc_errors,
        legend_stat_labels=legend_stat_labels,
        pdfs_plot_replicas=pdfs_plot_replicas,
    )


@figuregen
@check_pdf_normalize_to
@check_pdfs_noband
@check_scale("xscale", allow_none=True)
def plot_pdfs_kinetic_energy(
    pdfs,
    kinetic_xplotting_grids,
    xscale: (str, type(None)) = None,
    normalize_to: (int, str, type(None)) = None,
    ymin=None,
    ymax=None,
    pdfs_noband: (list, type(None)) = None,
    show_mc_errors: bool = True,
    legend_stat_labels: bool = True,
):
    """Band plotting of the "kinetic energy" of the PDF as a function of x for a given value of Q.
    The input of this function is similar to those of ``plot_pdfs``.
    """
    return plot_pdfs(
        pdfs,
        kinetic_xplotting_grids,
        xscale=xscale,
        normalize_to=normalize_to,
        ymin=ymin,
        ymax=ymax,
        pdfs_noband=pdfs_noband,
        show_mc_errors=show_mc_errors,
        legend_stat_labels=legend_stat_labels,
    )


@figuregen
@check_pdf_normalize_to
@check_scale('xscale', allow_none=True)
@_warn_any_pdf_not_montecarlo
def plot_pdfreplicas_kinetic_energy(
    pdfs,
    kinetic_xplotting_grids,
    xscale: (str, type(None)) = None,
    normalize_to: (int, str, type(None)) = None,
    ymin=None,
    ymax=None,
):
    """Plot the kinetic energy of the replicas of the specified PDFs.
    Otherise it works the same as ``plot_pdfs_kinetic_energy``.
    """
    return plot_pdfreplicas(
        pdfs,
        kinetic_xplotting_grids,
        xscale=xscale,
        normalize_to=normalize_to,
        ymin=ymin,
        ymax=ymax,
    )


@figuregen
@check_pdf_normalize_to
@check_pdfs_noband
@check_pdfs_plot_replicas
@check_scale("xscale", allow_none=True)
def plot_pdfs_mixed_kinetic_energy(
    pdfs,
    kinetic_xplotting_grids,
    xscale: (str, type(None)) = None,
    normalize_to: (int, str, type(None)) = None,
    ymin=None,
    ymax=None,
    pdfs_noband: (list, type(None)) = None,
    show_mc_errors: bool = True,
    legend_stat_labels: bool = True,
    pdfs_plot_replicas: (list, type(None)) = None,
):
    """Mixed band and replica plotting of the "kinetic energy" of the PDF as a
    function of x for a given value of Q.
    """
    return plot_pdfs_mixed(
        pdfs,
        kinetic_xplotting_grids,
        xscale=xscale,
        normalize_to=normalize_to,
        ymin=ymin,
        ymax=ymax,
        pdfs_noband=pdfs_noband,
        show_mc_errors=show_mc_errors,
        legend_stat_labels=legend_stat_labels,
        pdfs_plot_replicas=pdfs_plot_replicas,
    )


class FlavoursPlotter(AllFlavoursPlotter, BandPDFPlotter):
    def get_title(self, parton_name):
        return f'{self.pdfs[0]} Q={self.Q : .1f} GeV'

@figure
@check_scale('xscale', allow_none=True)
def plot_flavours(pdf, xplotting_grid, xscale:(str,type(None))=None,
                      normalize_to:(int,str,type(None))=None,ymin=None,ymax=None):
    """Plot the absolute central value and the uncertainty of all the flavours
    of a pdf as a function of x for a given value of Q.

    xscale: One of the matplotlib allowed scales. If undefined, it will be
    set based on the scale in xgrid, which should be used instead.

    """
    return FlavoursPlotter([pdf], [xplotting_grid], xscale, normalize_to=None, ymin= ymin, ymax=ymax)()

@figure
@check_pdf_normalize_to
@check_pdfs_noband
def plot_lumi1d(
    pdfs,
    pdfs_lumis,
    lumi_channel,
    sqrts: numbers.Real,
    y_cut: (numbers.Real, type(None)) = None,
    normalize_to=None,
    show_mc_errors: bool = True,
    ymin: (numbers.Real, type(None)) = None,
    ymax: (numbers.Real, type(None)) = None,
    pdfs_noband=None,
    scale="log",
    legend_stat_labels: bool=True,
):
    """Plot PDF luminosities at a given center of mass energy.
    sqrts is the center of mass energy (GeV).

    This action plots the luminosity (as computed by `lumigrid1d`) as a
    function of invariant mass for all PDFs for a single lumi channel.
    ``normalize_to`` works as for `plot_pdfs` and allows to plot a ratio to the
    central value of some of the PDFs. `ymin` and `ymax` can be used to set
    exact bounds for the scale. `y_cut` can be used to specify a rapidity cut
    over the integration range. `show_mc_errors` controls whether the 1σ error
    bands are shown in addition to the 68% confidence intervals for Monte Carlo
    PDFs. A list `pdfs_noband` can be passed to supress the error bands for
    certain PDFs and plot the central values only. `legend_stat_labels` controls
    whether to show detailed information on what kind of confidence interval
    is being plotted in the legend labels.
    """

    fig, ax = plt.subplots()
    if normalize_to is not None:
        norm = pdfs_lumis[normalize_to].grid_values.central_value()
        ylabel = f"Ratio to {pdfs[normalize_to]}"
    else:
        norm = 1
        ylabel = r"$\mathcal{L} (GeV^{-2})$"

    hatchit = plotutils.hatch_iter()
    pcycler = plotutils.color_iter()
    handles = []
    labels = []

    for pdf, lumigrid1d in zip(pdfs, pdfs_lumis):
        mx = lumigrid1d.m
        gv = lumigrid1d.grid_values

        cv = gv.central_value()

        err68down, err68up = gv.errorbar68()
        errstddown, errstdup = gv.errorbarstd()

        color = next(pcycler)
        hatch = next(hatchit)

        alpha = 0.5

        (central_line,) = ax.plot(mx, cv / norm, color=color)

        if pdfs_noband and pdf in pdfs_noband:
            handles.append(central_line)
            labels.append(pdf.label)
            continue

        ax.fill_between(
            mx, err68down / norm, err68up / norm, color=color, alpha=alpha, zorder=1
        )
        ax.fill_between(
            mx,
            err68down / norm,
            err68up / norm,
            facecolor="None",
            alpha=alpha,
            edgecolor=color,
            hatch=hatch,
            zorder=1,
        )

        if isinstance(gv, MCStats) and show_mc_errors:
            ax.plot(mx, errstddown / norm, linestyle="--", color=color)
            ax.plot(mx, errstdup / norm, linestyle="--", color=color)
            label_add = r"($68%$ c.l.+$1\sigma$)" if legend_stat_labels else ""
            outer = True
        else:
            label_add = r"($68\%$ c.l.)" if legend_stat_labels else ""
            outer = False

        handle = plotutils.HandlerSpec(
            color=color, alpha=alpha, hatch=hatch, outer=outer
        )
        handles.append(handle)
        labels.append(f"{pdf.label} {label_add}")

    ax.legend(
        handles,
        labels,
        handler_map={plotutils.HandlerSpec: plotutils.ComposedHandler()},
    )

    ax.set_ylabel(ylabel)
    ax.set_xlabel("$m_{X}$ (GeV)")
    ax.set_xlim(mx[0], mx[-1])
    ax.set_ylim(ymin, ymax)
    ax.set_xscale(scale)
    ax.grid(False)
    if y_cut==None:
        ax.set_title(
            f"${LUMI_CHANNELS[lumi_channel]}$ luminosity\n"
            f"$\\sqrt{{s}}={format_number(sqrts/1000)}$ TeV"
        )
    else:
        ax.set_title(
            f"${LUMI_CHANNELS[lumi_channel]}$ luminosity\n"
            f"$\\sqrt{{s}}={format_number(sqrts/1000)}$ TeV   "
            f"$\\|y|<{format_number(y_cut)}$"
        )

    return fig


@figure
@check_pdf_normalize_to
def plot_lumi1d_uncertainties(
    pdfs,
    pdfs_lumis,
    lumi_channel,
    sqrts: numbers.Real,
    y_cut: (numbers.Real, type(None)) = None,
    normalize_to=None,
    ymin: (numbers.Real, type(None)) = None,
    ymax: (numbers.Real, type(None)) = None,
    scale = "log",
):
    """Plot PDF luminosity uncertainties at a given center of mass energy.
    sqrts is the center of mass energy (GeV).

    If `normalize_to` is set, the values are normalized to the central value of
    the corresponding PDFs. `y_cut` can be used to specify a rapidity cut
    over the integration range.
    """

    fig, ax = plt.subplots()
    if normalize_to is not None:
        norm = pdfs_lumis[normalize_to].grid_values.central_value()
        ylabel = f"Ratio to {pdfs[normalize_to]}"
    else:
        norm = None
        ylabel = r"$\sigma\left(\mathcal{L} (GeV^{-2})\right)$"

    for pdf, lumigrid1d, color in zip(pdfs, pdfs_lumis, plotutils.color_iter()):
        mx = lumigrid1d.m
        gv = lumigrid1d.grid_values

        err = gv.std_error()

        if norm is not None:
            err /= norm
        ax.plot(mx, err, color=color, label=pdf.label)

    ax.legend()

    ax.set_ylabel(ylabel)
    ax.set_xlabel("$m_{X}$ (GeV)")
    ax.set_xlim(mx[0], mx[-1])
    ax.set_xscale(scale)
    ax.grid(False)
    if y_cut==None:
        ax.set_title(
            f"${LUMI_CHANNELS[lumi_channel]}$ luminosity uncertainty\n"
            f"$\\sqrt{{s}}={format_number(sqrts/1000)}$ TeV"
        )
    else:
        ax.set_title(
            f"${LUMI_CHANNELS[lumi_channel]}$ luminosity uncertainty\n"
            f"$\\sqrt{{s}}={format_number(sqrts/1000)}$ TeV   "
            f"$\\|y|<{format_number(y_cut)}$"
        )
    ax.set_ylim(ymin, ymax)
    current_ymin, _ = ax.get_ylim()
    ax.set_ylim(max(0, current_ymin), None)


    return fig



@figure
@check_pdf_normalize_to
def plot_lumi1d_replicas(
    pdfs,
    pdfs_lumis,
    lumi_channel,
    sqrts: numbers.Real,
    y_cut: (numbers.Real, type(None)) = None,
    normalize_to=None,
    ymin: (numbers.Real, type(None)) = None,
    ymax: (numbers.Real, type(None)) = None,
    scale="log",
):
    """This function is similar to `plot_lumi1d`, but instead of plotting
    the standard deviation and 68% c.i. it plots the luminosities for
    individual replicas.

    Plot PDF replica luminosities at a given center of mass energy.
    sqrts is the center of mass energy (GeV).

    This action plots the luminosity (as computed by `lumigrid1d`) as a
    function of invariant mass for all PDFs for a single lumi channel.
    ``normalize_to`` works as for `plot_pdfs` and allows to plot a ratio to the
    central value of some of the PDFs. `ymin` and `ymax` can be used to set
    exact bounds for the scale. `y_cut` can be used to specify a rapidity cut
    over the integration range. `show_mc_errors` controls whether the 1σ error
    bands are shown in addition to the 68% confidence intervals for Monte Carlo
    PDFs.
    """

    fig, ax = plt.subplots()
    if normalize_to is not None:
        norm = pdfs_lumis[normalize_to].grid_values.central_value()
        ylabel = f"Ratio to {pdfs[normalize_to]}"
    else:
        norm = 1
        ylabel = r"$\mathcal{L} (GeV^{-2})$"

    pcycler = plotutils.color_iter()

    lines = []
    labels = []
    for pdf, lumigrid1d in zip(pdfs, pdfs_lumis):
        mx = lumigrid1d.m
        gv = lumigrid1d.grid_values

        cv = gv.central_value()

        replicas = gv.data

        color = next(pcycler)

        ax.plot(mx, (replicas/norm).T, alpha=0.2, linewidth=0.5,
                color=color, zorder=1)
        line, = ax.plot(mx, cv/norm, color=color,
                linewidth=2)
        lines.append(line)
        labels.append(pdf.label)

    ax.set_ylabel(ylabel)
    ax.set_xlabel("$m_{X}$ (GeV)")
    ax.set_xlim(mx[0], mx[-1])
    ax.set_ylim(ymin, ymax)
    ax.set_xscale(scale)
    ax.legend(lines,labels)
    ax.grid(False)
    if y_cut==None:
        ax.set_title(
            f"${LUMI_CHANNELS[lumi_channel]}$ luminosity\n"
            f"$\\sqrt{{s}}={format_number(sqrts/1000)}$ TeV"
        )
    else:
        ax.set_title(
            f"${LUMI_CHANNELS[lumi_channel]}$ luminosity\n"
            f"$\\sqrt{{s}}={format_number(sqrts/1000)}$ TeV   "
            f"$\\|y|<{format_number(y_cut)}$"
        )

    return fig



#TODO: Move these to utils somewhere? Find better implementations?
def _reflect_matrl(mat, odd=False):
    """Reflect a matrix with positive values in the first axis to have the
    same balues for the nwgative axis. The first value is not reflected.

    If ``odd`` is set, the negative part will be multiplied by -1.

    """
    mat = np.asarray(mat)
    res = np.empty(shape=(mat.shape[0]*2-1, *mat.shape[1:]),dtype=mat.dtype)
    neglen = mat.shape[0]-1
    fact = -1 if odd else 1
    res[:neglen,...] = fact*mat[:0:-1,...]
    res[neglen:,...] = mat
    return res

def _reflect_matud(mat, odd=False):
    """Reflect a matrix with positive values in the second axis to have the
    same balues for the nwgative axis. The first value is not reflected.

    If ``odd`` is set, the negative part will be multiplied by -1.

    """
    mat = np.asarray(mat)
    res = np.empty(shape=(mat.shape[0], mat.shape[1]*2-1, *mat.shape[2:]),
        dtype=mat.dtype)
    neglen = mat.shape[1]-1
    fact = -1 if odd else 1
    res[:,:neglen,...] = fact*mat[:,:0:-1,...]
    res[:,neglen:,...] = mat
    return res


@figure
def plot_lumi2d(pdf, lumi_channel, lumigrid2d, sqrts,
                display_negative:bool=True):
    """Plot the absolute luminosity on a grid of invariant mass and
    rapidity for a given center of mass energy `sqrts`.
    The color scale is logarithmic.
    If `display_negative` is True, mark the negative values.

    The luminosity is calculated for positive rapidity, and reflected for
    negative rapidity for display purposes.

    """


    cmap = copy.copy(cm.viridis_r)
    cmap.set_bad("white", alpha=0)
    fig, ax = plt.subplots()
    gv = lumigrid2d.grid_values
    mat = gv.central_value()

    fig, ax = plt.subplots()

    mat = _reflect_matud(mat)
    y = _reflect_matrl(lumigrid2d.y, odd=True)
    masked_weights = np.ma.masked_invalid(mat, copy=False)

    #TODO: SymLogNorm is really the right thing to do here, but I can't be
    #bothered to make it work. Mostly the ticks around zero are completely
    #broken and looks like it takes a lot of fidlling wirh the mpl internals
    #to fix it.

    with np.errstate(invalid='ignore'):
        positive_mask = masked_weights>0
    linlim = np.nanpercentile(masked_weights[positive_mask],90)/1e5

    #norm = mcolors.SymLogNorm(linlim, vmin=None)

    norm = mcolors.LogNorm(vmin=linlim)
    with np.errstate(invalid='ignore'):
        masked_weights[masked_weights<linlim] = linlim

    mesh = ax.pcolormesh(y, lumigrid2d.m, masked_weights, cmap=cmap,
        shading='gouraud',
        linewidth=0,
        edgecolor='None',
        rasterized=True,
        norm=norm,
    )
    #Annoying code because mpl does the defaults horribly
    #loc = mticker.SymmetricalLogLocator(base=10, linthresh=linlim,)
    #loc.numticks = 5

    #fig.colorbar(mesh, ticks=loc)

    if display_negative:
        cmap_neg =  mcolors.ListedColormap(['red', (0,0,0,0)])
        neg_norm = mcolors.BoundaryNorm([0,0.5,1],2)
        ax.pcolormesh(y, lumigrid2d.m, positive_mask, cmap=cmap_neg,
            shading='gouraud',
            linewidth=0,
            edgecolor='None',
            rasterized=True,
            norm=neg_norm,
        )

    fig.colorbar(mesh, extend='min', label="Differential luminosity ($GeV^{-1}$)")
    ax.set_ylabel('$m_{X}$ (GeV)')
    ax.set_xlabel('y')
    ax.set_yscale('log')
    ax.grid(False)

    ax.set_title("$%s$ luminosity\n%s - "
             "$\\sqrt{s}=%.1f$ GeV" % (LUMI_CHANNELS[lumi_channel],
                     pdf.label, sqrts))

    return fig

@figure
def plot_lumi2d_uncertainty(pdf, lumi_channel, lumigrid2d, sqrts:numbers.Real):
    """
    Plot 2D luminosity unciertainty plot at a given center of mass energy.
    Porting code from https://github.com/scarrazza/lumi2d.

    The luminosity is calculated for positive rapidity, and reflected for
    negative rapidity for display purposes.
    """
    norm = mcolors.LogNorm(vmin=1, vmax=50)
    cmap = copy.copy(cm.viridis_r)
    cmap.set_bad("white", alpha=0)

    grid = lumigrid2d
    channel = lumi_channel

    gv = grid.grid_values
    mat = gv.std_error()/np.abs(gv.central_value())*100

    fig, ax = plt.subplots()

    mat = _reflect_matud(mat)
    y = _reflect_matrl(grid.y, odd=True)


    masked_weights = np.ma.masked_invalid(mat, copy=False)

    mesh = ax.pcolormesh(y, grid.m, masked_weights, norm=norm, cmap=cmap,
                         shading='gouraud',
                         linewidth=0,
                         edgecolor='None',
                         rasterized=True)

    # some extra options
    extup = np.nanmax(masked_weights) > 50
    extdown = np.nanmin(masked_weights) < 1

    #TODO: Wrap this somewhere
    if extup:
        if extdown:
            extend = 'both'
        else:
            extend = 'max'
    elif extdown:
        extend = 'min'
    else:
        extend = None

    fig.colorbar(mesh, label="Relative uncertainty (%)",
        ticks=[1,5,10,25,50], format='%.0f', extend=extend)
    ax.set_yscale('log')
    ax.set_title("Relative uncertainty for $%s$-luminosity\n%s - "
                 "$\\sqrt{s}=%.1f$ GeV" % (LUMI_CHANNELS[channel],
                         pdf.label, sqrts))
    ax.set_ylabel('$m_{X}$ (GeV)')
    ax.set_xlabel('y')
    ax.grid(False)

    return fig


class MixMatchPDFPlotter(BandPDFPlotter):
    """Special wrapper class to plot, in the same figure, PDF bands and PDF replicas
    depending on the type of PDF.
    Practical use: plot together the PDF central values with the NNPDF bands
    """
    def __init__(self, *args, pdfs_plot_replicas, **kwargs):
        self.pdfs_plot_replicas = pdfs_plot_replicas
        super().__init__(*args, **kwargs)

    def draw(self, pdf, grid, flstate):
        if pdf in self.pdfs_plot_replicas:
            labels = flstate.labels
            handles = flstate.handles
            ax = flstate.ax
            next_prop = next(ax._get_lines.prop_cycler)
            color = next_prop['color']
            stats = grid.select_flavour(flstate.flindex).grid_values
            gv = stats.data
            ax.plot(grid.xgrid, gv.T, alpha=0.2, linewidth=0.5,
                    color=color, zorder=1)
            cv_line = ax.plot(grid.xgrid[0:1], stats.central_value()[0:1], 
                    color=color, linewidth=2)
            handle = cv_line[0]
            labels.append(pdf.label)
            handles.append(handle)
            return gv
        return super().draw(pdf, grid, flstate)
