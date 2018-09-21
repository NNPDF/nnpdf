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

from validphys import plotutils
from validphys.core import MCStats
from validphys.gridvalues import LUMI_CHANNELS
from validphys.utils import scale_from_grid
from validphys.checks import check_pdf_normalize_to, check_scale, check_have_two_pdfs

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
            normvals = normalize_pdf.stats_class(
                            normalize_grid.grid_values).central_value()

            #Handle division by zero more quietly
            def fp_error(tp, flag):
                log.warn("Invalid values found computing normalization to %s: "
                 "Floating point error (%s).", normalize_pdf, tp)
                #Show warning only once
                np.seterr(all='ignore')

            newgrids = []
            with np.errstate(all='call'):
                np.seterrcall(fp_error)
                for grid in self._xplotting_grids:
                    newvalues = grid.grid_values/normvals
                    #newgrid is like the old grid but with updated values
                    newgrid = type(grid)(**{**grid._asdict(),
                                             'grid_values':newvalues})
                    newgrids.append(newgrid)

            return newgrids
        return self._xplotting_grids

    @property
    def normalize_pdf(self):
        if self.normalize_to is None:
            raise AttributeError("Need to set a normalize_to index.")
        return self.pdfs[self.normalize_to]

    def get_ylabel(self, parton_name):
        if self.normalize_to is not None:
            return "Ratio to {}".format(self.normalize_pdf.label)
        else:
            return '$x{}(x)$'.format(parton_name)

    def get_title(self, parton_name):
        return "$%s$ at %.1f GeV" % (parton_name, self.Q)

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
    et = pdf.ErrorType
    if et != 'replicas':
        log.warn("Plotting members of a non-Monte Carlo PDF set:"
        " %s with error type '%s'.", pdf.name, et)

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
        gv = grid.grid_values[:,flstate.flindex,:]
        ax.plot(grid.xgrid, gv.T, alpha=0.2, linewidth=0.5,
                color=color, zorder=1)
        stats = pdf.stats_class(gv)
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
    """Plot the replicas of the specifid PDFs. Otherise it works the same as
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
        return r"$\sigma/\sigma_{ref}$"

    def draw(self, pdf, grid, flstate):
        ax = flstate.ax
        flindex = flstate.flindex
        gv = grid.grid_values[:,flindex,:]
        stats = pdf.stats_class(gv)

        res = stats.std_error()

        ax.plot(grid.xgrid, res, label=pdf.label)

        return res

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

        plotutils.frame_center(ax, self.firstgrid.xgrid, np.concatenate(all_vals))
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

        gv = grid.grid_values[flindex,:]
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
    def setup_flavour(self, flstate):
        flstate.handles=[]
        flstate.labels=[]
        flstate.hatchit=plotutils.hatch_iter()

    def draw(self, pdf, grid, flstate):
        ax = flstate.ax
        flindex = flstate.flindex
        hatchit = flstate.hatchit
        labels = flstate.labels
        handles = flstate.handles
        stats = pdf.stats_class(grid.grid_values[:,flindex,:])
        pcycler = ax._get_lines.prop_cycler
        #This is ugly but can't think of anything better

        next_prop = next(pcycler)
        cv = stats.central_value()
        xgrid = grid.xgrid
        #Ignore spurious normalization warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            err68down, err68up = stats.errorbar68()

        color = next_prop['color']
        ax.plot(xgrid, cv, color=color)
        alpha = 0.5
        ax.fill_between(xgrid, err68up, err68down, color=color, alpha=alpha,
                        zorder=1)
        #http://stackoverflow.com/questions/5195466/matplotlib-does-not-display-hatching-when-rendering-to-pdf
        hatch = next(hatchit)
        ax.fill_between(xgrid, err68up, err68down, facecolor='None', alpha=alpha,
                        edgecolor=color,
                        hatch=hatch,
                        zorder=1)
        if isinstance(stats, MCStats):
            errorstdup, errorstddown = stats.errorbarstd()
            ax.plot(xgrid, errorstdup, linestyle='--', color=color)
            ax.plot(xgrid, errorstddown, linestyle='--', color=color)
            label  = rf"{pdf.label} ($68%$ c.l.+$1\sigma$)"
            outer = True
        else:
            outer = False
            label = rf"{pdf.label} ($68\%$ c.l.)"
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
@check_scale('xscale', allow_none=True)
def plot_pdfs(pdfs, xplotting_grids, xscale:(str,type(None))=None,
                      normalize_to:(int,str,type(None))=None,ymin=None,ymax=None):
    """Plot the central value and the uncertainty of a list of pdfs as a
    function of x for a given value of Q. If normalize_to is given, plot the
    ratios to the corresponding PDF. Otherwise, plot absolute values.
    See the help for ``xplotting_grid`` for information on how to set basis,
    flavours and x ranges. Yields one figure per PDF flavour.

    normalize_to:  Either the name of one of the PDFs or its corresponding
    index in the list, starting from one, or None to plot absolute values.

    xscale: One of the matplotlib allowed scales. If undefined, it will be
    set based on the scale in xgrid, which should be used instead.

    """
    yield from BandPDFPlotter(pdfs, xplotting_grids, xscale, normalize_to, ymin, ymax)


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
def plot_lumi1d(pdfs, pdfs_lumis,
                lumi_channel, sqrts:numbers.Real,
                normalize_to=None):
    """Plot PDF luminosities at a given center of mass energy.
    sqrts is the center of mass energy (GeV).
    """

    fig, ax = plt.subplots()
    if normalize_to is not None:
        norm = pdfs_lumis[normalize_to].grid_values.central_value()
        ylabel = f"Ratio to {pdfs[normalize_to]}"
    else:
        norm = 1
        ylabel = r"$L(GeV^{-2})$"

    for pdf, lumigrid1d in zip(pdfs, pdfs_lumis):
        mx = lumigrid1d.m
        gv = lumigrid1d.grid_values

        cv = gv.central_value()
        err = gv.std_error()

        ax.fill_between(mx, (cv -err)/norm, (cv+err)/norm, alpha=0.5)
        ax.plot(mx, cv/norm, label='%s' % pdf.label)
    ax.legend(loc='best')
    ax.set_ylabel(ylabel)
    ax.set_xlabel('$M_{X}$ (GeV)')
    ax.set_xscale('log')
    ax.grid(False)
    ax.set_title("$%s$ luminosity\n"
                 "$\\sqrt{s}=%.1f$ GeV" % (LUMI_CHANNELS[lumi_channel],
                                           sqrts))

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
    ax.set_ylabel('$M_{X}$ (GeV)')
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
    ax.set_ylabel('$M_{X}$ (GeV)')
    ax.set_xlabel('y')
    ax.grid(False)

    return fig