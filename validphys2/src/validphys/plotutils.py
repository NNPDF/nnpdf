# -*- coding: utf-8 -*-
"""
Basic utilities for plotting functions.
Created on Thu Apr 21 18:41:43 2016

@author: Zahari Kassabov
"""
import functools
import itertools
from collections import namedtuple
import logging

import scipy.stats as stats

import numpy as np
from matplotlib.figure import Figure
import matplotlib as mpl
import matplotlib.scale as mscale
import matplotlib.patches as mpatches
import matplotlib.collections as mcollections
from matplotlib  import transforms
from matplotlib.markers import MarkerStyle
from matplotlib import ticker

from reportengine.floatformatting import format_number

log = logging.getLogger(__name__)

def subplots(figsize=None,nrows=1, ncols =1,sharex=False, sharey=False):
    """
    Use matplotlib.figure.Figure() objects and not import pyplot anywhere. 
    The reason is that pyplot maintains a global state that makes it misbehave 
    in multithreaded applications such when executed under dask parallel mode.

    Parameters
    ----------
    figsize : 2-tuple of floats
            defaults is None
    
    nrows, ncols : int, default 1

    sharex, sharey : bool, default False

    Returns
    -------
    tuple
        fig, ax = (matplotlib.figure.Figure, fig.subplots)
    """
    fig = Figure(figsize=figsize)
    ax = fig.subplots(nrows=nrows, ncols = ncols, sharex=sharex, sharey=sharey)
    return fig, ax

def add_subplot(figsize=None,projection=None):
    """
    Use matplotlib.figure.Figure() objects and not import pyplot anywhere. 
    The reason is that pyplot maintains a global state that makes it misbehave 
    in multithreaded applications such when executed under dask parallel mode.
    
    Parameters
    ----------
    figsize : 2-tuple of floats
            default is None
    projections : The projection type of the subplot (Axes).
                default is None
    
    Returns
    -------
    tuple
        fig, ax = (matplotlib.figure.Figure, fig.add_subplot)
    """
    fig = Figure(figsize=figsize)    
    ax = fig.add_subplot(projection=projection)
    return fig, ax

def ax_or_gca(f):
    """A decorator. When applied to a function, the keyword argument  ``ax``
    will automatically be filled with the current axis, if it was None."""
    @functools.wraps(f)
    def _f(*args, **kwargs):
            if 'ax' not in kwargs or kwargs['ax'] is None:
                kwargs['ax'] = Figure().add_subplot(1, 1, 1)
            return f(*args, **kwargs)
    return _f

def ax_or_newfig(f):
    """A decorator. When applied to a function, the keyword argument  ``ax``
    will automatically be filled with the a new axis corresponding to an empty,
    if it was None."""
    @functools.wraps(f)
    def _f(*args, **kwargs):
        noax = 'ax' not in kwargs or kwargs['ax'] is None
        if noax:
                fig = Figure()
                kwargs['ax'] = fig.add_subplot(1,1,1)
        result = f(*args, **kwargs)
        if noax:
            kwargs['ax'].legend(loc = 'best')
        return result

    return _f

def frame_center(ax, x, values):
    """Set the `ylims` of the axis ``ax`` to appropriately display
    ``values``, which can be 1 or 2D and are assumed to be sampled uniformly
    in the coordinates of the plot (in the second dimension, for 2D arrays)."""
    values = np.atleast_2d(values)
    scale = mscale.scale_factory(ax.xaxis.get_scale(), ax.xaxis)
    t = scale.get_transform()
    tx = t.transform(x)
    absxmax, absxmin =  np.max(tx) , np.min(tx)

    l = absxmax - absxmin
    center = l/2 + absxmin

    dist = np.abs(tx-center)

    #Require some margin around the center
    close_region_mask = np.where(dist < l/10)
    close_vals = values[:, close_region_mask]
    close_anchor_max = np.percentile(np.percentile(close_vals, 95, axis=0), 95)
    close_anchor_min = np.percentile(np.percentile(close_vals, 5, axis=0), 5)
    close_anchor_min, close_anchor_max = expand_margin(close_anchor_min,
                                                       close_anchor_max, 1.4)

    #And to see at least 50% everywhere in t
    medium_region_mask = np.where(dist < l/2.5)
    medium_vals = values[:, medium_region_mask]
    medium_anchor_max = np.percentile(np.percentile(medium_vals, 75, axis=0),95)
    medium_anchor_min = np.percentile(np.percentile(medium_vals, 25, axis=0),5)

    medium_anchor_min, medium_anchor_max = expand_margin(medium_anchor_min,
                                                         medium_anchor_max,
                                                         1.1)

    view_max = max((close_anchor_max, medium_anchor_max))
    view_min = min((close_anchor_min, medium_anchor_min))

    amin, amax = ax.get_ylim()
    #Fix edge cases where the limits are nan or infinite
    #(e.g. when dividing by zero in the whole range)
    if not np.isfinite(view_min):
        view_min = amin
    if not np.isfinite(view_max):
        view_max = amax

    ax.set_ylim(max(view_min, amin), min(view_max, amax))


def expand_margin(a,b,proportion):
    """Return a pair of numbers that have the same mean as ``(a,b)`` and their
    distance is ``proportion`` times bigger."""
    halfdiff = (b-a)/2
    center = a + halfdiff
    expansion = halfdiff*proportion
    return center - expansion, center+expansion

def hatch_iter():
    """An infinite iterator that yields increasingly denser patterns of
    hatches suitable for passing as the ``hatch`` argument of matplotlib
    functions."""
    hatches = "/ \\ - + o O".split()
    i = 1
    while True:
        for hatch in hatches:
            yield hatch*i
        i+=1

def marker_iter_scatter():
    """Yield the possible matplotplib.markers.Markersyle instances with
    different fillsyles and markers. This can be passed to
    ``plt.scatter``.
    For ``plt.plot``, use ``marker_iter_scatter``.
    """
    while True:
        for fill in MarkerStyle.fillstyles:
            for shape in MarkerStyle.filled_markers:
                yield MarkerStyle(marker=shape, fillstyle=fill)

def marker_iter_plot():
    """Because of the mpl strange interface, markers work differently in plots
    and scatter. This is the same as `marker_iter_scatter`, but
    returns kwargs to be passed to ``plt.plot()``"""
    for ms in marker_iter_scatter():
        yield {'marker':ms.get_marker(),'fillstyle': ms.get_fillstyle()}


def color_iter():
    """Yield the colors in the cycle defined in the matplotlib style.  When the
    colores are exhausted a warning will be logged and the cycle will be
    repeated infinitely. Therefore this avoids the overflow error at runtime
    when using matplotlib's ``f'C{i}'`` color specification (equivalent to
    ``colors[i]``) when ``i>len(colors)`` """
    color_list = [prop['color'] for prop in mpl.rcParams['axes.prop_cycle']]
    yield from color_list
    log.warning("Color cycle exhausted. Will repeat colors.")
    yield from itertools.cycle(color_list)


def scalar_log_formatter():
    """Return a matplotlib formatter to display powers of 10 in a log rather
    than exponential notation.

    Returns
    -------
    formatter : ticker.FuncFormatter
        an object that can be passed to the ``set_major_formatter`` matplotlib
        functions.

    Examples
    --------
    >>> from matplotlib.figure import Figure
    >>> fig = Figure()
    >>> ax = fig.subplots()
    >>> ax.plot([0.01, 0.1, 1, 10, 100])
    >>> ax.set_yscale("log")
    >>> ax.yaxis.set_major_formatter(scalar_log_formatter())
    """
    # See https://stackoverflow.com/a/33213196
    def formatter(y, _pos):
        decimalplaces = int(np.maximum(-np.log10(y), 0))  # =0 for numbers >=1
        # Insert that number into a format string
        formatstring = f"{{:.{decimalplaces}f}}"
        # Return the formatted tick label
        return formatstring.format(y)

    return ticker.FuncFormatter(formatter)

HandlerSpec = namedtuple('HandelrSpec', ["color", "alpha", "hatch", "outer"])

class ComposedHandler:
    """Legend artist for PDF plots."""
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height

        patches = []
        if orig_handle.outer:
            wpad = width*0.1
            hpad = height*0.1
            edges = 'none'
            outer = mpatches.Rectangle([x0, y0], width, height,
                                   facecolor='none',
                                   linestyle= 'dashed',
                                   edgecolor = orig_handle.color,
                                   transform=handlebox.get_transform())
            handlebox.add_artist(outer)
            patches.append(outer)
        else:
            wpad = hpad = 0
            edges = 'black'

        patch = mpatches.Rectangle([x0+wpad, y0+hpad],
                                   width-2*wpad, height-2*hpad,
                                   facecolor=orig_handle.color,
                                   alpha = orig_handle.alpha,
                                   hatch=orig_handle.hatch,
                                   edgecolor=edges,
                                   transform=handlebox.get_transform())

        handlebox.add_artist(patch)
        patches.append(patch)
        return patches

def offset_xcentered(n, ax,*, offset_prop=0.05):
    """Yield ``n`` matplotlib transforms in such a way that the corresponding
    ``n`` transofrmed x values are centered around the middle. The offset
    between to consecutive points is ``offset_prop`` in units of the figure
    dpi scale."""
    first_offset = +(n//2)
    #http://matplotlib.org/users/transforms_tutorial.html
    for i in range(n):

            dx = offset_prop*(i-first_offset)
            offset = transforms.ScaledTranslation(dx, 0,
                                                  ax.figure.dpi_scale_trans)
            offset_transform = ax.transData + offset
            yield offset_transform


def centered_range(n, value=0, distance=1):
    """Generte a range of ``n`` points centered
    around ``value``, unifirmely sampled at
    intervals of ``distance``."""
    first_offset = +(n/2) - 0.5
    for i in range(n):
        yield distance*(i-first_offset) + value


def barplot(values, collabels, datalabels, orientation='auto'):
    """The barplot as matplotlib should have it.
    It resizes on overflow.
    ``values``  should be one or two dimensional and should contain the
    values for the barplot. ``collabels`` must have as many elements
    as ``values`` has columns (or total elements if it is one dimensional),
    and contains the labels for each column in the
    bar plot.  ``datalabels`` should have as many elements as values has
    rows, and contains the labels for the individual items to be
    compared. If ``orientation`` is ``"auto"``, the barplot will be
    horizontal or vertical depending on the number of items.
    Otherwise, the orientation can ve fixes as ``"horizontal"`` or
    ``"vertical"``.

    Parameters
    ----------
    values : array of dimensions *M×N* or *N*.
        The input data.
    collabels : Iterable[str] of dimensions N
        The labels for each of the bars.
    datalabels : Iterable[str] of dimensions M or 1
        The label for each of the datasets to be compared.
    orientation : {'auto', 'horizontal', 'vertical'}, 'optional'
        The orientation of the bars.

    Returns
    -------
    (fig, ax) : tuple
        a tuple of a matplotlib figure and an axis, like
        matplotlib.pyplot.subplots. The axis will have a ``_bar_orientation``
        attribute that will either be 'horizontal' or 'vertical' and will
        correspond to the actual orientaion of the plot.

    Examples
    --------
    >>> import numpy as np
    >>> from validphys.plotutils import barplot
    >>> vals = np.random.rand(2,5)
    >>> collabels = ["A", "B", "C", "D", "e"]
    >>> fig, ax = barplot(vals, collabels, ['First try', 'Second try'])
    >>> ax.legend()
    """
    values = np.atleast_2d(values)
    ntypes, l = values.shape
    lc = len(collabels)
    if lc != l:
        raise ValueError(f"Mismatch between the number of data points ({l}) and "
                         f"the number axis labels ({lc})")

    width = 2
    #The tick positions
    x = np.linspace(0, 1.1*width*ntypes*l-1, l)
    w,h = mpl.rcParams["figure.figsize"]

    #Rescale if we have too much data
    rescale = max(1, 1+(width*l*ntypes-15)*0.05)
    #Rescale if the labels are too long
    lbrescale = max(1, 1+0.04*(max(len(l) for l in collabels)-5))

    if orientation == 'auto':
        if l*ntypes > 20:
            orientation = 'horizontal'
        else:
            orientation = 'vertical'

    if orientation == 'vertical':
        fig = Figure(figsize=(w*rescale, h*lbrescale))
        ax = fig.subplots()
        barfunc = ax.bar
        infoaxis = ax.xaxis
        infolim = ax.set_xlim
        otheraxis = ax.yaxis
        rotation = 80
        def get_pos(val):
            if val >= 0:
                xytext = (0,5)
            else:
                xytext = (0,-5)
            horizontalalignment='center'
            return {'xytext': xytext,
                    'horizontalalignment':horizontalalignment}

        def xytext(x,y): return x,y
    elif orientation =='horizontal':
        fig = Figure(figsize=(w*lbrescale, h*rescale))
        ax = fig.subplots()
        barfunc = ax.barh
        infoaxis = ax.yaxis
        infolim = ax.set_ylim
        otheraxis = ax.xaxis
        rotation = 10
        def get_pos(val):
            if val >= 0:
                xytext = (5,0)
                horizontalalignment= 'left'
            else:
                xytext = (-5,0)
                horizontalalignment='right'
            return {'xytext': xytext,
                    'horizontalalignment':horizontalalignment}
        def xytext(x,y): return y,x

    else:
        raise ValueError("orientation must be one of ('auto', 'horizontal', 'vertical')")

    infoaxis.set_ticks(x)

    infoaxis.set_ticklabels(collabels,rotation=rotation)
    deltas = list(centered_range(ntypes, distance=width))
    for row, delta, datalabel in zip(values, deltas, datalabels):
        thisx = x+delta
        barfunc(thisx, row, width, label=datalabel)
        for xp, v in zip(thisx, row):
            #NaN coords cause error (https://github.com/NNPDF/nnpdf/issues/363)
            if np.all(np.isfinite([xp, v])):
                ax.annotate(f'{format_number(v,3)}', xy=xytext(xp,v),
                            textcoords='offset points',
                            size='small', wrap=True, **get_pos(v)
                            )
            else:
            #place label at zero for nan coordinate -> ensure `get_pos` is fed altered coords
                new_pos = [val if np.isfinite(val) else 0 for val in [xp, v]]
                ax.annotate(f'{format_number(v,3)}', xy=xytext(*new_pos),
                            textcoords='offset points',
                            size='small', wrap=True, **get_pos(new_pos[0])
                            )


    infolim(x[0]+deltas[0] - width/2, x[-1]+deltas[-1]+width/2)
    otheraxis.set_visible(False)
    ax.tick_params(length=0)
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.grid(False)
    fig.tight_layout()
    ax._bar_orientation = orientation
    return fig, ax

def plot_horizontal_errorbars(cvs, errors, categorylabels, datalabels=None,
                              xlim=None):
    """A plots with a list of horizontal errorbars oriented vertically.
    ``cvs`` and ``errors`` are the central values and errors both of shape
    `ndatasets x ncategories`, ``cateogorylabels`` are the labels of each
    element for which errorbars are drawn and ``datalabels`` are the labels of
    the different datasets that are compared.
    """
    w,h = mpl.rcParams["figure.figsize"]
    rescale = max(1, 1 + 0.1*(len(categorylabels) - 7))
    fig = Figure(figsize=(w*1.5, h*rescale))
    ax = fig.subplots()
    if datalabels is None:
        datalabels = itertools.repeat(None)
    y = np.arange(len(categorylabels))

    ax.yaxis.set_ticks(y)
    ax.yaxis.set_ticklabels(categorylabels)
    mi = marker_iter_plot()


    distance = 0.5/len(cvs)
    pos = centered_range(len(cvs), distance=distance)

    for cv, err, lb, markerspec, shift in zip(cvs, errors, datalabels, mi, pos):
        ax.errorbar(cv, y+shift, xerr=err, linestyle='none', label=lb,
                    **markerspec)
    if xlim is None:
        #Support both single error and up, low.
        errmat = np.atleast_2d(errors)
        low = cvs - errmat[0,:]
        high = cvs + errmat[-1,:]
        xlim =  expand_margin(np.nanpercentile(low, 15),
                              np.nanpercentile(high,85),
                              1.5)

    ax.set_xlim(xlim)
    ax.set_ylim(-0.5, len(categorylabels)-0.5)
    ax.invert_yaxis()
    ax.grid(axis='y')
    return fig, ax

@ax_or_gca
def kde_plot(a, height=0.05, ax=None, label=None, color=None, max_marks=100000):
    """Plot a Kernel Density Estimate of a 1D array, togther with individual
    occurrences .

    This plot provides a quick visualizaton of the distribution of one
    dimensional data in a more complete way than an histogram would.  It
    produces both a `Kernel Density Estimate (KDE)
    <https://en.wikipedia.org/wiki/Kernel_density_estimation>`_ and individual
    occurences of the data (`rug plot
    <https://en.wikipedia.org/wiki/Rug_plot>`_). The KDE uses a Gaussian Kernel
    with the Silverman rule to select the bandwidth (this is the optimal choice
    if the input data is Gaussian). The individual ocurrences are displayed as
    marks along the bottom axis. For performance reasons, and to avoid
    cluttering the plot, a maximum of ``max_marks`` marks are displayed; if the
    length of the data is bigger, a random sample of ``max_marks`` is taken.

    Parameters
    ----------
    a: vector
       1D array of observations.
    height: scalar, optional
       Height of marks in the rug plot as proportion of the axis height.
    ax: matplotlib axes, optional
       Axes to draw plot into; otherwise grabs current axes.
    label: string, optional
       The label for the legend (note that you have to generate the legend yourself).
    color: optional
       A matplotlib color specification, used for both the KDE and the rugplot. If not given,
       the next in the underlying axis cycle will be consumed and used.
    max_marks: integer, optional
       The maximum number of points that will be displayed individually.

    Returns
    -------
    ax: matplotlib axes
       The Axes object with the plot on it, allowing further customization.

    Example
    -------

    >>> import numpy as np
    >>> dist = np.random.normal(size=100)
    >>> ax = kde_plot(dist)

    """

    a = np.asarray(a).ravel()
    if color is None:
        next_prop = next(ax._get_lines.prop_cycler)
        color = next_prop["color"]
    kde_func = stats.gaussian_kde(a, bw_method="silverman")
    kde_x = np.linspace(*expand_margin(np.min(a), np.max(a), 1.3), 100)
    ax.plot(kde_x, kde_func(kde_x), label=label, color=color)
    if len(a) > max_marks:
        segment_data = np.random.choice(a, max_marks, replace=False)
    else:
        segment_data = a
    # Define segments by a pair of points (data, 0) and (data, height).
    # Bug in pylint
    # pylint: disable=too-many-function-args
    segments = np.c_[
        segment_data,
        np.zeros_like(segment_data),
        segment_data,
        np.full_like(segment_data, height),
    ].reshape(-1, 2, 2)
    rugs = mcollections.LineCollection(
        segments,
        # Make the x coordinate refer to the data but the y (the height
        # relative to the plot height.
        # https://matplotlib.org/tutorials/advanced/transforms_tutorial.html?highlight=transforms%20blended_transform_factory#blended-transformations
        transform=transforms.blended_transform_factory(ax.transData, ax.transAxes),
        color=color,
    )
    ax.add_collection(rugs)
    ax.set_ylim(ymin=0)
    return ax

@ax_or_gca
def spiderplot(xticks, vals, label, ax=None):
    """
    Makes a spider/radar plot.

    xticks: list of names of x tick labels, e.g. datasets
    vals: list of values to plot corresponding to each xtick
    label: label for values, e.g. fit name
    """
    N = len(xticks)

    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    # Add this on so that the plot line connects back to the start
    angles += angles[:1]
    vals = list(vals)
    vals += vals[:1]

    maxval = np.max(vals)

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    ax.set_ticks(angles[:-1])
    ax.set_xticklabels(xticks, size=8, zorder=6)

    # Draw ylabels

    ax.plot(angles, vals, linewidth=2, label=label, linestyle="solid", zorder=1)
    ax.fill(angles, vals, alpha=0.4, zorder=1)
    ax.grid(linewidth=1)
    ax.legend(fontsize=12)

    return ax
