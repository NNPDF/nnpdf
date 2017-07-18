# -*- coding: utf-8 -*-
"""
Basic utilities for plotting functions.
Created on Thu Apr 21 18:41:43 2016

@author: Zahari Kassabov
"""
import functools
import itertools
from collections import namedtuple

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.scale as mscale
import matplotlib.patches as mpatches
from matplotlib  import transforms
from matplotlib.markers import MarkerStyle

def ax_or_gca(f):
    @functools.wraps(f)
    def _f(*args, **kwargs):
            if 'ax' not in kwargs or kwargs['ax'] is None:
                kwargs['ax'] = plt.gca()
            return f(*args, **kwargs)
    return _f

def ax_or_newfig(f):
    @functools.wraps(f)
    def _f(*args, **kwargs):
        noax = 'ax' not in kwargs or kwargs['ax'] is None
        if noax:
                plt.figure()
                kwargs['ax'] = plt.gca()
        result = f(*args, **kwargs)
        if noax:
            plt.legend(loc = 'best')
        return result

    return _f

def frame_center(ax, x, values):
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
    ax.set_ylim(max(view_min, amin), min(view_max, amax))


def expand_margin(a,b,proportion):
    halfdiff = (b-a)/2
    center = a + halfdiff
    expansion = halfdiff*proportion
    return center - expansion, center+expansion

def hatch_iter():
    hatches = "/ \\ - + o 0".split()
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
    """The barplot as matplotlib should have it. It resizes on overflow.
    ``values``  should be one or two dimensional and should contain the
    values for the barplot. ``colllabels`` must have as many element
    s ``values`` has rows, and contains the labels for each column in the
    bar plot.  ``datalabels`` should have as many elements as values has
    columns, and contains the labels for the individual items to be
    compared. If ``orientation`` is ``"auto"``, the barplot will be
    horizontal or vertical depending on the number of items.
    Otherwise, the orientation can ve fixes as ``"horizontal"`` or
    ``"vertical"``

    Returns a tuple figure, axis like plt.subplots.
    """
    values = np.atleast_2d(values)
    ntypes, l = values.shape

    width = 2
    #The tick positions
    x = np.linspace(0, 1.1*width*ntypes*l-1, l)
    w,h = plt.rcParams["figure.figsize"]

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
        fig, ax = plt.subplots(figsize=(w*rescale, h*lbrescale))
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
        fig, ax = plt.subplots(figsize=(w*lbrescale, h*rescale))
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
        for xp,v in zip(thisx,row):

            ax.annotate(f'{v:.2f}', xy=xytext(xp,v),
                         textcoords='offset points',
                        size='small', wrap=True, **get_pos(v)
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

    return fig, ax

def plot_horizontal_errorbars(cvs, errors, categorylabels, datalabels=None):
    """A plots with a list of horizontal errorbars oriented vertically.
    ``cvs`` and ``errors`` are the central values and errors both of shape
    `ndatasets x ncategories`, ``cateogorylabels`` are the labels of each
    element for which errorbars are drawn and ``datalabels`` are the labels of
    the different datasets that are compared.
    """
    w,h = plt.rcParams["figure.figsize"]
    rescale = max(1, 1 + 0.1*(len(categorylabels) - 7))
    fig, ax = plt.subplots(figsize=(w, h*rescale))
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
    ax.set_xlim(*expand_margin(np.nanpercentile(cvs, 15),
                               np.nanpercentile(cvs, 85),
                               1.1))

    ax.set_ylim(-0.5, len(categorylabels)+0.5)
    ax.grid(axis='y')
    return fig, ax
