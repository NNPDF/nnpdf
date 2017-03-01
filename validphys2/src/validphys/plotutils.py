# -*- coding: utf-8 -*-
"""
Basic utilities for plotting functions.
Created on Thu Apr 21 18:41:43 2016

@author: Zahari Kassabov
"""
import functools
from collections import namedtuple

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.scale as mscale
import matplotlib.patches as mpatches
from matplotlib  import transforms


def setup_ax(ax):
    """Change properties that are not correctly handled with styles.
    Eventually this should be deprecated."""

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


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
    ``n`` transofrmed x values are centeres around the middle. The offset
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
