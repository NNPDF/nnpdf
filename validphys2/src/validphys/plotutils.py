# -*- coding: utf-8 -*-
"""
Basic utilities for plotting functions.
Created on Thu Apr 21 18:41:43 2016

@author: Zahari Kassabov
"""
import functools

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.scale as mscale


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

