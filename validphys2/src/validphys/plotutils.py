# -*- coding: utf-8 -*-
"""
Basic utilities for plotting functions.
Created on Thu Apr 21 18:41:43 2016

@author: Zahari Kassabov
"""
import functools

import matplotlib.pyplot as plt

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
