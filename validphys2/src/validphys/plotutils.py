# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 18:41:43 2016

@author: Zahari Kassabov
"""

def setup_ax(ax):
    """Change properties that are not correctly handled with styles.
    Eventually this should be deprecated."""

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
