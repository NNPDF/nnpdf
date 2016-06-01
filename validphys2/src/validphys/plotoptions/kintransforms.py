# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 12:48:55 2016

@author: Zahari Kassabov
"""

from validphys.plotoptions import utils
from math import sqrt

@utils.new_labels('$2k_1$', '$3k_1$', '$4k_3$')
def dummy_transform(k1,k2,k3):
    return k1*2, k2*3, k3*4

@utils.new_labels('$k_1$', '$\sqrt{k_2}$', '$k_3$')
def sqrt_scale(k1,k2,k3):
    return k1, sqrt(k2), k3