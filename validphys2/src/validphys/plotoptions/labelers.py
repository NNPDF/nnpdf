# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 20:44:21 2016

@author: Zahari Kassabov
"""

def high_xq(k1, k2, k3, **kwargs):
    return k1 > 1e-2 and k2 > 1000

def pt_ratio(k1, k2, k3 , **kwargs):
    return k1/k2
