# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 20:44:21 2016

@author: Zahari Kassabov
"""
from validphys.plotoptions.utils import label

@label(r"$I(x>10^{-2})\times I(Q > 1000 GeV)$")
def high_xq(k1, k2, k3, **kwargs):
    return k1 > 1e-2 and k2 > 1000

def pt_ratio(k1, k2, k3 , **kwargs):
    return k1/k2

def jet_eta(k1,k2,k3,**kwargs):
    return k1