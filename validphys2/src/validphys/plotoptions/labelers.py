# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 20:44:21 2016

@author: Zahari Kassabov
"""
from validphys.plotoptions.utils import label
from validphys.plotoptions.utils import bins

@label(r"$I(x>10^{-2})\times I(Q > 1000 GeV)$")
def high_xq(k1, k2, k3, **kwargs):
    return (k1 > 1e-2) & (k2 > 1000)

def pt_ratio(k1, k2, k3 , **kwargs):
    return k1/k2

def jet_eta(k1,k2,k3,**kwargs):
    return k1

def k2bins5(k1,k2,k3,**kwargs):
    qbin = bins(k2)
    qbin[:] = [int(x / 5) for x in qbin]
    return qbin

def k2bins6(k1,k2,k3,**kwargs):
    qbin = bins(k2)
    qbin[:] = [int(x / 6) for x in qbin]
    return qbin

def k2bins10(k1,k2,k3,**kwargs):
    qbin = bins(k2)
    qbin[:] = [int(x / 10) for x in qbin]
    return qbin