# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 12:48:55 2016

@author: Zahari Kassabov
"""

from validphys.plotoptions import utils
from numpy import sqrt, ceil, nditer

@utils.new_labels('$2k_1$', '$3k_1$', '$4k_3$')
def dummy_transform(k1,k2,k3):
    return k1*2, k2*3, k3*4

@utils.new_labels('$y$', '$M$ (GeV)', '$\sqrt(s)$')
def dyp_sqrt_scale(k1,k2,k3):
    return k1, sqrt(k2), k3

@utils.new_labels('$|y|$', '$p_T$ (GeV)', '$\sqrt(s)$')
def jet_sqrt_scale(k1,k2,k3):
    return k1, sqrt(k2), k3

@utils.new_labels('$p_T$ (GeV)', '$|\eta|$', '$\sqrt(s)$')
def new_eta_labels(k1,k2,k3):
    return k1, k2, k3

@utils.new_labels('$p_T$ (GeV)', '$M$', '$\sqrt(s)$')
def new_mass_labels(k1,k2,k3):
    return k1, k2, k3

@utils.new_labels('$x$', '$Q$ (GeV)', '$\sqrt{(s)}$')
def dis_sqrt_scale(k1,k2,k3):
    ecm = sqrt(k2/(k1*k3))
    return k1, sqrt(k2), ceil(ecm)

@utils.new_labels('$x$', '$Q$ (GeV)', '$\sqrt{(s)}$')
def nmc_process(k1,k2,k3):
    xBins = [0.0045, 0.008, 0.0125, 0.0175,
    0.025, 0.035, 0.05, 0.07, 0.09, 0.11,
    0.14, 0.18, 0.225, 0.275, 0.35, 0.5]
    for x in nditer(k1, op_flags=['readwrite']):
        x[...] = min(xBins, key=lambda y:abs(x-y))
    ecm = sqrt(k2/(k1*k3))
    return k1, sqrt(k2), ceil(ecm)