# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 20:44:21 2016

@author: Zahari Kassabov
"""
import numpy as np

from validphys.plotoptions.utils import bins, label


@label(r"$I(x>10^{-2})\times I(Q > 1000 GeV)$")
def high_xq(k1, k2, k3, **kwargs):
    return (k1 > 1e-2) & (k2 > 1000)


def pt_ratio(k1, k2, k3, **kwargs):
    return k1 / k2


def jet_eta(k1, k2, k3, **kwargs):
    return k1


def k2bins5(k1, k2, k3, **kwargs):
    qbin = bins(k2)
    qbin[:] = [int(x / 5) for x in qbin]
    return qbin


def k2bins6(k1, k2, k3, **kwargs):
    qbin = bins(k2)
    qbin[:] = [int(x / 6) for x in qbin]
    return qbin


def k2bins10(k1, k2, k3, **kwargs):
    qbin = bins(k2)
    qbin[:] = [int(x / 10) for x in qbin]
    return qbin


@label("$Q^2$ (GeV²)")
def two_Q2_bins(k1, k2, k3, **kwargs):
    min_, median, max_ = np.percentile(k2, (0, 50, 100))
    firstlabel = '[%.2f, %.2f)' % (min_, median)
    # Use dtype=object to avoid longer strings in secondlabel getting trimmed
    res = np.array([firstlabel] * len(k2), dtype=object)

    secondlabel = '[%.2f, %.2f]' % (median, max_)
    res[k2 >= median] = secondlabel
    return res
