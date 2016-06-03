# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 09:52:43 2016

@author: Zahari Kassabov


Transform the result (central value and error) for plotting.

The functions here receive as arguments the original cental value and
error (from e.g. the convolution or commondata) as well as all the labels
defined in the plotting file, as keyword arguments.
They are expected to return a new central value and a new error. Therefore
the signature is:


.. code-block:: python

    def xbinexp(cv, error, **labels):
        ...
        return newcv, newerror

Note that these functions will be called point by point
(rather than on lists of values).

"""

from validphys.plotoptions.utils import bins
import numpy

class MissingLabelError(KeyError):
    def __init__(self, key_error):
        msg = "A label is required to perform the operation: %s" % key_error.args[0]
        super().__init__(msg)


def qbinEMC(cv, error, **labels):
    q = labels['k2']
    qbin = numpy.sqrt(q)
    return (10**qbin)*cv, (10**qbin)*error

def qbinexp(cv, error, **labels):
    q = labels['k2']
    qbin = bins(q)
    return 10**qbin*cv, 10**qbin*error

def qbindis(cv, error, **labels):
    qbin = labels['k2']
    return (10**qbin)*cv, (10**qbin)*error

def qbinNMC(cv, error, **labels):
    NMCdict = { 0.0015:0, 0.0030:1, 0.0050:2, 0.0080:3, 0.0125:4, 0.0175:5, 0.025:6, 0.035:7, 0.050:8, 0.070:9, 0.090:10, 0.110:11, 0.140:12, 0.180:13, 0.225:14, 0.275:15, 0.350:16, 0.450:17, 0.550:18, 0.675:19 }
    qbin = NMCdict[labels['k1']]
    return 10**(qbin)*cv, 10**(qbin)*error

def qbinjets(cv, error, **labels):
    qbin = labels['k1']
    return 1000**(5-qbin)*cv, 1000**(5-qbin)*error

def qbindyp(cv, error, **labels):
    qbin = labels['k1']
    return 10000**(qbin)*cv, 10000**(qbin)*error
