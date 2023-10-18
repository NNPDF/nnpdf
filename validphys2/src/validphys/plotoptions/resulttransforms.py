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


"""

import numpy

from validphys.plotoptions.utils import bins


class MissingLabelError(KeyError):
    def __init__(self, key_error):
        msg = "A label is required to perform the operation: %s" % key_error.args[0]
        super().__init__(msg)


def half(cv, error, **labels):
    return cv / 2, error / 2


# TODO: Refactor these so we don't write the same code over and over?


def qbinEMC(cv, error, **labels):
    q = labels['k2']
    qbin = numpy.sqrt(q)
    k = float(10) ** qbin
    return k * cv, k * error


def qbinexp(cv, error, **labels):
    q = labels['k2']
    qbin = bins(q)
    k = float(10) ** qbin
    return k * cv, k * error


def qbindis(cv, error, **labels):
    q = labels['k1']
    qbin = bins(q)
    k = float(10) ** (10 - qbin)
    return k * cv, k * error


def qbinjets(cv, error, **labels):
    qbin = labels['k1']
    k = float(1000) ** (5 - qbin)
    return k * cv, k * error


def qbindyp(cv, error, **labels):
    qbin = labels['k1']
    k = float(10000) ** (qbin)
    return k * cv, k * error
