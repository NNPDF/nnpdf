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

class MissingLabelError(KeyError):
    def __init__(self, key_error):
        msg = "A label is required to perform the operation: %s" % key_error.args[0]
        super().__init__(msg)

def qbinexp(cv, error, **labels):
    try:
        qbin = labels['qbin']
    except KeyError as e:
        raise MissingLabelError(e)
    return 10**qbin*cv, 10**qbin*error

def qbindis(cv, error, **labels):
    qbin = labels['k2']
    return (10**qbin)*cv, (10**qbin)*error

def qbinjets(cv, error, **labels):
    qbin = labels['k1']
    return 1000**(5-qbin)*cv, 1000**(5-qbin)*error

def qbindyp(cv, error, **labels):
    qbin = labels['k1']
    return 10000**(qbin)*cv, 10000**(qbin)*error
