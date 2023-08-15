# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 12:49:39 2016

@author: Zahari Kassabov
"""
import collections
import inspect

import numpy as np

def bins(arr):
    """Return bins corresponding to unique values of ``arr`` sorted by value.

    .. code-block:: python

        bins([-3, 5, -3, -3, 0, 1,1,0])

        array([0, 3, 0, 0, 1, 2, 2, 1])

    """
    arr = np.atleast_1d(arr)
    return np.unique(arr, return_inverse=True)[1]

def new_labels(k1label, k2lebel, k3label):
    def closure(f):
        f.new_labels = k1label, k2lebel, k3label
        return f
    return closure

def label(label):
    def closure(f):
        f.label = label
        return f
    return closure

def get_subclasses(obj, base):
    """Return the classes in ``obj`` that are subclasses of ``base``"""
    predicate =  lambda x: inspect.isclass(x) and issubclass(x, base)
    return collections.OrderedDict(inspect.getmembers(obj, predicate))


def apply_to_all_columns(df, func):
    """Apply a function to all columns of a dataframe at the same time.
    The parameter names are the names of the column and the values are arrays
    containing the each column's values."""
    params = dict((col,df[col].values) for col in df.columns)
    #import ipdb; ipdb.set_trace()
    return func(**params)
