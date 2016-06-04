# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 12:49:39 2016

@author: Zahari Kassabov
"""
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

def apply_to_all_columns(df, func):
    params = dict((col,df[col].as_matrix()) for col in df.columns)
    return func(**params)
