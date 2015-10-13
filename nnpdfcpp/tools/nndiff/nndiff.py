# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 11:39:05 2015

@author: zah
"""

import sympy
import numpy as np
import numba

import nnets

DEFAULT_ARCH = (2,5,3,1)

def nnpdf_derivative(params, xvals, arch=None):
    if arch is None:
        arch = DEFAULT_ARCH
    n = nnets.build_network(*arch)
    n._set_params(params)
    x, lx = list(n.input_symbols)
    form = n.total_output_formula.subs(lx, sympy.log(x))
    diff = form.diff(x)
    sdiff = diff.subs({s:val for s, val in zip(n.parameter_symbols, 
                                               n.parameter_values)})
    df = sympy.lambdify(x, sdiff)
    #Compile func
    df = numba.vectorize(df)
    return df(np.array(xvals, copy=False))
    #else
    #return list(map(df, xvals))
