# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 15:41:00 2015

@author: zah
"""

import sympy

import nnets
DEFAULT_ARCH = (2,5,3,1)

DEFAULT_REAL_TYPE = 'double'


template = (
"""
#include <math.h>

{float_t} nnval({float_t} {x}, {float_t} * params)
{{
    {defs}
    return {valexpr};
}}

{float_t} nndiff({float_t} {x}, {float_t} * params)
{{
    {defs}
    return {diffexpr};
}}
"""
)

def defs_from_list(symbols, label, identation="    ", 
                   float_t=DEFAULT_REAL_TYPE):
    statement = "{float_t} {symbol} = {label}[{i}];"
    sep = "\n" + identation
    return sep.join(statement.format(symbol=symbol, i=i, label=label,
                                     float_t=float_t) 
                                     for i,symbol in enumerate(symbols))


def nnpdf_derivative_to_c(arch=None, float_t=DEFAULT_REAL_TYPE):
    if arch is None:
        arch = DEFAULT_ARCH
    n = nnets.build_network(*arch)
    x, lx = list(n.input_symbols)
    #Add prepo here
    form = n.total_output_formula.subs(lx, sympy.log(x))
    valexpr = sympy.ccode(form)    
    diff = form.diff(x)
    defs = defs_from_list(n.parameter_symbols, 'params')
    diffexpr = sympy.ccode(diff)
    return template.format(x=x,defs=defs, diffexpr=diffexpr, 
                           float_t=float_t,
                           valexpr=valexpr)
    
if __name__ == '__main__':
    source = nnpdf_derivative_to_c()
    with open('nndiff.c', 'w') as f:
        f.write(source)
