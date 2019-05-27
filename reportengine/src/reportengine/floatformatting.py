"""
floatformatting.py

Tools to format floating point number properly. This is more difficult than
it looks like.
"""
import decimal
import numbers
from typing import NamedTuple

import numpy as np

def significant_digits(value, digits):
    """Return a `Decimal` object with all the digits less signingicant than
    `digits` trimmed (that is, with floor rounding)."""
    cv = decimal.getcontext().copy()
    cv.prec = digits
    fval =  cv.create_decimal(value)
    return fval

def write_in_adequate_representation(n, minexp = -4, maxexp = 5):
    """Return a decimal string representatin of `n` if its most signigicative
    power of 10 is between ``minexp`` and ``maxexp``. Otherwise return a
    scientific reporesentation.
    Values of ``None``
    for either minexp or maxexp signifies that the value is unbounded"""
    dec = decimal.Decimal(n)
    if not dec.is_finite():
        return str(dec)
    if minexp is None:
        minexp = -float('inf')
    if maxexp is None:
        maxexp = float('inf')
    sigexp = dec.adjusted()
    lowexp = dec.as_tuple().exponent
    if sigexp < 0:
        nexp = lowexp
    else:
        nexp = sigexp


    if nexp < minexp or nexp > maxexp:
        return f'{dec:E}'

    return f'{dec:f}'

def format_number(n, digits=4, minexp=-4):
    """Return a string representation of n with at most ``digits``
    significative figures"""
    if isinstance(n, np.number):
        n = np.asscalar(n)
    sig = significant_digits(n, digits)
    return write_in_adequate_representation(sig, minexp, digits)

def format_value_error(value, error, error_digits=2, **kwargs):
    error = significant_digits(error, error_digits)
    try:
        value = decimal.Decimal(value).quantize(error)
    except decimal.InvalidOperation:
        pass
    return (write_in_adequate_representation(value, **kwargs),
            write_in_adequate_representation(error, **kwargs))

def format_error_value_columns(df, valcol, errcol, inplace=False, **kwargs):
    if not inplace:
        df = df.copy()
    func = lambda x: format_value_error(x[valcol], x[errcol], **kwargs)
    cols = df[[valcol, errcol]]
    #I couldn't find a less stupid way to do this
    df[valcol] = cols.apply(lambda x : func(x)[0], axis=1)
    df[errcol] = cols.apply(lambda x : func(x)[1], axis=1)

    if not inplace:
        return df

class ValueErrorTuple(NamedTuple):
    """A class to represent a value and its error. The only functionality it
    adds is the ``__str__`` method, to represent the quantity."""
    value: numbers.Real
    error: numbers.Real
    def __str__(self):
        valstr, errstr = format_value_error(self.value, self.error)
        return f'{valstr}±{errstr}'
