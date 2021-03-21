# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 21:12:41 2016

@author: Zahari Kassabov
"""
import numpy as np

from validobj import parse_input, ValidationError


def parse_yaml_inp(inp, spec):
    """Helper function to parse yaml using the `validobj` library and print
    useful error messages in case of a parsing error.

    https://validobj.readthedocs.io/en/latest/examples.html#yaml-line-numbers
    """
    try:
        return parse_input(inp, spec)
    except ValidationError as e:
        current_exc = e
        current_inp = inp
        error_text_lines = []
        while current_exc:
            if hasattr(current_exc, 'wrong_field'):
                wrong_field = current_exc.wrong_field
                # Mappings compping from ``round_trip_load`` have an
                # ``lc`` attribute that gives a tuple of
                # ``(line_number, column)`` for a given item in
                # the mapping.
                line = current_inp.lc.item(wrong_field)[0]
                error_text_lines.append(f"Problem processing key at line {line}:")
                current_inp = current_inp[wrong_field]
            elif hasattr(current_exc, 'wrong_index'):
                wrong_index = current_exc.wrong_index
                # Similarly lists allow to retrieve the line number for
                # a given item.
                line = current_inp.lc.item(wrong_index)[0]
                current_inp = current_inp[wrong_index]
                error_text_lines.append(f"Problem processing list item at line {line}:")
            elif hasattr(current_exc, 'unknown'):
                unknown_lines = []
                for u in current_exc.unknown:
                    unknown_lines.append((current_inp.lc.item(u)[0], u))
                unknown_lines.sort()
                for line, key in unknown_lines:
                    error_text_lines.append(
                        f"Unknown key {key!r} defined at line {line}:"
                    )
            error_text_lines.append(str(current_exc))
            current_exc = current_exc.__cause__
        raise ValidationError('\n'.join(error_text_lines)) from e

def split_by(it, crit):
    """Split ``it`` in two lists, the first is such that ``crit`` evaluates to
    True and the second such it doesn't. Crit can be either a function or an
    iterable (in this case the original ``it`` will be sliced if the length of
    ``crit`` is smaller)."""

    true, false = [], []
    if callable(crit):
        for ele in it:
            if crit(ele):
                true.append(ele)
            else:
                false.append(ele)
    elif hasattr(crit, '__iter__'):
        for keep, ele in zip(crit,it):
            if keep:
                true.append(ele)
            else:
                false.append(ele)
    else:
        raise TypeError("Crit must be  a function or a sequence")

    return true, false

#Copied from smpdf.utils
def split_ranges(a,cond=None,*, filter_falses=False):
    """Split ``a`` so that each range has the same
    value for ``cond`` . If ``filter_falses`` is true, only the ranges
    for which the
    condition is true will be returned."""
    if cond is None:
        cond = a
    cond = cond.astype(bool)
    d = np.r_[False, cond[1:]^cond[:-1]]
    split_at = np.argwhere(d)
    splits = np.split(a, np.ravel(split_at))
    if filter_falses:
        #Evaluate condition at split points
        it = iter(cond[np.r_[0, np.ravel(split_at)]])
        return [s for s in splits if next(it)]
    else:
        return splits


def sane_groupby_iter(df, by, *args, **kwargs):
    """Iterate groupby in such a way that  first value is always the tuple of
    the common values.

    As a concenience for plotting, if by is None, yield the empty string and
    the whole dataframe.
    """
    if by is None or not by:
        yield ('',), df
        return
    gb = df.groupby(by, *args,**kwargs)
    for same_vals, table in gb:
        if not isinstance(same_vals, tuple):
            same_vals = (same_vals,)
        yield same_vals, table

def common_prefix(*s):
    """Return the longest string that is a prefix to both s1 and s2"""
    small, big = min(s), max(s)
    for i, c in enumerate(small):
        if big[i] != c:
            return small[:i]
    return small

def scale_from_grid(grid):
    """Guess the appropriate matplotlib scale from a grid object.
    Returns ``'linear'`` if the scale of the grid object is linear,
    and otherwise ``' log'``."""
    return 'linear' if grid.scale == 'linear' else 'log'
