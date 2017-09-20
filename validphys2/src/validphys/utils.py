# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 21:12:41 2016

@author: Zahari Kassabov
"""
import numpy as np

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

def common_prefix(s1,s2):
    """Return the longest string that is a prefix to both s1 and s2"""
    small, big = sorted((s1,s2), key=len)
    i = 0
    for test in range(0, len(small)):
        if big[test] != small[test]:
            break
        i += 1
    return small[:i]
