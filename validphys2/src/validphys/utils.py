# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 21:12:41 2016

@author: Zahari Kassabov
"""

def split_by(it, crit):

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
