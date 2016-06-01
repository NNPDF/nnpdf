# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 12:49:39 2016

@author: Zahari Kassabov
"""

def new_labels(k1label, k2lebel, k3label):
    def closure(f):
        f.new_labels = k1label, k2lebel, k3label
        return f
    return closure
