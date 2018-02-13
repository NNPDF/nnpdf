#!/usr/bin/python

import numpy as np

with open("output/tables/default_theory0_experiments_corrmat.csv") as f:
    ncols = len(f.readline().split(','))

corrmat = np.loadtxt(open("output/tables/default_theory0_experiments_corrmat.csv"), delimiter="\t", skiprows=8, usecols=range(8,ncols))

print(corrmat[1,1])
