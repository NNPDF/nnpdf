#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""\ 
Using libnnpdf in python:
Before running this code do:
  compile the libnnpdf
  python setup.py install
  python example.py
"""

__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@cern.ch'

from nnpdf import *

## Loading NNPDF from NNPDF::LHAPDFSet
pdf = LHAPDFSet('NNPDF30_nlo_as_0118', PDFSet.ER_MC)

## Loading FKTable
filename = "../../nnpdfcpp/data/theory_2/fastkernel/FK_ATLASR04JETS36PB.dat"
fk  = FKTable(filename,[])

## Build Theoretical Predictions
thp = ThPredictions(pdf,fk)

## Retrieving CV predictions
import numpy as np
cvs = np.array([ thp.GetObsCV(el) for el in range(thp.GetNData) ])
ers = np.array([ thp.GetObsError(el) for el in range(thp.GetNData) ])

## Plot results and save to file
import matplotlib.pyplot as plt
plt.figure()
plt.errorbar(range(len(cvs)), cvs/cvs, yerr=ers/cvs,
             fmt='o', label=thp.GetSetName)
plt.legend(loc='best')
plt.title('NNPDF::ThPredictions with ' + thp.GetPDFName)
plt.ylabel('Prediction [pb/bin] / CV')
plt.xlabel('Data Point')
plt.savefig(thp.GetSetName + '_' + thp.GetPDFName + '.pdf')
