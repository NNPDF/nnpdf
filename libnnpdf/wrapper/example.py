
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Using libnnpdf in python
This requires the FK table 'FK_ATLASR04JETS2P76TEV.dat'
to be in the current directory
"""

__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.1'
__email__ = 'stefano.carrazza@cern.ch'

from nnpdf import (LHAPDFSet,
                   PDFSet,
                   FKTable,
                   ThPredictions,)

## Loading NNPDF from NNPDF::LHAPDFSet
pdf = LHAPDFSet('NNPDF30_nlo_as_0118', LHAPDFSet.erType_ER_MC)

## Loading FKTable
filename = 'FK_ATLASR04JETS2P76TEV.dat'
fk  = FKTable(filename, [])

## Build Theoretical Predictions
thp = ThPredictions(pdf,fk)

## Retrieving CV predictions
import numpy as np
cvs = np.array([ thp.GetObsCV(el) for el in range(thp.GetNData()) ])
ers = np.array([ thp.GetObsError(el) for el in range(thp.GetNData()) ])

for i in range(len(cvs)):
    print('N=%.2i CV=%.2e ERROR=%.2e' %(i, cvs[i], ers[i]))

## Plot results and save to file
import matplotlib.pyplot as plt
plt.figure()
plt.errorbar(range(len(cvs)), cvs/cvs, yerr=ers/cvs,
             fmt='o', label=thp.GetSetName())
plt.legend(loc='best')
plt.title('NNPDF::ThPredictions with ' + thp.GetPDFName())
plt.ylabel('Prediction / CV')
plt.xlabel('Data Point')
plt.savefig('%s_%s.pdf' % (thp.GetSetName(),thp.GetPDFName()))
