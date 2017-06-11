"""
gridvalues.py

Core functionality needed to obtain a set of values from
LHAPDF. The tools for representing these grids are in pdfgrids.py
(the validphys provider module), and the
basis transformations are in pdfbases.py
"""
import itertools

import numpy as np

from NNPDF import _lhapdfset, LHAPDFSet

from validphys.core import PDF

#Numpy is unhappy about downcasting double to float implicitly, so we have
#to manually make sure all input arrays correspond to NNPDF::real.
FLTYPE = np.int32
REALTYPE = np.double if _lhapdfset.REALDOUBLE else np.float32

#This represents some canonical ordering of all the relevant flavours
#we may want to query from LHAPDF
ALL_FLAVOURS = (-6,-5,-4,-3,-2,-1,21,1,2,3,4,5,6,22)

LUMI_CHANNELS = {
    'gg': r'gg',
    'gq': r'gq',
    'qqbar': r'q\bar{q}',
    'qq': r'qq',
    'udbar': r'u\bar{d}',
    'dubar': r'd\bar{u}',
}

QUARK_FLAVORS = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]



def grid_values(pdf:PDF, flmat, xmat, qmat):
     """Returns a 4-dimension array with the PDF values at the input parameters
     for each replica. The return value is indexed as follows:

     grid_values[replica][flavour][x][Q]

     This uses libnnpdf, and therefore follows the convention
     to throw away replica 0 for Monte Carlo ensembles
     (so index 0 corresponds to replica 1). The higher level function
     `central_and_error_grid_values` sets this straight.
     """
     flmat = np.atleast_1d(np.asanyarray(flmat, dtype=FLTYPE))
     xmat, qmat =  (np.atleast_1d(np.asarray(x, dtype=REALTYPE))
                           for x in (xmat,qmat))
     lpdf = pdf.load()
     return lpdf.grid_values(flmat, xmat, qmat)

#TODO: Investigate writting these in cython/cffi/numba/...

def evaluate_luminosity(pdf_set: LHAPDFSet, n: int, s:float, mx: float,
                        x1: float, x2: float, channel):
    """Returns PDF luminosity at specified values of mx, x1, x2, sqrts**2
    for a given channel.

    pdf_set: The interested PDF set
    s: The center of mass energy GeV.
    mx: The invariant mass bin GeV.
    x1 and x2: The the partonic x1 and x2.
    channel: The channel tag name from LUMI_CHANNELS.
    """


    res = 0
    if channel == 'gg':
        res = pdf_set.xfxQ(x1, mx, n, 21) * pdf_set.xfxQ(x2, mx, n, 21)
    elif channel == 'gq':
        for i in QUARK_FLAVORS:
            res += (pdf_set.xfxQ(x1, mx, n, i) * pdf_set.xfxQ(x2, mx, n, 21)
                    + pdf_set.xfxQ(x1, mx, n, 21) * pdf_set.xfxQ(x2, mx, n, i))
    elif channel == 'qqbar':
        for i in QUARK_FLAVORS:
            res += pdf_set.xfxQ(x1, mx, n, i) * pdf_set.xfxQ(x2, mx, n, -i)
    elif channel == 'qq':
        r1 = []
        r2 = []
        for i in QUARK_FLAVORS:
            r1.append(pdf_set.xfxQ(x1, mx, n, i))
            r2.append(pdf_set.xfxQ(x2, mx, n, i))

        res = sum(a*b for a,b in itertools.product(r1,r2))
    elif channel == 'udbar':
        res = (pdf_set.xfxQ(x1, mx, n, 2) * pdf_set.xfxQ(x2, mx, n, -1)
               + pdf_set.xfxQ(x1, mx, n, -1) * pdf_set.xfxQ(x2, mx, n, 2))
    elif channel == 'dubar':
        res = (pdf_set.xfxQ(x1, mx, n, 1) * pdf_set.xfxQ(x2, mx, n, -2)
               + pdf_set.xfxQ(x1, mx, n, -2) * pdf_set.xfxQ(x2, mx, n, 1))

    else:
        raise ValueError("Bad channel")

    return res/x1/x2/s

def uvalence_sum_rule_integrand(x, lpdf:LHAPDFSet, irep, Q):
    return (lpdf.xfxQ(x, Q=Q, n=irep, fl=2) - lpdf.xfxQ(x, Q=Q, n=irep, fl=-2))/x

def dvalence_sum_rule_integrand(x, lpdf:LHAPDFSet, irep, Q):
    return (lpdf.xfxQ(x, Q=Q, n=irep, fl=1) - lpdf.xfxQ(x, Q=Q, n=irep, fl=-1))/x

def svalence_sum_rule_integrand(x, lpdf:LHAPDFSet, irep, Q):
    return (lpdf.xfxQ(x, Q=Q, n=irep, fl=3) - lpdf.xfxQ(x, Q=Q, n=irep, fl=-3))/x

def momentum_sum_rule_integrand(x, lpdf:LHAPDFSet, irep, Q):
    return sum([lpdf.xfxQ(x, Q=10, n=irep, fl=fl) for fl in ALL_FLAVOURS])

SUM_RULES = {
    'uvalence': uvalence_sum_rule_integrand,
    'dvalence': dvalence_sum_rule_integrand,
    'svalence': svalence_sum_rule_integrand,
    'momentum': momentum_sum_rule_integrand,
}

SUM_RULES_EXPECTED= {
    'uvalence': 2,
    'dvalence': 1,
    'svalence': 0,
    'momentum': 1,
}
