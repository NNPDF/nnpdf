"""
sumrules.py

Module for the computation of sum rules

Note that this contains only the code for the computation of sum rules from
scratch using LHAPDF tables. The code reading the sum rule information output
from the fit is present in fitinfo.py
"""
import numbers
from collections import namedtuple

import numpy as np
import pandas as pd
import scipy.integrate as integrate

from NNPDF import LHAPDFSet
from reportengine.table import table
from reportengine.checks import check_positive

from validphys.core import PDF
from validphys.pdfbases import ALL_FLAVOURS


def _uvalence_sum_rule_integrand(x, lpdf:LHAPDFSet, irep, Q):
    return (lpdf.xfxQ(x, Q=Q, n=irep, fl=2) - lpdf.xfxQ(x, Q=Q, n=irep, fl=-2))/x

def _dvalence_sum_rule_integrand(x, lpdf:LHAPDFSet, irep, Q):
    return (lpdf.xfxQ(x, Q=Q, n=irep, fl=1) - lpdf.xfxQ(x, Q=Q, n=irep, fl=-1))/x

def _svalence_sum_rule_integrand(x, lpdf:LHAPDFSet, irep, Q):
    return (lpdf.xfxQ(x, Q=Q, n=irep, fl=3) - lpdf.xfxQ(x, Q=Q, n=irep, fl=-3))/x

def _momentum_sum_rule_integrand(x, lpdf:LHAPDFSet, irep, Q):
    return sum([lpdf.xfxQ(x, Q=Q, n=irep, fl=fl) for fl in ALL_FLAVOURS])


#NOTE: For the moment we rely on this order being the same as in the .sumrules
#file produced by nnfit.
SUM_RULES = {
    'momentum': _momentum_sum_rule_integrand,
    'uvalence': _uvalence_sum_rule_integrand,
    'dvalence': _dvalence_sum_rule_integrand,
    'svalence': _svalence_sum_rule_integrand,
}

SUM_RULES_EXPECTED = {
    'momentum': 1,
    'uvalence': 2,
    'dvalence': 1,
    'svalence': 0,
}

#Output result tuple
SumRulesGrid = namedtuple('SumRulesGrid', SUM_RULES)
def _sum_rules(lpdf, Q):
    """Compute a SumRulesGrid from the loaded PDF, at Q"""
    nmembers = lpdf.GetMembers()
    #TODO: Write this in something fast
    #If nothing else, at least allocate and store the result contiguously
    res = np.zeros((len(SUM_RULES), nmembers))
    integrands = SUM_RULES.values()
    def integral(f, a, b, irep):
        #We increase the limit to capture the log scale fluctuations
        return integrate.quad(f, a, b, args=(lpdf, irep, Q),
               limit=1000,
               epsabs=1e-4, epsrel=1e-4)[0]

    for irep in range(nmembers):
        for r, f in enumerate(integrands):
            res[r,irep] =  (integral(f, 0, 1e-5, irep)  +
               integral(f, 1e-5, 1e-3, irep) + integral(f, 1e-3, 1, irep))
    return SumRulesGrid(*res)


@check_positive('Q')
def sum_rules(pdf:PDF, Q:numbers.Real):
    """Compute the sum rules for each member (as defined by libnnpdf), at the
    energy scale ``Q``. Return a SumRulesGrid object with the list of valued "
    "for each sum rule.
    The integration is performed with absolute and relative tolerance of 1e-4."""
    lpdf = pdf.load()
    return _sum_rules(lpdf, Q)

@check_positive('Q')
def central_sum_rules(pdf:PDF, Q:numbers.Real):
    """Compute the sum rules for the central member, at the scale Q"""
    lpdf = pdf.load_t0()
    return _sum_rules(lpdf, Q)


@table
def sum_rules_table(sum_rules):
    """Return a table with the descriptive statistics of the sum rules,
    over members of the PDF."""
    #We don't  really want the count, which is going to be the same for all.
    #Hence the .iloc[1:,:].
    return pd.DataFrame(sum_rules._asdict()).describe().iloc[1:,:]

@table
def central_sum_rules_table(central_sum_rules):
    """Construct a table with the value of each sum rule for the central
    member"""
    return pd.DataFrame(central_sum_rules._asdict())


@table
def bad_replica_sumrules(pdf, sum_rules, threshold:numbers.Real=0.01):
    """Return a table with the sum rules for the replica where some sum rule is
    farther from the correct value than ``threshold`` (in absolute value).
    """
    ncomputed = len(sum_rules[0])
    #TODO: Get rid of this nonsense
    if pdf.ErrorType == 'replicas':
        x = np.arange(1, ncomputed + 1)
    else:
        x = np.arange(ncomputed)
    df = pd.DataFrame(sum_rules._asdict(), index=x)
    filt = ((df - pd.Series(SUM_RULES_EXPECTED)).abs() > threshold).any(axis=1)
    return df[filt]
