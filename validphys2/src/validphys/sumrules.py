"""
sumrules.py

Module for the computation of sum rules

Note that this contains only the code for the computation of sum rules from
scratch using LHAPDF tables. The code reading the sum rule information output
from the fit is present in fitinfo.py
"""
import numbers

import numpy as np
import pandas as pd
from scipy.integrate import quad

from reportengine.checks import check_positive
from reportengine.floatformatting import format_error_value_columns
from reportengine.table import table
from validphys.core import PDF
from validphys.pdfbases import parse_flarr


def _momentum_sum_rule_integrand(x, lpdf, Q):
    xqvals = lpdf.xfxQ(x, Q)
    return sum([xqvals[f] for f in lpdf.flavors()])


def _make_momentum_fraction_integrand(fldict):
    """Make a suitable integrand function, which takes x to be integrated over
    and a PDF member and Q that computes the momentum fraction based on ``fldict``.

    The keys of ``fldict`` are free form values corresponding to PDG parton ids
    (that end up being passed by :py:func:`validphys.pdfbases.parse_flarr` and
    then to LHAPDF) and the values are multipliers for each parton. The
    integrand is the sum of ``x*flavour(x)*multiplier`` for all the given
    entries.

    Parameters
    ----------
    fldict : Mapping[int, int]
        A map from PDG parton id to multipliers

    Returns
    -------
    f : Callable
        An integrand function.
    """
    # Do this outside to aid integration time
    fldict = {parse_flarr([k])[0]: v for k, v in fldict.items()}

    def f(x, lpdf, Q):
        return sum(multiplier * lpdf.xfxQ(x, Q)[flavour] for flavour, multiplier in fldict.items())

    return f


def _make_pdf_integrand(fldict):
    """Make a suitable integrand function, which takes x to be integrated over
    and a PDF member and Q that computes the integrand of the PDFs based on ``fldict``.

    The keys of ``fldict`` are free form values corresponfing to PDG parton ids
    (that end up being passed :py:func:`validphys.pdfbases.parse_flarr` and
    then to LHAPDF) and the values are multipliers for each parton. The
    integrand is the sum of ``x*flavour(x)*multiplier`` for all the given
    entries.

    Parameters
    ----------
    fldict : Mapping[int, int]
        A map from PDG parton id to multipliers

    Returns
    -------
    f : Callable
        An integrand function.
    """
    # Do this outsde to aid integration time
    fldict = {parse_flarr([k])[0]: v for k, v in fldict.items()}

    def f(x, lpdf, Q):
        return (
            sum(multiplier * lpdf.xfxQ(x, Q)[flavour] for flavour, multiplier in fldict.items()) / x
        )

    return f


KNOWN_SUM_RULES = {
    "momentum": _momentum_sum_rule_integrand,
    "uvalence": _make_pdf_integrand({"u": 1, "ubar": -1}),
    "dvalence": _make_pdf_integrand({"d": 1, "dbar": -1}),
    "svalence": _make_pdf_integrand({"s": 1, "sbar": -1}),
    "cvalence": _make_pdf_integrand({"c": 1, "cbar": -1}),
}

UNKNOWN_SUM_RULES = {
    "u momentum fraction": _make_momentum_fraction_integrand({"u": 1}),
    "ubar momentum fraction": _make_momentum_fraction_integrand({"ubar": 1}),
    "d momentum fraction": _make_momentum_fraction_integrand({"d": 1}),
    "dbar momentum fraction": _make_momentum_fraction_integrand({"dbar": 1}),
    "s momentum fraction": _make_momentum_fraction_integrand({"s": 1}),
    "sbar momentum fraction": _make_momentum_fraction_integrand({"sbar": 1}),
    "cp momentum fraction": _make_momentum_fraction_integrand({"c": 1, "cbar": 1}),
    "cm momentum fraction": _make_momentum_fraction_integrand({"c": 1, "cbar": -1}),
    "g momentum fraction": _make_momentum_fraction_integrand({"g": 1}),
    "T3": _make_pdf_integrand({"u": 1, "ubar": 1, "d": -1, "dbar": -1}),
    "T8": _make_pdf_integrand({"u": 1, "ubar": 1, "d": 1, "dbar": 1, "s": -2, "sbar": -2}),
}

POLARIZED_SUM_RULES = {
    "momentum": _momentum_sum_rule_integrand,
    "T3": _make_pdf_integrand({"u": 1, "ubar": 1, "d": -1, "dbar": -1}),
    "T8": _make_pdf_integrand({"u": 1, "ubar": 1, "d": 1, "dbar": 1, "s": -2, "sbar": -2}),
}


KNOWN_SUM_RULES_EXPECTED = {
    'momentum': 1,
    'uvalence': 2,
    'dvalence': 1,
    'svalence': 0,
    'cvalence': 0,
}


def _integral(rule_f, pdf_member, Q, config=None):
    """Integrate `rule_f` for a given `pdf_member` at a given energy
    separating the regions of integration. Uses quad.
    """
    if config is None:
        config = {"limit": 1000, "epsabs": 1e-4, "epsrel": 1e-4}
    res = 0.0
    lims = [(1e-9, 1e-5), (1e-5, 1e-3), (1e-3, 1)]
    for lim in lims:
        res += quad(rule_f, *lim, args=(pdf_member, Q), **config)[0]
    return res


def _sum_rules(rules_dict, lpdf, Q):
    """Compute a SumRulesGrid from the loaded PDF, at Q"""
    return {k: [_integral(r, m, Q) for m in lpdf.members] for k, r in rules_dict.items()}


@check_positive('Q')
def sum_rules(pdf: PDF, Q: numbers.Real):
    """Compute the momentum, uvalence, dvalence, svalence and cvalence sum rules for
    each member, at the energy scale ``Q``.
    Return a SumRulesGrid object with the list of values for each sum rule.
    The integration is performed with absolute and relative tolerance of 1e-4."""
    lpdf = pdf.load()
    return _sum_rules(KNOWN_SUM_RULES, lpdf, Q)


@check_positive('Q')
def polarized_sum_rules(pdf: PDF, Q: numbers.Real):
    """Compute the polarized sum rules. Return a SumRulesGrid object with the list of
    values for each sum rule. The integration is performed with absolute and relative
    tolerance of 1e-4."""
    lpdf = pdf.load()
    return _sum_rules(POLARIZED_SUM_RULES, lpdf, Q)


@check_positive('Q')
def central_sum_rules(pdf: PDF, Q: numbers.Real):
    """Compute the sum rules for the central member, at the scale Q"""
    lpdf = pdf.load_t0()
    return _sum_rules(KNOWN_SUM_RULES, lpdf, Q)


@check_positive('Q')
def unknown_sum_rules(pdf: PDF, Q: numbers.Real):
    """Compute the following integrals
    - u momentum fraction
    - ubar momentum fraction
    - d momentum fraction
    - dbar momentum fraction
    - s momentum fraction
    - sbar momentum fraction
    - cp momentum fraction
    - cm momentum fraction
    - g momentum fraction
    - T3
    - T8
    """
    lpdf = pdf.load()
    return _sum_rules(UNKNOWN_SUM_RULES, lpdf, Q)


def _simple_description(d):
    res = {}
    for k, arr in d.items():
        res[k] = d = {}
        d["mean"] = np.mean(arr)
        d["std"] = np.std(arr)
        d["min"] = np.min(arr)
        d["max"] = np.max(arr)

    return pd.DataFrame(res).T


def _err_mean_table(d, polarized=False):
    res = {}
    for k, arr in d.items():
        res[k] = d = {}
        d["mean"] = np.mean(arr)
        d["std"] = np.std(arr)
        if polarized:
            d["min"] = np.min(arr)
            d["max"] = np.max(arr)
    df = pd.DataFrame(res)
    return format_error_value_columns(df.T, "mean", "std")


@table
def sum_rules_table(sum_rules):
    """Return a table with the descriptive statistics of the sum rules,
    over members of the PDF."""
    return _simple_description(sum_rules)


@table
def polarized_sum_rules_table(polarized_sum_rules):
    """Return a table with the descriptive statistics of the polarized sum rules,
    over members of the PDF."""
    return _err_mean_table(polarized_sum_rules, polarized=True)


@table
def central_sum_rules_table(central_sum_rules):
    """Construct a table with the value of each sum rule for the central
    member"""
    return pd.DataFrame(central_sum_rules, index=["Central value"]).T


@table
def unknown_sum_rules_table(unknown_sum_rules):
    return _err_mean_table(unknown_sum_rules)


@table
def bad_replica_sumrules(pdf, sum_rules, threshold: numbers.Real = 0.01):
    """Return a table with the sum rules for the replica where some sum rule is
    farther from the correct value than ``threshold`` (in absolute value).
    """
    ncomputed = len(sum_rules[0])
    if pdf.error_type == "replicas":
        x = np.arange(1, ncomputed + 1)
    else:
        x = np.arange(ncomputed)
    df = pd.DataFrame(sum_rules._asdict(), index=x)
    filt = ((df - pd.Series(KNOWN_SUM_RULES_EXPECTED)).abs() > threshold).any(axis=1)
    return df[filt]
