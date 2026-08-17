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

# Limits of the partial integration when computing (Sum) Rules
LIMS = [(1e-9, 1e-5), (1e-5, 1e-3), (1e-3, 1)]
POL_LIMS = ((1e-4, 1e-3), (1e-3, 1))

# Mellin moments x**n * f(x) are dominated by large-x. Small-x follows LIMS (one slice per
# order of magnitude); the extra large-x slices are phenomenological, not numerical: they
# show which large-x region drives the moment, as needed for lattice comparisons.
MELLIN_LIMS = [
    (1e-9, 1e-5),
    (1e-5, 1e-3),
    (1e-3, 1e-1),
    (1e-1, 0.3),
    (0.3, 0.6),
    (0.6, 0.9),
    (0.9, 1.0),
]


def _momentum_sum_rule_integrand(x, lpdf, Q):
    xqvals = lpdf.xfxQ(x, Q)
    return sum([xqvals[f] for f in lpdf.flavors()])


def _make_mellin_moment_integrand(fldict, n: int):
    """Integrand for the n-th Mellin moment ``int x**n f(x) dx`` of the flavour combination
    ``fldict`` (PDG parton ids -> multipliers, see :py:func:`validphys.pdfbases.parse_flarr`).
    LHAPDF returns ``x*f(x)``, hence the ``x**(n-1)``. n=0 gives the PDF integral (valence
    sum rules), n=1 the momentum fraction.
    """
    # Do this outside to aid integration time
    fldict = {parse_flarr([k])[0]: v for k, v in fldict.items()}
    exponent = n - 1

    def f(x, lpdf, Q):
        xfvals = lpdf.xfxQ(x, Q)
        return x**exponent * sum(multiplier * xfvals[fl] for fl, multiplier in fldict.items())

    return f


KNOWN_SUM_RULES = {
    "momentum": _momentum_sum_rule_integrand,
    "uvalence": _make_mellin_moment_integrand({"u": 1, "ubar": -1}, n=0),
    "dvalence": _make_mellin_moment_integrand({"d": 1, "dbar": -1}, n=0),
    "svalence": _make_mellin_moment_integrand({"s": 1, "sbar": -1}, n=0),
    "cvalence": _make_mellin_moment_integrand({"c": 1, "cbar": -1}, n=0),
}

UNKNOWN_SUM_RULES = {
    "u momentum fraction": _make_mellin_moment_integrand({"u": 1}, n=1),
    "ubar momentum fraction": _make_mellin_moment_integrand({"ubar": 1}, n=1),
    "d momentum fraction": _make_mellin_moment_integrand({"d": 1}, n=1),
    "dbar momentum fraction": _make_mellin_moment_integrand({"dbar": 1}, n=1),
    "s momentum fraction": _make_mellin_moment_integrand({"s": 1}, n=1),
    "sbar momentum fraction": _make_mellin_moment_integrand({"sbar": 1}, n=1),
    "c momentum fraction": _make_mellin_moment_integrand({"c": 1}, n=1),
    "cbar momentum fraction": _make_mellin_moment_integrand({"cbar": 1}, n=1),
    "g momentum fraction": _make_mellin_moment_integrand({"g": 1}, n=1),
    "T3": _make_mellin_moment_integrand({"u": 1, "ubar": 1, "d": -1, "dbar": -1}, n=0),
    "T8": _make_mellin_moment_integrand(
        {"u": 1, "ubar": 1, "d": 1, "dbar": 1, "s": -2, "sbar": -2}, n=0
    ),
}

POLARIZED_SUM_RULES = {
    "singlet": _make_mellin_moment_integrand(
        {'u': 1, 'ubar': 1, 'd': 1, 'dbar': 1, 's': 1, 'sbar': 1}, n=0
    ),
    "g": _make_mellin_moment_integrand({"g": 1}, n=0),
    "momentum": _momentum_sum_rule_integrand,
    "T3": _make_mellin_moment_integrand({"u": 1, "ubar": 1, "d": -1, "dbar": -1}, n=0),
    "T8": _make_mellin_moment_integrand(
        {"u": 1, "ubar": 1, "d": 1, "dbar": 1, "s": -2, "sbar": -2}, n=0
    ),
    "xV": _make_mellin_moment_integrand({'u': 1, 'ubar': -1, 'd': 1, 'dbar': -1}, n=1),
    "xV3": _make_mellin_moment_integrand({'u': 1, 'ubar': -1, 'd': -1, 'dbar': 1}, n=1),
}

MELLIN_MOMENTS = {
    "m_5": _make_mellin_moment_integrand({"u": 1, "ubar": -1, "d": -1, "dbar": 1}, n=4),
    "m_4": _make_mellin_moment_integrand({"u": 1, "ubar": -1, "d": -1, "dbar": 1}, n=3),
    "m_3": _make_mellin_moment_integrand({"u": 1, "ubar": -1, "d": -1, "dbar": 1}, n=2),
    "m_2": _make_mellin_moment_integrand({"u": 1, "ubar": -1, "d": -1, "dbar": 1}, n=1),
    "gluon m_2": _make_mellin_moment_integrand({"g": 1}, n=1),
    "gluon m_3": _make_mellin_moment_integrand({"g": 1}, n=2),
    "gluon m_4": _make_mellin_moment_integrand({"g": 1}, n=3),
    "gluon m_5": _make_mellin_moment_integrand({"g": 1}, n=4),
}


KNOWN_SUM_RULES_EXPECTED = {
    'momentum': 1,
    'uvalence': 2,
    'dvalence': 1,
    'svalence': 0,
    'cvalence': 0,
}


def _integral(rule_f, pdf_member, Q, lim, config=None):
    """Integrate `rule_f` for a given `pdf_member` at a given energy
    for a given region of integration. Uses quad.
    """
    if config is None:
        config = {"limit": 1000, "epsabs": 1e-4, "epsrel": 1e-4}
    return quad(rule_f, *lim, args=(pdf_member, Q), **config)[0]


def _sum_rules(rules_dict, lpdf, Q, lims=LIMS):
    """Compute a SumRulesGrid from the loaded PDF, at Q"""
    return [
        {k: [_integral(r, m, Q, lim=l) for m in lpdf.members] for k, r in rules_dict.items()}
        for l in lims
    ]


def _combine_limits(res: list[dict]):
    """Sum the various limits together for all SR and return a dictionary."""
    return {k: np.sum([v[k] for v in res], axis=0) for k in res[0].keys()}


@check_positive('Q')
def partial_polarized_sum_rules(pdf: PDF, Q: numbers.Real, lims: tuple = POL_LIMS):
    """Compute the partial low- and large-x polarized sum rules. Return a SumRulesGrid
    object with the list of values for each sum rule. The integration is performed with
    absolute and relative tolerance of 1e-4."""
    lpdf = pdf.load()
    return _sum_rules(POLARIZED_SUM_RULES, lpdf, Q, lims=lims)


@check_positive('Q')
def sum_rules(pdf: PDF, Q: numbers.Real):
    """Compute the momentum, uvalence, dvalence, svalence and cvalence sum rules for
    each member, at the energy scale ``Q``.
    Return a SumRulesGrid object with the list of values for each sum rule.
    The integration is performed with absolute and relative tolerance of 1e-4."""
    lpdf = pdf.load()
    return _combine_limits(_sum_rules(KNOWN_SUM_RULES, lpdf, Q))


@check_positive('Q')
def polarized_sum_rules(partial_polarized_sum_rules):
    """Compute the full polarized sum rules. The integration is performed with absolute
    and relative tolerance of 1e-4."""
    return _combine_limits(partial_polarized_sum_rules)


@check_positive('Q')
def central_sum_rules(pdf: PDF, Q: numbers.Real):
    """Compute the sum rules for the central member, at the scale Q"""
    lpdf = pdf.load_t0()
    return _combine_limits(_sum_rules(KNOWN_SUM_RULES, lpdf, Q))


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
    return _combine_limits(_sum_rules(UNKNOWN_SUM_RULES, lpdf, Q))


@check_positive('Q')
def partial_mellin_moments(pdf: PDF, Q: numbers.Real, lims: list = MELLIN_LIMS):
    """Per-interval contributions to the Mellin moments defined in
    ``MELLIN_MOMENTS``. Returns a list (one entry per element of ``lims``) of dicts
    ``{moment_name: [value_per_member]}``, showing which x region drives each moment.
    """
    lpdf = pdf.load()
    return _sum_rules(MELLIN_MOMENTS, lpdf, Q, lims=lims)


@check_positive('Q')
def mellin_moments(partial_mellin_moments):
    """Sum the partial Mellin-moment integrals over all x sub-intervals."""
    return _combine_limits(partial_mellin_moments)


def _simple_description(data, pdf):
    """Return a table with simple descriptive statistics over the members of the PDF.
    The statistics used depend on the type of PDF (MC or Hessian)."""
    res = {}
    stats_class = pdf.stats_class
    for k, arr in data.items():
        res[k] = stats_dict = {}
        # Each entry in `data` is expected to be a vector with shape
        # (central_member + error_members,), hence we reshape to (-1, 1) before
        # passing to the stats class.
        stats = stats_class(arr.reshape(-1, 1))
        stats_dict["mean"] = stats.central_value().item()
        stats_dict["std"] = stats.std_error().item()
        stats_dict["min"] = np.min(arr)
        stats_dict["max"] = np.max(arr)

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
    df = df[["T3", "T8", "singlet", "g", "xV", "xV3"]] if polarized else df
    return format_error_value_columns(df.T, "mean", "std")


@table
def sum_rules_table(sum_rules, pdf):
    """Return a table with the descriptive statistics of the sum rules,
    over members of the PDF."""
    return _simple_description(sum_rules, pdf)


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
def mellin_moments_table(mellin_moments, pdf):
    """Mean ± std (over PDF members) of the higher Mellin moments."""
    return _simple_description(mellin_moments, pdf)


@table
def partial_mellin_moments_table(partial_mellin_moments, lims: list = MELLIN_LIMS):
    """Per-interval Mellin-moment contributions, averaged over PDF members.
    One row per moment, one column per ``lims`` interval.
    """
    rows = {}
    for interval_dict, lim in zip(partial_mellin_moments, lims):
        col = f"({lim[0]:.0e}, {lim[1]:.0e})"
        rows[col] = {k: float(np.mean(v)) for k, v in interval_dict.items()}
    return pd.DataFrame(rows)


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
