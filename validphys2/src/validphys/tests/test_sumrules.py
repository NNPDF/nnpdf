"""
    Test for the sum rules
"""
import pandas as pd
import pytest

from validphys.api import API
from validphys.tableloader import sane_load

from .conftest import HESSIAN_PDF, PDF
from .test_regressions import make_table_comp

Q = 10  # GeV


def parse_sumrules(filename):
    """Parse the regression"""
    return sane_load(filename, header=0)


@pytest.mark.parametrize("pdf_name", [PDF, HESSIAN_PDF])
def test_central_value(pdf_name):
    """Check that the central value is the same as the mean"""
    all_sr = API.sum_rules_table(pdf=pdf_name, Q=Q)
    cv_sr = API.central_sum_rules_table(pdf=pdf_name, Q=Q)
    pd.testing.assert_series_equal(cv_sr.squeeze(), all_sr["mean"], atol=1e-5, check_names=False)


def _regression_sum_rules(pdf_name):
    known_sumrules = API.sum_rules_table(pdf=pdf_name, Q=Q)
    central_val = API.central_sum_rules_table(pdf=pdf_name, Q=Q)
    # First concatenate the central value to the sum rules
    ret = pd.concat([known_sumrules, central_val], axis=1)
    # The unknown sum rules need to be treated a bit because
    # 1) They don't contain min or max
    # 2) They are of type object?
    unk_sumrules = API.unknown_sum_rules_table(pdf=pdf_name, Q=Q)
    for col in ret.columns:
        if col in unk_sumrules.columns:
            unk_sumrules[col] = unk_sumrules[col].astype(float)
        if col not in unk_sumrules.columns:
            unk_sumrules[col] = -1.0
    # all together now:
    return pd.concat([ret, unk_sumrules])


@make_table_comp(parse_sumrules, tolerance=1e-5)
def test_sum_rules_MC():
    return _regression_sum_rules(PDF)


@make_table_comp(parse_sumrules, tolerance=1e-5)
def test_sum_rules_hessian():
    return _regression_sum_rules(HESSIAN_PDF)
