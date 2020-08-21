import pytest

from validphys.loader import FallbackLoader as Loader
from validphys.filters import (
    Rule,
    RuleProcessingError,
    default_filter_settings_input,
    default_filter_rules_input,
    PerturbativeOrder,
    BadPerturbativeOrder,
)
from validphys.tests.conftest import THEORYID

bad_rules = [
    {'dataset': 'NMC'},
    {'rule': 'x < 0.1'},
    {'dataset': 'NOT_EXISTING', 'rule': 'x < 0.1'},
    {'dataset': 'NMC', 'rule': 'x < 0.1', 'PTO': 'bogus'},
    {'dataset': 'NMC', 'rule': 'x < 0.1', 'PTO': {'bog': 'us'}},
    {'dataset': 'NMC', 'rule': 'x < 0.1', 'local_variables': 'bogus'},
    {'dataset': 'NMC', 'rule': 'bogus syntax'},
    {'dataset': 'NMC', 'rule': 'unknown_variable > 10'},
    {'dataset': 'NMC', 'local_variables': {'z': 'bogus syntax'}, 'rule': 'z > 10'},
    {
        'dataset': 'NMC',
        'local_variables': {'z': 'unknown_variable + 1'},
        'rule': 'z > 10',
    },
    {'dataset': 'NMC', 'local_variables': {'z': 'v+1', 'v': '10'}, 'rule': 'z > 10'},
]

# Note: Don't change the order here. In this way it tests all cases.
good_rules = [
    {'process_type': 'DIS_ALL', 'PTO': 'N3LO', 'rule': 'x < 1e-2'},
    {'process_type': 'DIS_ALL', 'IC': 'False', 'rule': 'x < 1e-2'},
    {'process_type': 'JET', 'rule': 'p_T2 < 10'},
]


def mkrule(inp):
    l = Loader()
    th = l.check_theoryID(THEORYID)
    desc = th.get_description()
    defaults = default_filter_settings_input()
    return Rule(initial_data=inp, defaults=defaults, theory_parameters=desc)


def test_PTO():
    assert 2 in PerturbativeOrder("NNLO")
    assert 2 in PerturbativeOrder("N2LO")
    assert 2 not in PerturbativeOrder("NNLO!")
    assert 2 in PerturbativeOrder("NNLO+")
    assert 2 not in PerturbativeOrder("NNLO-")
    with pytest.raises(BadPerturbativeOrder):
        PerturbativeOrder("NBogus+")


def test_bad_rules():
    for rule_inp in bad_rules:
        with pytest.raises(RuleProcessingError):
            mkrule(rule_inp)


def test_default_rules():
    l = Loader()
    dsnames = ['NMC', 'LHCBWZMU8TEV']
    for dsname in dsnames:
        ds = l.check_dataset(dsname, cuts='internal', theoryid=THEORYID)
        assert ds.cuts.load() is not None


def test_good_rules():
    l = Loader()
    rules = [mkrule(inp) for inp in good_rules]
    dsnames = ['ATLAS1JET11', 'NMC']
    for dsname in dsnames:
        ds = l.check_dataset(dsname, cuts='internal', rules=rules, theoryid=THEORYID)
        assert ds.cuts.load() is not None
