import numpy as np
import pytest

from validphys.api import API
from validphys.filters import (
    BadPerturbativeOrder,
    PerturbativeOrder,
    Rule,
    RuleProcessingError,
    default_filter_settings_input,
)
from validphys.loader import FallbackLoader as Loader
from validphys.tests.conftest import PDF
from validphys.tests.conftest import THEORYID as THEORYID

bad_rules = [
    {'dataset': "NMC_NC_NOTFIXED_DW_EM-F2"},
    {'rule': 'x < 0.1'},
    {'dataset': 'NOT_EXISTING', 'rule': 'x < 0.1'},
    {'dataset': "NMC_NC_NOTFIXED_DW_EM-F2", 'rule': 'x < 0.1', 'PTO': 'bogus'},
    {'dataset': "NMC_NC_NOTFIXED_DW_EM-F2", 'rule': 'x < 0.1', 'PTO': {'bog': 'us'}},
    {'dataset': "NMC_NC_NOTFIXED_DW_EM-F2", 'rule': 'x < 0.1', 'local_variables': 'bogus'},
    {'dataset': "NMC_NC_NOTFIXED_DW_EM-F2", 'rule': 'bogus syntax'},
    {'dataset': "NMC_NC_NOTFIXED_DW_EM-F2", 'rule': 'unknown_variable > 10'},
    {
        'dataset': "NMC_NC_NOTFIXED_DW_EM-F2",
        'local_variables': {'z': 'bogus syntax'},
        'rule': 'z > 10',
    },
    {
        'dataset': "NMC_NC_NOTFIXED_DW_EM-F2",
        'local_variables': {'z': 'unknown_variable + 1'},
        'rule': 'z > 10',
    },
    {
        'dataset': "NMC_NC_NOTFIXED_DW_EM-F2",
        'local_variables': {'z': 'v+1', 'v': '10'},
        'rule': 'z > 10',
    },
]

# Note: Don't change the order here. In this way it tests all cases.
good_rules = [
    {'process_type': 'DIS_ALL', 'PTO': 'N3LO', 'rule': 'x < 1e-2'},
    {'process_type': 'DIS_ALL', 'IC': 'False', 'rule': 'x < 1e-2'},
    {'process_type': 'JET', 'rule': 'pT < 3.16'},
]


def mkrule(inp):
    l = Loader()
    th = l.check_theoryID(THEORYID)
    desc = th.get_description()
    defaults = default_filter_settings_input()
    return Rule(initial_data=inp, defaults=defaults, theory_parameters=desc)


def test_rule_caching():
    rule_list_1, *rule_list_2 = good_rules
    rule_list_1 = [rule_list_1]

    cut_list = []
    for rule_list in (rule_list_1, rule_list_2):
        cut_list.append(
            API.cuts(
                dataset_input={"dataset": "NMC_NC_NOTFIXED_DW_EM-F2"},
                use_cuts="internal",
                theoryid=THEORYID,
                filter_rules=rule_list,
            )
        )
    assert not cut_list[0] == cut_list[1]


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
    dsnames = ['NMC_NC_NOTFIXED_EM-F2', 'LHCB_Z0_8TEV_MUON_Y']
    variants = ["legacy_dw", None]
    for dsname, v in zip(dsnames, variants):
        ds = l.check_dataset(dsname, cuts='internal', theoryid=THEORYID, variant=v)
        assert ds.cuts.load() is not None


def test_good_rules():
    l = Loader()
    rules = [mkrule(inp) for inp in good_rules]
    dsnames = ['ATLAS_1JET_8TEV_R06_PTY', 'NMC_NC_NOTFIXED_EM-F2']
    variants = ["legacy", "legacy_dw"]
    for dsname, v in zip(dsnames, variants):
        ds = l.check_dataset(
            dsname, cuts='internal', rules=tuple(rules), theoryid=THEORYID, variant=v
        )
        assert ds.cuts.load() is not None


def test_added_rules():
    inp = {
        "theoryid": THEORYID,
        "pdf": PDF,
        "use_cuts": "internal",
        "dataset_inputs": [{"dataset": "ATLAS_1JET_8TEV_R06_PTY", "variant": "legacy"}],
        "filter_rules": [],
        "dataspecs": [
            {"speclabel": "Original", "added_filter_rules": None},
            {
                "speclabel": "fewer data",
                "added_filter_rules": [
                    {"dataset": "ATLAS_1JET_8TEV_R06_PTY", "rule": "pT < 1000", "reason": "pt cut"}
                ],
            },
            {
                "speclabel": "empty data",
                "added_filter_rules": [
                    {"dataset": "ATLAS_1JET_8TEV_R06_PTY", "rule": "y < 0", "reason": "empty data"}
                ],
            },
        ],
    }
    tb = API.dataspecs_chi2_table(**inp)
    assert tb["empty data"]["ndata"].iloc[0] == 0
    assert np.isnan(tb["empty data"].iloc[1, 1])
    assert tb["empty data"]["ndata"].iloc[0] == 0
    assert np.all(tb[1:]["fewer data"] != tb[1:]["Original"])


def test_drop_internal_rules(data_internal_cuts_config, test_dataset="CMS_Z0J_8TEV_PT-Y"):
    """Check that the key drop_internal_rules work as expected:
    - Drops all cuts for a given dataset
    - It is applied before added_filter_rules
    """
    assert test_dataset in [
        i["dataset"] for i in data_internal_cuts_config["dataset_inputs"]
    ], "If you updated the test DATA, please update this test as well"

    def test_fun(**config):
        """Use some internal validphy function which will for sure use cuts and separate
        the results for the test dataset.
        """
        # Get data and predictions separated by dataset (drop grouping)
        ret = API.group_result_central_table_no_table(**config).droplevel(0)
        # Now separate the test dataset from the rest
        df_test = ret.loc[test_dataset]
        df_rest = ret.drop(index=test_dataset)
        return df_test, df_rest

    # Use internal cuts
    def_test, def_all = test_fun(**data_internal_cuts_config)

    # Drop all rules for the test dataset only
    drop_test, drop_all = test_fun(**data_internal_cuts_config, drop_internal_rules=[test_dataset])

    assert len(drop_test) > len(def_test), "Cuts have not been dropped!"
    assert len(drop_all) == len(def_all), "Drop cuts have affected other datasets!"

    # Add a new rule for this dataset while dropping all previous rules
    new_rule = {"dataset": test_dataset, "rule": "pT >= 80"}
    add_test, add_all = test_fun(
        **data_internal_cuts_config,
        added_filter_rules=[new_rule],
        drop_internal_rules=[test_dataset]
    )
    assert len(new_rule) < len(drop_test), "New rule has not been added after dropping the cuts!"
