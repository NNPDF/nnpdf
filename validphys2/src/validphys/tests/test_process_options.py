"""
Tests for the kinematic maps (x-Q2, z-Q2) that can be attached to a process in ``process_options``
"""

import dataclasses
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from nnpdf_data.process_options import DIS, PROCESSES, _Process, _Vars


def _xq2(kin_info):
    return kin_info[_Vars.x], kin_info[_Vars.Q2]


def _zq2(kin_info):
    return kin_info[_Vars.z], kin_info[_Vars.Q2]


def _metadata(kinematic_coverage):
    """Minimal stand-in for the metadata of an observable"""
    return SimpleNamespace(
        name="FAKE_DATASET",
        plotting=SimpleNamespace(extra_labels=None),
        kinematic_coverage=kinematic_coverage,
        cm_energy=100.0,
    )


def _process(**kwargs):
    return _Process("FAKE", "fake process", accepted_variables=("x", "z", "Q2"), **kwargs)


KIN_DF = pd.DataFrame({"a": [0.1, 0.2], "b": [0.3, 0.4], "c": [4.0, 9.0]})


def test_existing_processes_are_xq2_only():
    """Processes defined before the z-Q2 map existed must keep behaving as before"""
    for name, process in PROCESSES.items():
        if process.xq2map_function is None:
            assert process.kinematic_maps == (), name
        else:
            assert process.kinematic_maps == ("xq2",), name
        assert process.zq2map_function is None, name


def test_positivity_processes_have_no_maps():
    for name in ("POS_XPDF", "POS_DIS", "POS_POL", "POS_DY", "INTEG"):
        assert PROCESSES[name].kinematic_maps == ()


def test_replace_keeps_the_maps():
    assert dataclasses.replace(DIS, name="DIS_OTHER").kinematic_maps == ("xq2",)


@pytest.mark.parametrize(
    "functions, expected",
    [
        ({"xq2map_function": _xq2}, ("xq2",)),
        ({"zq2map_function": _zq2}, ("zq2",)),
        ({"xq2map_function": _xq2, "zq2map_function": _zq2}, ("xq2", "zq2")),
        ({}, ()),
    ],
)
def test_maps_are_inferred_from_the_functions(functions, expected):
    assert _process(**functions).kinematic_maps == expected


def test_explicit_maps():
    process = _process(xq2map_function=_xq2, zq2map_function=_zq2, kinematic_maps=("xq2", "zq2"))
    assert process.kinematic_maps == ("xq2", "zq2")
    # A list is also accepted but it is stored as a tuple
    process = _process(zq2map_function=_zq2, kinematic_maps=["zq2"])
    assert process.kinematic_maps == ("zq2",)


@pytest.mark.parametrize(
    "kwargs",
    [
        # declared but there is no function
        {"xq2map_function": _xq2, "kinematic_maps": ("xq2", "zq2")},
        {"kinematic_maps": ("zq2",)},
        # function but not declared
        {"xq2map_function": _xq2, "zq2map_function": _zq2, "kinematic_maps": ("xq2",)},
        {"zq2map_function": _zq2, "kinematic_maps": ()},
    ],
)
def test_inconsistent_maps_are_rejected(kwargs):
    with pytest.raises(ValueError, match="inconsistent"):
        _process(**kwargs)


def test_unknown_map_is_rejected():
    with pytest.raises(ValueError, match="Unknown kinematic maps"):
        _process(xq2map_function=_xq2, kinematic_maps=("xq2", "yq2"))


def test_adding_a_function_requires_declaring_it():
    with pytest.raises(ValueError, match="inconsistent"):
        dataclasses.replace(DIS, name="DIS_OTHER", zq2map_function=_zq2)
    both = dataclasses.replace(
        DIS, name="DIS_OTHER", zq2map_function=_zq2, kinematic_maps=("xq2", "zq2")
    )
    assert both.kinematic_maps == ("xq2", "zq2")


def test_zq2map_not_implemented():
    with pytest.raises(NotImplementedError, match="zq2map is not implemented for DIS"):
        DIS.zq2map(KIN_DF, _metadata(["x", "Q2", "y"]))


def test_xq2map_not_implemented():
    process = _process(zq2map_function=_zq2)
    with pytest.raises(NotImplementedError, match="xq2map is not implemented for FAKE"):
        process.xq2map(KIN_DF, _metadata(["x", "z", "Q2"]))


def test_both_maps_use_the_same_kinematics():
    process = _process(xq2map_function=_xq2, zq2map_function=_zq2)
    metadata = _metadata(["x", "z", "Q2"])
    x, q2_x = process.xq2map(KIN_DF, metadata)
    z, q2_z = process.zq2map(KIN_DF, metadata)
    np.testing.assert_allclose(x, KIN_DF["a"])
    np.testing.assert_allclose(z, KIN_DF["b"])
    np.testing.assert_allclose(q2_x, KIN_DF["c"])
    np.testing.assert_allclose(q2_z, KIN_DF["c"])


def test_zq2map_checks_the_accepted_variables():
    process = _process(zq2map_function=_zq2)
    with pytest.raises(NotImplementedError, match="kinematic variables are not supported"):
        process.zq2map(KIN_DF, _metadata(["x", "not_a_variable", "Q2"]))


def test_zq2map_missing_variable():
    """A variable required by the map that is not in the kinematics is reported"""
    process = _process(zq2map_function=_zq2)
    with pytest.raises(NotImplementedError, match="Error trying to compute zq2map"):
        process.zq2map(KIN_DF, _metadata(["x", "Q2", "extra_0"]))
