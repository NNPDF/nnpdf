"""
Tests for the kinematic maps (x-Q2, z-Q2) that can be attached to a process in ``process_options``
"""

import dataclasses
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from nnpdf_data.process_options import _HADRON_MASSES, DIS, PROCESSES, _Process, _Vars


def _xq2(kin_info):
    return kin_info[_Vars.x], kin_info[_Vars.Q2]


def _zq2(kin_info):
    return kin_info[_Vars.z], kin_info[_Vars.Q2]


def _metadata(kinematic_coverage, process="FAKE-PROCESS", cm_energy=100.0):
    """Minimal stand-in for the metadata of an observable"""
    return SimpleNamespace(
        name="FAKE_DATASET",
        process=process,
        plotting=SimpleNamespace(extra_labels=None),
        kinematic_coverage=kinematic_coverage,
        cm_energy=cm_energy,
    )


def _process(**kwargs):
    return _Process("FAKE", "fake process", accepted_variables=("x", "z", "Q2"), **kwargs)


KIN_DF = pd.DataFrame({"a": [0.1, 0.2], "b": [0.3, 0.4], "c": [4.0, 9.0]})

# Processes that were added together with the z-Q2 map
FRAGMENTATION_PROCESSES = ("SIDIS", "SIA")


def test_existing_processes_are_xq2_only():
    """Processes defined before the z-Q2 map existed must keep behaving as before"""
    for name, process in PROCESSES.items():
        if name in FRAGMENTATION_PROCESSES:
            continue
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


def test_sidis_has_both_maps():
    assert PROCESSES["SIDIS"].kinematic_maps == ("xq2", "zq2")


@pytest.mark.parametrize("scale", ["Q2", "Q"])
def test_sidis_maps(scale):
    """Both the maps can be computed from x, z and the scale (either Q or Q2)"""
    sidis = PROCESSES["SIDIS"]
    q2 = KIN_DF["c"]
    kin_df = KIN_DF.assign(c=q2 if scale == "Q2" else np.sqrt(q2))
    metadata = _metadata(["x", "z", scale])
    x, q2_x = sidis.xq2map(kin_df, metadata)
    z, q2_z = sidis.zq2map(kin_df, metadata)
    np.testing.assert_allclose(x, KIN_DF["a"])
    np.testing.assert_allclose(z, KIN_DF["b"])
    np.testing.assert_allclose(q2_x, q2)
    np.testing.assert_allclose(q2_z, q2)


def test_sidis_xq2map_is_the_dis_one():
    x_sidis, q2_sidis = PROCESSES["SIDIS"].xq2map(KIN_DF, _metadata(["x", "z", "Q2"]))
    x_dis, q2_dis = DIS.xq2map(KIN_DF, _metadata(["x", "y", "Q2"]))
    np.testing.assert_allclose(x_sidis, x_dis)
    np.testing.assert_allclose(q2_sidis, q2_dis)


def test_sidis_maps_are_independent():
    """Each map only needs the variables it uses, the other one fails if they are missing"""
    sidis = PROCESSES["SIDIS"]
    # No z: only the x-Q2 map can be computed
    metadata = _metadata(["x", "y", "Q2"])
    sidis.xq2map(KIN_DF, metadata)
    with pytest.raises(NotImplementedError, match="Error trying to compute zq2map"):
        sidis.zq2map(KIN_DF, metadata)
    # No x: only the z-Q2 map can be computed
    metadata = _metadata(["z", "y", "Q2"])
    sidis.zq2map(KIN_DF, metadata)
    with pytest.raises(NotImplementedError, match="Error trying to compute xq2map"):
        sidis.xq2map(KIN_DF, metadata)


# The variables in which SIA can be measured, besides z: the momentum of the hadron,
# given in different ways, from which z is computed using the mass of the hadron
SIA_MOMENTUM_VARIABLES = ("xp", "ph", "xi")


def _sia_zq2map(variable, values, hadron="PI", cm_energy=10.0):
    """Apply the z-Q2 map of the SIA process to the given values of one of its variables"""
    kin_df = pd.DataFrame({variable: values, "pad_1": 0.0, "pad_2": 0.0})
    metadata = _metadata([variable, "extra_0", "extra_1"], f"SIA-{hadron}", cm_energy)
    return PROCESSES["SIA"].zq2map(kin_df, metadata)


def test_sia_has_only_zq2_map():
    assert PROCESSES["SIA"].kinematic_maps == ("zq2",)


def test_sia_z_map():
    """When measured as a function of z the map is the identity, at Q2 = sqrts^2"""
    z_values = np.array([0.1, 0.5, 0.9])
    # The hadron is not needed in this case
    z, q2 = _sia_zq2map("z", z_values, hadron="UNKNOWN")
    np.testing.assert_allclose(z, z_values)
    np.testing.assert_allclose(q2, 10.0**2)
    assert q2.shape == z.shape


@pytest.mark.parametrize("hadron", list(_HADRON_MASSES))
def test_sia_momentum_maps(hadron):
    """z follows from the momentum of the hadron, which can be given as xp, ph or xi.
    The three variables must agree with the same closed form for the same points."""
    mass, sqrts = _HADRON_MASSES[hadron], 10.0
    xp = np.array([0.05, 0.3, 0.9])
    expected_z = np.sqrt(xp**2 + 4 * mass**2 / sqrts**2)
    variables = {"xp": xp, "ph": xp * sqrts / 2, "xi": np.log(1 / xp)}
    for variable in SIA_MOMENTUM_VARIABLES:
        z, q2 = _sia_zq2map(variable, variables[variable], hadron, sqrts)
        np.testing.assert_allclose(z, expected_z, err_msg=variable)
        np.testing.assert_allclose(q2, sqrts**2, err_msg=variable)
        assert q2.shape == z.shape


def test_sia_sqrts_as_a_variable():
    """The scale is taken from sqrts when it is one of the kinematic variables"""
    kin_df = pd.DataFrame({"z": [0.1, 0.4], "sqrts": [10.0, 20.0], "pad": 0.0})
    metadata = _metadata(["z", "sqrts", "extra_0"], "SIA-PI", cm_energy=None)
    _, q2 = PROCESSES["SIA"].zq2map(kin_df, metadata)
    np.testing.assert_allclose(q2, [100.0, 400.0])


@pytest.mark.parametrize("variable", SIA_MOMENTUM_VARIABLES)
def test_sia_unknown_hadron(variable):
    with pytest.raises(NotImplementedError, match="Error trying to compute zq2map") as error:
        _sia_zq2map(variable, [0.1], hadron="UNKNOWN")
    assert "Unknown hadron" in str(error.value.__cause__)


@pytest.mark.parametrize("variable", ["z", *SIA_MOMENTUM_VARIABLES])
def test_sia_needs_the_center_of_mass_energy(variable):
    with pytest.raises(NotImplementedError, match="Error trying to compute zq2map") as error:
        _sia_zq2map(variable, [0.1], cm_energy=None)
    assert "center of mass energy" in str(error.value.__cause__)


@pytest.mark.parametrize("variables", [[], ["xp", "ph"], ["z", "xp", "xi"]])
def test_sia_needs_exactly_one_variable(variables):
    """The observable is given as a function of only one of z, xp, ph and xi"""
    coverage = [*variables, *(f"extra_{i}" for i in range(3 - len(variables)))]
    kin_df = pd.DataFrame(np.full((2, 3), 0.1))
    with pytest.raises(NotImplementedError, match="Error trying to compute zq2map") as error:
        PROCESSES["SIA"].zq2map(kin_df, _metadata(coverage, "SIA-PI"))
    assert "Exactly one of" in str(error.value.__cause__)
