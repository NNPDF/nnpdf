import logging
import pathlib
import shutil
import subprocess as sp

from evolven3fit import eko_utils, utils
import numpy as np
import pytest

from eko import EKO, runner
from validphys.api import API
from validphys.utils import yaml_safe

REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regressions")
log = logging.getLogger(__name__)


def assert_sorted(arr, title):
    """Assert that the array is sorted"""
    if not np.all(np.diff(arr) >= 0):
        raise ValueError(f"The values of {title} are not sorted!")


def check_lhapdf_info(info_path):
    """Check the LHAPDF info file is correct"""
    info = yaml_safe.load(info_path.open("r", encoding="utf-8"))

    alphas_qs = info["AlphaS_Qs"]
    alphas = info["AlphaS_Vals"]

    np.testing.assert_equal(info["QMin"], alphas_qs[0])
    np.testing.assert_equal(len(alphas), len(alphas_qs))
    assert isinstance(info["Particle"], int)
    np.testing.assert_allclose(info["AlphaS_OrderQCD"], 2)  # Test theories are NNLO
    assert info["ErrorType"] == "replicas"

    # Check that the values of Q are sorted
    assert_sorted(alphas_qs, "Qs in the .info file")

    for flavor in info["Flavors"]:
        assert isinstance(flavor, int)

    assert info["NumFlavors"] <= 6

    return info


def check_lhapdf_dat(dat_path, info):
    """Check that the dat file itself makes sense
    and that it is consistent with the info file
    """
    blocks = dat_path.read_text().split("---")[:-1]

    for i, block in enumerate(blocks[1:]):
        sp_block = block.strip().split("\n", 3)
        x = np.fromstring(sp_block[0], sep=" ")
        q = np.fromstring(sp_block[1], sep=" ")
        flavs = np.fromstring(sp_block[2], dtype=int, sep=" ")

        np.testing.assert_array_equal(info["Flavors"], flavs)

        np.testing.assert_equal(x[0], info["XMin"])
        np.testing.assert_equal(x[-1], info["XMax"])
        assert_sorted(x, "x-grid")

        if i == 0:
            np.testing.assert_equal(q[0], info["QMin"])
        assert_sorted(q, "q-grid")

    # Use allclose here to avoid failing because of a different in the 7th place
    np.testing.assert_allclose(q[-1], info["QMax"])


def test_generate_q2grid():
    """Tests the creation of the default grids is as expected"""

    # Test if the correct errors are given
    with pytest.raises(TypeError):
        grid = utils.generate_q2grid(None, None, None, {}, None, None)

    with pytest.raises(ValueError):
        grid = utils.generate_q2grid(
            1, 1, 1, {"mc": 1.502, "mb": 4.936, "kcThr": 1, "kbThr": 1}, 1, 1
        )

    # Test if the grid contains the threshold values and the boundaries
    matched_grid = utils.generate_q2grid(
        1.65,
        1.0,
        100,
        {"mc": 1.502, "mb": 4.936, "kcThr": 1, "kbThr": 1},
        total_points=10,
        total_points_ic=3,
    )

    for n in {1.502, 1.65, 4.936}:
        assert any(np.allclose(n**2, x) for x in matched_grid)

    np.testing.assert_allclose((1.0) ** 2, matched_grid[0])
    np.testing.assert_allclose((100.0) ** 2, matched_grid[-1])

    # Test the legacy40 grid
    legacy40 = utils.generate_q2grid(1, 1, 1, {}, 1, 1, legacy40=True)
    np.testing.assert_allclose((1.65) ** 2, legacy40[0])
    np.testing.assert_allclose((1e5) ** 2, legacy40[-1])


def test_utils():
    # Testing read_runcard
    runcard = utils.read_runcard(REGRESSION_FOLDER)
    assert isinstance(runcard["description"], str)
    assert runcard["datacuts"]["t0pdfset"] == "NNPDF40_nnlo_as_01180"
    # Testing get_theoryID_from_runcard
    ID = utils.get_theoryID_from_runcard(REGRESSION_FOLDER)
    assert ID == 40_000_000


def test_eko_utils(tmp_path, nnpdf_theory_card):
    # Testing construct eko cards
    q_ini = 1.0
    q_fin = 1e5
    x_grid = [1.0e-3, 0.1, 1.0]
    pto = 2
    comments = "Test"
    t_card, op_card = eko_utils.construct_eko_cards(
        nnpdf_theory_card,
        x_grid,
        op_card_dict={"configs": {"interpolation_polynomial_degree": 2}},
        theory_card_dict={"Comments": comments},
    )
    t_card_dict = t_card.raw
    op_card_dict = op_card.raw
    assert (
        t_card_dict["order"][0] == pto + 1
    )  # This is due to a different convention in eko orders due to QED
    print(op_card_dict["mugrid"])
    np.testing.assert_allclose(op_card_dict["xgrid"], x_grid)
    # We should find two entries for each threshold energy,
    # one for nf=3(4) and another one for nf=4(5)
    for n in {(1.51, 3), (1.51, 4), (4.92, 4), (4.92, 5)}:
        assert any(np.allclose(n, x) for x in op_card_dict["mugrid"])

    # Testing if the endpoints are correct
    np.testing.assert_allclose(op_card_dict["mugrid"][0], (q_ini, 3))
    np.testing.assert_allclose(op_card_dict["mugrid"][-1], (q_fin, 5))

    # Testing computation of eko
    save_path = tmp_path / "ekotest.tar"
    runner.solve(t_card, op_card, save_path)
    eko_op = EKO.read(save_path)
    np.testing.assert_allclose(eko_op.operator_card.raw["xgrid"], x_grid)
    np.testing.assert_allclose(list(eko_op.operator_card.raw["mugrid"]), op_card_dict["mugrid"])


# Keep old theory 399 and 398 EKOs here to avoid having to download many GB of data
# the only _actual_ eko test is given by the fit-bot
@pytest.mark.parametrize(
    "fitname", ["Basic_runcard_3replicas_lowprec_399", "Basic_runcard_qed_3replicas_lowprec_398"]
)
def test_perform_evolution(tmp_path, fitname):
    """Test that evolven3fit is able to utilize the current eko in the respective theory.
    In addition checks that the generated .info files are correct
    """
    fit = API.fit(fit=fitname)
    _ = API.theoryid(theoryid=fit.as_input()['theory']['theoryid'])
    # Move the fit to a temporary folder
    tmp_fit = tmp_path / fitname
    shutil.copytree(fit.path, tmp_fit)
    # Clear the .log and .dat files
    (tmp_fit / "evolven3fit_new.log").unlink()
    tmp_nnfit = tmp_fit / "nnfit"
    tmp_info = tmp_nnfit / f"{fitname}.info"
    tmp_info.unlink()
    for datpath in tmp_nnfit.glob("replica_*/*.dat"):
        datpath.unlink()

    # And re-evolve the fit
    sp.run(["evolven3fit", "evolve", fitname], cwd=tmp_path, check=True)

    # check that everything worked!
    info = check_lhapdf_info(tmp_info)
    for datpath in tmp_nnfit.glob("replica_*/*.dat"):
        check_lhapdf_dat(datpath, info)
