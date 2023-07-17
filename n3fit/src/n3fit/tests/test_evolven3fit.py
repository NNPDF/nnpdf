import logging
import pathlib
import shutil
import subprocess as sp

from evolven3fit_new import eko_utils, utils
import numpy as np
from numpy.testing import assert_allclose
import pytest

from eko import EKO, runner
from reportengine.compat import yaml
from validphys.api import API
from validphys.pdfbases import PIDS_DICT

REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regressions")
log = logging.getLogger(__name__)


def assert_sorted(arr, title):
    """Assert that the array is sorted"""
    if not np.all(np.diff(arr) >= 0):
        raise ValueError(f"The values of {title} are not sorted!")


def check_lhapdf_info(info_path):
    """Check the LHAPDF info file is correct"""
    info = yaml.load(info_path.open("r", encoding="utf-8"))

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
    """Tests the creation of the default grids for different values of nf
    and whether the matched grid is generating points in the desired locations
    """
    # nf 3, q0 = 1.0
    grid = utils.generate_q2grid(None, None, None, {}, 3)
    assert grid[0] == 1.0**2
    # nf 4, q0 = 1.65
    grid = utils.generate_q2grid(None, None, None, {}, 4)
    assert grid[0] == 1.65**2

    for nf in [1, 2, 5, 6]:
        with pytest.raises(NotImplementedError):
            grid = utils.generate_q2grid(None, None, None, {}, nf)

    with pytest.raises(ValueError):
        grid = utils.generate_q2grid(None, None, None, {})

    matched_grid = utils.generate_q2grid(1.65, 1.0e5, 100, {4.92: 2.0, 100: 1.0})
    t1 = 4.92 * 2.0
    t2 = 100.0 * 1.0

    assert_allclose((1.65) ** 2, matched_grid[0])
    assert_allclose((1.0e5) ** 2, matched_grid[-1])
    assert t1**2 in matched_grid
    assert t2**2 in matched_grid


def test_utils():
    # Testing the fake LHAPDF class
    q20 = 1.65**2
    x_grid = np.geomspace(1.0e-7, 1.0, 30)
    fake_grids = [[x * (1.0 - x) for x in x_grid] for _ in PIDS_DICT.keys()]
    pdf_grid = dict([(pid, v) for pid, v in zip(range(len(PIDS_DICT)), fake_grids)])
    my_PDF = utils.LhapdfLike(pdf_grid, q20, x_grid)
    assert my_PDF.hasFlavor(6)
    assert not my_PDF.hasFlavor(0)
    for pid in PIDS_DICT:
        for x in x_grid:
            assert_allclose(my_PDF.xfxQ2(pid, x, q20), x * (1.0 - x))
    # Testing read_runcard
    runcard = utils.read_runcard(REGRESSION_FOLDER)
    assert runcard["description"] == "n3fit regression test"
    assert runcard["datacuts"]["t0pdfset"] == "NNPDF40_nnlo_as_01180"
    # Testing get_theoryID_from_runcard
    ID = utils.get_theoryID_from_runcard(REGRESSION_FOLDER)
    assert ID == 162


def test_eko_utils(tmp_path):
    # Testing construct eko cards
    theoryID = 162
    q_fin = 100
    q_points = 5
    x_grid = [1.0e-3, 0.1, 1.0]
    pto = 2
    comments = "Test"
    t_card, op_card = eko_utils.construct_eko_cards(
        theoryID,
        q_fin,
        q_points,
        x_grid,
        op_card_dict={"configs": {"interpolation_polynomial_degree": 2}},
        theory_card_dict={"Comments": comments},
    )
    t_card_dict = t_card.raw
    op_card_dict = op_card.raw
    assert (
        t_card_dict["order"][0] == pto + 1
    )  # This is due to a different convention in eko orders due to QED
    assert_allclose(op_card_dict["xgrid"], x_grid)
    # In theory 162 the charm threshold is at 1.51
    # and we should find two entries, one for nf=3 and another one for nf=4
    assert_allclose(op_card_dict["mugrid"][0], (1.51, 3))
    assert_allclose(op_card_dict["mugrid"][1], (1.51, 4))
    # Then (with the number of points we chosen it will happen in position 2,3
    # we will find the bottom threshold at two different nf
    assert_allclose(op_card_dict["mugrid"][2], (4.92, 4))
    assert_allclose(op_card_dict["mugrid"][3], (4.92, 5))
    assert_allclose(op_card_dict["mugrid"][-1], (q_fin, 5))
    # Testing computation of eko
    save_path = tmp_path / "ekotest.tar"
    runner.solve(t_card, op_card, save_path)
    eko_op = EKO.read(save_path)
    assert_allclose(eko_op.operator_card.raw["xgrid"], x_grid)
    assert_allclose(list(eko_op.operator_card.raw["mugrid"]), op_card_dict["mugrid"])


@pytest.mark.parametrize("fitname", ["Basic_runcard_3replicas_lowprec_399", "Basic_runcard_qed_3replicas_lowprec_398"])
def test_perform_evolution(tmp_path, fitname):
    """Test that evolven3fit_new is able to utilize the current eko in the respective theory.
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
    sp.run(["evolven3fit_new", "evolve", fitname], cwd=tmp_path, check=True)

    # check that everything worked!
    info = check_lhapdf_info(tmp_info)
    for datpath in tmp_nnfit.glob("replica_*/*.dat"):
        check_lhapdf_dat(datpath, info)

