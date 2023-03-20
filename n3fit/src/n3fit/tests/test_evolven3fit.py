import pathlib
import logging
from numpy.testing import assert_allclose
import numpy as np
from validphys.pdfbases import PIDS_DICT
from evolven3fit_new import utils, eko_utils
from eko import EKO

REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regressions")
log = logging.getLogger(__name__)


def check_consecutive_members(grid, value):
    """Check if the first occurrence of value in grid is followed by value again"""
    return np.allclose(grid[list(grid).index(value) + 1], value)


def test_utils():
    # Testing the default grid
    grid = utils.generate_q2grid(1.65, None, None, {})
    assert_allclose(1.65**2, grid[0])
    assert len(grid) == 50
    # We expect the bottom mass to be repeated twice because it is intended once in 4 flavors and once in 5 flavors.
    assert check_consecutive_members(grid, 4.92**2)
    # Testing if the points of the matching are correctly repeated twice
    matched_grid = utils.generate_q2grid(1.65, 1.0e5, 100, {4.92: 2.0, 100: 1.0})
    assert len(matched_grid) == 100
    assert_allclose((1.0e5) ** 2, matched_grid[-1])
    assert check_consecutive_members(matched_grid, (4.92 * 2.0) ** 2)
    assert check_consecutive_members(matched_grid, (100.0 * 1.0) ** 2)
    # Testing the fake LHAPDF class
    q20 = 1.65**2
    x_grid = np.geomspace(1.0e-7, 1.0, 30)
    fake_grids = [[x * (1.0 - x) for x in x_grid] for pid in PIDS_DICT.keys()]
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
    assert_allclose(op_card_dict["rotations"]["xgrid"], x_grid)
    assert_allclose(op_card_dict["_mugrid"][0], 1.65)
    assert_allclose(op_card_dict["_mugrid"][-1], q_fin)
    # In this case there are not enough points to have twice the bottom matching scale
    assert_allclose(op_card_dict["_mugrid"][1], 4.92)
    # Testing construct_eko_for_fit
    save_path = tmp_path / "ekotest.tar"
    eko_utils.construct_eko_for_fit(t_card, op_card, save_path=save_path)
    eko_op = EKO.read(save_path)
    assert_allclose(eko_op.operator_card.raw["rotations"]["xgrid"], x_grid)
    assert_allclose(list(eko_op.operator_card.raw["_mugrid"]), op_card_dict["_mugrid"])
