from evolven3fit import utils
from numpy.testing import assert_allclose
import numpy as np
import pytest 

pids_dict = {
        -6: "TBAR",
        -5: "BBAR",
        -4: "CBAR",
        -3: "SBAR",
        -2: "UBAR",
        -1: "DBAR",
        21: "GLUON",
        1: "D",
        2: "U",
        3: "S",
        4: "C",
        5: "B",
        6: "T",
        22: "PHT",
    }

pids_order = [
        "TBAR",
        "BBAR",
        "CBAR",
        "SBAR",
        "UBAR",
        "DBAR",
        "GLUON",
        "D",
        "U",
        "S",
        "C",
        "B",
        "T",
        "PHT",
    ]

def test_utils():
    grid = utils.generate_q2grid(1.65, None, None, {})
    assert_allclose(1.65**2, grid[0])
    assert len(grid) == 50
    assert_allclose(4.92**2, grid[11])
    assert_allclose(4.92**2, grid[12])
    matched_grid = utils.generate_q2grid(1.65, 1.0e5, 100, {4.92:2.0, 100:1.0} )
    assert len(matched_grid) ==  100
    assert_allclose((1.0e5)**2, matched_grid[-1])
    assert_allclose(matched_grid[list(matched_grid).index((4.92*2.0)**2)+1], (4.92*2.0)**2)
    assert_allclose(matched_grid[list(matched_grid).index((100*1.0)**2)+1], (100*1.0)**2)
    q20 = 1.65**2
    x_grid = np.geomspace(1.0e-7, 1.0, 30)
    pdf_grid = dict([(pid,v) for pid,v in zip(range(len(pids_order)), [[x*(1.-x) for x in x_grid] for pid in pids_order] )])
    my_PDF = utils.LhapdfLike(pdf_grid, q20, x_grid)
    assert my_PDF.hasFlavor(6)
    assert not my_PDF.hasFlavor(0)
    with pytest.raises(ValueError):
        my_PDF.xfxQ2(21, 0.01, 5.0**2)
    for pid in pids_dict.keys():
        for x in x_grid:
            assert_allclose(my_PDF.xfxQ2(pid, x, q20), x*(1.-x) )
