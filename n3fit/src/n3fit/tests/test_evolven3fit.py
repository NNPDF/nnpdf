from evolven3fit import utils
from numpy.testing import assert_allclose
import numpy as np
import pytest 
from validphys.pdfbases import PIDS_DICT


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
    pdf_grid = dict([(pid,v) for pid,v in zip(range(len(PIDS_DICT)), [[x*(1.-x) for x in x_grid] for pid in PIDS_DICT.keys()] )])
    my_PDF = utils.LhapdfLike(pdf_grid, q20, x_grid)
    assert my_PDF.hasFlavor(6)
    assert not my_PDF.hasFlavor(0)
    with pytest.raises(ValueError):
        my_PDF.xfxQ2(21, 0.01, 5.0**2)
    for pid in PIDS_DICT.keys():
        for x in x_grid:
            assert_allclose(my_PDF.xfxQ2(pid, x, q20), x*(1.-x) )
