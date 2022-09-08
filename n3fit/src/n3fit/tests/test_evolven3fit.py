import evolven3fit
from numpy.testing import assert_allclose

def test_utils():
    grid = evolven3fit.utils.generate_q2grid(1.65, None, None)
    assert_allclose(1.65**2, grid[0])