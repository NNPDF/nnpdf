"""
    Test hyperoptimization features
"""

from numpy.testing import assert_approx_equal
from n3fit.hyper_optimization.rewards import HyperLoss

def test_rewards():
    """ Ensure that rewards continue doing what they are supposed to do """
    losses = [0.0, 1.0, 2.0]
    loss_average = HyperLoss("average")
    assert_approx_equal(loss_average.compute(losses), 1.0)

    loss_best_worst = HyperLoss("best_worst")
    assert_approx_equal(loss_best_worst.compute(losses), 2.0)

    loss_std = HyperLoss("std")
    assert_approx_equal(loss_std.compute(losses), 0.816496580927726)
