"""
    Test hyperoptimization features
"""

import numpy as np
from numpy.testing import assert_approx_equal
from n3fit.hyper_optimization.rewards import HyperLoss

def test_rewards():
    """ Ensure that rewards continue doing what they are supposed to do """
    losses = np.array([[0.0, 2.0], [1.0, 3.0], [2.0, 4.0]])
    loss_average = HyperLoss(replica_statistic="average", fold_statistic="average")
    assert_approx_equal(loss_average.compute(losses), 2.0)

    loss_best_worst = HyperLoss(replica_statistic="average", fold_statistic="best_worst")
    assert_approx_equal(loss_best_worst.compute(losses), 3.0)

    loss_std = HyperLoss(replica_statistic="average", fold_statistic="std")
    assert_approx_equal(loss_std.compute(losses), 0.816496580927726)

    loss_best_worst_best_worst = HyperLoss(replica_statistic="best_worst", fold_statistic="best_worst")
    assert_approx_equal(loss_best_worst_best_worst.compute(losses), 4.0)
