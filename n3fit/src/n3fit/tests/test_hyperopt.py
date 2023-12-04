"""
    Test hyperoptimization features
"""

import numpy as np
from numpy.testing import assert_approx_equal

from n3fit.hyper_optimization.rewards import HyperLoss


def test_rewards():
    """Ensure that rewards continue doing what they are supposed to do"""
    losses = np.array([1.0, 2.0, 3.0])
    loss_average = HyperLoss(fold_statistic="average")
    assert_approx_equal(loss_average.reduce_over_folds(losses), 2.0)

    loss_std = HyperLoss(fold_statistic="std")
    assert_approx_equal(loss_std.reduce_over_folds(losses), 0.816496580927726)

    loss_best_worst_best_worst = HyperLoss(fold_statistic="best_worst")
    assert_approx_equal(loss_best_worst_best_worst.reduce_over_folds(losses), 3.0)
