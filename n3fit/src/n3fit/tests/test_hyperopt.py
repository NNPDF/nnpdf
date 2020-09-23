"""
    Test hyperoptimization features
"""

from numpy.testing import assert_approx_equal
from n3fit.hyper_optimization import rewards

def test_rewards():
    """ Ensure that rewards continue doing what they are supposed to do """
    losses = [0.0, 1.0, 2.0]
    assert_approx_equal(rewards.average(losses), 1.0)
    assert_approx_equal(rewards.best_worst(losses), 2.0)
    assert_approx_equal(rewards.std(losses), 0.816496580927726)
