"""
test_weights.py
"""
import numpy as np

def test_weights_have_same_commondata(weighted_data):
    _, exps = weighted_data
    normal, weighted = exps
    normalds, weightedds = normal.datasets[0].load(), weighted.datasets[0].load()
    assert normalds.GetSys(0, 0).mult == weightedds.GetSys(0, 0).mult
    assert normalds.GetSys(0, 0).add == weightedds.GetSys(0, 0).add


def test_chi2_arithmetic(weighted_chi2data):
    normal, weighted = weighted_chi2data
    assert np.allclose(weighted[0].data/normal[0].data, 100)
