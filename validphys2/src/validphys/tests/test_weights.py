"""
test_weights.py
"""
import numpy as np

from validphys.api import API

def test_weights_have_same_commondata(weighted_data_witht0_config):
    exps = API.experiments(**weighted_data_witht0_config)
    normal, weighted = exps
    normalds, weightedds = normal.datasets[0].load(), weighted.datasets[0].load()
    assert normalds.GetSys(0, 0).mult == weightedds.GetSys(0, 0).mult
    assert normalds.GetSys(0, 0).add == weightedds.GetSys(0, 0).add


def test_chi2_arithmetic(weighted_data_witht0_config):
    normal, weighted = API.experiments_chi2(**weighted_data_witht0_config)
    assert np.allclose(weighted[0].data/normal[0].data, 100)
