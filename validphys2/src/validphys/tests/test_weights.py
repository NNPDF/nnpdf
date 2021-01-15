"""
test_weights.py
"""
import numpy as np

from validphys.api import API

def test_weights_have_same_commondata(weighted_data_witht0_config):
    data = API.data(**weighted_data_witht0_config)
    normal, weighted = data.datasets
    normalds, weightedds = normal.load(), weighted.load()
    assert normalds.GetSys(0, 0).mult == weightedds.GetSys(0, 0).mult
    assert normalds.GetSys(0, 0).add == weightedds.GetSys(0, 0).add


def test_chi2_arithmetic(weighted_data_witht0_config):
    ((normal, weighted,),) = API.groups_datasets_chi2_data(
        **weighted_data_witht0_config
    )
    assert np.allclose(weighted[0].data / normal[0].data, 100)


def test_disable_weights(weighted_data_witht0_config):
    ((normal, weighted,),) = API.groups_datasets_chi2_data(
        **weighted_data_witht0_config, use_weights_in_covmat=False
    )
    assert np.allclose(weighted[0].data / normal[0].data, 1)
    unweighted = API.groups_chi2(
        **weighted_data_witht0_config, use_weights_in_covmat=False
    )[0].central_result
    weighted = API.groups_chi2(
        **weighted_data_witht0_config, use_weights_in_covmat=True
    )[0].central_result
    assert np.allclose(weighted / unweighted, (100 + 1) / (1 + 1))
