"""
test_weights.py
"""
import numpy as np

from validphys.api import API

def test_weights_have_same_commondata(weighted_data_witht0_config):
    data = API.data(**weighted_data_witht0_config)
    normal, weighted = data.datasets
    normalds, weightedds = normal.load_commondata(), weighted.load_commondata()
    assert (
        normalds.systematics_table["MULT"].iloc[0][0]
        == weightedds.systematics_table["MULT"].iloc[0][0]
    )
    assert (
        normalds.systematics_table["ADD"].iloc[0][0]
        == weightedds.systematics_table["ADD"].iloc[0][0]
    )


def test_chi2_arithmetic(weighted_data_witht0_internal_cuts_config):
    ((normal, weighted,),) = API.groups_datasets_chi2_data(
        **weighted_data_witht0_internal_cuts_config
    )
    assert np.allclose(weighted[0].data / normal[0].data, 100)


def test_disable_weights(weighted_data_witht0_internal_cuts_config):
    ((normal, weighted,),) = API.groups_datasets_chi2_data(
        **weighted_data_witht0_internal_cuts_config, use_weights_in_covmat=False
    )
    assert np.allclose(weighted[0].data / normal[0].data, 1)
    unweighted = API.groups_chi2(
        **weighted_data_witht0_internal_cuts_config, use_weights_in_covmat=False
    )[0].central_result
    weighted = API.groups_chi2(
        **weighted_data_witht0_internal_cuts_config, use_weights_in_covmat=True
    )[0].central_result
    assert np.allclose(weighted / unweighted, (100 + 1) / (1 + 1))

def test_python_weights(weighted_data_witht0_config):
    """Test python implementation of weighted covmats is constent with
    libnnpdf and that ``use_weights_in_covmat`` is working correctly in
    python interface.

    """
    weighted_data_witht0_config = dict(weighted_data_witht0_config)
    weighted_data_witht0_config["use_cuts"] = "internal"
    cov = API.dataset_inputs_covariance_matrix(
        **weighted_data_witht0_config
    )
    py_cov = API.dataset_inputs_t0_covmat_from_systematics(
        **weighted_data_witht0_config
    )

    np.testing.assert_allclose(cov, py_cov, rtol=1e-05, atol=1e-08)

    # now test without weights - assumes that libnnpdf tests pass.
    unweighted = API.dataset_inputs_covariance_matrix(
        **weighted_data_witht0_config, use_weights_in_covmat=False,
    )
    # use t0 here
    py_unweighted = API.dataset_inputs_t0_covmat_from_systematics(
        **weighted_data_witht0_config, use_weights_in_covmat=False,
    )
    np.testing.assert_allclose(py_unweighted, unweighted, rtol=1e-05, atol=1e-08)
