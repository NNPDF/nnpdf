"""
test_totalchi2.py

test that the action which calculates the total chi2 produces sensible
results for both MC and hessian pdfs
"""
import numpy as np

from validphys.api import API


def test_hessian_total_chi2(hessian_data_internal_cuts_config):
    """testing total chi2 for hessian pdf

    In particular check that the sum across experiments is handled correctly

    and that calculating the total chi2 from the flat list of datasets gives
    the same answer as using ``total_chi2_data``

    """
    member_chi2, cent_chi2, ndata = API.total_chi2_data(**hessian_data_internal_cuts_config)

    # this is only true for hessian PDF
    np.testing.assert_allclose(member_chi2.central_value(), cent_chi2)

    exps_chi2_data = API.experiments_chi2_data(**hessian_data_internal_cuts_config)
    exps_member_chi2, exps_cent_chi2, exps_ndata = list(zip(*exps_chi2_data))

    assert np.sum(exps_ndata) == ndata

    np.testing.assert_allclose(np.sum(exps_cent_chi2), cent_chi2)

    exps_chi2_error_mem = [stats_obj.error_members() for stats_obj in exps_member_chi2]
    np.testing.assert_allclose(
        np.sum(exps_chi2_error_mem, axis=0), member_chi2.error_members()
    )

    dsinp_mem_chi2, dsinp_cent_chi2, dsinp_ndata = API.dataset_inputs_abs_chi2_data(
        **hessian_data_internal_cuts_config
    )

    np.testing.assert_allclose(dsinp_mem_chi2.data, member_chi2.data)

    np.testing.assert_allclose(dsinp_cent_chi2, cent_chi2)

    assert dsinp_ndata == ndata


def test_mc_total_chi2(data_internal_cuts_config):
    """Testing total chi2 for mc pdf

    In particular check that the sum across experiments is handled correctly

    and that calculating the total chi2 from the flat list of datasets gives
    the same answer as using ``total_chi2_data``

    """
    member_chi2, cent_chi2, ndata = API.total_chi2_data(**data_internal_cuts_config)

    exps_chi2_data = API.experiments_chi2_data(**data_internal_cuts_config)
    exps_member_chi2, exps_cent_chi2, exps_ndata = list(zip(*exps_chi2_data))

    assert np.sum(exps_ndata) == ndata

    np.testing.assert_allclose(np.sum(exps_cent_chi2), cent_chi2)

    exps_chi2_error_mem = [stats_obj.error_members() for stats_obj in exps_member_chi2]
    np.testing.assert_allclose(
        np.sum(exps_chi2_error_mem, axis=0), member_chi2.error_members()
    )

    dsinp_mem_chi2, dsinp_cent_chi2, dsinp_ndata = API.dataset_inputs_abs_chi2_data(
        **data_internal_cuts_config
    )

    np.testing.assert_allclose(dsinp_mem_chi2.data, member_chi2.data)

    np.testing.assert_allclose(dsinp_cent_chi2, cent_chi2)

    assert dsinp_ndata == ndata
