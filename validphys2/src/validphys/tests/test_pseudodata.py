"""
Test to ensure the validphys.pseudodata.get_pseudodata action
correctly obtains the appropriate pseudodata for an n3fit fit.

A fit has been generated called pseudodata_test_fit_n3fit
which has the pseudodata saved as training and validation splits.
This is used to benchmark the correctness of the pseudodata
recreation.
"""
import pandas as pd
import pytest

from validphys.api import API
from validphys.tests.conftest import FIT, PSEUDODATA_FIT


def test_read_fit_pseudodata():
    fit_pseudodata = API.read_fit_pseudodata(fit=PSEUDODATA_FIT)

    nrep = API.num_fitted_replicas(fit=PSEUDODATA_FIT)
    assert nrep == len(fit_pseudodata)

    for data, tr_idx, val_idx in fit_pseudodata:
        assert set(tr_idx).isdisjoint(set(val_idx))
        assert set(tr_idx).union(val_idx) == set(data.index)


def test_read_pdf_pseudodata():
    pdf_pseudodata = API.read_pdf_pseudodata(fit=PSEUDODATA_FIT)

    pdf = API.pdf(pdf=PSEUDODATA_FIT)
    # -1 because we ignore replica 0
    assert len(pdf) - 1 == len(pdf_pseudodata)

    for data, tr_idx, val_idx in pdf_pseudodata:
        assert set(tr_idx).isdisjoint(set(val_idx))
        assert set(tr_idx).union(val_idx) == set(data.index)


def test_recreate_fit_pseudodata():
    config = {}
    config["use_t0_sampling"] = False
    config["separate_multiplicative"]= True
    config["output_path"] = None
    config["theory_covmat_flag"] = False
    config["use_user_uncertainties"] = None 
    config["use_scalevar_uncertainties"] = None
    config["use_thcovmat_in_sampling"]=None
    config["dataset_t0_predictions"] = None
    
    fit_pseudodata = API.recreate_fit_pseudodata(fit=PSEUDODATA_FIT, **config)
    nrep = API.num_fitted_replicas(fit=PSEUDODATA_FIT)
    assert nrep == len(fit_pseudodata)

    for data, tr_idx, val_idx in fit_pseudodata:
        assert set(tr_idx).isdisjoint(set(val_idx))
        assert set(tr_idx).union(val_idx) == set(data.index)


def test_recreate_pdf_pseudodata():
    config = {}
    config["use_t0_sampling"] = False
    config["separate_multiplicative"]= True
    config["output_path"] = None
    config["theory_covmat_flag"] = False
    config["use_user_uncertainties"] = None 
    config["use_scalevar_uncertainties"] = None
    config["use_thcovmat_in_sampling"]=None
    config["dataset_t0_predictions"] = None

    pdf_pseudodata = API.recreate_pdf_pseudodata(fit=PSEUDODATA_FIT, **config)

    pdf = API.pdf(pdf=PSEUDODATA_FIT)
    # -1 because we ignore replica 0
    assert len(pdf) - 1 == len(pdf_pseudodata)

    for data, tr_idx, val_idx in pdf_pseudodata:
        assert set(tr_idx).isdisjoint(set(val_idx))
        assert set(tr_idx).union(val_idx) == set(data.index)


def test_no_savepseudodata():
    for func in (API.read_fit_pseudodata, API.read_pdf_pseudodata):
        with pytest.raises(FileNotFoundError):
            # Check a FileNotFoundError is raised
            # if the input fit wasn't generated
            # with the savepseudodata flag set to true
            func(fit=FIT)


def test_read_matches_recreate():
    config = {}
    config["use_t0_sampling"] = False
    config["separate_multiplicative"]= True
    config["output_path"] = None
    config["theory_covmat_flag"] = False
    config["use_user_uncertainties"] = None 
    config["use_scalevar_uncertainties"] = None
    config["use_thcovmat_in_sampling"]=None
    config["dataset_t0_predictions"] = None

    reads = API.read_fit_pseudodata(fit=PSEUDODATA_FIT)
    recreates = API.recreate_fit_pseudodata(fit=PSEUDODATA_FIT, **config)
    for read, recreate in zip(reads, recreates):
        # We ignore the absolute ordering of the dataframes and just check
        # that they contain identical elements.
        pd.testing.assert_frame_equal(
            read.pseudodata, recreate.pseudodata, check_like=True
        )
        pd.testing.assert_index_equal(read.tr_idx, recreate.tr_idx, check_order=False)
        pd.testing.assert_index_equal(read.val_idx, recreate.val_idx, check_order=False)
