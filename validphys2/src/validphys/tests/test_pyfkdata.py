import pytest
import pandas as pd
import numpy as np
from numpy.testing import assert_allclose

from validphys.api import API
from validphys.loader import Loader
from validphys.results import ThPredictionsResult, PositivityResult
from validphys.fkparser import load_fktable
from validphys.convolution import predictions, central_predictions, linear_predictions
from validphys.tests.conftest import PDF, HESSIAN_PDF, THEORYID, POSITIVITIES


def test_basic_loading():
    l = Loader()
    # Test both with and without cfactors, and load both DIS and hadronic
    for cfac in ((), ("QCD",)):
        fk = l.check_fktable(setname="ATLASTTBARTOT", theoryID=THEORYID, cfac=cfac)
        res = load_fktable(fk)
        assert res.ndata == 3
        assert isinstance(res.sigma, pd.DataFrame)
    fk = l.check_fktable(setname="H1HERAF2B", theoryID=THEORYID, cfac=())
    res = load_fktable(fk)
    assert res.ndata == 12
    assert isinstance(res.sigma, pd.DataFrame)

    # Check if cfactors for datasets having one entry are correctly parsed
    fk = l.check_fktable(setname="CMSTTBARTOT7TEV", theoryID=THEORYID, cfac=("QCD",))
    res = load_fktable(fk)
    assert res.ndata == 1


def test_cuts():
    l = Loader()
    ds = l.check_dataset("ATLASTTBARTOT", theoryid=THEORYID, cfac=("QCD",))
    table = load_fktable(ds.fkspecs[0])
    # Check explicit cuts
    newtable = table.with_cuts([0, 1])
    assert set(newtable.sigma.index.get_level_values(0)) == {0, 1}
    assert newtable.ndata == 2
    assert newtable.metadata["GridInfo"].ndata == ds.commondata.ndata
    # Check empty cuts
    assert newtable.with_cuts(None) is newtable
    # Check loaded cuts
    ds = l.check_dataset("H1HERAF2B", theoryid=THEORYID)
    table = load_fktable(ds.fkspecs[0])
    newtable = table.with_cuts(ds.cuts)
    assert len(newtable.sigma.index.get_level_values(0).unique()) == len(ds.cuts.load())


@pytest.mark.parametrize("pdf_name", [PDF, HESSIAN_PDF])
def test_predictions(pdf_name):
    """Test that the ThPredictionsResult class do not break the raw predictions
    coming from the convolution module and that they are compatible with the API result"""
    l = Loader()
    pdf = l.check_pdf(pdf_name)
    datasets = [
        {"name": "ATLASTTBARTOT", "cfac": ("QCD",)},  # cfactors
        {"name": "H1HERAF2B"},  # DIS, op: NULL
        {"name": "D0ZRAP"},  # op: RATIO
        {"name": "D0WEASY"},  # op: ASY
        {"name": "CMSWCHARMTOT"},  # op: ADD
        {"name": "ATLASWPT31PB"},  # op: SMN
        {"name": "DYE906R"},  # op: COM <----
        {"name": "DYE906_D"},  # op: SMT <----
    ]
    for daset in datasets:
        ds = l.check_dataset(**daset, theoryid=THEORYID)
        preds = predictions(ds, pdf)
        core_predictions = ThPredictionsResult.from_convolution(pdf, ds)
        assert_allclose(preds.values, core_predictions._rawdata, atol=1e-8, rtol=1e-3)
        # Now check that the stats class does the right thing with the data
        cv_predictions = central_predictions(ds, pdf).squeeze()
        stats_predictions = pdf.stats_class(preds.T)
        # rtol to 1e-2 due to DYE906R and DYE906_D for MC sets
        # TODO: check whether this tolerance can be decreased when using double precision
        assert_allclose(cv_predictions, stats_predictions.central_value(), rtol=1e-2)


@pytest.mark.parametrize("pdf_name", [PDF, HESSIAN_PDF])
def test_positivity(pdf_name):
    """Test that the PositivityResult is sensible and like test_predictions
    that no internal step modifies the result and that they are compatible
    with the API result
    """
    l = Loader()
    pdf = l.check_pdf(pdf_name)
    for posset in POSITIVITIES:
        # Use the loader to load the positivity dataset
        ps = l.check_posset(setname=posset, theoryID=THEORYID, postlambda=1e6)
        preds = predictions(ps, pdf)
        core_predictions = PositivityResult.from_convolution(pdf, ps)
        assert_allclose(preds.values, core_predictions.rawdata.T)
        # Now do the same with the API
        api_predictions = API.positivity_predictions_data_result(
            theoryid=THEORYID, pdf=pdf_name, posdataset={"dataset": posset, "maxlambda": 1e6}
        )
        assert_allclose(preds.values, api_predictions.rawdata.T)
        # And now check that the results are correct for any kind of PDF
        cv_predictions = central_predictions(ps, pdf).squeeze()
        assert_allclose(cv_predictions, api_predictions.central_value, atol=1e-3)


def test_extended_predictions():
    """Test the python predictions dataframe stasts with MC sets"""
    l = Loader()
    pdf = l.check_pdf(PDF)
    had = l.check_dataset("ATLASTTBARTOT", theoryid=THEORYID, cfac=("QCD",))
    dis = l.check_dataset("H1HERAF2B", theoryid=THEORYID)
    dis_all = predictions(dis, pdf).T
    dis_central = central_predictions(dis, pdf).T
    assert np.allclose(dis_all.mean().values, dis_central.values)
    dis_linear = linear_predictions(dis, pdf).T
    assert np.all(dis_all == dis_linear)
    had_all = predictions(had, pdf).T
    had_central = central_predictions(had, pdf).T
    had_linear = linear_predictions(had, pdf).T
    assert np.allclose(had_linear.mean().values, had_central)
    assert not np.allclose(had_all.mean().values, had_central)
    assert np.all((had_linear - had_all).std() < had_all.std())
