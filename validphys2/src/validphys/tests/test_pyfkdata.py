import numpy as np
from numpy.testing import assert_allclose
import pandas as pd
import pytest

from validphys.api import API
from validphys.convolution import central_predictions, linear_predictions, predictions
from validphys.fkparser import load_fktable
from validphys.loader import FallbackLoader as Loader
from validphys.results import PositivityResult, ThPredictionsResult
from validphys.tests.conftest import HESSIAN_PDF, PDF, POSITIVITIES, THEORYID

# Use specific datasets (not necessarily from conftest) to check loading also with cfactors
DS1 = "CMS_Z0J_8TEV_PT-Y"  # hadronic, needs NRM cfactor
CFAC_DS1 = "NRM"
NDATA_DS1 = 50
DS2 = "HERA_NC_318GEV_EP-SIGMARED"
NDATA_DS2 = 485


def test_basic_loading():
    """Test the loading of an old theory using directly the legacy name"""
    l = Loader()
    # Test both with and without cfactors, and load both DIS and hadronic
    for cfac in ((), (CFAC_DS1,)):
        ds = l.check_dataset(DS1, theoryid=THEORYID, cfac=cfac)
        res = load_fktable(ds.fkspecs[0])
        assert res.ndata == NDATA_DS1
        assert isinstance(res.sigma, pd.DataFrame)
    ds = l.check_dataset(DS2, theoryid=THEORYID)
    res = load_fktable(ds.fkspecs[0])
    assert res.ndata == NDATA_DS2
    assert isinstance(res.sigma, pd.DataFrame)


def test_cuts():
    l = Loader()
    ds = l.check_dataset(DS1, theoryid=THEORYID, cfac=(CFAC_DS1,), variant="legacy_10")
    table = load_fktable(ds.fkspecs[0])
    # Check explicit cuts
    newtable = table.with_cuts([0, 1])
    assert set(newtable.sigma.index.get_level_values(0)) == {0, 1}
    assert newtable.ndata == 2
    assert table.ndata == ds.commondata.ndata
    # Check empty cuts
    assert newtable.with_cuts(None) is newtable
    # Check loaded cuts
    ds = l.check_dataset(DS2, theoryid=THEORYID, variant="legacy")
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
        {"name": DS1, "cfac": (CFAC_DS1,)},  # cfactors
        {"name": DS2},  # DIS, op: NULL
        {"name": "D0_Z0_1P96TEV_ZRAP"},  # op: RATIO
        {"name": "D0_WPWM_1P96TEV_ASY"},  # op: ASY
        {"name": "CMS_SINGLETOP_7TEV_TCHANNEL-XSEC"},  # op: ADD
        {"name": "DYE906_Z0_120GEV_DW_PDXSECRATIO"},  # op: COM
        # Not used in any dataset:
        #         {"name": "DYE906_D"},  # op: SMT
        #         {"name": "ATLASWPT31PB"},  # op: SMN
    ]
    for daset in datasets:
        daset["variant"] = "legacy"
        ds = l.check_dataset(**daset, theoryid=THEORYID)
        preds = predictions(ds, pdf)
        core_predictions = ThPredictionsResult.from_convolution(pdf, ds)
        # Uses rawdata since we want to check all members for which we computed the convolution
        assert_allclose(preds.values, core_predictions.rawdata, atol=1e-8, rtol=1e-3)
        # Now check that the stats class does the right thing with the data
        cv_predictions = central_predictions(ds, pdf).squeeze()
        stats_predictions = pdf.stats_class(preds.T)
        # rtol to 1e-2 due to DYE906R and DYE906_D for MC sets
        # TODO: check whether this tolerance can be decreased when using double precision
        assert_allclose(cv_predictions, stats_predictions.central_value(), rtol=1e-2)
        assert_allclose(
            core_predictions.error_members, stats_predictions.error_members().T, rtol=1e-3
        )
        assert_allclose(
            core_predictions.central_value, stats_predictions.central_value(), rtol=1e-2
        )


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
        ps = l.check_posset(setname=posset, theoryID=THEORYID, postlambda=1e6, rules=())
        preds = predictions(ps, pdf)
        core_predictions = PositivityResult.from_convolution(pdf, ps)
        assert_allclose(preds.values, core_predictions.rawdata)
        # Now do the same with the API
        api_predictions = API.positivity_predictions_data_result(
            theoryid=THEORYID,
            use_cuts="internal",
            pdf=pdf_name,
            posdataset={"dataset": posset, "maxlambda": 1e6},
        )
        assert_allclose(preds.values, api_predictions.rawdata)
        # And now check that the results are correct for any kind of PDF
        cv_predictions = central_predictions(ps, pdf).squeeze()
        assert_allclose(cv_predictions, api_predictions.central_value, atol=1e-3)


def test_extended_predictions():
    """Test the python predictions dataframe stats with MC sets"""
    l = Loader()
    pdf = l.check_pdf(PDF)
    had = l.check_dataset(DS1, theoryid=THEORYID, cfac=(CFAC_DS1,), variant="legacy_10")
    dis = l.check_dataset(DS2, theoryid=THEORYID)
    # For DIS we expect the central value to be exactly the same as the average
    dis_all = predictions(dis, pdf).T
    dis_central = central_predictions(dis, pdf).T
    np.testing.assert_allclose(dis_all.mean().values.flat, dis_central.values.flat)
    # For hadronic predictions instead we expect some small differences
    had_all = predictions(had, pdf).T
    had_central = central_predictions(had, pdf).T
    assert not np.allclose(had_all.mean().values, had_central)
    assert np.all(np.abs(had_all.mean() - had_central) < had_all.std())
