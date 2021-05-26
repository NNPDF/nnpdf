import pandas as pd
import numpy as np
from numpy.testing import assert_allclose

from validphys.loader import Loader
from validphys.results import ThPredictionsResult
from validphys.fkparser import load_fktable
from validphys.convolution import predictions, central_predictions, linear_predictions
from validphys.tests.conftest import PDF, THEORYID


def test_basic_loading():
    l = Loader()
    # Test both with and without cfactors, and load both DIS and hadronic
    for cfac in ((), ('QCD',)):
        fk = l.check_fktable(setname='ATLASTTBARTOT', theoryID=THEORYID, cfac=cfac)
        res = load_fktable(fk)
        assert res.ndata == 3
        assert isinstance(res.sigma, pd.DataFrame)
    fk = l.check_fktable(setname='H1HERAF2B', theoryID=THEORYID, cfac=())
    res = load_fktable(fk)
    assert res.ndata == 12
    assert isinstance(res.sigma, pd.DataFrame)

    # Check if cfactors for datasets having one entry are correctly parsed
    fk = l.check_fktable(setname='CMSTTBARTOT7TEV', theoryID=THEORYID, cfac=('QCD',))
    res = load_fktable(fk)
    assert res.ndata == 1


def test_cuts():
    l = Loader()
    ds = l.check_dataset('ATLASTTBARTOT', theoryid=THEORYID, cfac=('QCD',))
    table = load_fktable(ds.fkspecs[0])
    # Check explicit cuts
    newtable = table.with_cuts([0, 1])
    assert set(newtable.sigma.index.get_level_values(0)) == {0, 1}
    assert newtable.ndata == 2
    assert newtable.metadata['GridInfo'].ndata == ds.commondata.ndata
    # Check empty cuts
    assert newtable.with_cuts(None) is newtable
    # Check loaded cuts
    ds = l.check_dataset('H1HERAF2B', theoryid=THEORYID)
    table = load_fktable(ds.fkspecs[0])
    newtable = table.with_cuts(ds.cuts)
    assert len(newtable.sigma.index.get_level_values(0).unique()) == len(ds.cuts.load())


def test_predictions():
    l = Loader()
    pdf = l.check_pdf(PDF)
    dss = [
        l.check_dataset(
            'ATLASTTBARTOT', theoryid=THEORYID, cfac=('QCD',)
        ),  # cfactors
        l.check_dataset('H1HERAF2B', theoryid=THEORYID),  # DIS, op: NULL
        l.check_dataset('D0ZRAP', theoryid=THEORYID),  # op: RATIO
        l.check_dataset('D0WEASY', theoryid=THEORYID),  # op: ASY
        l.check_dataset('CMSWCHARMTOT', theoryid=THEORYID),  # op: ADD
        l.check_dataset('ATLASWPT31PB', theoryid=THEORYID),  # op: SMN
        l.check_dataset('DYE906R', theoryid=THEORYID), # op: COM
        l.check_dataset('DYE906_D', theoryid=THEORYID), # op: SMT
    ]
    for ds in dss:
        preds = predictions(ds, pdf)
        cppres = ThPredictionsResult.from_convolution(pdf, ds)
        # Change the atol and rtol from 1e-8 and 1e-7 since DYE906R
        # fails with the default setting
        assert_allclose(preds.values, cppres._rawdata, atol=1e-8, rtol=1e-3)

def test_extended_predictions():
    l = Loader()
    pdf = l.check_pdf(PDF)
    had = l.check_dataset('ATLASTTBARTOT', theoryid=THEORYID, cfac=('QCD',))
    dis = l.check_dataset('H1HERAF2B', theoryid=THEORYID)
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
