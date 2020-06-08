import pandas as pd
import numpy as np

from validphys.loader import Loader
from validphys.results import ThPredictionsResult
from validphys.fkparser import load_fktable
from validphys.convolution import predictions, central_predictions, linear_predictions


def test_basic_loading():
    l = Loader()
    # Test both with and without cfactors, and load both DIS and hadronic
    for cfac in ((), ('QCD',)):
        fk = l.check_fktable(setname='ATLASTTBARTOT', theoryID=162, cfac=cfac)
        res = load_fktable(fk)
        assert res.ndata == 3
        assert isinstance(res.sigma, pd.DataFrame)
    fk = l.check_fktable(setname='H1HERAF2B', theoryID=162, cfac=())
    res = load_fktable(fk)
    assert res.ndata == 12
    assert isinstance(res.sigma, pd.DataFrame)


def test_cuts():
    l = Loader()
    ds = l.check_dataset('ATLASTTBARTOT', theoryid=162, cfac=('QCD',))
    table = load_fktable(ds.fkspecs[0])
    # Check explicit cuts
    newtable = table.with_cuts([0, 1])
    assert set(newtable.sigma.index.get_level_values(0)) == {0, 1}
    assert newtable.ndata == 2
    assert newtable.metadata['GridInfo'].ndata == ds.commondata.ndata
    # Check empty cuts
    assert newtable.with_cuts(None) is newtable
    # Check loaded cuts
    ds = l.check_dataset('H1HERAF2B', theoryid=162)
    table = load_fktable(ds.fkspecs[0])
    newtable = table.with_cuts(ds.cuts)
    assert len(newtable.sigma.index.get_level_values(0).unique()) == len(ds.cuts.load())


def test_preditions():
    l = Loader()
    pdf = l.check_pdf('NNPDF31_nnlo_as_0118')
    dss = [
        l.check_dataset(
            'ATLASTTBARTOT', theoryid=162, cfac=('QCD',), cuts=None
        ),  # Had, cfactors, no cuts
        l.check_dataset('H1HERAF2B', theoryid=162),  # DIS, op: NULL
        l.check_dataset('D0ZRAP', theoryid=162),  # op: RATIO
        l.check_dataset('D0WEASY', theoryid=162),  # op: ASY
        l.check_dataset('CMSWCHARMTOT', theoryid=162, cuts=None),  # op: ADD
        l.check_dataset('ATLASWPT31PB', theoryid=162),  # op: SMN
    ]
    for ds in dss:
        preds = predictions(ds, pdf)
        cppres = ThPredictionsResult.from_convolution(pdf, ds)
        assert np.allclose(preds.values, cppres._rawdata)

def test_extended_predictions():
    l = Loader()
    pdf = l.check_pdf('NNPDF31_nnlo_as_0118')
    had = l.check_dataset('ATLASTTBARTOT', theoryid=162, cfac=('QCD',))
    dis = l.check_dataset('H1HERAF2B', theoryid=162)
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
