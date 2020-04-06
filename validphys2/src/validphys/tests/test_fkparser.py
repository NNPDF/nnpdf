import pandas as pd

from validphys.fkparser import load_fktable
from validphys.loader import Loader


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
