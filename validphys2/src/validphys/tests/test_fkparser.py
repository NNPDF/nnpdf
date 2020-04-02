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

