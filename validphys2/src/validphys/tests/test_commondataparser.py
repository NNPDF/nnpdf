import pandas as pd

from validphys.commondataparser import load_commondata
from validphys.loader import Loader

def test_basic_commondata_loading():
    l = Loader()
    cd = l.check_commondata(setname='H1HERAF2B')
    res = load_commondata(cd)
    assert res.ndata == 12
    assert isinstance(res.data, pd.DataFrame)