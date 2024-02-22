import pandas as pd
import pytest

from validphys.api import API
from validphys.commondataparser import load_commondata
from validphys.loader import FallbackLoader as Loader
from validphys.tests.conftest import FIT, THEORYID


def test_basic_commondata_loading():
    l = Loader()
    cd = l.check_commondata(setname="H1HERAF2B")
    res = load_commondata(cd)
    # Test commondata loading
    assert res.ndata == 12
    assert isinstance(res.commondata_table, pd.DataFrame)
    # Test systype loading
    assert res.nsys == 25
    assert isinstance(res.systype_table, pd.DataFrame)

    # Test a dataset with no systematics
    emptysyscd = l.check_posset(theoryID=THEORYID, setname='POSDYCBD', postlambda=1e-10)
    emptysysres = load_commondata(emptysyscd.commondata)
    assert emptysysres.nsys == 0
    assert emptysysres.systype_table.empty is True


def test_commondata_with_cuts():
    l = Loader()
    setname = "NMC"

    cd = l.check_commondata(setname=setname)
    loaded_cd = load_commondata(cd)

    fit_cuts = l.check_fit_cuts(fit=FIT, commondata=cd)
    internal_cuts = l.check_internal_cuts(cd, API.rules(theoryid=THEORYID, use_cuts="internal"))

    loaded_cd_fit_cuts = loaded_cd.with_cuts(fit_cuts)
    # We must do these - 1 subtractions due to the fact that cuts indexing
    # starts at 0 while commondata indexing starts at 1
    assert all(loaded_cd_fit_cuts.commondata_table.index - 1 == fit_cuts.load())
    assert all(loaded_cd_fit_cuts.additive_errors.index - 1 == fit_cuts.load())
    assert all(loaded_cd_fit_cuts.multiplicative_errors.index - 1 == fit_cuts.load())

    loaded_cd_internal_cuts = loaded_cd.with_cuts(internal_cuts)
    assert all(loaded_cd_internal_cuts.commondata_table.index - 1 == internal_cuts.load())

    loaded_cd_nocuts = loaded_cd.with_cuts(None)
    assert all(loaded_cd_nocuts.commondata_table.index == range(1, cd.ndata + 1))

    preloaded_fit_cuts = fit_cuts.load()
    loaded_cd_preloaded_cuts = loaded_cd.with_cuts(fit_cuts)
    assert all(loaded_cd_preloaded_cuts.commondata_table.index - 1 == preloaded_fit_cuts)

    assert all(loaded_cd.with_cuts([1, 2, 3]).commondata_table.index - 1 == [1, 2, 3])

    # Check that giving cuts for another dataset raises the correct ValueError exception
    cd_bad = l.check_commondata(setname="NMCPD")
    bad_cuts = l.check_fit_cuts(fit=FIT, commondata=cd_bad)
    with pytest.raises(ValueError):
        loaded_cd.with_cuts(bad_cuts)
