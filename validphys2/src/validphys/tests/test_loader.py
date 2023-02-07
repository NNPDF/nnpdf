"""
test_loader.py

Test loading utilities.
"""
import numpy as np
import pytest
from hypothesis.strategies import sampled_from, sets, composite
from hypothesis import given, settings

from validphys.core import Cuts, CommonDataSpec
from validphys.loader import FallbackLoader, rebuild_commondata_without_cuts, FitNotFound
from validphys.plotoptions import kitable, get_info
from validphys.tests.conftest import FIT

l = FallbackLoader()
#The sorted is to appease hypothesis
dss = sorted(l.available_datasets - {"PDFEVOLTEST"})

class MockCuts():
    def __init__(self, arr):
        self.arr = arr
    def load(self):
        return self.arr


@composite
def commondata_and_cuts(draw):
    cd = l.check_commondata(draw(sampled_from(dss)))
    ndata = cd.metadata.ndata
    # Get a cut mask with at least one selected datapoint
    masks = sets(sampled_from(range(ndata)), min_size=1)
    mask = sorted(draw(masks))
    return cd, mask


@given(arg=commondata_and_cuts())
@settings(deadline=None)
def test_rebuild_commondata_without_cuts(tmp_path_factory, arg):
    # We need to create a new directory for each call of the test
    # otherwise we get files mixed together
    tmp = tmp_path_factory.mktemp("test_loader")

    cd, cuts = arg
    lcd = cd.load()
    cutspec = None
    if cuts:
        cutpath = tmp / "cuts.txt"
        np.savetxt(cutpath, np.asarray(cuts, dtype=int), fmt="%u")
        cutspec = Cuts(cd, cutpath)
        lcd = lcd.with_cuts(cuts)
    lcd.export(tmp)
    # We have to reconstruct the name here...
    with_cuts = tmp / f"DATA_{cd.name}.dat"
    newpath = tmp / "commondata.dat"
    rebuild_commondata_without_cuts(with_cuts, cutspec, cd.datafile, newpath)
    newcd = CommonDataSpec(newpath, cd.sysfile, cd.plotfiles)
    # Note this one is without cuts
    t1 = kitable(cd, get_info(cd))
    t2 = kitable(newcd, get_info(newcd))
    assert (t1 == t2).all
    lncd = newcd.load()
    if cuts:
        assert np.allclose(lncd.get_cv()[cuts], lcd.get_cv())
        nocuts = np.ones(cd.ndata, dtype=bool)
        nocuts[cuts] = False
        assert (lncd.get_cv()[nocuts] == 0).all()

@given(inp=commondata_and_cuts())
def test_kitable_with_cuts(inp):
    cd, cuts = inp
    info = get_info(cd, cuts=cuts)
    tb = kitable(cd, info, cuts=MockCuts(cuts))
    assert len(tb) == len(cuts)

def test_load_fit():
    assert l.check_fit(FIT)
    with pytest.raises(FitNotFound):
        l.check_fit(f"{FIT}/")
