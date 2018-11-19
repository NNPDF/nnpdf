"""
test_loader.py

Test loading utilities.
"""
import numpy as np
from hypothesis.strategies import sampled_from, sets, composite
from hypothesis import given

from validphys.core import Cuts, CommonDataSpec
from validphys.loader import Loader, rebuild_commondata_without_cuts
from validphys.plotoptions import kitable, get_info

l = Loader()
#The sorted is to appease hypothesis
dss = sorted(l.available_datasets - {'PDFEVOLTEST'})

@composite
def commodata_and_cuts(draw):
    cd = l.check_commondata(draw(sampled_from(dss)))
    ndata = cd.metadata.ndata
    #TODO: Maybe upgrade to this
    #https://github.com/HypothesisWorks/hypothesis/issues/1115
    mask = sorted(draw(sets(sampled_from(range(ndata)))))
    return cd, mask

ipath = 0

@given(arg=commodata_and_cuts())
def test_rebuild_commondata_without_cuts(tmp, arg):
    #We have to do this because otherwise files get mixed together. Note that
    #tmp does fire once for all hypothesis runs.
    global ipath
    tmp = tmp/f'{ipath}'
    tmp.mkdir()
    ipath += 1

    cd, cuts = arg
    lcd = cd.load()
    cutspec = None
    if cuts:
        cutpath = tmp/'cuts.txt'
        np.savetxt(cutpath, np.asarray(cuts, dtype=int), fmt='%u')
        cutspec = Cuts(cd.name, cutpath)
        lcd = type(lcd)(lcd, cuts)
    lcd.Export(str(tmp))
    #We have to reconstruct the name here...
    with_cuts = tmp/f'DATA_{cd.name}.dat'
    newpath = tmp/'commondata.dat'
    rebuild_commondata_without_cuts(with_cuts, cutspec, cd.datafile, newpath)
    newcd = CommonDataSpec(newpath, cd.sysfile, cd.plotfiles)
    #Note this one is without cuts
    t1 = kitable(cd, get_info(cd))
    t2 = kitable(newcd, get_info(newcd))
    assert (t1==t2).all
    lncd = newcd.load()
    if cuts:
        assert np.allclose(lncd.get_cv()[cuts],  lcd.get_cv())
        nocuts = np.ones(cd.ndata, dtype=bool)
        nocuts[cuts] = False
        assert (lncd.get_cv()[nocuts] == 0).all()
