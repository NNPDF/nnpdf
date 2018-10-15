from validphys.loader import FallbackLoader as Loader
from validphys.fitdata import print_systype_overlap
from validphys.core import ExperimentSpec

def test_print_systype_overlap():
    l = Loader()
    ds = l.check_dataset(name='ATLASWZRAP11', theoryid=162, cuts=None)
    e = ExperimentSpec('e2', [ds])
    ebis = ExperimentSpec('e2bis', [ds])
    match = print_systype_overlap([e,ebis])
    assert isinstance(match, tuple)
    match2 = print_systype_overlap([e])
    assert isinstance(match2, str)
    ds2 = l.check_dataset(name='ATLAS1JET11', theoryid=162, cuts=None)
    ediffbis = ExperimentSpec('e3bis', [ds2])
    match3 = print_systype_overlap([e, ediffbis])
    assert isinstance(match3, tuple)
    match4 = print_systype_overlap([])
    assert isinstance(match4, str)
