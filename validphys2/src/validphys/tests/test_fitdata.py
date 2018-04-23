from validphys.loader import FallbackLoader as Loader
from validphys.fitdata import print_systype_overlap
from validphys.core import ExperimentSpec

def test_print_systype_overlap():
    l = Loader()
    ds = l.check_dataset(name='ATLASWZRAP11', theoryid=162, use_cuts=False)
    e = ExperimentSpec('e2', [ds])
    ebis = ExperimentSpec('e2bis', [ds])
    match = print_systype_overlap([e,ebis])
    assert isinstance(match, tuple)
