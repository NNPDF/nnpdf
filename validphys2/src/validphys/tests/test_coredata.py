from validphys.tests.conftest import SINGLE_DATASET
from validphys.loader import Loader
import numpy as np

def test_with_MULT_sys():
    """
    test that replacement of MULT columns with themselves
    works
    """
    l = Loader()
    observable = SINGLE_DATASET["dataset"]
    commondata = l.check_commondata(observable).load_commondata_instance()
    mult_sys = commondata.commondata_table["MULT"].copy()

    new_commondata = commondata.with_MULT_sys(mult_sys=mult_sys)

    assert( np.sum((commondata.commondata_table["MULT"] - new_commondata.commondata_table["MULT"]).to_numpy()) <= 1e-12 )

