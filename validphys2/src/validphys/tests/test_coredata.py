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

def test_with_ADD_sys():
    """
    test that replacement of ADD columns with themselves
    works
    """
    l = Loader()
    observable = SINGLE_DATASET["dataset"]
    commondata = l.check_commondata(observable).load_commondata_instance()
    add_sys = commondata.commondata_table["ADD"].copy()

    new_commondata = commondata.with_ADD_sys(add_sys=add_sys)

    assert( np.sum((commondata.commondata_table["ADD"] - new_commondata.commondata_table["ADD"]).to_numpy()) <= 1e-12 )


def test_multiplicative_errors_rescale():
    """
    tests that rescaling both CORR and UNCORR MULT systematics
    by 1 leaves the MULT sys unchanged
    """
    l = Loader()
    observable = SINGLE_DATASET["dataset"]
    commondata = l.check_commondata(observable).load_commondata_instance()

    mult_table = commondata.commondata_table["MULT"]
    rescaled_mult_table = commondata.multiplicative_errors_rescale(CORR=True,UNCORR=True,sys_rescaling_factor=1)
    
    assert( np.sum((rescaled_mult_table - mult_table).to_numpy()) <= 1e-12 )

def test_additive_errors_rescale():
    """
    tests that rescaling both CORR and UNCORR ADD systematics
    by 1 leaves the ADD sys unchanged
    """
    l = Loader()
    observable = SINGLE_DATASET["dataset"]
    commondata = l.check_commondata(observable).load_commondata_instance()

    add_table = commondata.commondata_table["ADD"]
    rescaled_add_table = commondata.additive_errors_rescale(CORR=True,UNCORR=True,sys_rescaling_factor=1)
    
    assert( np.sum((rescaled_add_table - add_table).to_numpy()) <= 1e-12 )