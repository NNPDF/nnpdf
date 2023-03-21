import numpy as np
from numpy.testing import assert_allclose
from validphys.loader import Loader
from validphys.tests.conftest import SINGLE_DATASET
from validphys.inconsistent_ct import InconsistentCommonData

l = Loader()
cd = l.check_commondata(SINGLE_DATASET['dataset']).load_commondata_instance()

inconsys_cd = InconsistentCommonData(setname=cd.setname, ndata=cd.ndata, 
                        commondataproc=cd.commondataproc, 
                        nkin=cd.nkin, nsys=cd.nsys, 
                        commondata_table = cd.commondata_table, 
                        systype_table = cd.systype_table)



def test_with_MULT_sys():
    """
    test if MULT commondata_table is
    replaced correctly by
    dataclasses.replace(self, commondata_table = new_table)
    """

    mult_sys_tab = 3 * cd.commondata_table["MULT"].to_numpy()
    
    inc_mult_sys_tab = inconsys_cd.with_MULT_sys(mult_sys_tab).commondata_table["MULT"].to_numpy()
    
    assert_allclose(mult_sys_tab, inc_mult_sys_tab, rtol=1e-5, atol=1e-5)


def test_with_ADD_sys():
    """
    test if ADD commondata_table is
    replaced correctly by
    dataclasses.replace(self, commondata_table = new_table)
    """

    mult_sys_tab = 3 * cd.commondata_table["ADD"].to_numpy()
    
    inc_mult_sys_tab = inconsys_cd.with_ADD_sys(mult_sys_tab).commondata_table["ADD"].to_numpy()
    
    assert_allclose(mult_sys_tab, inc_mult_sys_tab, rtol=1e-5, atol=1e-5)


def test_rescale_sys_CORR_MULT():
    """
    Check whether rescaling of
    CORR MULT uncertainties works
    as expected
    """

    rescaling_factor = 2.
    type_err = "MULT"
    new_icd = inconsys_cd.with_MULT_sys(inconsys_cd.rescale_sys(type_err=type_err, CORR=True, 
                                UNCORR=False, sys_rescaling_factor=rescaling_factor))
    

    err_table = cd.systematics_table.loc[:,[type_err]].copy()
    # get indices of CORR sys
    systype_corr = cd.systype_table[(cd.systype_table["type"] == type_err) 
                        & (~cd.systype_table["name"].isin(["UNCORR","THEORYUNCORR"]))]

    tab2 = rescaling_factor * err_table.iloc[:,systype_corr.index-1].to_numpy()

    tab1 = new_icd.systematics_table.loc[:,[type_err]].iloc[:,systype_corr.index -1]
    
    assert_allclose(tab1,tab2, rtol=1e-5, atol=1e-5)


def test_rescale_sys_CORR_ADD():
    """
    Check whether rescaling of
    CORR ADD uncertainties works
    as expected
    """

    rescaling_factor = 2.
    type_err = "ADD"
    new_icd = inconsys_cd.with_ADD_sys(inconsys_cd.rescale_sys(type_err=type_err, CORR=True, 
                                UNCORR=False, sys_rescaling_factor=rescaling_factor))


    err_table = cd.systematics_table.loc[:,[type_err]].copy()
    # get indices of CORR sys
    systype_corr = cd.systype_table[(cd.systype_table["type"] == type_err) 
                        & (~cd.systype_table["name"].isin(["UNCORR","THEORYUNCORR"]))]

    tab2 = rescaling_factor * err_table.iloc[:,systype_corr.index-1].to_numpy()

    tab1 = new_icd.systematics_table.loc[:,[type_err]].iloc[:,systype_corr.index -1]
    
    assert_allclose(tab1,tab2, rtol=1e-5, atol=1e-5)


def test_process_commondata():
    """
    Check whether process_commondata
    leaves the commondata instance
    unchanged when told to do so.
    """
    
    new_icd = inconsys_cd.process_commondata(
                    ADD=False,MULT=False,
                    CORR=False,UNCORR=False,
                    inconsistent_datasets=[SINGLE_DATASET['dataset']],
                    sys_rescaling_factor=1
                    )
    tab1 = new_icd.commondata_table.drop(['process'],axis=1).to_numpy()
    tab2 = inconsys_cd.commondata_table.drop(['process'], axis=1).to_numpy()
    
    assert_allclose(tab1,tab2, rtol=1e-5, atol=1e-5)


def test_process_commondata_CORR_MULT():
    """
    Check whether rescaling of
    CORR MULT uncertainties works
    as expected with process_commondata
    method
    """
    type_err = "MULT"
    rescaling_factor = 2.
    new_icd = inconsys_cd.process_commondata(
                    ADD=False,MULT=True,
                    CORR=True,UNCORR=False,
                    inconsistent_datasets=[SINGLE_DATASET['dataset']],
                    sys_rescaling_factor=rescaling_factor
                    )

    err_table = cd.systematics_table.loc[:,[type_err]].copy()
    # get indices of CORR sys
    systype_corr = cd.systype_table[(cd.systype_table["type"] == type_err) 
                        & (~cd.systype_table["name"].isin(["UNCORR","THEORYUNCORR"]))]

    tab2 = rescaling_factor * err_table.iloc[:,systype_corr.index-1].to_numpy()

    tab1 = new_icd.systematics_table.loc[:,[type_err]].iloc[:,systype_corr.index -1]
    
    assert_allclose(tab1,tab2, rtol=1e-5, atol=1e-5)


def test_process_commondata_CORR_ADD():
    """
    Check whether rescaling of
    CORR ADD uncertainties works
    as expected with process_commondata
    method
    """
    type_err = "ADD"
    rescaling_factor = 2.
    new_icd = inconsys_cd.process_commondata(
                    ADD=True,MULT=False,
                    CORR=True,UNCORR=False,
                    inconsistent_datasets=[SINGLE_DATASET['dataset']],
                    sys_rescaling_factor=rescaling_factor
                    )

    err_table = cd.systematics_table.loc[:,[type_err]].copy()
    # get indices of CORR sys
    systype_corr = cd.systype_table[(cd.systype_table["type"] == type_err) 
                        & (~cd.systype_table["name"].isin(["UNCORR","THEORYUNCORR"]))]

    tab2 = rescaling_factor * err_table.iloc[:,systype_corr.index-1].to_numpy()

    tab1 = new_icd.systematics_table.loc[:,[type_err]].iloc[:,systype_corr.index -1]
    
    assert_allclose(tab1,tab2, rtol=1e-5, atol=1e-5)


