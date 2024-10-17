"""
Module to test the InconsistentCommonData class.
"""

from numpy.testing import assert_allclose

from validphys.tests.conftest import SINGLE_DATASET
from validphys.closuretest.inconsistent_closuretest.inconsistent_ct import InconsistentCommonData


def load_cd():
    """
    Load a commondata instance and an inconsistent commondata instance.
    """
    # avoid circular import
    from validphys.api import API

    cd = API.commondata(**{"dataset_input": {**SINGLE_DATASET}}).load()

    inconsys_cd = InconsistentCommonData(
        setname=cd.setname,
        ndata=cd.ndata,
        commondataproc=cd.commondataproc,
        nkin=cd.nkin,
        nsys=cd.nsys,
        commondata_table=cd.commondata_table,
        systype_table=cd.systype_table,
    )
    return cd, inconsys_cd


def test_with_MULT_sys():
    """
    test if MULT commondata_table is
    replaced correctly by
    dataclasses.replace(self, commondata_table = new_table)
    """
    cd, inconsys_cd = load_cd()

    mult_sys_tab = 3 * cd.commondata_table["MULT"].to_numpy()

    inc_mult_sys_tab = inconsys_cd.with_MULT_sys(mult_sys_tab).commondata_table["MULT"].to_numpy()

    assert_allclose(mult_sys_tab, inc_mult_sys_tab)


def test_with_ADD_sys():
    """
    test if ADD commondata_table is
    replaced correctly by
    dataclasses.replace(self, commondata_table = new_table)
    """
    cd, inconsys_cd = load_cd()
    mult_sys_tab = 3 * cd.commondata_table["ADD"].to_numpy()

    inc_mult_sys_tab = inconsys_cd.with_ADD_sys(mult_sys_tab).commondata_table["ADD"].to_numpy()

    assert_allclose(mult_sys_tab, inc_mult_sys_tab)


def test_rescale_sys_CORR_MULT():
    """
    Check whether rescaling of
    CORR MULT uncertainties works
    as expected
    """
    cd, inconsys_cd = load_cd()

    rescaling_factor = 2.0
    treatment_err = "MULT"
    new_icd = inconsys_cd.with_MULT_sys(
        inconsys_cd.rescale_sys(
            treatment_err=treatment_err,
            CORR=True,
            UNCORR=False,
            SPECIAL=False,
            sys_rescaling_factor=rescaling_factor,
        )
    )

    # get indices of CORR sys
    systype_corr = cd.systype_table[
        (cd.systype_table["treatment"] == treatment_err)
        & (~cd.systype_table["name"].isin(["UNCORR", "THEORYUNCORR"]))
    ]

    tab2 = rescaling_factor * cd.systematics_table.iloc[:, systype_corr.index - 1].to_numpy()

    tab1 = new_icd.systematics_table.iloc[:, systype_corr.index - 1]

    assert_allclose(tab1, tab2)


def test_rescale_sys_CORR_ADD():
    """
    Check whether rescaling of
    CORR ADD uncertainties works
    as expected
    """
    cd, inconsys_cd = load_cd()

    rescaling_factor = 2.0
    treatment_err = "ADD"
    new_icd = inconsys_cd.with_ADD_sys(
        inconsys_cd.rescale_sys(
            treatment_err,
            CORR=True,
            UNCORR=False,
            SPECIAL=False,
            sys_rescaling_factor=rescaling_factor,
        )
    )

    # get indices of CORR sys
    systype_corr = cd.systype_table[
        (cd.systype_table["treatment"] == treatment_err)
        & (~cd.systype_table["name"].isin(["UNCORR", "THEORYUNCORR"]))
    ]

    tab2 = rescaling_factor * cd.systematics_table.iloc[:, systype_corr.index - 1].to_numpy()

    tab1 = new_icd.systematics_table.iloc[:, systype_corr.index - 1]

    assert_allclose(tab1, tab2)


def test_process_commondata():
    """
    Check whether process_commondata
    leaves the commondata instance
    unchanged when told to do so.
    """
    cd, inconsys_cd = load_cd()
    new_icd = inconsys_cd.process_commondata(
        ADD=False,
        MULT=False,
        CORR=False,
        UNCORR=False,
        SPECIAL=False,
        inconsistent_datasets=[SINGLE_DATASET['dataset']],
        sys_rescaling_factor=1,
    )
    tab1 = new_icd.commondata_table.drop(['process'], axis=1).to_numpy()
    tab2 = inconsys_cd.commondata_table.drop(['process'], axis=1).to_numpy()

    assert_allclose(tab1, tab2)


def test_process_commondata_CORR_MULT():
    """
    Check whether rescaling of
    CORR MULT uncertainties works
    as expected with process_commondata
    method
    """
    cd, inconsys_cd = load_cd()
    treatment_err = "MULT"
    rescaling_factor = 2.0
    new_icd = inconsys_cd.process_commondata(
        ADD=False,
        MULT=True,
        CORR=True,
        UNCORR=False,
        SPECIAL=False,
        inconsistent_datasets=[SINGLE_DATASET['dataset']],
        sys_rescaling_factor=rescaling_factor,
    )

    # get indices of CORR sys
    systype_corr = cd.systype_table[
        (cd.systype_table["treatment"] == treatment_err)
        & (~cd.systype_table["name"].isin(["UNCORR", "THEORYUNCORR"]))
    ]

    tab2 = rescaling_factor * cd.systematics_table.iloc[:, systype_corr.index - 1].to_numpy()

    tab1 = new_icd.systematics_table.iloc[:, systype_corr.index - 1]

    assert_allclose(tab1, tab2)


def test_process_commondata_CORR_ADD():
    """
    Check whether rescaling of
    CORR ADD uncertainties works
    as expected with process_commondata
    method
    """
    cd, inconsys_cd = load_cd()
    treatment_err = "ADD"
    rescaling_factor = 2.0
    new_icd = inconsys_cd.process_commondata(
        ADD=True,
        MULT=False,
        CORR=True,
        UNCORR=False,
        SPECIAL=False,
        inconsistent_datasets=[SINGLE_DATASET['dataset']],
        sys_rescaling_factor=rescaling_factor,
    )

    # get indices of CORR sys
    systype_corr = cd.systype_table[
        (cd.systype_table["treatment"] == treatment_err)
        & (~cd.systype_table["name"].isin(["UNCORR", "THEORYUNCORR"]))
    ]

    tab2 = rescaling_factor * cd.systematics_table.iloc[:, systype_corr.index - 1].to_numpy()

    tab1 = new_icd.systematics_table.iloc[:, systype_corr.index - 1]

    assert_allclose(tab1, tab2)
