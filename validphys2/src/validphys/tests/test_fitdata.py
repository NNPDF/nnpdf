from validphys.api import API
from validphys.fitdata import print_systype_overlap, print_different_cuts
from validphys.tests.conftest import FIT_3REPLICAS, FIT_3REPLICAS_DCUTS


def test_print_different_cuts():
    """Checks the print_different_cuts functions
    using two fits with a different choice of q2min and w2min in the runcard
    One of the datasets (SLACP) gets 0 points in in the most restrictive case
    The different cuts are:
    q2min: 3.49 - 13.49
    w2min: 12.5 - 22.5
    """
    fit_1 = API.fit(fit=FIT_3REPLICAS)
    fit_2 = API.fit(fit=FIT_3REPLICAS_DCUTS)
    fits = [fit_1, fit_2]
    testi = API.test_for_same_cuts(fits=[FIT_3REPLICAS, FIT_3REPLICAS_DCUTS], use_cuts="fromfit")
    res = print_different_cuts(fits, testi)
    assert "121 out of 260" in res
    assert "59 out of 260" in res
    assert "33 out of 211" in res
    assert "0 out of 211" in res


def test_print_systype_overlap():
    """Test that print_systype_overlap does expected thing
    for some simple examples. We can't use the API directly
    here because we want to create fictional groups where
    overlaps do exist.

    Note that the first input of py:func:`print_systype_overlap` is
    ``groups_commondata`` which is a list of lists, the outer list usually
    contains an inner list for each ``metadata_group``. Each inner list contains
    a ``CommonDataSpec`` for each dataset which is part of that group. In this
    test we create fake groups and ensure the output of the function is correct.

    The second input is ``group_dataset_inputs_by_metadata`` which is a list
    containing a dictionary for each ``metadata_group``. The function gets
    ``group_name`` from each dictionary and uses to label each group, but the
    actual value is unimportant for these tests.

    """
    cd_1 = API.commondata(dataset_input={"dataset":"ATLASWZRAP11"})
    cd_2 = API.commondata(dataset_input={"dataset": "ATLAS1JET11"})
    cd_3 = API.commondata(dataset_input={"dataset": "NMC"})

    # group names don't affect results, set arbitrarily.
    group_1 = {"group_name": "group_1"}
    group_2 = {"group_name": "group_2"}

    # each group contains same dataset, so systypes will overlap
    match = print_systype_overlap([[cd_1],[cd_1]], [group_1, group_2])
    assert isinstance(match, tuple)
    # single group passed so systype won't overlap
    match2 = print_systype_overlap([[cd_1]], [group_1])
    assert isinstance(match2, str)
    # cd in each group are different but share a systematic so overlap.
    match3 = print_systype_overlap([[cd_1],[cd_2]], [group_1, group_2])
    assert isinstance(match3, tuple)
    # test non-overlapping groups
    match4 = print_systype_overlap([[cd_1, cd_2], [cd_3]], [group_1, group_2])
    assert isinstance(match4, str)
    # no groups, no overlap
    match5 = print_systype_overlap([], [])
    assert isinstance(match5, str)
