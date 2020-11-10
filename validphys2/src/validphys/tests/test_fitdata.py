from validphys.api import API
from validphys.fitdata import print_systype_overlap

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
