from validphys.api import API
from validphys.fitdata import print_systype_overlap

def test_print_systype_overlap():
    """Test that print_systype_overlap does expected thing
    for some simple examples. We can't use the API directly
    here because we want to create fictional groups where
    overlaps do exist.

    """
    dsinp_1 = {"dataset":"ATLASWZRAP11"}
    dsinp_2 = {"dataset": "ATLAS1JET11"}

    cd_1 = API.commondata(dataset_input=dsinp_1)
    cd_2 = API.commondata(dataset_input=dsinp_2)

    group = {"group_name": "e2"}
    group_bis = {"group_name": "e2bis"}

    group_diffbis = {"group_name": "e3bis"}

    match = print_systype_overlap([[cd_1],[cd_1]], [group, group_bis])
    assert isinstance(match, tuple)
    match2 = print_systype_overlap([[cd_1]], [group])
    assert isinstance(match2, str)
    match3 = print_systype_overlap([[cd_1],[cd_2]], [group, group_diffbis])
    assert isinstance(match3, tuple)
    match4 = print_systype_overlap([], [])
    assert isinstance(match4, str)
