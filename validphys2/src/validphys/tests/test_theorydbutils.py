import pytest

from validphys.loader import Loader
from validphys.api import API
from validphys.theorydbutils import fetch_theory, TheoryNotFoundInDatabase, fetch_all


L = Loader()
DBPATH = L.theorydb_file


def test_fetch_theory():
    with pytest.raises(TheoryNotFoundInDatabase):
        fetch_theory(DBPATH, 999999)
    th53 = fetch_theory(DBPATH, 53)
    assert th53['ID'] == 53


def test_fetch_all():
    res = fetch_all(DBPATH)
    assert res.loc[53, 'Comments'] == 'NNPDF3.1 NNLO central'


def test_vp_theoryinfo_tables():
    all_tables = API.all_theory_info_table()
    assert all_tables.loc[53, 'Comments'] == 'NNPDF3.1 NNLO central'
    tb53 = API.theory_info_table(theory_db_id=53)
    assert tb53.loc['Comments'].iloc[0] == 'NNPDF3.1 NNLO central'
