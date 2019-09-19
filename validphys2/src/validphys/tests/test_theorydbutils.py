import pytest

from validphys.loader import Loader
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
