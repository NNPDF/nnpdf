import pytest
from validobj import ValidationError

from nnpdf_data.theorydbutils import TheoryNotFoundInDatabase, fetch_all, fetch_theory
from validphys.api import API
from validphys.loader import Loader
from validphys.utils import yaml_safe

L = Loader()
DBPATH = L.theorydb_folder

# TODO: take it directly from conftest
THEORYID = 40_000_000
TEST_COMMENT = 'NNLO QCD theory. Equivalent to 4.0 baseline, DIS use nFONLL.'


def test_fetch_theory():
    with pytest.raises(TheoryNotFoundInDatabase):
        fetch_theory(DBPATH, 999999)
    th40 = fetch_theory(DBPATH, THEORYID)
    assert th40['ID'] == THEORYID


def test_fetch_all():
    res = fetch_all(DBPATH)
    assert res.loc[THEORYID, 'Comments'] == TEST_COMMENT


def test_vp_theoryinfo_tables():
    all_tables = API.all_theory_info_table()
    assert all_tables.loc[THEORYID, 'Comments'] == TEST_COMMENT
    tb40 = API.theory_info_table(theory_db_id=THEORYID)
    assert tb40.loc['Comments'].iloc[0] == TEST_COMMENT


def _dump_and_check_error(tdict, tmp, bad_number=999):
    """Dump theory dict to a file and load expecting an error"""
    tdict["ID"] = bad_number
    ofile = tmp / f"{bad_number}.yaml"
    yaml_safe.dump(tdict, ofile.open("w"))
    with pytest.raises(ValidationError):
        fetch_theory(tmp, bad_number)


def test_fetch_with_errors(tmp):
    """Test some of the errors that theories can raise"""
    # Get a good theory and copy it before doing evil to it
    th_good = fetch_theory(DBPATH, THEORYID)
    # Wrong PTO
    bad_pto = dict(th_good)
    bad_pto["PTO"] = 7
    _dump_and_check_error(bad_pto, tmp)
    # Wrong ModEv
    bad_ev = dict(th_good)
    bad_ev["ModEv"] = "WRONG"
    _dump_and_check_error(bad_ev, tmp)
