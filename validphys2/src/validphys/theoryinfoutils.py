"""Module for querying the theory database file and representing the data
as a panda dataframe.
"""
from pathlib import Path
# Keep the sqlite3 dependence in one location
import sqlite3

import pandas as pd

class TheoryNotFoundInDatabase(Exception): pass

def make_query(query: str, dbpath: Path):
    """base level database query function"""
    #Note this still requires a string and not a path
    conn = sqlite3.connect(str(dbpath))
    with conn:
        cursor = conn.cursor()
        res = cursor.execute(query)
    return res

def fetch_theory(theory_database: Path, theoryID: int):
    """ Looks in the theory database and returns a dictionary of theory info for the
    theory number specified by `theoryID`.
    """
    #int casting is intentional to avoid malformed querys.
    query = f"SELECT * FROM TheoryIndex WHERE ID={int(theoryID)};"
    res = make_query(query, theory_database)
    val = res.fetchone()
    if not val:
        raise TheoryNotFoundInDatabase(f"ID {theoryID} not found in database.")
    return dict([(k[0], v) for k, v in zip(res.description, val)])

def fetch_all(theory_database: Path):
    """Looks in the theory database and returns a dataframe with theory info
    for all theories
    """
    query = "SELECT * FROM TheoryIndex;"
    res = make_query(query, theory_database)
    val = res.fetchall()
    cols = [k[0] for k in res.description]
    df = pd.DataFrame(val, columns=cols)
    return df.set_index(['ID'])
