# -*- coding: utf-8 -*-
"""
theorydbutils.py

low level utilities for querying the theory database file and representing the
data as a python object.
"""
from pathlib import Path

# Keep the sqlite3 dependence in one location
import sqlite3

import pandas as pd


class TheoryNotFoundInDatabase(Exception):
    pass


def make_query(query: str, dbpath: Path):
    """Base level function which executes a `query` given the path to a sqlite3
    database

    Parameters
    ----------
    query: str
        valid sqlite3 database query
    dbpath: Path
        pathlib.Path object pointing to sqlite3 database

    Returns
    -------
    res: sqlite3.Cursor
        an sqlite cursor with relevant data saved

    Example
    -------
    >>> from pathlib import Path as p
    >>> from validphys.theorydbutils import make_query
    >>> query = "SELECT * FROM TheoryIndex WHERE ID=53;"
    >>> dbpath = p("./validphys2/src/validphys/datafiles/theory.db")
    >>> res = make_query(query, dbpath)
    >>> val = res.fetchone()
    >>> theory_info_dict = {k[0]: v for k, v in zip(res.description, val)}

    See Also
    --------
    sqlite3: https://docs.python.org/3/library/sqlite3.html
    """
    # python 3.6 still requires path as string, but 3.7 onwards can take
    # pathlib.Path. This could be changed if we stop supporting 3.6
    conn = sqlite3.connect(str(dbpath))
    with conn:
        cursor = conn.cursor()
        res = cursor.execute(query)
    return res


def fetch_theory(theory_database: Path, theoryID: int):
    """Looks in the theory database and returns a dictionary of theory info for the
    theory number specified by `theoryID`.

    Parameters
    ----------
    theory_database: Path
        pathlib.Path pointing to theorydb.db
    theoryID: int
        numeric identifier of theory to query info

    Returns
    -------
    theory_info_dict: dict
        dictionary filled with relevant entry from theory database
    """
    # int casting is intentional to avoid malformed querys.
    query = f"SELECT * FROM TheoryIndex WHERE ID={int(theoryID)};"
    res = make_query(query, theory_database)
    val = res.fetchone()
    if not val:
        raise TheoryNotFoundInDatabase(f"ID {theoryID} not found in database.")
    return {k[0]: v for k, v in zip(res.description, val)}


def fetch_all(theory_database: Path):
    """Looks in the theory database and returns a dataframe with theory info
    for all theories

    Parameters
    ----------
    theory_database: Path
        pathlib.Path pointing to theorydb.db

    Returns
    -------
    theory_info_dataframe: pd.Dataframe
        dataframe filled with all entries in theorydb file
    """
    query = "SELECT * FROM TheoryIndex;"
    res = make_query(query, theory_database)
    val = res.fetchall()
    cols = [k[0] for k in res.description]
    df = pd.DataFrame(val, columns=cols)
    return df.set_index(['ID'])
