# -*- coding: utf-8 -*-
"""
theoryinfo.py

Actions for displaying theory info for one or more theories.
"""
from pandas import DataFrame

from reportengine.table import table

from validphys.theorydbutils import fetch_theory, fetch_all

@table
def all_theory_info_table(theory_database):
    """Produces a DataFrame with all theory info and saves it

    Returns
    -------
    all_theory_info_table: pd.Dataframe
        dataframe filled with all entries in theorydb file

    Example
    -------
    >>> from validphys.api import API
    >>> df = API.all_theory_info_table()
    >>> df['Comments'].iloc[:5]
    ID
    1                 3.0 LO benchmark
    2                3.0 NLO benchmark
    3               3.0 NNLO benchmark
    4     3.0 NLO - Q0=1.3 For IC Test
    5    3.0 NNLO - Q0=1.3 For IC Test
    Name: Comments, dtype: object
    """
    return fetch_all(theory_database)

@table
def theory_info_table(theory_database, theoryID):
    """fetches theory info for given `theoryID` constructs DataFrame from it

    Parameters
    ----------
    theoryID: int
        numeric identifier of theory to be queried. Can be specified at the
        runcard level.

        NOTE: this parameter differs from the runcard key `theoryid`
        (case sensitive) and doesn't require the user to have the theory locally

    Returns
    -------
    theory_info_table: pd.Dataframe
        dataframe filled with theory info for specified `theoryID`

    Example
    -------
    >>> from validphys.api import API
    >>> df = API.theory_info_table(theoryID=53)
    >>> df.loc['Comments']
    Info for theory 53    NNPDF3.1 NNLO central
    Name: Comments, dtype: object
    """
    res_dict = fetch_theory(theory_database, theoryID)
    res_df = DataFrame(
        list(res_dict.values()),
        index=res_dict.keys(),
        columns=[f'Info for theory {theoryID}'],
    )
    return res_df
