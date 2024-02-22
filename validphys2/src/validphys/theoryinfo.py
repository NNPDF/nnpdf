# -*- coding: utf-8 -*-
"""
theoryinfo.py

Actions for displaying theory info for one or more theories.
"""
from pandas import DataFrame

from reportengine.table import table
from validphys.theorydbutils import fetch_all, fetch_theory


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
def theory_info_table(theory_database, theory_db_id):
    """fetches theory info for given `theory_db_id` constructs DataFrame from it

    Parameters
    ----------
    theory_db_id: int
        numeric identifier of theory to be queried. Can be specified at the
        runcard level.

    Returns
    -------
    theory_info_table: pd.Dataframe
        dataframe filled with theory info for specified `theory_db_id`

    Example
    -------
    >>> from validphys.api import API
    >>> df = API.theory_info_table(theory_db_id=53)
    >>> df.loc['Comments']
    Info for theory 53    NNPDF3.1 NNLO central
    Name: Comments, dtype: object
    """
    res_dict = fetch_theory(theory_database, theory_db_id)
    res_df = DataFrame(
        list(res_dict.values()), index=res_dict.keys(), columns=[f'Info for theory {theory_db_id}']
    )
    return res_df
