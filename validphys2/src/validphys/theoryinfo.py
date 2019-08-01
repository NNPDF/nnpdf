from pandas import DataFrame

from reportengine.table import table

from validphys.theorydbutils import fetch_theory, fetch_all

@table
def all_theory_info_table(theory_database):
    """Produces a DataFrame with all theory info and saves it"""
    return fetch_all(theory_database)

@table
def theory_info_table(theory_database, theoryID):
    """fetches theory info for given `theoryID` constructs DataFrame from it
    """
    res_dict = fetch_theory(theory_database, theoryID)
    res_df = DataFrame(
        list(res_dict.values()),
        index=res_dict.keys(),
        columns=[f'Info for theory {theoryID}'],
    )
    return res_df
