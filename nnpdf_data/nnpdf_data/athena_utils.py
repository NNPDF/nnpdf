import pandas as pd
from pathlib import Path

def read_excel(path_xlsx: Path, beams: tuple) -> pd.DataFrame:
    """Parse the xlsx file containing all the information regarding
    the projections and returns the ones corresponding to the chosen
    beam energies.

    Parameters
    ----------
    path_xlsx : Path
        path to the xlsx file
    beams: tuple
        tuple specifying the beam energies of the lepton & proton

    Returns
    -------
    pd.DataFrame
        returns a panda table corresponding to the chosen beams
    """
    xdf = pd.read_excel(path_xlsx)
    el, ep = beams
    df_beam = xdf[(xdf["El"] == el) & (xdf["Eh"] == ep)]
    return df_beam[df_beam["y"] <= 1]
