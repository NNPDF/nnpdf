# -*- coding: utf-8 -*-
"""
theorydbutils.py

low level utilities for querying the theory database file and representing the
data as a python object.
"""
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path

import pandas as pd

from validphys.utils import parse_yaml_inp


@dataclass(frozen=True)
class _TheoryCard:
    ID: int
    PTO: int
    FNS: str
    DAMP: int
    IC: int
    ModEv: str
    XIR: float
    XIF: float
    NfFF: int
    MaxNfAs: int
    MaxNfPdf: int
    Q0: float
    alphas: float
    Qref: float
    QED: int
    alphaqed: float
    Qedref: float
    SxRes: int
    SxOrd: str
    HQ: str
    mc: float
    Qmc: float
    kcThr: float
    mb: float
    Qmb: float
    kbThr: float
    mt: float
    Qmt: float
    ktThr: float
    CKM: list[float]
    MZ: float
    MW: float
    GF: float
    SIN2TW: float
    TMC: int
    MP: float
    Comments: str
    global_nx: int
    EScaleVar: int


class TheoryNotFoundInDatabase(Exception):
    pass


@lru_cache
def parse_theory_card(theory_card):
    """Read the theory card using validobj parsing
    Returns the theory as a dictionary
    """
    if theory_card.exists():
        return asdict(parse_yaml_inp(theory_card, _TheoryCard))
    raise TheoryNotFoundInDatabase(f"Theory card {theory_card} not found")


def fetch_theory(theory_database: Path, theoryID: int):
    """Looks in the theory card folder and returns a dictionary of theory info for the
    theory number specified by `theoryID`.

    Parameters
    ----------
    theory_database: Path
        pathlib.Path pointing to the folder with the theory cards
    theoryID: int
        numeric identifier of theory to query info

    Returns
    -------
    theory_info_dict: dict
        dictionary filled with relevant entry from theory database

    Example
    ------
    >>> from validphys.datafiles import theory_cards
    >>> from validphys.theorydbutils import fetch_theory
    >>> theory = fetch_theory(theory_cards, 700)
    """
    filepath = theory_database / f"{theoryID}.yaml"
    tdict = parse_theory_card(filepath)
    if tdict["ID"] != int(theoryID):
        raise ValueError(f"The theory ID in {filepath} doesn't correspond with its ID entry")
    return tdict


def fetch_all(theory_database: Path):
    """Looks in the theory database and returns a dataframe with theory info
    for all theories

    Parameters
    ----------
    theory_database: Path
        pathlib.Path pointing to the folder with theory cards

    Returns
    -------
    theory_info_dataframe: pd.Dataframe
        dataframe filled with all entries in theorydb file

    Example
    ------
    >>> from validphys.datafiles import theory_cards
    >>> from validphys.theorydbutils import fetch_all
    >>> theory_df = fetch_all(theory_cards)
    """
    theories = []
    for theory_path in theory_database.glob("*.yaml"):
        theories.append(parse_theory_card(theory_path))
    df = pd.DataFrame(theories)
    return df.set_index(['ID']).sort_index()
