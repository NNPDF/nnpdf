# -*- coding: utf-8 -*-
"""
theorydbutils.py

low level utilities for querying the theory database file and representing the
data as a python object.
"""
from functools import lru_cache
from pathlib import Path

import pandas as pd

from .theory import TheoryCard
from .utils import parse_yaml_inp


class TheoryNotFoundInDatabase(Exception):
    pass


@lru_cache
def parse_theory_card(theory_card):
    """Read the theory card using validobj parsing
    Returns the theory as a dictionary
    """
    if theory_card.exists():
        tcard = parse_yaml_inp(theory_card, TheoryCard)
        return tcard.asdict()
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
    >>> from nnpdf_data import theory_cards
    >>> from nnpdf_data.theorydbutils import fetch_theory
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
    >>> from nnpdf_data.theorydbutils import fetch_all
    >>> theory_df = fetch_all(theory_cards)
    """
    theories = []
    for theory_path in theory_database.glob("*.yaml"):
        theories.append(parse_theory_card(theory_path))
    df = pd.DataFrame(theories)
    return df.set_index(['ID']).sort_index()
