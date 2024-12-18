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
def _parse_theory_card(theory_card):
    """Read the theory card using validobj parsing
    Returns the theory as a dictionary
    """
    tcard = parse_yaml_inp(theory_card, TheoryCard)
    return tcard.asdict()


@lru_cache
def _get_available_theory_cards(path):
    """Since theoryfile names may contain underscores, the theoryfile name cannot uniquely be
    determined from the theoryID. Therefore we create a mapping between the theoryid as python
    int and the corresponding file.
    """
    available_theories = {}
    for theoryfile in path.glob("*.yaml"):
        tmp_id = int(theoryfile.stem)
        if tmp_id in available_theories:
            another = available_theories[tmp_id]
            raise ValueError(f"Two theory files with same id: {theoryfile} and {another}")
        available_theories[tmp_id] = theoryfile
    return available_theories


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

    available_theories = _get_available_theory_cards(theory_database)
    theoryID = int(theoryID)
    try:
        theoryfile = available_theories[theoryID]
    except KeyError as e:
        raise TheoryNotFoundInDatabase(f"Theorycard for theory not found: {e}")

    tdict = _parse_theory_card(theoryfile)
    if tdict["ID"] != theoryID:
        raise ValueError(f"The theory ID in {theoryfile} doesn't correspond with its ID entry")
    return tdict


def fetch_all(theory_database: Path):
    """Looks in the theory database and returns a dataframe with theory info
    for all theories that can be found in the folder ``theory_database``.
    Every ``.yaml`` file found in the folder will be considered a possible theory.

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
    >>> from nnpdf_data import theory_cards
    >>> from nnpdf_data.theorydbutils import fetch_all
    >>> theory_df = fetch_all(theory_cards)
    """
    available_theories = _get_available_theory_cards(theory_database)
    theories = [fetch_theory(theory_database, thid) for thid in available_theories]
    df = pd.DataFrame(theories)
    return df.set_index(['ID']).sort_index()
