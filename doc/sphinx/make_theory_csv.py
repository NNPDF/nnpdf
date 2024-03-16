#!/usr/bin/env python
"""
    Generate the theory csv using validphys functions
"""
from argparse import ArgumentParser
from pathlib import Path

from validphys.datafiles import theory_cards
from validphys.theorydbutils import fetch_all

if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument("csvpath", type=Path, help="Path to write the csv to")

    args = parser.parse_args()

    theory_df = fetch_all(theory_cards)
    theory_df.to_csv(args.csvpath)
