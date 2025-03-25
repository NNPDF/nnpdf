#!/usr/bin/env python
"""
    Generate the theory csv using validphys functions
"""
from argparse import ArgumentParser
from pathlib import Path

from nnpdf_data import theory_cards
from nnpdf_data.theorydbutils import fetch_all

if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument("csvpath", type=Path, help="Path to write the csv to")
    parser.add_argument("--zero-only", action="store_true", help="If enabled, only the _000 entry is printed out. Old theories are ignored")

    args = parser.parse_args()

    theory_df = fetch_all(theory_cards)

    # Enforce the following order in the table:
    order = ["PTO", "QED", "Comments", "IC", "Q0", "ModEv"]
    for i, c in enumerate(order):
        theory_df.insert(i, c, theory_df.pop(c))

    if args.zero_only:
        theory_df = theory_df[theory_df.index.astype(str).str.endswith("000") & (theory_df.index.astype(str).str.len() > 5)]
        # Drop entries with useless info...
        theory_df = theory_df[~theory_df["Comments"].str.startswith("Same as", na=False)]

    theory_df.to_csv(args.csvpath)
