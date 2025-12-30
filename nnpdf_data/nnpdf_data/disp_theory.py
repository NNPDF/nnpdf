#!/usr/bin/env python

from pathlib import Path

from tabulate import tabulate
from yaml import safe_load

theory_cards = Path("theory_cards")

if __name__ == "__main__":

    theory_desc = {}
    for theory in theory_cards.glob("*.yaml"):
        tdict = safe_load(theory.read_text())
        tid = tdict["ID"]
        tdesc = tdict["Comments"]
        theory_desc[tid] = tdesc

    head = ["ID", "Description"]
    print(tabulate(sorted(theory_desc.items()), headers=head))
