#!/usr/bin/env python
"""
A script which utilises the validphys loader method `check_theoryinfo` to allow
user to quickly print table to terminal.

By default the table is just printed however the table can also be saved to
a CSV file with the name theory_<theoryid>_info.csv by using the command line
option `-d`, which stands for `--dumptable`.
"""

__authors__ = 'Michael Wilson, Zahari Kassabov'

import argparse
import logging
import pathlib
import sys

from reportengine import colors
from reportengine.table import savetable

from validphys.loader import Loader
from validphys.theorydbutils import TheoryNotFoundInDatabase
from validphys.theoryinfo import theory_info_table

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())

DBPATH = Loader().theorydb_file

def main():
    parser = argparse.ArgumentParser(description='Script to check theory info')
    parser.add_argument(
        'theoryid',
        help=(
            "Numeric identifier of theory to look up info of"
        ),
        type=int
    )
    parser.add_argument(
        '--dumptable',
        '-d',
        help=(
            "Boolen flag which causes table to be dumped to CSV"
            " file in current working directory, with name "
            "theory_<theoryid>_info.csv"
        ),
        action='store_true',
    )
    args = parser.parse_args()

    try:
        df = theory_info_table(DBPATH ,args.theoryid)
    except TheoryNotFoundInDatabase as e:
        log.error(e)
        sys.exit(1)
    print(df)
    if args.dumptable:
        outpath = pathlib.Path(f"theory_{args.theoryid}_info.csv")
        if outpath.exists():
            log.error(
                f"The file `theory_{args.theoryid}_info.csv` already exists in "
                "your current working directory."
            )
            sys.exit(1)
        log.info(f"Saving info table to theory_{args.theoryid}_info.csv")
        savetable(df, outpath)


if __name__ == "__main__":
    main()
