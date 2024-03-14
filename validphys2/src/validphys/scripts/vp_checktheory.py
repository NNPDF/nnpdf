#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A script which utilises the validphys loader method `check_theoryinfo` to allow
user to quickly print table to terminal.

Example
-------

$ vp-checktheory 53

will print the info of theory 53 to terminal

$ vp-checktheory 53 --dumptable

in addition to printing the theory info table will save the table to csv in
the current working directory with the name theory_<theoryid>_info.csv

$ vp-checktheory --fit NNPDF31_nlo_as_0118

will parse the `theoryid` from supplied fit, downloading the fit if required.

Notes
-----
User can either specify `theoryid` as a positional argument or `--fit FIT` but
not both at the same time.

"""

__authors__ = 'Michael Wilson, Zahari Kassabov'

import argparse
import logging
import pathlib
import sys

from reportengine import colors
from reportengine.table import savetable

from validphys.loader import FallbackLoader
from validphys.theorydbutils import TheoryNotFoundInDatabase
from validphys.theoryinfo import theory_info_table

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())

loader_log = logging.getLogger('validphys.loader')
loader_log.setLevel(logging.INFO)
loader_log.addHandler(colors.ColorHandler())


LOADER = FallbackLoader()
DBPATH = LOADER.theorydb_folder


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        'theoryid',
        nargs='?',
        default=None,
        help=(
            "Numeric identifier of theory to look up info of"
        ),
        type=int
    )
    group.add_argument(
        '--fit',
        type=str,
        help=(
            "Name of a fit from which to parse `theoryid` from, instead of "
            "supplying theoryid on command line"
        ),
        default=None
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
    # get theoryid from command line or --fit
    if args.fit is not None:
        fitspec = LOADER.check_fit(args.fit)
        theoryid = fitspec.as_input()['theory']['theoryid']
        log.info(f"Fit used theory {theoryid}")
    else:
        theoryid = args.theoryid

    try:
        df = theory_info_table(DBPATH, theoryid)
    except TheoryNotFoundInDatabase as e:
        log.error(e)
        sys.exit(1)
    print(df)
    if args.dumptable:
        outpath = pathlib.Path(f"theory_{theoryid}_info.csv")
        if outpath.exists():
            log.error(
                f"The file `theory_{theoryid}_info.csv` already exists in "
                "your current working directory."
            )
            sys.exit(1)
        log.info(f"Saving info table to theory_{theoryid}_info.csv")
        savetable(df, outpath)


if __name__ == "__main__":
    main()
