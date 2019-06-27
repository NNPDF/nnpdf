#!/usr/bin/env python
"""
A script which utilises the validphys loader method `check_theoryinfo` to allow
user to quickly print table to terminal, leveraging the panda dataframe pretty
printing.

By default the table is just printed however the table can also be saved by
specifying an output folder with the command line option `-o`.
"""

__authors__ = 'Michael Wilson, Zahari Kassabov'

import argparse
import logging
import pathlib

from pandas import DataFrame

from reportengine import colors
from reportengine.table import savetable

from validphys.loader import Loader

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())

def theory_info_table(theoryid):
    res_dict = Loader().check_theoryinfo(theoryid)
    res_df = DataFrame(
        list(res_dict.values()),
        index=res_dict.keys(),
        columns=[f'Info for theory {theoryid}'])
    return res_df

def main():
    parser = argparse.ArgumentParser(description='Script to check theory info')
    parser.add_argument(
        'theoryid', help='Numeric identifier of theory to look up info of', type=int)
    parser.add_argument('-o','--output',
                        help="Output folder for theory info table",
                        default=None)
    args = parser.parse_args()
    df = theory_info_table(args.theoryid)

    if args.output is not None:
        outpath = pathlib.Path(args.output)
        log.info(
            f"Saving info table to {outpath.absolute()}/theory_{args.theoryid}"
            "_info.csv")
        if outpath.is_dir():
            log.warning(
                f"Saving info table to {outpath.absolute()}/theory_"
                f"{args.theoryid}_info.csv already exists, attempting to override")
        try:
            outpath.mkdir(exist_ok=True)
        except FileExistsError:
            log.error(
                "A file already exists with the same name as the output folder "
                "please delete this file or choose a different name"
            )
            raise
        savetable(df, outpath/f"theory_{args.theoryid}_info.csv")

    print(df)
