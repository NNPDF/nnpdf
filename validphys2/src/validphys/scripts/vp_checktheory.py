#!/usr/bin/env python
"""
A script which utilises the validphys loader method `check_theoryinfo` to allow
user to quickly print table to terminal, leveraging the panda dataframe pretty
printing
"""

__authors__ = 'Michael Wilson, Zahari Kassabov'

import argparse

from pandas import DataFrame

from validphys.loader import Loader


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
    args = parser.parse_args()
    df = theory_info_table(args.theoryid)
    print(df)
