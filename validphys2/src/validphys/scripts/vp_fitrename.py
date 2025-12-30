#!/usr/bin/env python
"""
        vp-fitrename - command line tool to rename fits

        vp-fitrename allows for command line renaming of fits.
        To do so, call vp-fitrename with a path to the original fit and
        the requested new fit. The optional flags allow for preserving
        a copy of the original fit and also to change fits in the NNPDF
        results directory.
"""

__authors__ = 'Shayan Iranipour, Zahari Kassabov, Michael Wilson'

import argparse
import pathlib
import shutil
import sys
import tempfile
import logging

from reportengine import colors

from validphys.renametools import change_name
from validphys.loader import Loader


#Taking command line arguments
def process_args():
    parser = argparse.ArgumentParser(description='Script to rename fits')
    parser.add_argument('initial', help='Name of the fit to be changed')
    parser.add_argument('final', help='Desired new name of fit')
    parser.add_argument(
        '-r',
        '--result_path',
        action='store_true',
        help='Use to change name of a fit in results path')
    parser.add_argument(
        '-c',
        '--copy',
        action='store_true',
        help='Use to create a copy of the original fit')
    args = parser.parse_args()
    return args


def main():
    log = logging.getLogger()
    log.setLevel(logging.DEBUG)
    log.addHandler(colors.ColorHandler())

    args = process_args()
    if args.final[-1] == '/':
        args.final = args.final[:-1]

    initial_dir = pathlib.Path(args.initial)
    initial_fit_name = initial_dir.name
    if args.result_path:
        if len(initial_dir.parts) != 1:
            log.error("Enter a fit name and not a path with the -r option")
            sys.exit(1)
        fitpath = Loader().resultspath / initial_fit_name
    else:
        fitpath = initial_dir
    if not fitpath.is_dir():
        log.error(f"Could not find fit. Path '{fitpath.absolute()}' is not a directory.")
        sys.exit(1)
    if not (fitpath/'filter.yml').exists():
        log.error(f"Path {fitpath.absolute()} does not appear to be a fit. "
                  "File 'filter.yml' not found in the directory")
        sys.exit(1)

    dest = fitpath.with_name(args.final)
    if dest.exists():
        log.error(f"Destination path {dest.absolute()} already exists.")
        sys.exit(1)
    with tempfile.TemporaryDirectory(dir=fitpath.parent) as tmp:
        tmp = pathlib.Path(tmp)
        copied_fit = tmp/initial_fit_name
        shutil.copytree(fitpath, copied_fit, symlinks=True)
        newpath = change_name(copied_fit, args.final)
        newpath.rename(dest)
        if args.copy:
            log.info("Renaming completed with copy")
        else:
            shutil.rmtree(fitpath)
            log.info("Renaming completed")
