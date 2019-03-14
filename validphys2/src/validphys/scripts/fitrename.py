#!/usr/bin/env python
"""
        fitrename - command line tool to rename fits

        fitrename allows for command line renaming of fits.
        To do so, call fitrename with a path to the original fit and
        the requested new fit. The optional flags allow for preserving
        a copy of the original fit and also to change fits in the NNPDF
        results directory.
"""

__authors__ = 'Shayan Iranipour, Zahari Kassabov, Michael Wilson'

import argparse
import os
import pathlib
import shutil
import sys
import tempfile
import logging

import NNPDF as nnpath
from reportengine import colors

log = logging.getLogger()
log.setLevel(logging.DEBUG)
log.addHandler(colors.ColorHandler())

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


def rename_pdf(pdf_folder, initial_fit_name, final_name):
    for item in os.listdir(pdf_folder):
        p = pdf_folder/item
        if p.is_symlink():
            replica = p.resolve().parent.name
            pointer = f'../../nnfit/{replica}/{final_name}.dat'
            p.unlink()
            p.symlink_to(pointer)
        newname = p.name.replace(initial_fit_name, final_name)
        p.rename(p.with_name(newname))
    pdf_folder.rename(pdf_folder.with_name(final_name))


def rename_nnfit(nnfit_path, initial_fit_name, final_name):
    info_file = nnfit_path/f'{initial_fit_name}.info'
    info_file.rename(info_file.with_name(f'{final_name}.info'))
    #Some older fits have the PDF here
    pdf_folder = nnfit_path / initial_fit_name
    if pdf_folder.is_dir():
        rename_pdf(pdf_folder, initial_fit_name, final_name)
    #Change replica names
    for item in nnfit_path.glob('replica*'):
        if item.is_dir():
            files = item.glob(initial_fit_name + '*')
            for i in files:
                newname = i.name.replace(initial_fit_name, final_name)
                i.rename(item/newname)


def rename_postfit(postfit_path, initial_fit_name, final_name):
    pdf_folder = postfit_path / initial_fit_name
    rename_pdf(pdf_folder, initial_fit_name, final_name)
    os.system(f'sed -i -e "s/{initial_fit_name}/{final_name}/g" {postfit_path/"postfit.log"}')

def change_name(initial_path, final_name):
    """Function that takes initial fit name and final fit name
    and performs the renaming"""
    initial_fit_name = initial_path.name
    nnfit = initial_path/'nnfit'
    if nnfit.exists():
        rename_nnfit(nnfit, initial_fit_name, final_name)
    postfit = initial_path/'postfit'
    if postfit.exists():
        rename_postfit(postfit, initial_fit_name, final_name)
    newpath = initial_path.with_name(final_name)
    initial_path.rename(newpath)
    return newpath

def main():
    args = process_args()
    if args.final[-1] == '/':
        args.final = args.final[:-1]

    initial_dir = pathlib.Path(args.initial)
    initial_fit_name = initial_dir.name
    if args.result_path:
        if len(initial_dir.parts) != 1:
            log.error("Enter a fit name and not a path with the -r option")
            sys.exit(1)
        fitpath = pathlib.Path(nnpath.get_results_path())/initial_fit_name
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
