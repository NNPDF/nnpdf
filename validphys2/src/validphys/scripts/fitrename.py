#!/usr/bin/env python
"""
        fitrename - tool to rename fits

        fitrename is a script which renames fits from the command line.
        The path of the initial fit and requested final fit name are provided
        as command line arguments. Optional flags can be used to request a
        copy of the original fit or to rename a fit that is in the 
        NNPDF results path.
"""

__authors__ = 'Shayan Iranipour, Zahari Kassabov, Michael Wilson'

import argparse
import os
import pathlib
import shutil
import tempfile
import sys

from reportengine import colors
import NNPDF as nnpath

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

class FitExistsError(FileExistsError): pass
class FitNotFoundError(FileNotFoundError): pass
class InputError(SyntaxError): pass

def main():
    args = process_args()
    initial_dir = pathlib.Path(args.initial)
    initial_fit_name = initial_dir.name
    if args.result_path:
        try:
            if len(initial_dir.parts) != 1:
                raise InputError("Enter a fit name and not a path with the -r option")
        except InputError as e:
            print(colors.color_exception(e.__class__, e, e.__traceback__))
            sys.exit(1)
        fitpath = pathlib.Path(nnpath.get_results_path())/initial_fit_name
    else:
        fitpath = initial_dir
    try:
        if not fitpath.is_dir():
            raise FitNotFoundError(f"Could not find fit. Path '{fitpath.absolute()}' is not a directory.")
    except FitNotFoundError as e:
        print(colors.color_exception(e.__class__, e, e.__traceback__))
        sys.exit(1)
    try:
        if not (fitpath/'filter.yml').exists():
            raise FitNotFoundError(f"Path {fitpath.absolute()} does not appear to be a fit. "
                    "File 'filter.yml' not found in the directory")
    except FitNotFoundError as e:
        print(colors.color_exception(e.__class__, e, e.__traceback__))
        sys.exit(1)

    dest = fitpath.with_name(args.final)
    if dest.exists():
        raise FitExistsError(f"Destination path {dest.absolute()} already exists.")
    if args.copy:
        with tempfile.TemporaryDirectory(dir=fitpath.parent) as tmp:
            tmp = pathlib.Path(tmp)
            copied_fit = tmp/initial_fit_name
            shutil.copytree(fitpath, copied_fit, symlinks=True)
            newpath = change_name(copied_fit, args.final)
            newpath.rename(dest)

    else:
        change_name(fitpath, args.final)
