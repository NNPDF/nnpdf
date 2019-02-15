#!/usr/bin/env python
import NNPDF as nnpath 
import argparse
import os
import pathlib, shutil

#Taking command line arguments
parser = argparse.ArgumentParser(description = 'Script to rename fits')
parser.add_argument('initial', help='Name of the fit to be changed')
parser.add_argument('final', help='Desired new name of fit')
parser.add_argument('-r', '--result_path', action='store_true', help='Use to change name of a fit in results path')
parser.add_argument('-c', '--copy', action='store_true', help='Use to create a copy of the original fit')
args = parser.parse_args()

initial_dir = os.path.abspath(args.initial)
if args.initial[-1] == '/':
    args.initial = args.initial[:-1] 
    pass 

initial_fit_name = os.path.basename(args.initial)
def change_name(initial, final):
    """Function that takes initial fit name and final fit name
    and performs the renaming"""
    os.chdir(os.path.abspath(initial))
    nnfit = pathlib.Path('nnfit') 
    if nnfit.exists():
        os.chdir('nnfit')
        #Change .info filename
        info_file = pathlib.Path(f'{initial_fit_name}.info')
        info_file.rename(f'{final}.info')

        #Change replica names
        for item in os.listdir():
            q = pathlib.Path(item)
            if q.is_dir():
                os.chdir(item)
                p = pathlib.Path('.')
                files = list(p.glob(initial_fit_name + '*'))
                for i in files:
                    i.rename(final + i.suffix)
                os.chdir('..')
        os.chdir('..')

    postfit = pathlib.Path('postfit')
    if postfit.exists():
        os.chdir('postfit')
        folder = pathlib.Path(initial)
        folder.rename(final)
        os.system('sed -i -e "s/{}/{}/g" postfit.log'.format(initial, final))

        #Change symlinks
        os.chdir(final)
        for item in os.listdir():
            p = pathlib.Path(item)
            if p.is_symlink():
                replica = p.resolve().parent.name
                pointer = f'../../nnfit/{replica}/{final}.dat' 
                p.unlink()
                p.symlink_to(pointer)
                p.rename(p.name.replace(initial, final))
            else:
                p.rename(p.name.replace(initial, final))
        os.chdir('../../..')

    pass

def main():
    if args.result_path:
        os.chdir(nnpath.get_results_path())
    else:
        os.chdir(initial_dir + '/..')
    if args.copy:
        shutil.copytree(initial_fit_name, initial_fit_name + 'copy', symlinks=True)

    copied_fit = pathlib.Path(initial_fit_name + 'copy')
    fit = pathlib.Path(initial_fit_name)
    if copied_fit.exists():
        change_name(initial_fit_name, args.final)
        fit.rename(args.final)
        copied_fit.rename(initial_fit_name)
    else:   
        change_name(initial_fit_name, args.final)
        fit.rename(args.final)
