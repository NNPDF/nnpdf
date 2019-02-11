#!/usr/bin/env python
import NNPDF as nnpath 
import argparse
import os
import pathlib

#Taking command line arguments
parser = argparse.ArgumentParser(description = 'Script to rename fits')
group = parser.add_mutually_exclusive_group()
group.add_argument('-r', '--result_path', action='store_true', help='Use to change name of a fit in results path')
parser.add_argument('initial', help='Name of the fit to be changed')
parser.add_argument('final', help='Desired new name of fit')
args = parser.parse_args()

initial_dir = os.path.abspath(args.initial)
if args.initial[-1] == '/':
    args.initial = args.initial[:-1] 
    pass 

initial_fit_name = os.path.basename(args.initial)
def change_name(initial, final):
    """Function that takes initial fit name and final fit name
    and performs the renaming"""
    os.chdir(os.path.abspath(initial) + '/nnfit')

    #Change .info filename
    os.system('mv {}.info {}.info'.format(initial, final))
    #Change replica names
    for item in os.listdir():
        try:
            os.chdir(item)
            p = pathlib.Path('.')
            files = list(p.glob(initial_fit_name + '*'))
            for i in files:
                i.rename(final + i.suffix)
            os.chdir('..')
        except:
            pass
    os.chdir('../postfit')
    os.system('mv {} {}'.format(initial, final)) 
    os.system('sed -i -e "s/{}/{}/g" postfit.log'.format(initial, final))
    os.chdir('../../')
    pass

def main():
    if args.result_path:
        os.chdir(nnpath.get_results_path())
        change_name(initial_fit_name, args.final)
    else:
        os.chdir(initial_dir + '/..')
        change_name(initial_fit_name, args.final)

    os.system('mv {} {}'.format(initial_fit_name, args.final))
