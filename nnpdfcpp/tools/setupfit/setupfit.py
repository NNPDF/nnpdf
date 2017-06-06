#!/usr/bin/env python
"""
setupfit.py
Given a fit on which filter is run, and a destination file, copy all the
files neccesary to run nnfit. The inputs are the result folder produced by
filter, and the destination file.

Zahari Kassabov
"""
__version__ = '0.1'
import pathlib
import argparse
import yaml
import shutil

to_copy = [
'{nnpdfcpp}/bin/nnfit',
'{resultpath}',
'{nnpdfcpp}/config/{fitname}.yml',
'{nnpdfcpp}/data/theory_{theoryid}',
'{nnpdfcpp}/data/theory.db',
'{nnpdfcpp}/data/commondata', #Needed???
]



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('filterres', help="Output path of the completed "
                    "filter was executed. It is assumed to be embedded in the "
                    "nnpdfsrc path labyrinth.")

    parser.add_argument('dest', help="Destination folder where the relevant "
    "nnpdfcpp directory structure is to be reproduced.")

    args = parser.parse_args()
    resultpath = pathlib.Path(args.filterres).absolute()
    nnpdfcpp = resultpath.absolute().parents[1]
    fitname = resultpath.name

    config_yaml = resultpath / 'filter.yml'
    with config_yaml.open() as f:
        config_data = yaml.load(f)
        theoryid = config_data['theory']['theoryid']

    destroot = pathlib.Path(args.dest)/'nnpdfcpp'
    destroot.mkdir()

    for loc in to_copy:
        inp = pathlib.Path(loc.format(**locals()))
        print(inp)
        rel_inp = inp.relative_to(nnpdfcpp)
        dest = destroot / rel_inp
        dest.parent.mkdir(parents=True, exist_ok=True)
        #TODO: remove str() for Python 3.6
        if inp.is_dir():
            shutil.copytree(str(inp), str(dest))
        else:
            shutil.copy2(str(inp),str(dest))

if __name__ == '__main__':
    main()

