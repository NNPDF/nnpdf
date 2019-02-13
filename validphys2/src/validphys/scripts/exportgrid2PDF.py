#!/usr/bin/env python3
""" Wrapper script which allows one to run a validphys report where the `pdfs` label is filled
with fake PDFs which are created from export grids from a fit you have ran. You don't need to have
ran postfit on the fit, and this is primarily to produce validphys reports for different epochs
along the training.
"""
__author__ = "Michael Wilson"

import os
from pathlib import Path
import argparse
from math import sqrt

from logging import warning

import pandas as pd
import numpy as np

from reportengine.compat import yaml

import validphys.lhaindex as lhaindex
from validphys.lhio import write_replica, generate_replica0
from validphys.core import PDF
import NNPDF

NNPATH = Path(NNPDF.get_results_path())

FLID = {
    'TBAR':-6, 'BBAR':-5, 'CBAR':-4, 'SBAR':-3, 'UBAR':-2, 'DBAR':-1, 'GLUON':21,
    'D':1, 'U':2, 'S':3, 'C':4, 'B':5, 'T':6, 'PHT': 13}

def exportgridtopanda(gridfile):
    """Given an export file, creates a subgrids dataframe which can be written by vp to lhapdf"""
    with open(gridfile, "r") as f:
        filtermap = yaml.safe_load(f)
        # assume export grid has just one Q and then fill other as 10^5 since lhapdf needs 2 Qs
    qvals = [sqrt(filtermap['q20']), 1E5]
    valmat = np.array(filtermap['pdfgrid'])
    pdfvals = np.concatenate((valmat[:, 1:-2], valmat[:, 1:-2]), axis=1).flatten()
    #pdfvals = np.repeat(np.array(filtermap['pdfgrid']).flatten(), 2)
    xvals = filtermap['xgrid']
    flvals = [FLID[l] for l in filtermap['labels']]
    return pd.DataFrame(
        pdfvals, index=pd.MultiIndex.from_product(
            ([sqrt(filtermap['q20'])], xvals, qvals, flvals[1:-2])))


def move_info_file(infopath, lhapdf_path, nreps):
    """ Given the path to an info file, copies the file to new location. Modifying
    the number of replicas field (accounting for replica zero)

    So if your fit has replica 1 and 2 then give nrep=2 and the info file written will
    have the line NumMembers: 3 since lhapdf expects the zeroth replica to be counted.
    """
    ipath = Path(infopath)
    lhapdf_path = Path(lhapdf_path)
    #infoname = ipath.name
    if os.path.isdir(lhapdf_path):
        infoname = lhapdf_path.name + ".info"
        lhapdf_file = lhapdf_path / infoname
    else:
        lhapdf_file = lhapdf_path
    with open(lhapdf_file, "w+") as outfile:
        with open(infopath, "r") as infile:
            for line in infile:
                if "NumMembers" in line:
                    # Add one to account for replica zero
                    outfile.write(f"NumMembers: {nreps+1}\n")
                else:
                    outfile.write(f"{line}")

def clean_temporary(fitpath):
    """Given a fitpath, deletes the temporary folder for lhapdf objects
    based on assumption that there exists only fitpath/temp_lhapdf/<directory>/<file>
    or fitpath/temp_lhapdf/<file>
    """
    if os.path.isdir(fitpath/"temp_lhapdf"):
        for d in os.listdir(fitpath/"temp_lhapdf"):
            if os.path.isdir(fitpath/"temp_lhapdf"/d):
                for f in os.listdir(fitpath/"temp_lhapdf"/d):
                    os.remove(fitpath/"temp_lhapdf"/d/f)
                os.rmdir(fitpath/"temp_lhapdf"/d)
            else:
                os.remove(fitpath/"temp_lhapdf"/d)
        os.rmdir(fitpath/"temp_lhapdf")
    else:
        os.remove(fitpath/"temp_lhapdf")

def create_pdf(fitpath, gen):
    if isinstance(fitpath, str):
        fitpath = Path(fitpath)
    outpath = Path(fitpath/"temp_lhapdf"/f"{fitpath.name}_{gen}")
    # len() here will count the info file so take off 1
    nreps = len(os.listdir(fitpath/"nnfit")) - 1
    # count replicas
    i = 1
    os.mkdir(outpath)
    for f in os.listdir(fitpath/"nnfit"):
        if f.endswith("info"):
            # move_info_file adds 1 to nreps to account for rep zero
            move_info_file(fitpath/"nnfit"/f, outpath, nreps)
        else:
            a = exportgridtopanda(fitpath/f"nnfit/{f}/{fitpath.name}_{gen}.exportgrid")
            #get original rep number from folder name to store in header
            rep = [int(s) for s in f.split("_") if s.isdigit()][0]
            preamble = f"PdfType: replica\nFormat: lhagrid1\nFromMCReplica: {rep}\n"
            write_replica(i, outpath, str.encode(preamble), a)
            i += 1
    # symlink into the first lha_path folder
    os.symlink(f"{outpath.absolute()}", f"{lhaindex.get_lha_paths()[0]}/{outpath.name}")
    # Generate replica 0
    pdf = PDF(f"{fitpath.name}_{gen}")
    generate_replica0(pdf)

def create_fit(fitpath, gen):
    if isinstance(fitpath, str):
        fitpath = Path(fitpath)
    os.symlink(f"{fitpath.absolute()}", f"{NNPATH}/{fitpath.name}_{gen}")

def create_runcard(runcardpath, fitpath, names):
    """Given a runcard and a list of PDFs, create a temporary runcard with pdfs filled from list"""
    with open(runcardpath, "r") as f:
        filtermap = yaml.safe_load(f)
    outpath = Path.joinpath(fitpath, "temp_lhapdf", "runcard.yaml")
    filtermap["pdfs"] = names
    filtermap["fits"] = names
    with open(outpath, "w+") as f:
        yaml.safe_dump(filtermap, stream=f)

def main():
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fitpath', type=str, help="Path to a completed fit with valid .exportgrid "
                        "files, which will be used to create the PDFs used in report")
    parser.add_argument('runcardpath', type=str, help="Path to base validphys runcard "
                        "the `pdfs` field will be changed to be the "
                        "pdfs created in wrapper script")
    parser.add_argument('--interval', default='0', type=str,
                        help="generation interval either a single number "
                        " i.e 10 or a linspace style specification, colon seperated "
                        " i.e 0:100:1000 which would be generations 0 to 1000 in gaps of 100.")
    args, vpargs = parser.parse_known_args()
    vpargstring = ""
    for arg in vpargs:
        vpargstring += f"{arg} "

    fitpath = Path(args.fitpath)
    fitname = fitpath.name
    genspec_input = args.interval.split(":")
    genspec = list(map(int, genspec_input))
    # create list of generations from start:step:stop
    # accounting for range having final entry start-step
    try:
        generations = list(range(genspec[0], genspec[2]+genspec[1], genspec[1]))
    except IndexError:
        generations = genspec
    if (not os.path.isdir(fitpath)
            or not os.path.exists(fitpath/"filter.yml")
            or not os.path.exists(fitpath/ f"nnfit/replica_1/")
            or not all(
                [os.path.exists(
                    fitpath/f"nnfit/replica_1/{fitname}_{i}.exportgrid") for i in generations])):
        raise RuntimeError("first argument should point at a fit directory with at least 1 "
                           "completed replica with `<fitname>_<i>.exportgrid` files with i "
                           "specified by the interval given as 2nd argument (you gave "
                           f"{':'.join(genspec_input)}).")

    if len(generations) > 5:
        warning("More than 5 generations being plotted, this could look really bad.")

    if os.path.exists(fitpath/"temp_lhapdf"):
        warning("temporary folder already exists, overriding content")
        clean_temporary(fitpath)

    for i in generations:
        try:
            os.unlink(f"{lhaindex.get_lha_paths()[0]}/{fitname}_{i}")
            warning(f"pdf with name {fitname}_{i} already exists in LHAPDF folder, "
                    f"attempting to delete.")
        except FileNotFoundError:
            pass

        try:
            os.unlink(f"{NNPATH}/{fitname}_{i}")
            warning(f"fit with name {fitname}_{i} already exists in results folder, "
                    f"attempting to delete.")
        except FileNotFoundError:
            pass


    os.mkdir(fitpath/"temp_lhapdf")
    for i in generations:
        create_pdf(fitpath, i)
        create_fit(fitpath, i)
    names = [fitname + f"_{i}" for i in generations]
    create_runcard(Path(f"{args.runcardpath}"), fitpath, names)
    os.system(f'validphys {fitpath/"temp_lhapdf/runcard.yaml"} {vpargstring}')
    for i in generations:
        os.unlink(f"{lhaindex.get_lha_paths()[0]}/{fitname}_{i}")
        os.unlink(f"{NNPATH}/{fitname}_{i}")
    clean_temporary(fitpath)

if __name__ == "__main__":
    main()
