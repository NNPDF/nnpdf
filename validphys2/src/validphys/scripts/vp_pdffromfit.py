#!/usr/bin/env python

import argparse
import logging
import os
import pathlib
import subprocess
import shutil
import tarfile
import tempfile

import NNPDF

from lhapdf import paths

from reportengine.compat import yaml

from validphys.loader import FallbackLoader as Loader
from validphys.scripts.vp_fitrename import change_name

log = logging.getLogger()

# Taking command line arguments
def process_args():
    parser = argparse.ArgumentParser(description='Script to obtain an LHAPDF grid from a fit')
    parser.add_argument('Fit', help='Path to fit')
    parser.add_argument('PDF', help='Name of the desired PDF set output')
    parser.add_argument(
        '--reference',
        help='The reference to be added to the PDF .info file, usually an arXiv reference.')
    parser.add_argument(
        '-n',
        '--nnpdf_path',
        action='store_true',
        help='Use a PDF fit stored in the NNPDF results directory.')
    parser.add_argument(
        '-l',
        '--lhapdf_path',
        action='store_true',
        help='Place the output LHAPDF in the LHAPDF directory.')
    parser.add_argument(
        '-c',
        '--compress',
        action='store_true',
        help='Compress the resulting fit.')
    args = parser.parse_args()
    return args

def fixup_ref(pdf_path: pathlib.Path, reference: str):
    pdf_name = pdf_path.name
    infopath = pdf_path/f"postfit/{pdf_name}/{pdf_name}.info"
    with open(infopath) as f:
        y = yaml.YAML()
        res = y.load(infopath)
        res["Reference"] = reference
    with open(infopath, "w") as f:
        y.dump(res, f)

def postfit_path(path: pathlib.Path) -> pathlib.Path:
    pdf_name = path.name
    return pathlib.Path(path/f"postfit/{pdf_name}")

def compress(lhapdf_path: pathlib.Path):
    output = lhapdf_path.name + ".tar.gz"
    with tarfile.open(output, "w:gz", dereference=True) as tar:
        tar.add(str(lhapdf_path), arcname=os.path.basename(str(lhapdf_path)))

def main():
    args = process_args()
    fit_path, pdf_name = pathlib.Path(args.Fit).resolve(), args.PDF
    reference = args.reference

    with tempfile.TemporaryDirectory(dir=fit_path.parent) as tmp:
        tmp = pathlib.Path(tmp)
        copied_fit = tmp/fit_path.name
        shutil.copytree(fit_path, copied_fit)

        if reference is not None:
            fixup_ref(copied_fit, reference)
            log.info(f"Reference {reference} added to info file")

        new_path = change_name(copied_fit, pdf_name)
        lhapdf_path = postfit_path(new_path)

        if args.lhapdf_path:
            lhapdf_path = lhapdf_path.rename(pathlib.Path(paths()[-1])/pdf_name)
            log.info(f"PDF generated and placed in {paths()[-1]}")
        else:
            lhapdf_path = lhapdf_path.rename(new_path.parent.with_name(pdf_name))
            log.info(f"PDF generated and placed in {os.getcwd()}")

        if args.compress:
            log.info("Compressing output")
            compress(lhapdf_path)

    return 0
