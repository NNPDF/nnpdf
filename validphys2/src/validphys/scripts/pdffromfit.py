#!/usr/bin/env python

import argparse
import pathlib
import subprocess
import shutil
import tempfile

import NNPDF

from lhapdf import paths

from reportengine.compat import yaml

from validphys.loader import FallbackLoader as Loader
from validphys.scripts.fitrename import change_name

# Taking command line arguments
def process_args():
    parser = argparse.ArgumentParser(description='Script to obtain an LHAPDF grid from a fit')
    parser.add_argument('Fit', help='Path to fit')
    parser.add_argument('PDF', help='Name of the desired PDF set output')
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
    args = parser.parse_args()
    return args

def fixup_ref(new):
    l = Loader()
    p = l.check_pdf(new)
    fit = l.check_fit(new)
    desc = fit.as_input()["description"]
    infopath = pathlib.Path(p.infopath)
    with open(infopath) as f:
        y = yaml.YAML()
        res = y.load(infopath)
        res["SetDesc"] = desc
        res["Reference"] = "arxiv:1802.03398"
    with open(infopath, "w") as f:
        y.dump(res, f)


def postfit_path(path: pathlib.Path) -> pathlib.Path:
    pdf_name = path.name
    return pathlib.Path(path.resolve()/f"postfit/{pdf_name}")

def compress(new):
    fixup_ref(new)
    l = Loader()
    p = l.check_pdf(new)
    dst = pathlib.Path(p.infopath).parent
    subprocess.run(
        ["tar", "--dereference", "-czvf", f"res/{new}.tar.gz", "-C", str(dst.parent), new]
        , check = True
    )
    # shutil.make_archive(f"res/{new}", "gztar", root_dir=dst.parent, base_dir=new)


def main():
    args = process_args()
    fit_path, pdf_name = pathlib.Path(args.Fit).resolve(), args.PDF

    with tempfile.TemporaryDirectory(dir=fit_path.parent) as tmp:
        tmp = pathlib.Path(tmp)
        copied_fit = tmp/fit_path.name
        shutil.copytree(fit_path, copied_fit)
        new_path = change_name(copied_fit, pdf_name)
        lhapdf_path = postfit_path(new_path)
        lhapdf_path.rename(new_path.parent.with_name(pdf_name))

    return 0
