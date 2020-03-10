#!/usr/bin/env python
"""
    vp-pdffromfit - command line tool to obtain an LHAPDF from a NNPDF fit

    To obtain the PDF from an fit, simply run
    vp-pdffromfit <path-to-fit> <PDF name>. Optional flags allow for the
    resulting pdf to be placed in the LHAPDF directory, as well as modifying
    various fields of the info file. In addition, it is possible to compress
    the resulting PDF also using tar archiving.
"""

import argparse
import logging
import os
import pathlib
import subprocess
import shutil
import sys
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
    parser = argparse.ArgumentParser(
        description="Script to obtain an LHAPDF grid from a fit"
    )
    parser.add_argument("Fit", help="Path to fit")
    parser.add_argument("PDF", help="Name of the desired PDF set output")
    parser.add_argument(
        "--author",
        action="append",
        help='''The author to be added to the PDF .info file.
                Apply this argument multiple times for multiple authors,
                quotation marks can be used for an author name containing several words,
                e.g "The NNPDF collaboration"''',
    )
    parser.add_argument(
        "--description",
        help="""The set description to be added to the PDF .info file.
                Quotations should be used for this field.""",
    )
    parser.add_argument(
        "--data-version", help="The data version to be added to the PDF .info file."
    )
    parser.add_argument(
        "--index", help="The set index to be added to the PDF .info file."
    )
    parser.add_argument(
        "--reference",
        help="The reference to be added to the PDF .info file, usually an arXiv reference.",
    )
    parser.add_argument(
        "-l",
        "--lhapdf_path",
        action="store_true",
        help="Place the output LHAPDF in the LHAPDF directory.",
    )
    parser.add_argument(
        "-c", "--compress", action="store_true", help="Compress the resulting fit."
    )
    args = parser.parse_args()
    return args


def fixup_ref(pdf_path: pathlib.Path, field_dict):
    """Function to handle alterations of the info file.
    The argparser namespace is read in as a dictionary which
    is then used to write to the resulting output file.

    If the user did not provide a field then we revert to the
    pre existing field.
    """
    pdf_name = pdf_path.name
    infopath = pdf_path / f"postfit/{pdf_name}/{pdf_name}.info"
    with open(infopath) as f:
        y = yaml.YAML()
        res = y.load(infopath)
        # If a field entry is not provided, then we revert to the existing
        # field in pre-existing info file.
        res["Authors"] = (
            field_dict.get("author")
            if field_dict.get("author") is not None
            else res["Authors"]
        )
        res["SetDesc"] = (
            field_dict.get("description")
            if field_dict.get("description") is not None
            else res["SetDesc"]
        )
        res["DataVersion"] = (
            field_dict.get("data-version")
            if field_dict.get("data-version") is not None
            else res["DataVersion"]
        )
        res["SetIndex"] = (
            field_dict.get("index")
            if field_dict.get("index") is not None
            else res["SetIndex"]
        )
        res["Reference"] = (
            field_dict.get("reference")
            if field_dict.get("reference") is not None
            else res["Reference"]
        )

    with open(infopath, "w") as f:
        y.default_flow_style = True
        y.dump(res, f)


def postfit_path(path: pathlib.Path) -> pathlib.Path:
    pdf_name = path.name
    return pathlib.Path(path / f"postfit/{pdf_name}")


def compress(lhapdf_path: pathlib.Path):
    """ Function to compress the resulting PDF. Dereferences are handled
    in order to account for possible symbolic linking of grids.
    """
    output = lhapdf_path.name + ".tar.gz"
    with tarfile.open(output, "w:gz", dereference=True) as tar:
        tar.add(str(lhapdf_path), arcname=os.path.basename(str(lhapdf_path)))


def main():
    args = process_args()
    fit_path, pdf_name = pathlib.Path(args.Fit).resolve(), args.PDF

    if not fit_path.is_dir():
        log.error(f"Could not find fit. Path '{fit_path.absolute()}' is not a directory.")
        sys.exit(1)
    if not (fit_path/'filter.yml').exists():
        log.error(f"Path {fit_path.absolute()} does not appear to be a fit. "
                  "File 'filter.yml' not found in the directory")
        sys.exit(1)

    with tempfile.TemporaryDirectory(dir=fit_path.parent) as tmp:
        tmp = pathlib.Path(tmp)
        copied_fit = tmp / fit_path.name
        shutil.copytree(fit_path, copied_fit)

        fixup_ref(copied_fit, vars(args))

        new_path = change_name(copied_fit, pdf_name)
        lhapdf_path = postfit_path(new_path)

        if args.lhapdf_path:
            dest_path = pathlib.Path(paths()[-1]) / pdf_name
        else:
            dest_path = new_path.parent.with_name(pdf_name)

        if dest_path.exists():
            log.error(f"Destination path {dest_path.absolute()} already exists.")
            sys.exit(1)

        lhapdf_path = lhapdf_path.rename(dest_path)
        log.info(f"PDF generated and placed in {lhapdf_path.parent}")

        if args.compress:
            log.info("Compressing output")
            compress(lhapdf_path)
