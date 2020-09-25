#!/usr/bin/env python
"""
    vp-pdfrename - command line tool to rename LHAPDFs

    To obtain the PDF from an fit, simply run
    vp-pdfrename <path-to-fit> <PDF name>. Optional flags allow for the
    resulting pdf to be placed in the LHAPDF directory, as well as modifying
    various fields of the info file. You can use this script to subsample the
    number of replicas, provided the fit uses Monte Carlo replicas. In addition,
    it is possible to compress the resulting PDF also using tar archiving.
"""

import argparse
import logging
import os
import pathlib
import random
import shutil
import sys
import tarfile
import tempfile

import lhapdf

from reportengine import colors
from reportengine.compat import yaml

from validphys.renametools import rename_pdf

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())


def check_none_or_gt_one(value):
    if value is None:
        return value
    try:
        ivalue = int(value)
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            f"{value} cannot be interpreted as an integer."
        ) from e
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"{value} is an invalid positive int value.")
    return ivalue


# Taking command line arguments
def process_args():
    parser = argparse.ArgumentParser(description="Script to rename LHAPDFs")
    parser.add_argument("Source", help="Path to source PDF")
    parser.add_argument("Destination", help="Name of the desired PDF set output")
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
        "--data-version",
        type=int,
        help="The data version to be added to the PDF .info file.",
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
        "--replicas",
        type=check_none_or_gt_one,
        default=None,
        help=(
            "Number of replicas to keep, replicas will be randomly subsampled "
            "from source PDF. The number of kept replicas should be less than "
            "the number of replicas of the source PDF. Also the source PDF must "
            "use MC replicas."
        ),
    )
    parser.add_argument(
        "-c", "--compress", action="store_true", help="Compress the resulting PDF."
    )
    args = parser.parse_args()
    if args.replicas is not None and not args.lhapdf_path:
        parser.error(
            "Subsampling replicas requires --lhapdf_path, in order to generate replica 0"
        )
    return args


def fixup_ref(pdf_path: pathlib.Path, field_dict):
    """Function to handle alterations of the info file.
    The argparser namespace is read in as a dictionary which
    is then used to write to the resulting output file.

    If the user did not provide a field then we revert to the
    pre existing field.
    """
    pdf_name = pdf_path.name
    infopath = pdf_path / f"{pdf_name}.info"

    with open(infopath) as f:
        y = yaml.YAML()
        res = y.load(f)

    # If a field entry is not provided, then we revert to the existing
    # field in pre-existing info file.
    if field_dict["author"]:  # Note: bool(None) is False
        res["Authors"] = field_dict["author"]

    if field_dict["description"]:
        res["SetDesc"] = field_dict["description"]

    if field_dict["data_version"]:
        res["DataVersion"] = field_dict["data_version"]

    if field_dict["index"]:
        res["SetIndex"] = field_dict["index"]

    if field_dict["reference"]:
        res["Reference"] = field_dict["reference"]

    if field_dict["replicas"]:
        # add one for replica zero (lhapdf convention) and one if replicas is 1
        # see ``subsample_replicas`` for reason why
        res["NumMembers"] = (
            field_dict["replicas"] + 1 + int(field_dict["replicas"] == 1)
        )

    with open(infopath, "w") as f:
        y.default_flow_style = True
        y.dump(res, f)


def subsample_replicas(copied_fit, n_reps):
    """
    Randomly keep ``n_reps`` of the replicas from the fit.

    This function
    does not assume that the copied_fit is in the lhapdf path, so does
    not attempt to regenerate replica 0.

    It also doesn't fixup the .info file with the new number of replicas.

    """
    replicas = [
        replica
        for replica in os.listdir(copied_fit)
        if (replica.endswith(".dat") and not replica.endswith("0000.dat"))
    ]

    if len(replicas) < n_reps:
        log.error(
            f"Too many replicas requested: {n_reps}, source PDF has {len(replicas)}."
        )
        sys.exit(1)

    random.shuffle(replicas)

    for replica in replicas[:-n_reps]:
        os.remove(copied_fit / replica)

    for i, replica in enumerate(replicas[-n_reps:], start=1):
        new_name = replica[:-8] + f"{i}".zfill(4) + ".dat"
        os.rename(copied_fit / replica, copied_fit / new_name)

    # lhapdf requires 2 replicas, so if n_reps is 1 we must copy replica 1
    if n_reps == 1:
        log.warning(
            "PDFs in the LHAPDF format are required to have 2 replicas, copying "
            "replica 1 to replica 2"
        )
        shutil.copyfile(copied_fit / new_name, copied_fit / f"{replica[:-8]}0002.dat")


def compress(lhapdf_path: pathlib.Path):
    """ Function to compress the resulting PDF. Dereferences are handled
    in order to account for possible symbolic linking of grids.
    """
    output = lhapdf_path.name + ".tar.gz"
    with tarfile.open(output, "w:gz", dereference=True) as tar:
        tar.add(str(lhapdf_path), arcname=os.path.basename(str(lhapdf_path)))


def main():
    args = process_args()
    # We need to use the os method to avoid resolving the symlink
    # when we use pathlib.Path.resolve()
    source_path = pathlib.Path(os.path.abspath(args.Source))
    pdf_name = args.Destination

    if args.lhapdf_path:
        dest_path = pathlib.Path(lhapdf.paths()[-1]) / pdf_name
    else:
        dest_path = source_path.with_name(pdf_name)

    if dest_path.exists():
        log.error(f"Destination path {dest_path.absolute()} already exists.")
        sys.exit(1)

    if not source_path.is_dir():
        log.error(
            f"Could not find fit. Path '{source_path.absolute()}' is not a directory."
        )
        sys.exit(1)

    if args.replicas is not None:
        infopath = source_path / f"{source_path.name}.info"

        with open(infopath) as f:
            # y = yaml.YAML()
            # res = y.load(f)
            info_peek = yaml.safe_load(f)
        if info_peek["ErrorType"] != "replicas":
            log.error(
                "%s does not have ErrorType: replicas, and so replicas "
                "cannot be subsampled.",
                source_path.name
            )
            sys.exit(1)

    with tempfile.TemporaryDirectory(dir=dest_path.parent) as tmp:
        tmp = pathlib.Path(tmp)
        copied_fit = tmp / source_path.name
        shutil.copytree(source_path, copied_fit)
        if args.replicas is not None:
            subsample_replicas(copied_fit, args.replicas)

        fixup_ref(copied_fit, vars(args))

        rename_pdf(copied_fit, source_path.name, pdf_name)

        lhapdf_path = copied_fit.with_name(pdf_name)
        lhapdf_path.rename(dest_path)
        if args.replicas is not None:
            # these imports significantly slowdown script.
            from validphys.lhio import generate_replica0
            from validphys.loader import Loader

            l = Loader()
            pdf = l.check_pdf(pdf_name)
            generate_replica0(pdf)

        log.info(f"PDF generated and placed in {dest_path.parent}")

    if args.compress:
        from validphys.renametools import Spinner

        log.info("Compressing output")
        with Spinner():
            compress(dest_path)
