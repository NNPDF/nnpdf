"""
    Rename commands for use with vp-tools. Functions for renaming
    PDFs and fits are included with documentation within the argument
    help.
"""
import os
import pathlib
import shutil
import sys
import tarfile
import tempfile
import logging

import click

import lhapdf

import NNPDF as nnpath

from reportengine import colors
from reportengine.compat import yaml

from validphys.renametools import change_name, rename_pdf

log = logging.getLogger()
log.setLevel(logging.DEBUG)
log.addHandler(colors.ColorHandler())


@click.group()
def cli():
    """Utility for renaming fits or pdfs
    """
    pass


@cli.command()
@click.argument("src")
@click.argument("dst")
@click.option(
    "-r",
    "--result_path",
    "result_path",
    help="Flag to change the name of a fit in results path",
    is_flag=True,
)
@click.option(
    "-c",
    "--copy",
    "copy",
    help="Flag to create a copy of the original fit",
    is_flag=True,
)
def fit(src, dst, result_path, copy):
    """Rename a fit from SRC to DST
    """
    initial_dir = pathlib.Path(src)
    initial_fit_name = initial_dir.name
    if result_path:
        if len(initial_dir.parts) != 1:
            log.error("Enter a fit name and not a path with the -r option")
            sys.exit(1)
        fitpath = pathlib.Path(nnpath.get_results_path()) / initial_fit_name
    else:
        fitpath = initial_dir
    if not fitpath.is_dir():
        log.error(
            f"Could not find fit. Path '{fitpath.absolute()}' is not a directory."
        )
        sys.exit(1)
    if not (fitpath / "filter.yml").exists():
        log.error(
            f"Path {fitpath.absolute()} does not appear to be a fit. "
            "File 'filter.yml' not found in the directory"
        )
        sys.exit(1)

    dest = fitpath.with_name(dst)
    if dest.exists():
        log.error(f"Destination path {dest.absolute()} already exists.")
        sys.exit(1)
    with tempfile.TemporaryDirectory(dir=fitpath.parent) as tmp:
        tmp = pathlib.Path(tmp)
        copied_fit = tmp / initial_fit_name
        shutil.copytree(fitpath, copied_fit, symlinks=True)
        newpath = change_name(copied_fit, dst)
        newpath.rename(dest)
        if copy:
            log.info("Renaming completed with copy")
        else:
            shutil.rmtree(fitpath)
            log.info("Renaming completed")


@cli.command()
@click.argument("src")
@click.argument("dst")
@click.option(
    "--author",
    "author",
    help='''The author to be added to the
            PDF .info file. Apply this argument multiple times for multiple authors,
            quotation marks can be used for an author name containing several words, e.g
            "The NNPDF collaboration"''',
    type=str,
    multiple=True,
)
@click.option(
    "--description",
    "description",
    help="""The set description to be added to the PDF .info file.
            Quotations should be used for this field.""",
    type=str,
)
@click.option(
    "--data-version",
    "data_version",
    help="The data version to be added to the PDF .info file.",
    type=str,
)
@click.option(
    "--index",
    "index",
    help="The set index to be added to the PDF .info file.",
    type=str,
)
@click.option(
    "--reference",
    "reference",
    help="The reference to be added to the PDF .info file, usually an arXiv reference.",
    type=str,
)
@click.option(
    "--lhapdf-path",
    "lhapdf_path",
    help="Place the output LHAPDF in the LHAPDF directory.",
    is_flag=True,
)
@click.option("--compress", "compress", help="Compress the resulting PDF.", is_flag=True)
def pdf(src, dst, author, description, data_version, index, reference, lhapdf_path, compress):
    """Rename a LHAPDF from SRC to DST
    """
    # We need to use the os method to avoid resolving the symlink
    # when we use pathlib.Path.resolve()
    field_dict = dict(author=author, description=description, data_version=data_version, index=index, reference=reference)
    source_path = pathlib.Path(os.path.abspath(src))
    pdf_name = dst

    if lhapdf_path:
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

    with tempfile.TemporaryDirectory(dir=dest_path.parent) as tmp:
        tmp = pathlib.Path(tmp)
        copied_fit = tmp / source_path.name
        shutil.copytree(source_path, copied_fit)

        fixup_ref(copied_fit, field_dict)

        rename_pdf(copied_fit, source_path.name, pdf_name)

        lhapdf_path = copied_fit.with_name(pdf_name)
        lhapdf_path.rename(dest_path)
        log.info(f"PDF generated and placed in {dest_path.parent}")

    if compress:
        from validphys.renametools import Spinner

        log.info("Compressing output")
        with Spinner():
            compress_pdf(dest_path)


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

    with open(infopath, "w") as f:
        y.default_flow_style = True
        y.dump(res, f)


def compress_pdf(lhapdf_path: pathlib.Path):
    """ Function to compress the resulting PDF. Dereferences are handled
    in order to account for possible symbolic linking of grids.
    """
    output = lhapdf_path.name + ".tar.gz"
    with tarfile.open(output, "w:gz", dereference=True) as tar:
        tar.add(str(lhapdf_path), arcname=os.path.basename(str(lhapdf_path)))
