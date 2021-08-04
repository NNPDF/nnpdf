"""
mc2hessian.py

This module containts the functionality to compute reduced set using the `mc2hessian` 
(Appendix of `1505.06736 <http://arxiv.org/abs/1505.06736>`_) algorithm.
"""

import logging
import pathlib
import shutil

import numpy as np

from validphys import lhaindex
from validphys.lhio import hessian_from_lincomb
from validphys.pdfgrids import xplotting_grid

log = logging.getLogger(__name__)


def mc2hessian(pdf, Neig, output_path, Q, gridname, db=None, installgrids=False):
    """Produces a Hessian PDF by transfroming a Monte Carlo PDF set.

    Noteworthy args:
        Neig (float): Number of basis eigenvectors
        Q (float): Energy scale
        gridname (str): Name of the output Hessian PDF
        installgrids (bool, optional): If True, the Hessian PDF is installed in the LHAPDF 
            directory. Defaults to False.
    """
    gridpaths = []
    result = create_mc2hessian(pdf, Q=Q, Neig=Neig, output_path=output_path, name=gridname, db=db)
    gridpaths.append(result)
    if installgrids:
        install_grids(gridname, output_path)


def create_mc2hessian(pdf, Q, Neig, output_path, name=None, db=None):
    output_path = output_path
    X = get_X(pdf, Q, reshape=True)
    vec = compress_X(X, Neig)
    norm = _pdf_normalization(pdf)
    return hessian_from_lincomb(pdf, vec / norm, folder=output_path, set_name=name, db=db)


def install_grids(gridname, output_path):
    lhafolder = pathlib.Path(lhaindex.get_lha_datapath())
    dest = lhafolder / gridname
    if lhaindex.isinstalled(gridname):
        log.warning(
            "Target directory for new PDF, %s, already exists. " "Overwriting contents.", gridname
        )
        if dest.is_dir():
            shutil.rmtree(str(dest))
        else:
            dest.unlink()
    source = output_path / gridname
    shutil.copytree(source, dest)
    log.info(f"Hessian PDF {gridname} succesfully written to the folder {dest}")


def gridname(pdf, Neig, mc2hname=None):
    """If no custom `mc2hname' is specified, the name of the Hessian PDF is automatically generated.
    """
    if mc2hname is None:
        grid_name = f"{pdf.name}_hessian_{Neig}"
    else:
        grid_name = mc2hname
    return grid_name


def get_X(pdf, Q, reshape=False, xgrid=None):
    if xgrid is None:
        xgrid = make_xgrid()
    # Allow tuples that can be saved in cache
    elif isinstance(xgrid, tuple):
        xgrid = make_xgrid(*xgrid)
    pdfGrid = xplotting_grid(pdf, Q, xgrid=xgrid)
    pdfGrid_values = pdfGrid.grid_values
    replicas = pdfGrid_values
    mean = pdfGrid_values.mean(axis=0)
    Xt = replicas - mean
    if reshape:
        Xt = Xt.reshape(Xt.shape[0], Xt.shape[1] * Xt.shape[2])
    return Xt.T


def make_xgrid(xminlog=1e-5, xminlin=1e-1, xmax=1, nplog=50, nplin=50):
    """Provides the points in x to sample the PDF. `logspace` and `linspace`
    will be called with the respsctive parameters."""
    return np.append(
        np.logspace(np.log10(xminlog), np.log10(xminlin), num=nplog, endpoint=False),
        np.linspace(xminlin, xmax, num=nplin, endpoint=False),
    )


def compress_X(X, neig):
    _U, _S, V = np.linalg.svd(X, full_matrices=False)
    vec = V[:neig, :].T
    return vec


def _pdf_normalization(pdf):
    """Extract the quantity by which we have to divide the eigenvectors to
    get the correct errors, depending on the `ErrorType` of `pdf`."""
    nrep = len(pdf) - 1
    if pdf.ErrorType == "replicas":
        norm = np.sqrt(nrep - 1)
    elif pdf.ErrorType in ("hessian", "symmhessian"):
        norm = 1
    else:
        raise NotImplementedError(
            "This PDF error type is not supported." "PDF error: %s" % pdf.ErrorType
        )
    return norm
