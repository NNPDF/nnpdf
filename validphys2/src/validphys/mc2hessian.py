"""
mc2hessian.py

This module containts the functionality to compute reduced set using the `mc2hessian` algorithm 
(See section 2.1 of of `1505.06736 <https://arxiv.org/pdf/1602.00005.pdf#subsection.2.1>`_).
"""

import logging

import numpy as np

from validphys.lhio import hessian_from_lincomb, install_mc2hessian_grids
from validphys.pdfgrids import xplotting_grid

log = logging.getLogger(__name__)


def gridname(pdf, Neig, mc2hname=None):
    """If no custom `mc2hname' is specified, the name of the Hessian PDF is automatically generated.
    """
    if mc2hname is None:
        grid_name = f"{pdf.name}_hessian_{Neig}"
    else:
        grid_name = mc2hname
    return grid_name


def mc2hessian_xgrid(xmin=1e-5, xminlin=1e-1, xmax=1, nplog=50, nplin=50):
    """Provides the points in x to sample the PDF. `logspace` and `linspace`
    will be called with the respsctive parameters."""
    return np.append(
        np.geomspace(xmin, xminlin, num=nplog, endpoint=False),
        np.linspace(xminlin, xmax, num=nplin, endpoint=False),
    )


def mc2hessian(
    pdf, Q, mc2hessian_xgrid, Neig: int, output_path, gridname, installgrid: bool = True
):
    """Produces a Hessian PDF by transfroming a Monte Carlo PDF set.

    Noteworthy args:
        Neig (float): Number of basis eigenvectors
        Q (float): Energy scale
        gridname (str): Name of the output Hessian PDF
        installgrids (bool, optional): If True, the Hessian PDF is installed in the LHAPDF 
            directory. Defaults to False.
    """
    gridpaths = []
    result = _create_mc2hessian(
        pdf, Q=Q, xgrid=mc2hessian_xgrid, Neig=Neig, output_path=output_path, name=gridname
    )
    gridpaths.append(result)
    if installgrid:
        install_mc2hessian_grids(gridname, output_path)


def _create_mc2hessian(pdf, Q, xgrid, Neig, output_path, name=None):
    X = _get_X(pdf, Q, xgrid, reshape=True)
    vec = _compress_X(X, Neig)
    norm = _pdf_normalization(pdf)
    return hessian_from_lincomb(pdf, vec / norm, folder=output_path, set_name=name)


def _get_X(pdf, Q, xgrid, reshape=False):
    pdf_grid = xplotting_grid(pdf, Q, xgrid=xgrid)
    pdf_grid_values = pdf_grid.grid_values
    replicas = pdf_grid_values
    mean = pdf_grid_values.mean(axis=0)
    Xt = replicas - mean
    if reshape:
        Xt = Xt.reshape(Xt.shape[0], Xt.shape[1] * Xt.shape[2])
    return Xt.T


def _compress_X(X, neig):
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
