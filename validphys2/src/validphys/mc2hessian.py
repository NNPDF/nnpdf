"""
mc2hessian.py

This module containts the functionality to compute reduced set using the `mc2hessian` algorithm 
(See section 2.1 of of `1505.06736 <https://arxiv.org/pdf/1602.00005.pdf#subsection.2.1>`_).
"""

import logging

import numpy as np

from validphys.lhio import hessian_from_lincomb
from validphys.pdfgrids import xplotting_grid
from validphys.lhio import install_mc2hessian_grids


log = logging.getLogger(__name__)


def gridname(pdf, Neig, mc2hname=None):
    """If no custom `mc2hname' is specified, the name of the Hessian PDF is automatically generated.
    """
    if mc2hname is None:
        grid_name = f"{pdf.name}_hessian_{Neig}"
    else:
        grid_name = mc2hname
    return grid_name


def mc2hessian(pdf, Neig, output_path, Q, gridname, installgrids=False):
    """Produces a Hessian PDF by transfroming a Monte Carlo PDF set.

    Noteworthy args:
        Neig (float): Number of basis eigenvectors
        Q (float): Energy scale
        gridname (str): Name of the output Hessian PDF
        installgrids (bool, optional): If True, the Hessian PDF is installed in the LHAPDF 
            directory. Defaults to False.
    """
    gridpaths = []
    result = _create_mc2hessian(pdf, Q=Q, Neig=Neig, output_path=output_path, name=gridname)
    gridpaths.append(result)
    if installgrids:
        install_mc2hessian_grids(gridname, output_path)


def _create_mc2hessian(pdf, Q, Neig, output_path, name=None):
    X = _get_X(pdf, Q, reshape=True)
    vec = _compress_X(X, Neig)
    norm = _pdf_normalization(pdf)
    return hessian_from_lincomb(pdf, vec / norm, folder=output_path, set_name=name)


def _get_X(pdf, Q, reshape=False, xgrid=None):
    if xgrid is None:
        xgrid = _make_xgrid()
    # Allow tuples that can be saved in cache
    elif isinstance(xgrid, tuple):
        xgrid = _make_xgrid(*xgrid)
    pdf_grid = xplotting_grid(pdf, Q, xgrid=xgrid)
    pdf_grid_values = pdf_grid.grid_values
    replicas = pdf_grid_values
    mean = pdf_grid_values.mean(axis=0)
    Xt = replicas - mean
    if reshape:
        Xt = Xt.reshape(Xt.shape[0], Xt.shape[1] * Xt.shape[2])
    return Xt.T


def _make_xgrid(xminlog=1e-5, xminlin=1e-1, xmax=1, nplog=50, nplin=50):
    """Provides the points in x to sample the PDF. `logspace` and `linspace`
    will be called with the respsctive parameters."""
    return np.append(
        np.logspace(np.log10(xminlog), np.log10(xminlin), num=nplog, endpoint=False),
        np.linspace(xminlin, xmax, num=nplin, endpoint=False),
    )


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
