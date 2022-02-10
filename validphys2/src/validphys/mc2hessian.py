"""
mc2hessian.py

This module containts the functionality to compute reduced set using the `mc2hessian` algorithm
(See section 2.1 of of `1602.00005 <https://arxiv.org/pdf/1602.00005.pdf#subsection.2.1>`_).
"""
import logging
import numbers
import pathlib
import shutil

import numpy as np

from reportengine.checks import check, check_positive, make_argcheck

from validphys import lhaindex
from validphys.lhio import hessian_from_lincomb
from validphys.pdfgrids import xplotting_grid

from validphys.checks import check_pdf_is_montecarlo

log = logging.getLogger(__name__)


def gridname(pdf, Neig, mc2hname: (str, type(None)) = None):
    """If no custom `mc2hname' is specified, the name of the Hessian PDF is automatically generated.
    """
    if mc2hname is None:
        grid_name = f"{pdf.name}_hessian_{Neig}"
    else:
        grid_name = mc2hname
    return grid_name


@make_argcheck
def _check_xminmax(xmin, xminlin, xmax):
    check(0 < xmin < xminlin < xmax <= 1, "Expecting 0 < xmin < xminlin < xmax <= 1")


@_check_xminmax
@check_positive("nplog")
@check_positive("nplin")
def mc2hessian_xgrid(
    xmin: float = 1e-5,
    xminlin: float = 1e-1,
    xmax: numbers.Real = 1,
    nplog: int = 50,
    nplin: int = 50,
):
    """Provides the points in x to sample the PDF. `logspace` and `linspace`
    will be called with the respsctive parameters.

    Generates a grid with ``nplog`` logarithmically spaced points between
    ``xmin`` and ``xminlin`` followed by ``nplin`` linearly spaced points
    between ``xminlin`` and ``xmax``
    """
    return np.append(
        np.geomspace(xmin, xminlin, num=nplog, endpoint=False),
        np.linspace(xminlin, xmax, num=nplin, endpoint=False),
    )


@check_pdf_is_montecarlo
def mc2hessian(
    pdf, Q, Neig: int, mc2hessian_xgrid, output_path, gridname, installgrid: bool = False
):
    """Produces a Hessian PDF by transfroming a Monte Carlo PDF set.

    Parameters
    -----------
    pdf : validphys.core.PDF
        An existng validphys PDF object which will be converted into a Hessian PDF set
    Q : float
        Energy scale at which the Monte Carlo PDF is sampled
    Neig : int
        Number of basis eigenvectors in the Hessian PDF set
    mc2hessian_xgrid : numpy.ndarray
        The points in x at which to sample the Monte Carlo PDF set
    output path : pathlib.PosixPath
        The validphys output path where the PDF will be written
    gridname : str
        Name of the Hessian PDF set
    installgrid : bool, optional, default=``False``
        Whether to copyt the Hessian grid to the LHAPDF path
    """
    result_path = _create_mc2hessian(
        pdf, Q=Q, xgrid=mc2hessian_xgrid, Neig=Neig, output_path=output_path, name=gridname
    )
    if installgrid:
        lhafolder = pathlib.Path(lhaindex.get_lha_datapath())
        dest = lhafolder / gridname
        if lhaindex.isinstalled(gridname):
            log.warning(
                "Target directory for new PDF, %s, already exists. " "Removing contents.",
                dest,
            )
            if dest.is_dir():
                shutil.rmtree(str(dest))
            else:
                dest.unlink()
        shutil.copytree(result_path, dest)
        log.info("Hessian PDF set installed at %s", dest)


def _create_mc2hessian(pdf, Q, xgrid, Neig, output_path, name=None):
    X = _get_X(pdf, Q, xgrid, reshape=True)
    vec = _compress_X(X, Neig)
    norm = _pdf_normalization(pdf)
    return hessian_from_lincomb(pdf, vec / norm, folder=output_path, set_name=name)


def _get_X(pdf, Q, xgrid, reshape=False):
    pdf_grid = xplotting_grid(pdf, Q, xgrid=xgrid)
    pdf_grid_values = pdf_grid.stats_gv
    replicas = pdf_grid_values.error_members()
    mean = pdf_grid_values.central_value()
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
    get the correct errors, depending on the `error_type` of `pdf`."""
    nrep = len(pdf) - 1
    if pdf.error_type == "replicas":
        norm = np.sqrt(nrep - 1)
    elif pdf.error_type in ("hessian", "symmhessian"):
        norm = 1
    else:
        raise NotImplementedError(
            "This PDF error type is not supported." "PDF error: %s" % pdf.error_type
        )
    return norm
