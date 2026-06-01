"""
PDF interpolation backend selector.

The active backend is controlled by the ``NNPDF_PDF_BACKEND`` environment
variable (default: ``"lhapdf"``).  Set it to ``"neopdf"`` to use the
NeoPDF interpolation library instead.

Both backends expose the same public interface (``LHAPDFSet`` / ``NeoPDFSet``)
so they are interchangeable throughout validphys and n3fit.
"""

import os

_BACKEND_ENV_VAR = "NNPDF_PDF_BACKEND"
_VALID_BACKENDS = ("lhapdf", "neopdf")


class InvalidPDFBackend(Exception):
    pass


def _active_backend() -> str:
    backend = os.environ.get(_BACKEND_ENV_VAR, "lhapdf").lower()

    if backend not in _VALID_BACKENDS:
        raise InvalidPDFBackend(
            f"Unknown PDF backend {backend!r} in {_BACKEND_ENV_VAR}. "
            f"Valid options are: {_VALID_BACKENDS}"
        )
    return backend


def make_pdfset(name: str, error_type: str, *, backend: str | None = None):
    """Return a PDF set object for *name* using the requested backend.

    Parameters
    ----------
    name:
        LHAPDF set name (both backends resolve sets by this name).
    error_type:
        One of ``replicas``, ``symmhessian``, ``hessian``, ``t0``
    backend:
        ``lhapdf`` or ``neopdf``.  When *None* the value of the
        ``NNPDF_PDF_BACKEND`` environment variable is used (default
        ``lhapdf``).

    Returns
    -------
    ``LHAPDFSet`` or ``NeoPDFSet`` depending on *backend*.
    """
    active = backend if backend is not None else _active_backend()

    if active == "neopdf":
        from validphys.neopdfset import NeoPDFSet

        return NeoPDFSet(name, error_type)

    from validphys.lhapdfset import LHAPDFSet

    return LHAPDFSet(name, error_type)
