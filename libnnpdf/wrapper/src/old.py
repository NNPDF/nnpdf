import warnings

warnings.warn(
    f"Module {__name__} is deprecated. Use NNPDF directly instead",
    #DeprecationWarning,
    stacklevel=2,
    )
from .nnpdf import *
