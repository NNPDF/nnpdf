"""
Input validation and classes for the n3fit runcard

The parsers (declared at the end of the file as parse_) should be used.
"""

from dataclasses import dataclass
import logging
from typing import Optional

from validobj import ValidationError, parse_input
from validobj.custom import Parser

# Map the ``debug_options::log_level`` string from the runcard to python logging levels
LOG_LEVELS = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
}


@Parser
def ValidLogLevel(log_level: str) -> str:
    """Check whether the log level received is valid and, if it is, conver it to lowercase."""
    _valid_log_levels = list(LOG_LEVELS.keys())
    if log_level.lower() not in _valid_log_levels:
        raise ValidationError(
            f"Invalid debug_options::log_level '{log_level}'. "
            f"Available values are {_valid_log_levels}."
        )
    return log_level.lower()


@Parser
def ValidPrintEach(printeach: int) -> int:
    """Check whether the print_each value is an integer and greater than 1."""
    if printeach < 1 or not isinstance(printeach, int):
        raise ValidationError(f"Print each must be an integer > 0, received {printeach}")
    return printeach


@dataclass(frozen=True)
class DebugOptions:
    """Parsed ``debug_options`` namespace of the n3fit runcard.

    All debug-related options should live in this object, which
    should be passed around to the internal fitting routines.

    Attributes
    ----------
    log_level: str
        level of the n3fit logger
    printeach: int
        print the training/validation stats every ``printeach`` epochs
    timer: bool
        enable the per-epoch timing callback independently of ``debug``
    print_logs: bool
        force the per-epoch stats to be printed to stdout
    """

    log_level: Optional[ValidLogLevel] = "info"
    printeach: Optional[ValidPrintEach] = 100
    timer: Optional[bool] = False
    print_logs: Optional[bool] = False

    @property
    def print_summary(self) -> bool:
        """Whether the network summaries should be printed. Print only for info and debug.
        Note that hyperopt always supresses the network summaries.
        """
        return self.log_level in ("debug", "info")

    @property
    def logger_level(self):
        return LOG_LEVELS[self.log_level]


def parse_debug_options(debug_options=None):
    if debug_options is None:
        return DebugOptions()
    return parse_input(debug_options, DebugOptions)
