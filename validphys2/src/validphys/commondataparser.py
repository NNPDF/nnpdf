"""
This module implements parsers for commondata  and systype files into useful
datastructures, contained in the :py:mod:`validphys.coredata` module, which are
not backed by C++ managed memory, and so they can be easily pickled and
interfaces with common Python libraries.  The integration of these objects into
the codebase is currently work in progress, and at the moment this module
serves as a proof of concept.
"""
import io
import functools
import tarfile
import dataclasses

import numpy as np
import pandas as pd

from validphys.coredata import CommonData 

class BadCommonDataError(Exception):
    """Exception raised when a commondata file cannot be parsed correctly"""

@dataclasses.dataclass(frozen=True)
class CommonDataInfo:
    """Class containing the basic properties of a commondata file."""
    setname: str
    ndata: int
    proc: str
    nsys: int

def load_dataset(datafile):
    """Reads commondata file for dataset_name and returns a panda DataFrame with:
        entry   process kin1    kin2    kin3    data    stat    \
            sys.add.0   sys.mult.0 .... sys.add.N   sys.mult.N
    """
    # read raw commondata file
    table = pd.read_csv(datafile, sep=r'\s+', skiprows=1, header=None)

    # remove NaNs
    # TODO: replace commondata files with bad formatting
    table.dropna(axis='columns', inplace=True)

    # build header
    header = ['entry', 'process', 'kin1', 'kin2', 'kin3', 'data', 'stat']
    for i in range((table.shape[1]-len(header))//2):
        header += [f'sys.add.{i+1}', f'sys.mult.{i+1}']
    table.columns = header
    table.set_index('entry', inplace=True)
    return table