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
import re

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

def load_commondata(spec):
    """
    Load the data corresponding to a CommonDataSpec object.

    Reads commondata file for dataset_name and returns a pandas DataFrame with:
    entry   process kin1    kin2    kin3    data    stat    \
            sys.add.0   sys.mult.0 .... sys.add.N   sys.mult.N
    """
    commondatafile = spec.datafile

    # Getting set name from commondata file name
    setname = re.search('DATA_(.*).dat', str(commondatafile)).group(1)
    tabledata = parse_commondata(commondatafile, setname)

    return tabledata

def parse_commondata(f, setname):
    
    """Parse a commondata file into a CommonData. Raise a BadCommondATAError
    if problems are encountered.
    Parameters
    ----------
    f : file
        Open file-like object. 
    Returns
    -------
    commondata : CommonData
        An object containing the data and information from the commondata file.
    """
    table = pd.read_csv(f, sep=r'\s+', skiprows=1, header=None)

    # remove NaNs
    # TODO: replace commondata files with bad formatting
    table.dropna(axis='columns', inplace=True)

    # build header
    header = ['entry', 'process', 'kin1', 'kin2', 'kin3', 'data', 'stat']
    nsys  = (table.shape[1]-len(header))//2
    for i in range(nsys):
        header += [f'sys.add.{i+1}', f'sys.mult.{i+1}']
    table.columns = header
    table.set_index('entry', inplace=True)

    # Populate CommonData object
    return CommonData(
                    setname = setname,
                    ndata = len(table),
                    commondataproc = table["process"][1],
                    nkin = 3 ,
                    nsys = nsys,
                    data = table)
