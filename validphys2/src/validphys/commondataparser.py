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

from collections import namedtuple

from validphys.coredata import CommonData, SystypeData

class BadCommonDataError(Exception):
    """Exception raised when a commondata file cannot be parsed correctly"""
class BadSystypeError(Exception):
    """Exception raised when a systype file cannot be parsed correctly"""

CommondataTables = namedtuple(
    "CommondataTables", ("commondata_table", "systype_table")
)
def load_commondata(spec):
    """
    Load the data corresponding to a CommonDataSpec object.
    
    Returns an instance of the namedtuple CommondataTables,
    with:
    
    commondata_table being a CommonData instance with data arranged like
    entry   process kin1    kin2    kin3    data    stat    \
            sys.add.0   sys.mult.0 .... sys.add.N   sys.mult.N ;

    systype_table being a SystypeData instance with data arranged like
    sys_index   treatment   description.
    """
    commondatafile = spec.datafile

    # Getting set name from commondata file name
    setname = re.search('DATA_(.*).dat', str(commondatafile)).group(1)
    commondata = parse_commondata(commondatafile, setname)

    systypefile = spec.sysfile
    systypedata = parse_systype(systypefile, setname)
     
    return CommondataTables(
        commondata_table=commondata, systype_table=systypedata
    )

def parse_commondata(f, setname):
    
    """Parse a commondata file into a CommonData. Raise a BadCommondDataError
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
    try:
        table = pd.read_csv(f, sep=r'\s+', skiprows=1, header=None)
    except Exception as e:
        raise BadCommonDataError(f"Could not read file {f}. Please
    check there is a valid COMMONDATA file at this location.") from e
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


def parse_systype(f, setname):
    
    """Parse a systype file into a SystypeData. Raise a BadSystypeDataError
    if problems are encountered. 
    Parameters
    ----------
    f : file
        Open file-like object. 
    Returns
    -------
    systypes : SystypeData
        An object containing the data and information from the systype file.
    """
    try: 
        table = pd.read_csv(f, sep=r'\s+', skiprows=1, header=None)
    except Exception as e:
        raise BadSystypeError(f"Could not read file {f}. Please check
    there is a valid SYSTYPES file at this location.") from e
    # build header
    header = ["sys_index", "treatment", "description"]
    table.columns = header

    # Populate SystypeData object
    return SystypeData(
                    setname = setname,
                    nsys = len(table),
                    systypes = table)