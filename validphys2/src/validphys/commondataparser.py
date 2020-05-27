"""
This module implements parsers for commondata  and systype files into useful
datastructures, contained in the :py:mod:`validphys.coredata` module, which are
not backed by C++ managed memory, and so they can be easily pickled and
interfaces with common Python libraries.  The integration of these objects into
the codebase is currently work in progress, and at the moment this module
serves as a proof of concept.
"""
from collections import namedtuple

import numpy as np
import pandas as pd

from validphys.coredata import CommonData, SystypeData

CommondataInfo = namedtuple(
    "CommondataInfo", ("commondata", "systypes")
)


def load_commondata(spec):
    """
    Load the data corresponding to a CommonDataSpec object.
    Returns an instance of the namedtuple CommondataInfo,
    with:
    commondata_table being a CommonData instance with data arranged like
    entry   process kin1    kin2    kin3    data    stat    \
            sys.add.0   sys.mult.0 .... sys.add.N   sys.mult.N ;
    systype_table being a SystypeData instance with data arranged like
    sys_index   treatment   description.
    """
    commondatafile = spec.datafile

    # Getting set name from commondata file name
    setname = commondatafile.name[:-4] # removing the .dat suffix
    commondata = parse_commondata(commondatafile, setname)

    systypefile = spec.sysfile
    systypedata = parse_systype(systypefile, setname)

    return CommondataInfo(commondata=commondata, systypes=systypedata)


def parse_commondata(f, setname):
    """Parse a commondata file into a CommonData.

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
    # Remove NaNs
    # TODO: replace commondata files with bad formatting
    table.dropna(axis="columns", inplace=True)

    # Build header
    header = ['entry', 'process', 'kin1', 'kin2', 'kin3', 'data', 'stat']
    nsys  = (table.shape[1] - len(header)) // 2
    for i in range(nsys):
        header += [f"sys.add.{i+1}", f"sys.mult.{i+1}"]
    table.columns = header
    table.set_index("entry", inplace=True)

    # Populate CommonData object
    return CommonData(
        setname=setname,
        ndata=len(table),
        commondataproc=table["process"][0],
        nkin=3,
        nsys=nsys,
        commondata_table=table
    )


def parse_systype(f, setname):
    """Parse a systype file into a SystypeData.

    Parameters
    ----------
    f : file
        Open file-like object.

    Returns
    -------
    systypes : SystypeData
        An object containing the data and information from the systype file.
    """
    table = pd.read_csv(f, sep=r'\s+', skiprows=1, header=None)
    table.dropna(axis='columns', inplace=True)
    # Build header
    header = ["sys_index", "treatment", "description"]
    table.columns = header
    table.set_index("sys_index", inplace=True)

    # Populate SystypeData object
    return SystypeData(setname=setname, nsys=len(table), systype_table=table)
