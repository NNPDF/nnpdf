"""
This module implements parsers for commondata  and systype files into useful
datastructures, contained in the :py:mod:`validphys.coredata` module, which are
not backed by C++ managed memory, and so they can be easily pickled and
interfaces with common Python libraries.  The integration of these objects into
the codebase is currently work in progress, and at the moment this module
serves as a proof of concept.
"""
import numpy as np
import pandas as pd

from validphys.coredata import CommonData

def load_commondata(spec):
    """
    Load the data corresponding to a CommonDataSpec object.
    Returns an instance of CommonData,
    with:
    commondata_table being a pandas dataframe with data arranged like
    entry   process kin1    kin2    kin3    data    stat    \
            sys.add.0   sys.mult.0 .... sys.add.N   sys.mult.N ;
    systype_table being a pandas dataframe with data arranged like
    sys_index   treatment   description.
    """
    commondatafile = spec.datafile
    # Getting set name from commondata file name
    setname = commondatafile.name[:-4] # removing the .dat suffix
    systypefile = spec.sysfile

    commondata = parse_commondata(commondatafile, systypefile, setname)

    return commondata


def parse_commondata(commondatafile, systypefile, setname):
    """Parse a commondata file  and a systype file into a CommonData.

    Parameters
    ----------
    commondatafile : file
        Open file-like object.
    systypefile : file
        Open file-like object.

    Returns
    -------
    commondata : CommonData
        An object containing the data and information from the commondata 
        and systype files.
    """
    # First parse commondata file
    commondatatable = pd.read_csv(commondatafile, sep=r'\s+', skiprows=1, header=None)
    # Remove NaNs
    # TODO: replace commondata files with bad formatting
    commondatatable.dropna(axis="columns", inplace=True)
    # Build header
    commondataheader = ['entry', 'process', 'kin1', 'kin2', 'kin3', 'data', 'stat']
    nsys  = (commondatatable.shape[1] - len(commondataheader)) // 2
    for i in range(nsys):
        commondataheader += [f"sys.add.{i+1}", f"sys.mult.{i+1}"]
    commondataable.columns = commondataheader
    commmondatatable.set_index("entry", inplace=True)

    # Now parse systyle file
    systypetable = pd.read_csv(systypefile, sep=r'\s+', skiprows=1, header=None)
    systypetable.dropna(axis='columns', inplace=True)
    # Build header
    systypeheader = ["sys_index", "type", "name"]
    systypetable.columns = systypeheader
    systypetable.set_index("sys_index", inplace=True)

    # Populate CommonData object
    return CommonData(
        setname=setname,
        ndata=len(table),
        commondataproc=table["process"][1],
        nkin=3,
        nsys=nsys,
        commondata_table=commondatatable
        systype_table=systypetable
    )

