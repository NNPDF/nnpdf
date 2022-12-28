"""
This module implements parsers for commondata  and systype files into useful
datastructures, contained in the :py:mod:`validphys.coredata` module, which are
not backed by C++ managed memory, and so they can be easily pickled and
interfaces with common Python libraries.  The integration of these objects into
the codebase is currently work in progress, and at the moment this module
serves as a proof of concept.
"""
from operator import attrgetter

import pandas as pd

from validphys.core import peek_commondata_metadata
from validphys.coredata import CommonData
import logging

log = logging.getLogger(__name__)

def load_commondata(spec):
    """
    Load the data corresponding to a CommonDataSpec object.
    Returns an instance of CommonData
    """
    commondatafile = spec.datafile
    setname = spec.name
    systypefile = spec.sysfile

    commondata = parse_commondata(commondatafile, systypefile, setname)

    return commondata


def parse_commondata(commondatafile, systypefile, setname):
    """Parse a commondata file  and a systype file into a CommonData.

    Parameters
    ----------
    commondatafile : file or path to file
    systypefile : file or path to file

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
    # Build header
    commondataheader = ['entry', 'process', 'kin1', 'kin2', 'kin3', 'data', 'stat']
    nsys  = (commondatatable.shape[1] - len(commondataheader)) // 2

    commondataheader += ["ADD", "MULT"] * nsys
    commondatatable.columns = commondataheader
    commondatatable.set_index("entry", inplace=True)
    ndata = len(commondatatable)
    commondataproc = commondatatable["process"][1]
    # Check for consistency with commondata metadata
    cdmetadata =  peek_commondata_metadata(commondatafile)
    if (setname, nsys, ndata) != attrgetter('name', 'nsys', 'ndata')(cdmetadata):
        raise ValueError("Commondata table information does not match metadata")

    # Now parse the systype file
    systypetable = parse_systypes(systypefile)

    # Populate CommonData object
    return CommonData(
        setname=setname,
        ndata=ndata,
        commondataproc=commondataproc,
        nkin=3,
        nsys=nsys,
        commondata_table=commondatatable,
        systype_table=systypetable
    )

def parse_systypes(systypefile):
    """Parses a systype file and returns a pandas dataframe.
    """
    systypeheader = ["sys_index", "type", "name"]
    try:
        systypetable = pd.read_csv(
            systypefile, sep=r"\s+", names=systypeheader, skiprows=1, header=None
        )
        systypetable.dropna(axis='columns', inplace=True)
    # Some datasets e.g. CMSWCHARMRAT have no systematics
    except pd.errors.EmptyDataError:
        systypetable = pd.DataFrame(columns=systypeheader)

    systypetable.set_index("sys_index", inplace=True)

    return systypetable



def write_commondata_data(commondata, buffer):
    """
    write commondata table to buffer, this can be a memory map, 
    compressed archive or strings (using for instance StringIO)
    
    
    Parameters
    ----------
    
    commondata : validphys.coredata.CommonData
    
    buffer : memory map, compressed archive or strings
            example: StringIO object
            
            
    Example
    -------
    >>> from validphys.loader import Loader
    >>> from io import StringIO
    
    >>> l = Loader()
    >>> cd = l.check_commondata("NMC").load_commondata_instance()
    >>> sio = StringIO()
    >>> write_commondata_data(cd,sio)
    >>> print(sio.getvalue())
    
    """
    header = f"{commondata.setname} {commondata.nsys} {commondata.ndata}\n"
    buffer.write(header)
    commondata.commondata_table.to_csv(buffer, sep="\t", header=None)

def write_commondata_to_file(commondata,path):
    """
    write commondata table to file
    """
    with open(path,"w") as file:
        write_commondata_data(commondata,file)

def write_systype_data(commondata, buffer):
    """
    write systype table to buffer, this can be a memory map, 
    compressed archive or strings (using for instance StringIO)
    
    
    Parameters
    ----------
    
    commondata : validphys.coredata.CommonData
    
    buffer : memory map, compressed archive or strings
            example: StringIO object
            
            
    Example
    -------
    >>> from validphys.loader import Loader
    >>> from io import StringIO
    
    >>> l = Loader()
    >>> cd = l.check_commondata("NMC").load_commondata_instance()
    >>> sio = StringIO()
    >>> write_systype_data(cd,sio)
    >>> print(sio.getvalue())
    
    """
    header = f"{commondata.nsys}\n"
    buffer.write(header)
    commondata.systype_table.to_csv(buffer, sep="\t", header=None)

def write_systype_to_file(commondata,path):
    """
    write systype table to file
    """
    with open(path,"w") as file:
        write_systype_data(commondata,file)
