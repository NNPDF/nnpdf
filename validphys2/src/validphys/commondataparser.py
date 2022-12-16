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


def write_commondata_table_to_string(commondata,sio,table="commondata_table"):
    """
    GENERAL DESCRIPTION:
    
    Given a validphys.coredata.CommonData instance and a StringIO,
    update value of StringIO to contents of specified CommonData table.
    
    Parameters
    ----------

    commondata: validphys.coredata.CommonData
    
    sio: StringIO
    
    table: str, default 'commondata_table'

    Example
    -------
    
    >>> from validphys.loader import Loader
    >>> from io import StringIO

    >>> l = Loader()
    >>> obs = "NMC"
    >>> commondata = l.check_commondata(obs).load_commondata_instance()

    >>> sio = StringIO()
    >>> print(sio.getvalue())
    >>> write_commondata_table_to_string(commondata,sio,table="systype_table")
    >>> sio.getvalue()


    """
    
    data_frame = getattr(commondata,table)
    data_frame.to_csv(sio, sep = "\t", header = None)
    
    

def write_systypes_to_file(commondata,path):
    """
    GENERAL DESCRIPTION:
    
    write a systype file to a specified path 
    (e.g. path = path_to_fit_folder/filter/name_observable/systypes/SYSTYPE_name_observable_DEFAULT.dat)

    Parameters
    ----------

    commondata: validphys.coredata.CommonData
    
    path: str
         
    
    """
    from io import StringIO
    sio = StringIO()
    
    write_commondata_table_to_string(commondata,sio,table="systype_table")
    header = f"{commondata.nsys}\n"
    
    with open(path,"w") as file:
        file.write(header)
        file.write(sio.getvalue())

        
        
def write_commondata_to_file(commondata,path):
    """
    GENERAL DESCRIPTION:
    
    write a commondata table as file to a specified path
    
    Parameters
    ----------

    commondata: validphys.coredata.CommonData
    
    path: str
         
         
    """
    from io import StringIO
    
    sio = StringIO()
    write_commondata_table_to_string(commondata,sio)
    
    header = f"{commondata.setname} {commondata.nsys} {commondata.ndata}\n"
    
    with open(path,"w") as file:
        file.write(header)
        file.write(sio.getvalue())
