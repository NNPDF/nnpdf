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

def write_commondata(commondata_list, filter_path):
    """
    GENERAL DESCRIPTION:

    writes dataset data in filter folder using the commondata file format

    Parameters
    ----------

    commondata_list : list
                     commondata object (cuts should be already applied)

    path : str
         path to filter folder

    
    """

    for commondata_instance in commondata_list:
        # path
        path = str(filter_path) + f'/{commondata_instance.setname}'
        path_data = str(path) + f"/DATA_{commondata_instance.setname}.dat"
        commondata_tab = commondata_instance.commondata_table
        header = f"{commondata_instance.setname} {commondata_instance.nsys} {commondata_instance.ndata}\n"
        
        #==== write DATA =====#
        with open(path_data, "w+") as f:
            f.write(header)
            commondata_tab.to_csv(f, sep="\t", header=None)



def make_systype_dir(path):
    """
    GENERAL DESCRIPTION:

    creates directory named systypes 

    Parameters
    ----------

    path : str
         path to systypes filter/dataset_name/systypes folder
    

    """
    if path.exists():
        log.warning(f"systypes folder exists: {path} Overwriting contents")
    else:
        path.mkdir(exist_ok=True)


def write_systype(commondata_list, filter_path):
    """
    GENERAL DESCRIPTION:

    writes systype data in filter folder using the systype file format  
    
    Parameters
    ----------

    commondata_list : list
                    commondata object (cuts should be already applied)

    filter_path : str
                 path to filter folder

    
    """

    for commondata_instance in commondata_list: 
        path = filter_path / commondata_instance.setname / 'systypes'
        make_systype_dir(path)
        path_data = str(path) + f"/SYSTYPES_{commondata_instance.setname}_DEFAULT.dat"
        systype_tab = commondata_instance.systype_table
        header = f"{len(systype_tab.index)}\n"
        
        #==== write DATA =====#
        with open(path_data, "w+") as f:
            f.write(header)
            systype_tab.to_csv(f, sep="\t", header=None)