"""
Note: this module will be removed after the next tag, don't use anything from here
"""

import dataclasses
import logging
from operator import attrgetter

import pandas as pd

from nnpdf_data.coredata import CommonData

log = logging.getLogger(__name__)

log.warning(
    "You are loading deprecated functionality that use the old commondata parser. This is no longer supported and will be removed in the near future"
)


### Old commondata:
### All code below this line is deprecated and will be removed
def load_commondata_old(commondatafile, systypefile, setname):
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
    commondatatable = pd.read_csv(commondatafile, sep=r"\s+", skiprows=1, header=None)
    # Remove NaNs
    # TODO: replace commondata files with bad formatting
    # Build header
    commondataheader = ["entry", "process", "kin1", "kin2", "kin3", "data", "stat"]
    nsys = (commondatatable.shape[1] - len(commondataheader)) // 2

    commondataheader += ["ADD", "MULT"] * nsys
    commondatatable.columns = commondataheader
    commondatatable.set_index("entry", inplace=True)
    ndata = len(commondatatable)
    commondataproc = commondatatable["process"][1]
    # Check for consistency with commondata metadata
    cdmetadata = peek_commondata_metadata(commondatafile)
    if (nsys, ndata) != attrgetter("nsys", "ndata")(cdmetadata):
        raise ValueError(f"Commondata table information does not match metadata for {setname}")

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
        systype_table=systypetable,
        legacy=True,
    )


def parse_systypes(systypefile):
    """Parses a systype file and returns a pandas dataframe."""
    systypeheader = ["sys_index", "treatment", "name"]
    try:
        systypetable = pd.read_csv(
            systypefile, sep=r"\s+", names=systypeheader, skiprows=1, header=None
        )
        systypetable.dropna(axis="columns", inplace=True)
    # Some datasets e.g. CMSWCHARMRAT have no systematics
    except pd.errors.EmptyDataError:
        systypetable = pd.DataFrame(columns=systypeheader)

    systypetable.set_index("sys_index", inplace=True)

    return systypetable


@dataclasses.dataclass(frozen=True)
class CommonDataMetadata:
    """Contains metadata information about the data being read"""

    name: str
    nsys: int
    ndata: int
    process_type: str


def peek_commondata_metadata(commondatafilename):
    """Read some of the properties of the commondata object as a CommonData Metadata"""
    with open(commondatafilename) as f:
        try:
            l = f.readline()
            name, nsys_str, ndata_str = l.split()
            l = f.readline()
            process_type_str = l.split()[1]
        except Exception:
            log.error(f"Error processing {commondatafilename}")
            raise

    return CommonDataMetadata(
        name, int(nsys_str), int(ndata_str), get_kinlabel_key(process_type_str)
    )


def get_plot_kinlabels(commondata):
    """Return the LaTex kinematic labels for a given Commondata"""
    key = commondata.process_type

    # TODO: the keys in KINLABEL_LATEX need to be updated for the new commondata
    return KINLABEL_LATEX.get(key, key)


def get_kinlabel_key(process_label):
    """
    Since there is no 1:1 correspondence between latex keys and the old libNNPDF names
    we match the longest key such that the proc label starts with it.
    """
    l = process_label
    try:
        if process_label == "EWK_RAP_ASY":
            # TODO this function is disappearing in this PR
            l = "EWK_RAP"
        return next(k for k in sorted(KINLABEL_LATEX, key=len, reverse=True) if l.startswith(k))
    except StopIteration as e:
        raise ValueError(
            "Could not find a set of kinematic "
            "variables matching  the process %s Check the "
            "labels defined in commondata.cc. " % (l)
        ) from e
