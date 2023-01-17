"""
This module implements parsers for commondata  and systype files into useful
datastructures, contained in the :py:mod:`validphys.coredata` module, which are
not backed by C++ managed memory, and so they can be easily pickled and
interfaces with common Python libraries. 

The validphys commondata structure is an instance of :py:class:`validphys.coredata.CommonData`
"""
import dataclasses
from operator import attrgetter
import logging

import pandas as pd

from validphys.coredata import CommonData

log = logging.getLogger(__name__)

KINLABEL_LATEX = {
    "DIJET": ("\\eta", "$\\m_{1,2} (GeV)", "$\\sqrt{s} (GeV)"),
    "DIS": ("$x$", "$Q^2 (GeV^2)$", "$y$"),
    "DYP": ("$y$", "$M^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "EWJ_JPT": ("$p_T (GeV)$", "$M^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "EWJ_JRAP": ("$\\eta/y$", "$M^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "EWJ_MLL": ("$M_{ll} (GeV)$", "$M_{ll}^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "EWJ_PT": ("$p_T (GeV)$", "$M^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "EWJ_PTRAP": ("$\\eta/y$", "$p_T^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "EWJ_RAP": ("$\\eta/y$", "$M^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "EWK_MLL": ("$M_{ll} (GeV)$", "$M_{ll}^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "EWK_PT": ("$p_T$ (GeV)", "$M^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "EWK_PTRAP": ("$\\eta/y$", "$p_T^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "EWK_RAP": ("$\\eta/y$", "$M^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "HIG_RAP": ("$y$", "$M_H^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "HQP_MQQ": ("$M^{QQ} (GeV)$", "$\\mu^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "HQP_PTQ": ("$p_T^Q (GeV)$", "$\\mu^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "HQP_PTQQ": ("$p_T^{QQ} (GeV)$", "$\\mu^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "HQP_YQ": ("$y^Q$", "$\\mu^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "HQP_YQQ": ("$y^{QQ}$", "$\\mu^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "INC": ("$0$", "$\\mu^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "JET": ("$\\eta$", "$p_T^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "PHT": ("$\\eta_\\gamma$", "$E_{T,\\gamma}^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "SIA": ("$z$", "$Q^2 (GeV^2)$", "$y$"),
}


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
    if (setname, nsys, ndata) != attrgetter("name", "nsys", "ndata")(cdmetadata):
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
        systype_table=systypetable,
    )


def parse_systypes(systypefile):
    """Parses a systype file and returns a pandas dataframe."""
    systypeheader = ["sys_index", "type", "name"]
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
    """Read some of the properties of the commondata object as a CommonData Metadata
    """
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

    return KINLABEL_LATEX[key]


def get_kinlabel_key(process_label):
    """
    Since there is no 1:1 correspondence between latex keys and GetProc,
    we match the longest key such that the proc label starts with it.
    """
    l = process_label
    try:
        return next(k for k in sorted(KINLABEL_LATEX, key=len, reverse=True) if l.startswith(k))
    except StopIteration as e:
        raise ValueError(
            "Could not find a set of kinematic "
            "variables matching  the process %s Check the "
            "labels defined in commondata.cc. " % (l)
        ) from e
