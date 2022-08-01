"""
This module implements parsers for commondata and its associated metadata and uncertainties files
into useful structures that can be fed to the main :py:class:`validphys.coredata.CommonData` class.

In this module a few auxiliary dataclasses that hold special information:
    - TheoryMeta: contains the necessary information to read the (new style) fktables
    - ReferenceMeta: literature references for the dataset

The CommonMetaData defines how the CommonData file is to be loaded,
by modifying the CommonMetaData using one of the loaded Variants one can change the resulting
:py:class:`validphys.coredata.CommonData` object.
"""
from copy import copy
from operator import attrgetter
from dataclasses import dataclass, field
from pathlib import Path
import typing

import pandas as pd
from validobj.custom import Parser
from validobj import ValidationError, parse_input

from validphys import convolution
from validphys.utils import parse_yaml_inp
from validphys.core import peek_commondata_metadata
from validphys.coredata import CommonData

# Auxiliary parser for common types or sanity checks
@Parser
def ValidPath(path_str: str) -> Path:
    """Parse strings into paths"""
    try:
        return Path(path_str)
    except Exception as e:
        raise ValidationError(f"{path_str} is not a valid path") from e


@Parser
def ValidOperation(op_str: str) -> str:
    """Ensures that the operation defined in the commondata file is implemented in validphys"""
    ret = op_str.upper()
    if ret not in convolution.OP:
        raise ValidationError(f"The operation '{op_str}' is not implemented in validphys")
    return ret


# Auxiliary objects
@dataclass
class TheoryMeta:
    """Contains the necessary information to load the associated fktables"""

    FK_tables: list
    operation: ValidOperation
    conversion_factor: float = 1.0
    apfelcomb: dict = None

    @classmethod
    def parser(cls, meta: dict):
        return parse_input(meta, cls)


@dataclass
class ReferenceMeta:
    """Holds literature information for the dataset"""

    url: str
    version: int = None
    tables: typing.List[int] = None

    @classmethod
    def parser(cls, meta: dict):
        return parse_input(meta, cls)


@dataclass
class Variant:
    """Defines the keys of the CommonMetaData that can be overwritten"""

    data_uncertainties: typing.List[ValidPath]


# Define parsers for the more complicated structures
ValidTheory = Parser(TheoryMeta.parser)
ValidReference = Parser(ReferenceMeta.parser)


@Parser
def ValidVariants(variant_dict: dict) -> dict:
    """Variants of a dataset are allowed to overwrite a subset of the keys of a dataset
    (those defined in the Variant dataclass).
    This wrapper class runs over the dictionary of variant and parses them into valid Variants
    """
    return {k: parse_input(v, Variant) for k, v in variant_dict.items()}


@dataclass
class CommonMetaData:
    setname: str
    ndata: int
    observable: dict
    kinematics: dict
    kinematic_coverage: dict
    data_central: ValidPath
    data_uncertainties: typing.List[ValidPath]
    dataset_label: str
    plot_x: str
    figure_by: typing.List[str]
    theory: ValidTheory
    nnpdf_metadata: dict
    version: int
    version_comment: str = ""
    arXiv: ValidReference = None
    iNSPIRE: ValidReference = None
    hepdata: ValidReference = None
    variants: ValidVariants = None
    _enabled_variant: typing.List[str] = field(default=None, repr=False)
    _default: dict = field(default=None, repr=False)

    def enable_variant(self, variant):
        """Enable a variant for this class by giving its name.
        Note that more than one variant can be enabled at once, but the last one will take priority
        """
        if self.variants is None:
            raise ValueError(f"There are no variants defined for {self.setname}")
        if variant not in self.variants:
            raise ValueError(f"The variant {variant} is not defined for {self.setname}")
        # If there were not enabled variants, save the current state of the class
        if self._enabled_variant is None:
            # A shallow copy is enough because the variants update attributes of this class
            self._default = copy(self.__dict__)
            self._enabled_variant = []

        self.__dict__.update(self.variants[variant].__dict__.items())
        self._enabled_variant.append(variant)

    def disable_variants(self):
        """Get the CommonMetaData back to its original state"""
        if self._default is not None:
            self.__dict__.update(self._default)


def parse_commondata_metadata(metadata_file):
    """Transitional function, it will be part of the load_commondata function below"""
    return parse_yaml_inp(metadata_file, CommonMetaData)


##### TODO: the functions below have not been touched yet so they are targetting old commondata
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
