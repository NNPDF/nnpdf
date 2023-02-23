"""
This module implements parsers for commondata and its associated metadata and uncertainties files
into useful structures that can be fed to the main :py:class:`validphys.coredata.CommonData` class.

A CommonData file is completely defined by a name (which defines the folder in which the information is)
and a variant (which defines which files inside of the folder will be read).


In this module a few auxiliary dataclasses that hold special information:
    - TheoryMeta: contains the necessary information to read the (new style) fktables
    - KinematicsMeta: containins metadata about the kinematics
    - ReferenceMeta: literature references for the dataset

The CommonMetaData defines how the CommonData file is to be loaded,
by modifying the CommonMetaData using one of the loaded Variants one can change the resulting
:py:class:`validphys.coredata.CommonData` object.
"""
import dataclasses
import logging
from operator import attrgetter
from pathlib import Path
from typing import Optional, Dict

import pandas as pd
from reportengine.compat import yaml
from validobj.custom import Parser
from validobj import ValidationError, parse_input

from validphys.utils import parse_yaml_inp
from validphys.coredata import CommonData

EXT = "pineappl.lz4"
_INDEX_NAME = "entry"

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


@Parser
def ValidPath(path_str: str) -> Path:
    """Parse strings into paths"""
    return Path(path_str)


@Parser
def ValidOperation(op_str: str) -> str:
    """Ensures that the operation defined in the commondata file is implemented in validphys"""
    ret = op_str.upper()
    # TODO: move accepted operations to this module so that the convolution receives an operation to apply
    # instead of an operation to understand
    from validphys.convolution import OP

    if ret not in OP:
        raise ValidationError(f"The operation '{op_str}' is not implemented in validphys")
    return ret


@dataclasses.dataclass
class ValidApfelComb:
    """Some of the grids might have been converted from apfelcomb and introduce hacks.
    These are the allowed hacks:
        - repetition_flag:
            list of fktables which might need to be repeated
            necessary to apply c-factors in compound observables
        - normalization:
            mapping with the single fktables which need to be normalized and the factor
            note that when they are global factors they are promoted to conversion_factor
        - shifts:
            mapping with the single fktables and their respective shifts
            necessary to create "gaps" that some cfactors or datasets might expect
    """

    repetition_flag: Optional[list[str]] = None
    normalization: Optional[dict] = None
    shifts: Optional[dict] = None


@dataclasses.dataclass
class TheoryMeta:
    """Contains the necessary information to load the associated fktables"""

    FK_tables: list[list]
    operation: ValidOperation
    conversion_factor: float = 1.0
    comment: Optional[str] = None
    apfelcomb: Optional[ValidApfelComb] = None
    # The following options are transitional so that the yamldb can be used from the theory
    appl: Optional[bool] = False
    target_dataset: Optional[str] = None

    def fktables_to_paths(self, grids_folder):
        """Given a source for pineappl grids, constructs the lists of fktables
        to be loaded"""
        ret = []
        for operand in self.FK_tables:
            ret.append([grids_folder / f"{m}.{EXT}" for m in operand])
        return ret

    @classmethod
    def parser(cls, yaml_file):
        """The yaml databases in the server use "operands" instead of "FK_tables" """
        if not yaml_file.exists():
            raise FileNotFoundError(yaml_file)
        meta = yaml.safe_load(yaml_file.read_text())
        # Make sure the operations are upper-cased for compound-compatibility
        meta["operation"] = "NULL" if meta["operation"] is None else meta["operation"].upper()
        if "operands" in meta:
            meta["FK_tables"] = meta.pop("operands")
        return parse_input(meta, cls)


@dataclasses.dataclass
class ValidKinematics:
    """Contains all metadata for the kinematics of the dataset"""

    file: ValidPath
    variables: dict


@dataclasses.dataclass
class ValidReference:
    """Holds literature information for the dataset"""

    url: str
    version: int = 0
    tables: list[int] = dataclasses.field(default_factory=list)


@dataclasses.dataclass
class Variant:
    """Defines the keys of the CommonMetaData that can be overwritten"""

    data_uncertainties: list[ValidPath]


@Parser
def ValidVariants(variant_dict: dict) -> Dict[str, Variant]:
    """Variants of a dataset are allowed to overwrite a subset of the keys of a dataset
    (those defined in the Variant dataclass).
    This wrapper class runs over the dictionary of variant and parses them into valid Variants
    """
    return {k: parse_input(v, Variant) for k, v in variant_dict.items()}


# TODO: obviously the folder here is only for development purposes once the
# whole thing is finished the data will be installed in the right path as given by the profile
_folder_data = Path(__file__).parent / "../../../buildmaster"


@dataclasses.dataclass
class CommonMetaData:
    setname: str
    ndata: int
    observable: dict
    kinematics: ValidKinematics
    kinematic_coverage: dict
    data_central: ValidPath
    data_uncertainties: list[ValidPath]
    dataset_label: str
    plot_x: str
    figure_by: list[str]
    theory: TheoryMeta
    nnpdf_metadata: dict
    version: int
    version_comment: str = ""
    arXiv: Optional[ValidReference] = None
    iNSPIRE: Optional[ValidReference] = None
    hepdata: Optional[ValidReference] = None
    variants: Optional[ValidVariants] = dataclasses.field(default_factory=dict)
    _folder: Optional[Path] = None

    def apply_variant(self, variant_name):
        """Return a new instance of this class with the variant applied

        This class also defines how the variant is applied to the commondata
        """
        try:
            variant = self.variants[variant_name]
        except KeyError as e:
            raise ValueError(f"The requested variant does not exist in {ret.name}") from e

        return dataclasses.replace(self, data_uncertainties=variant.data_uncertainties)

    @property
    def folder(self):
        if self._folder is None:
            self._folder = _folder_data / self.setname
        return self._folder


@dataclasses.dataclass
class ValidKinematics:
    """Contains all metadata for the kinematics of the dataset"""

    file: ValidPath
    variables: dict


@dataclasses.dataclass
class ValidReference:
    """Holds literature information for the dataset"""

    url: str
    version: int = 0
    tables: list[int] = dataclasses.field(default_factory=list)


@dataclasses.dataclass
class Variant:
    """Defines the keys of the CommonMetaData that can be overwritten"""

    data_uncertainties: list[ValidPath]


@Parser
def ValidVariants(variant_dict: dict) -> Dict[str, Variant]:
    """Variants of a dataset are allowed to overwrite a subset of the keys of a dataset
    (those defined in the Variant dataclass).
    This wrapper class runs over the dictionary of variant and parses them into valid Variants
    """
    return {k: parse_input(v, Variant) for k, v in variant_dict.items()}


# TODO:
# These three parsers could just as well be methods of the CommonMetaData class
def _parse_data(metadata):
    """Given the metadata defining the commondata,
    returns a dataframe with the right central data loaded

    Parameters
    ----------
    metadata: CommonMetaData
        instance of CommonMetaData defining exactly the commondata to be loaded

    Returns
    -------
    pd.DataFrame
        a dataframe containing the data
    """
    data_file = metadata.folder / metadata.data_central
    datayaml = yaml.safe_load(data_file.read_text(encoding="utf-8"))
    data_df = pd.DataFrame(
        datayaml["data_central"], index=range(1, metadata.ndata + 1), columns=["data"]
    )
    data_df.index.name = _INDEX_NAME
    return data_df


def _parse_uncertainties(metadata):
    """Given the metadata defining the commondata,
    returns a dataframe with all appropiate uncertainties

    Parameters
    ----------
    metadata: CommonMetaData
        instance of CommonMetaData defining exactly the commondata to be loaded

    Returns
    -------
    pd.DataFrame
        a dataframe containing the uncertainties
    """
    uncertainity_files = [metadata.folder / i for i in metadata.data_uncertainties]
    all_df = []
    for ufile in uncertainity_files:
        uncyaml = yaml.safe_load(ufile.read_text())

        mindex = pd.MultiIndex.from_tuples(
            [(k, v["treatment"], v["type"]) for k, v in uncyaml["definition"].items()],
            names=["name", "treatment", "type"],
        )
        # I'm guessing there will be a better way of doing this than calling  dataframe twice for the same thing?
        final_df = pd.DataFrame(
            pd.DataFrame(uncyaml["bins"]).values,
            columns=mindex,
            index=range(1, metadata.ndata + 1),
        )
        final_df.index.name = _INDEX_NAME
        all_df.append(final_df)
    return pd.concat(all_df, axis=1)


def _parse_kinematics(metadata):
    """Given the metadata defining the commondata,
    returns a dataframe with the kinematic information

    Parameters
    ----------
    metadata: CommonMetaData
        instance of CommonMetaData defining exactly the commondata to be loaded

    Returns
    -------
    pd.DataFrame
        a dataframe containing the kinematics
    """
    kinematics_file = metadata.folder / metadata.kinematics.file
    kinyaml = yaml.safe_load(kinematics_file.read_text())
    kin_dict = {i + 1: pd.DataFrame(d).stack() for i, d in enumerate(kinyaml["bins"])}
    return pd.concat(kin_dict, axis=1, names=[_INDEX_NAME]).swaplevel(0, 1).T


def parse_commondata_folder(commondata_folder, variants=[]):
    """In the current iteration of the commondata, each of the commondata
    (i.e., an observable from a data publication) correspond to one single folder.
    This function returns a commondata object constructed by parsing the metadata.

    A commondata object is defined by the folder from where it is
    being read and the variants to be enabled.
    Once a variant is selected, it can no longer be changed

    Note that this function reproduces `parse_commondata` below, which parses the
    _old_ file format
    """
    # Parse the metadata by iterating over variants
    commondata = commondata_folder.name
    metadata_file = commondata_folder / "metadata.yaml"
    metadata = parse_yaml_inp(metadata_file, CommonMetaData)
    # TODO for debugging purposes
    # keep a reference folder for people trying to test this out
    metadata._folder = commondata_folder
    for variant in variants:
        metadata = metadata.apply_variant(variant)

    # Now parse the data
    data_df = _parse_data(metadata)
    # the uncertainties
    uncertainties_df = _parse_uncertainties(metadata)
    # and the kinematics
    kin_df = _parse_kinematics(metadata)

    # Once we have loaded all uncertainty files, let's check how many sys we have
    nsys = len([i for i in uncertainties_df.columns.get_level_values(0) if "syst" in i])

    # Backwards-compatibility
    # Finally, create the commondata by merging the dataframes in the old commondata_table

    # For the kinematis, forget all the interesting information
    procname = metadata.nnpdf_metadata["nnpdf31_process"]
    kinames = ["kin1", "kin2", "kin3"]
    kin_df.columns = kinames
    kin_df["process"] = procname
    kin_df = kin_df[["process"] + kinames]

    # For the uncertainties, create a simplified version to concatenate
    # and save the systype information
    new_columns = []
    systypes = {"type": [], "name": []}
    for col in uncertainties_df.columns:
        if col[0].startswith("stat"):
            new_columns.append("stat")
        else:
            # if it is syst add the ADD/MULT information
            new_columns.append(col[1])
            systypes["type"].append(col[1])
            systypes["name"].append(col[2])

    uncertainties_df.columns = new_columns

    commondata_table = pd.concat([kin_df, data_df, uncertainties_df], axis=1)
    systype_table = pd.DataFrame(systypes, index=range(1, nsys + 1))
    systype_table.index.name = "sys_index"

    return CommonData(
        setname=commondata,
        ndata=metadata.ndata,
        commondataproc=procname,
        nkin=3,
        nsys=nsys,
        commondata_table=commondata_table,
        systype_table=systype_table,
    )


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

    return KINLABEL_LATEX[key]


def get_kinlabel_key(process_label):
    """
    Since there is no 1:1 correspondence between latex keys and the old libNNPDF names
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
