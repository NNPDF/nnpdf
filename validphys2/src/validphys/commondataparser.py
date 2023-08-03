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
from typing import Optional, Dict, Any

import pandas as pd
from reportengine.compat import yaml
from validobj.custom import Parser
from validobj import ValidationError, parse_input

from validphys.utils import parse_yaml_inp
from validphys.coredata import CommonData, KIN_NAMES

EXT = "pineappl.lz4"
_INDEX_NAME = "entry"

log = logging.getLogger(__name__)

# TODO: obviously the folder here is only for development purposes once the
# whole thing is finished the data will be installed in the right path as given by the profile
_folder_data = "/Users/markcostantini/codes/nnpdfgit/nnpdf/buildmaster"


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
def ValidOperation(op_str: Optional[str]) -> str:
    """Ensures that the operation defined in the commondata file is implemented in validphys"""
    if op_str is None:
        op_str = "NONE"
    ret = op_str.upper()
    # TODO: move accepted operations to this module so that the convolution receives an operation to apply
    # instead of an operation to understand
    from validphys.convolution import OP

    if ret not in OP:
        raise ValidationError(f"The operation '{op_str}' is not implemented in validphys")
    return str(ret)


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
    operation: ValidOperation = "NULL"
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
    """Contains all metadata for the kinematics of the dataset
    The minimum number of variables is 3, any variable beyond 3 is ignored.
    """

    file: ValidPath
    variables: dict

    @property
    def keys(self):
        return [*self.variables][:3]


@dataclasses.dataclass
class ValidReference:
    """Holds literature information for the dataset"""

    url: str
    version: Optional[int] = None
    tables: list[int] = dataclasses.field(default_factory=list)


@dataclasses.dataclass
class Variant:
    """Defines the keys of the CommonMetaData that can be overwritten"""

    data_uncertainties: list[ValidPath]


ValidVariants = Dict[str, Variant]


@dataclasses.dataclass
class ObservableMetaData:
    observable_name: str
    observable: dict
    ndata: int
    # Plotting settings
    dataset_label: str
    plot_x: str
    figure_by: list[str]
    kinematic_coverage: dict
    # Data itself
    kinematics: ValidKinematics
    data_central: ValidPath
    data_uncertainties: list[ValidPath]
    # Optional data
    theory: Optional[TheoryMeta] = None
    tables: Optional[list] = dataclasses.field(default_factory=list)
    npoints: Optional[list] = dataclasses.field(default_factory=list)
    variants: Optional[ValidVariants] = dataclasses.field(default_factory=dict)
    _parent: Optional[
        Any
    ] = None  # Note that an observable without a parent will fail in many different ways

    def apply_variant(self, variant_name):
        """Return a new instance of this class with the variant applied

        This class also defines how the variant is applied to the commondata
        """
        try:
            variant = self.variants[variant_name]
        except KeyError as e:
            raise ValueError(f"The requested variant does not exist {self.observable_name}") from e

        return dataclasses.replace(self, data_uncertainties=variant.data_uncertainties)

    @property
    def path_data_central(self):
        return self._parent.folder / self.data_central

    @property
    def paths_uncertainties(self):
        return [self._parent.folder / i for i in self.data_uncertainties]

    @property
    def path_kinematics(self):
        return self._parent.folder / self.kinematics.file

    # Properties inherited from parent
    @property
    def nnpdf_metadata(self):
        return self._parent.nnpdf_metadata


@dataclasses.dataclass
class SetMetaData:
    """Metadata of the whole set"""

    setname: str
    version: int
    version_comment: str
    nnpdf_metadata: dict
    implemented_observables: list[ObservableMetaData]
    arXiv: Optional[ValidReference] = None
    iNSPIRE: Optional[ValidReference] = None
    hepdata: Optional[ValidReference] = None

    @property
    def folder(self):
        # TODO currently _folder_data is the buildmaster folder as set above
        # but it should come from share/NNPDF/whatever
        return _folder_data + "/" + self.setname

    def select_observable(self, obs_name_raw):
        """Check whether the observable is implemented and return said observable"""
        # TODO: should we check that we don't have two observables with the same name?
        obs_name = obs_name_raw.lower().strip()
        for observable in self.implemented_observables:
            if observable.observable_name.lower().strip() == obs_name:
                observable._parent = (
                    self  # Not very happy with this but not sure how to do in a better way?
                )
                return observable
        return ValueError(f"The selected observable {obs_name} does not exist in {self.setname}")


def _parse_data(metadata):
    """Given the metadata defining the commondata,
    returns a dataframe with the right central data loaded

    Parameters
    ----------
    metadata: ObservableMetaData
        instance of ObservableMetaData defining the data to be loaded

    Returns
    -------
    pd.DataFrame
        a dataframe containing the data
    """
    data_file = metadata.path_data_central
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
    metadata: ObservableMetaData
        instance of ObservableMetaData defining the uncertainties to be loaded

    Returns
    -------
    pd.DataFrame
        a dataframe containing the uncertainties
    """
    all_df = []
    for ufile in metadata.paths_uncertainties:
        uncyaml = yaml.safe_load(ufile.read_text())

        mindex = pd.MultiIndex.from_tuples(
            [(k, v["treatment"], v["type"]) for k, v in uncyaml["definitions"].items()],
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
    metadata: ObservableMetaData
        instance of ObservableMetaData defining the kinematics to be loaded

    Returns
    -------
    pd.DataFrame
        a dataframe containing the kinematics
    """
    kinematics_file = metadata.path_kinematics
    kinyaml = yaml.safe_load(kinematics_file.read_text())

    kin_dict = {}
    for i, dbin in enumerate(kinyaml["bins"]):
        bin_index = i + 1
        # TODO: for now we are dropping min/max information since it didn't exist in the past
        # unless the point doesn't have a mid value, in that case we need to generate it!
        for d in dbin.values():
            if d["mid"] is None:
                d["mid"] = 0.5 * (d["max"] + d["min"])
            d["min"] = None
            d["max"] = None
        kin_dict[bin_index] = pd.DataFrame(dbin).stack()

    return pd.concat(kin_dict, axis=1, names=[_INDEX_NAME]).swaplevel(0, 1).T


def parse_commondata_new(dataset_fullname, variants=[]):
    """In the current iteration of the commondata, each of the commondata
    (i.e., an observable from a data publication) correspond to one single observable
    inside a folder which is named as "<experiment>_<process>_<energy>_<extra>"
    The observable is defined by a last suffix of the form "_<obs>" so that the full name
    of the dataset is always:

        "<experiment>_<process>_<energy>{_<extra>}_<obs>"

    where <extra> is optional.

    This function right now works under the assumotion that the folder/observable
    is separated in the last _ so that:
        folder_name = <experiment>_<process>_<energy>{_<extra>}
    but note that this convention is still not fully defined.

    This function returns a commondata object constructed by parsing the metadata.

    Once a variant is selected, it can no longer be changed

    Note that this function reproduces `parse_commondata` below, which parses the
    _old_ file format
    """
    # Look at the folder & observable
    setfolder, observable_name = dataset_fullname.rsplit("_", 1)
    commondata_folder = _folder_data + "/" + setfolder

    metadata_file = commondata_folder + "/" + "metadata.yaml"

    # Parse the metadata by iterating over variants
    commondata = setfolder# commondata_folder.name
    # Load the entire set
    set_metadata = parse_yaml_inp(metadata_file, SetMetaData)

    # Select one observable
    metadata = set_metadata.select_observable(observable_name)

    for variant in variants:
        metadata = metadata.apply_variant(variant)

    # Now parse the data
    data_df = _parse_data(metadata)
    # the uncertainties
    uncertainties_df = _parse_uncertainties(metadata)
    # and the kinematics
    kin_df = _parse_kinematics(metadata)

    # Once we have loaded all uncertainty files, let's check how many sys we have
    nsys = len(
        [i for i in uncertainties_df.columns.get_level_values(0) if not i.startswith("stat")]
    )

    # Backwards-compatibility
    # Finally, create the commondata by merging the dataframes in the old commondata_table

    # For the kinematis, forget all the interesting information
    procname = metadata.nnpdf_metadata["nnpdf31_process"]
    kin_df = kin_df[metadata.kinematics.keys]
    kin_df.columns = KIN_NAMES
    kin_df["process"] = procname

    kin_df = kin_df[["process"] + KIN_NAMES]

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

    # TODO: Legacy compatibility
    # 1. Add a stat column if it doesn't exist
    # 2. Transform multiplicatie uncertainties into % as it was done in the older version

    if "stat" not in commondata_table:
        commondata_table["stat"] = 0.0

    if "MULT" in commondata_table:
        commondata_table["MULT"] = commondata_table["MULT"].multiply(
            100 / commondata_table["data"], axis="index"
        )

    return CommonData(
        setname=commondata,
        ndata=metadata.ndata,
        commondataproc=procname,
        nkin=3,
        nsys=nsys,
        commondata_table=commondata_table,
        systype_table=systype_table,
        legacy=False,
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
        legacy=True,
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
