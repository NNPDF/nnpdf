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
from pathlib import Path
from typing import Optional, Dict
import dataclasses

import pandas as pd
from reportengine.compat import yaml
from validobj.custom import Parser
from validobj import ValidationError, parse_input

from validphys.utils import parse_yaml_inp
from validphys.coredata import CommonData

EXT = "pineappl.lz4"


# Auxiliary parser for common types or sanity checks
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


# Auxiliary objects


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
    def parser(cls, meta: dict):
        """The yaml databases in the server use "operands" instead of "FK_tables" """
        if "operands" in meta:
            meta["FK_tables"] = meta.pop("operands")
        return parse_input(meta, cls)


ValidTheory = Parser(TheoryMeta.parser)


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
    theory: ValidTheory
    nnpdf_metadata: dict
    version: int
    version_comment: str = ""
    arXiv: Optional[ValidReference] = None
    iNSPIRE: Optional[ValidReference] = None
    hepdata: Optional[ValidReference] = None
    variants: Optional[ValidVariants] = dataclasses.field(default_factory=dict)

    def apply_variant(self, variant_name):
        """Return a new instance of this class with the variant applied

        This class also defines how the variant is applied to the commondata
        """
        try:
            variant = self.variants[variant_name]
        except KeyError as e:
            raise ValueError(f"The requested variant does not exist in {ret.name}") from e

        return dataclasses.replace(self, data_uncertainties=variant.data_uncertainties)


# TODO: will be moved to coredata.py and will substitute CommonData.py
# and obviously the folder here is only for development purposes
_folder_data = Path(__file__).parent / "../../../buildmaster"


@dataclasses.dataclass
class _CommonData:
    """
    Data, kinematics and uncertainties contained in Commondata files.
    A CommonData is only defined by its name and the enabled variant(s).

    The information from CommonData is provided by the following properties
        - metadata: all metadata information for the dataset
        - data: central data for the dataset
        - uncertainties: uncertainties of the dataset
        - kinematics: kinematic information
    """

    name: str
    variants: list[str] = dataclasses.field(default_factory=list)

    _metadata: Optional[CommonMetaData] = None

    def enable_variant(self, new_variant):
        """Function to enable a new variant for the CommonData
        This function checks whether the variant is accepted and, if it is,
        it gets added to the ``variants`` property of this class so that
        when the metadata is accessed it contains the variant information.

        The attribute ``metadata`` gets then replace with a new instance of ``CommonMetaData``
        """
        if new_variant is None:
            # Reset the metadata to the unvaried state
            self.variants = []
            self._metadata = None
            return

        if new_variant not in self.metadata.variants:
            raise ValueError(f"The requested variant does not exist in {self.name}")

        self.variants.append(new_variant)
        self._metadata = None

    # The files where the CommonData information is retreived from as defined by the metadata
    @property
    def metadata_file(self):
        return _folder_data / self.name / "metadata.yaml"

    @property
    def data_file(self):
        return _folder_data / self.name / self.metadata.data_central

    @property
    def uncertainity_files(self):
        return [_folder_data / self.name / i for i in self.metadata.data_uncertainties]

    @property
    def kinematics_file(self):
        return _folder_data / self.name / self.metadata.kinematics.file

    @property
    def ndata(self):
        return self.metadata.ndata

    @property
    def metadata(self):
        if self._metadata is None:
            # Read the metadata for the first time
            ret = parse_yaml_inp(self.metadata_file, CommonMetaData)

            # Iterate over variants to update the information in the original metadata
            for variant_name in self.variants:
                ret = ret.apply_variant(variant_name)

            self._metadata = ret
        return self._metadata

    @property
    def data(self):
        """Pandas DataFrame containing the central data"""
        datayaml = yaml.safe_load(self.data_file.read_text(encoding="utf-8"))
        data_df = pd.DataFrame(
            datayaml["data_central"], index=range(1, self.ndata + 1), columns=["data"]
        )
        data_df.index.name = "index"
        return data_df

    @property
    def uncertainties(self):
        """Pandas DataFrame containing all uncertainties"""
        # TODO: uncertainties are complicated enough that they _might_ need their own class
        all_df = []
        for ufile in self.uncertainity_files:
            uncyaml = yaml.safe_load(ufile.read_text())

            mindex = pd.MultiIndex.from_tuples(
                [(k, v["treatment"], v["type"]) for k, v in uncyaml["definition"].items()],
                names=["name", "treatment", "type"],
            )
            # I'm guessing there will be a better way of doing this than calling  dataframe twice for the same thing?
            final_df = pd.DataFrame(
                pd.DataFrame(uncyaml["bins"]).values,
                columns=mindex,
                index=range(1, self.ndata + 1),
            )
            final_df.index.name = "index"
            all_df.append(final_df)
        return pd.concat(all_df, axis=1)

    @property
    def kinematics(self):
        """Pandas DataFrame containing kinematic information"""
        kinyaml = yaml.safe_load(self.kinematics_file.read_text())

        kin_dict = {i + 1: pd.DataFrame(d).stack() for i, d in enumerate(kinyaml["bins"])}
        kin_df = pd.concat(kin_dict, axis=1, names=["index"]).swaplevel(0, 1).T
        return kin_df


def parse_commondata_folder(commondata, variants=[]):
    """Given a commondata folder, parse the entire content into the appropiate objects
    A list of variants can be also given
    """
    cmdata = _CommonData(commondata, variants=variants)
    return cmdata


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
    # TODO: the check for consistency would introduce a circular import
    # this entire function is to be removed _anyway_
    from validphys.core import peek_commondata_metadata
    from operator import attrgetter

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
