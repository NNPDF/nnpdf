"""
This module implements parsers for commondata and its associated metadata and uncertainties files
into useful structures that can be fed to the main :py:class:`validphys.coredata.CommonData` class.

A CommonData file is completely defined by a dataset name
(which defines the folder in which the information is)
and observable name (which defines the specific data, fktables and plotting settings to read).

<experiment>_<process>_<energy>{_<extras>}_<observable>

Where the folder name is ``<experiment>_<process>_<energy>{_<extras>}``

The definition of all information for a given dataset (and all its observable) is in the
``metadata.yaml`` file and its ``implemented_observables``.


This module defines a number of parsers using the ``validobj`` library.

The full ``metadata.yaml`` is read as a ``SetMetaData`` object
which contains a list of ``ObservableMetaData``.
These ``ObservableMetaData`` are the "datasets" of NNPDF for all intents and purposes.
The parent ``SetMetaData`` collects some shared variables such as the version of the dataset,
arxiv, inspire or hepdata ids, the folder in which the data is, etc.

The main class in this module is thus ``ObservableMetaData`` which holds _all_ information
about the particular dataset-observable that we are interested in (and a reference to its parent).

Inside the ``ObservableMetaData`` we can find:
    - ``TheoryMeta``: contains the necessary information to read the (new style) fktables
    - ``KinematicsMeta``: containins metadata about the kinematics
    - ``PlottingOptions``: plotting style and information for validphys
    - ``Variant``: variant to be used 

The CommonMetaData defines how the CommonData file is to be loaded,
by modifying the CommonMetaData using one of the loaded Variants one can change the resulting
:py:class:`validphys.coredata.CommonData` object.
"""
import dataclasses
from functools import cached_property
import logging
from operator import attrgetter
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
from validobj import ValidationError, parse_input
from validobj.custom import Parser

from reportengine.compat import yaml
from validphys.coredata import KIN_NAMES, CommonData
from validphys.plotoptions.plottingoptions import PlottingOptions
from validphys.utils import parse_yaml_inp

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


@dataclasses.dataclass
class Variant:
    """The new commondata format allow the usage of variants
    A variant can overwrite a number of keys, as defined by this dataclass
    """

    data_uncertainties: list[ValidPath]


ValidVariants = Dict[str, Variant]


### Theory metadata
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
            this flag is left here for compatibility purposes but has been moved to TheoryMeta
    """

    repetition_flag: Optional[list[str]] = None
    normalization: Optional[dict] = None
    shifts: Optional[dict] = None


@dataclasses.dataclass
class TheoryMeta:
    """Contains the necessary information to load the associated fktables

    The theory metadata must always contain a key ``FK_tables`` which defines
    the fktables to be loaded.
    The ``FK_tables`` is organized as a double list such that:

    The inner list is concatenated
    In practice these are different fktables that might refer to the same observable but
    that are divided in subgrids for practical reasons.
    The outer list instead are the operands for whatever operation needs to be computed
    in order to match the experimental data.

    In addition there are other flags that can affect how the fktables are read or used:
    - operation: defines the operation to apply to the outer list
    - shifts: mapping with the single fktables and their respective shifts
              useful to create "gaps" so that the fktables and the respective experimental data
              are ordered in the same way (for instance, when some points are missing from a grid)

    Example
    -------
    >>> from validphys.commondataparser import TheoryMeta
    ... from validobj import parse_input
    ... from reportengine.compat import yaml
    ... theory_raw = '''
    ... FK_tables:
    ...   - - fk1
    ...   - - fk2
    ...     - fk3
    ... operation: ratio
    ... apfelcomb:
    ...   repetition_flag:
    ...     - fk3
    ... '''
    ... theory = yaml.safe_load(theory_raw)
    ... parse_input(theory, TheoryMeta)
    TheoryMeta(FK_tables=[['fk1'], ['fk2', 'fk3']], operation='RATIO', shifts = None, conversion_factor=1.0, comment=None, apfelcomb=ValidApfelComb(repetition_flag=['fk3'], normalization=None))

    """

    FK_tables: list[list]
    operation: ValidOperation = "NULL"
    conversion_factor: float = 1.0
    comment: Optional[str] = None
    shifts: Optional[dict] = None
    apfelcomb: Optional[ValidApfelComb] = None
    # The following options are transitional so that the yamldb can be used from the theory
    appl: Optional[bool] = False
    target_dataset: Optional[str] = None

    def __post_init__(self):
        """If a ``shifts`` flag is found in the apfelcomb object, move it outside"""
        if self.apfelcomb is not None:
            if self.apfelcomb.shifts is not None and self.shifts is None:
                self.shifts = self.apfelcomb.shifts
                self.apfelcomb.shifts = None

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


###


### Kinematic data
@dataclasses.dataclass
class ValidVariable:
    """Defines the variables"""

    label: str
    description: str = ""
    units: str = ""

    def full_label(self):
        if self.units:
            return f"{self.label} ({self.units})"
        return self.label

    def apply_label(self, value):
        """Return a string formatted as label = value (units)"""
        tmp = f"{self.label} = {value}"
        if self.units:
            tmp += f" ({self.units})"
        return tmp


@dataclasses.dataclass
class ValidKinematics:
    """Contains the metadata necessary to load the kinematics of the dataset.
    The variables should be a dictionary with the key naming the variable
    and the content complying with the ``ValidVariable`` spec.

    Only the kinematics defined by the key ``kinematic_coverage`` will be loaded,
    which must be three.

    Three shall be the number of the counting and the number of the counting shall be three.
    Four shalt thou not count, neither shalt thou count two,
    excepting that thou then proceedeth to three.
    Once the number three, being the number of the counting, be reached,
    then the kinematics be loaded in the direction of thine validobject.
    """

    file: ValidPath
    variables: Dict[str, ValidVariable]

    def get_label(self, var):
        """For the given variable, return the label as
        label (unit)
        """
        return self.variables[var].full_label()

    def apply_label(self, var, value):
        """For a given value for a given variable, return the labels
        as label = value (unit)
        If the variable is not include in the list of variables, returns None
        as the variable could've been transformed by a kinematic transformation
        """
        if var not in self.variables:
            return None
        return self.variables[var].apply_label(value)


###


### Observable and dataset definitions
@dataclasses.dataclass
class ObservableMetaData:
    observable_name: str
    observable: dict
    ndata: int
    # Data itself
    kinematics: ValidKinematics
    data_central: ValidPath
    data_uncertainties: list[ValidPath]
    # Plotting options
    plotting: PlottingOptions
    process_type: str
    kinematic_coverage: list[str]
    # Optional data
    theory: Optional[TheoryMeta] = None
    tables: Optional[list] = dataclasses.field(default_factory=list)
    npoints: Optional[list] = dataclasses.field(default_factory=list)
    variants: Optional[ValidVariants] = dataclasses.field(default_factory=dict)
    _parent: Optional[
        Any
    ] = None  # Note that an observable without a parent will fail in many different ways

    def __post_init__(self):
        """Checks to be run after reading the metadata file"""
        # Check that plotting.plot_x is being filled
        if self.plotting.plot_x is None:
            ermsg = "No variable selected as x-axis in the plot for {self.name}. Please add plotting::plot_x."
            if self.plotting.x is not None:
                ermsg += "If you are using `plotting:x` please change it to `plotting::plot_x`"
            raise ValidationError(ermsg)

        # Ensure that all variables in the kinematic coverage exist
        for var in self.kinematic_coverage:
            if var not in self.kinematics.variables:
                raise ValidationError(
                    f"Variable {var} is in `kinematic_coverage` but not included in `kinematics` for {self.name}"
                )

        if len(self.kinematic_coverage) > 3:
            raise ValidationError(
                "Only a maximum of 3 variables can be used for `kinematic_coverage`"
            )

        # Since vp will rely on the kinematics being 3 variables,
        # fill the extra with whatever can be found in the kinematics dictionary
        if len(self.kinematic_coverage) < 3:
            unused = list(set(self.kinematics.variables) - set(self.kinematic_coverage))
            self.kinematic_coverage += unused[3 - len(self.kinematic_coverage) :]

        self.process_type = self.process_type.upper()

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

    @property
    def setname(self):
        return self._parent.setname

    @property
    def experiment(self):
        return self.setname.split("_")[0]

    @property
    def process(self):
        return self.setname.split("_")[1]

    @property
    def cm_energy(self):
        return self.setname.split("_")[2]

    @property
    def name(self):
        return f"{self.setname}_{self.observable_name}"

    @property
    def kinlabels(self):
        """Return the kinematic labels in the same order as they are set
        in ``kinematic_coverage`` (which in turns follow the key kinematic_coverage
        """
        return [self.kinematics.get_label(i) for i in self.kinematic_coverage]

    @cached_property
    def plotting_options(self):
        """Return the PlottingOptions metadata

        Fill in missing information that can be learnt from the other variables (xlabel/ylabel)
        or that is shared by the whole dataset.
        """
        if self.plotting.nnpdf31_process is None:
            self.plotting.nnpdf31_process = self.nnpdf_metadata["nnpdf31_process"]

        if self.plotting.experiment is None:
            self.plotting.experiment = self.nnpdf_metadata["experiment"]

        ## Swap variables by the k_idx
        # Internally validphys takes the x/y to be "k1" "k2" or "k3"
        # Therefore, for the time being, swap the actual keys by k1/k2/k3
        used_idx = []
        x_idx = self.kinematic_coverage.index(self.plotting.plot_x)
        used_idx.append(x_idx)
        self.plotting.x = f"k{x_idx + 1}"
        if self.plotting.x_label is None:
            self.plotting.x_label = self.kinematics.get_label(self.plotting.plot_x)

        # Swap the `figure_by` and `line_by` variables by k1/k2/k3
        # unless this is something coming from the "extra labels"
        if self.plotting.figure_by is not None:
            new_fig_by = []
            for var in self.plotting.figure_by:
                if var in self.kinematic_coverage:
                    fig_idx = self.kinematic_coverage.index(var)
                    used_idx.append(fig_idx)
                    new_fig_by.append(f"k{fig_idx + 1}")
                elif self.plotting.extra_labels is not None and var in self.plotting.extra_labels:
                    new_fig_by.append(var)
                else:
                    raise ValueError(f"Cannot find {var} in the kinematic coverage or extra labels")

            self.plotting.figure_by = new_fig_by

        if self.plotting.line_by is not None:
            new_line_by = []
            for var in self.plotting.figure_by:
                line_idx = self.kinematic_coverage.index(var)
                used_idx.append(line_idx)
                new_line_by.append(f"k{line_idx + 1}")
            self.plotting.line_by = new_line_by

        return self.plotting


@dataclasses.dataclass
class ValidReference:
    """Holds literature information for the dataset"""

    url: str
    version: Optional[int] = None
    tables: list[int] = dataclasses.field(default_factory=list)


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
    _folder: Optional[Path] = None

    @property
    def folder(self):
        # TODO: at the moment the folder is set manually by the parser of the metadata
        # since the new commondata is still not installed (or declared in the profile)
        return self._folder
        # return _folder_data / self.setname

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


###


### Parsers
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
            pd.DataFrame(uncyaml["bins"]).values, columns=mindex, index=range(1, metadata.ndata + 1)
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


def parse_new_metadata(metadata_file, observable_name, variants=[]):
    """Given a metadata file in the new format and the specific observable to be read
    load and parse the metadata and select the observable.
    If any variants are selected, apply them.
    """
    # Note: we are re-loading many times the same yaml file, possibly a good target for lru_cache
    set_metadata = parse_yaml_inp(metadata_file, SetMetaData)
    set_metadata._folder = metadata_file.parent

    # Select one observable from the entire metadata
    metadata = set_metadata.select_observable(observable_name)

    # And apply variants
    for variant in variants:
        metadata = metadata.apply_variant(variant)

    return metadata


def parse_commondata_new(metadata):
    """

    TODO: update this docstring since now the parse_commondata_new takes the information from
    the metadata, and the name -> split is done outside

    In the current iteration of the commondata, each of the commondata
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

    procname = metadata.process_type  # nnpdf_metadata["nnpdf31_process"]
    kin_df = kin_df[metadata.kinematic_coverage]
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

    # TODO: this will be removed because the old ones will be loaded with the new names
    # but during the implementation this is useful for the cuts (filters.py, __call__)
    names_file = metadata.path_kinematics.parent.parent / "dataset_names.yml"
    names_dict = yaml.YAML().load(names_file)
    if names_dict is not None:
        legacy_name = names_dict.get(metadata.name)
    else:
        legacy_name = metadata.name

    return CommonData(
        setname=metadata.name,
        ndata=metadata.ndata,
        commondataproc=procname,
        nkin=3,
        nsys=nsys,
        commondata_table=commondata_table,
        systype_table=systype_table,
        legacy=False,
        legacy_name=legacy_name,
        kin_variables=metadata.kinematic_coverage,
    )


###


def load_commondata(spec):
    """
    Load the data corresponding to a CommonDataSpec object.
    Returns an instance of CommonData
    """
    if spec.legacy:
        commondatafile = spec.datafile
        setname = spec.name
        systypefile = spec.sysfile

        commondata = parse_commondata(commondatafile, systypefile, setname)
    else:
        commondata = parse_commondata_new(spec.metadata)

    return commondata


### Old commondata:
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

    # TODO: the keys in KINLABEL_LATEX need to be updated for the new commondata
    return KINLABEL_LATEX.get(key, key)


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
