"""
This module implements parsers for commondata and its associated metadata and uncertainties files
into useful structures that can be fed to the main :py:class:`nnpdf_data.coredata.CommonData` class.

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
    - ``Variant``: variant to be used
    - ``PlottingOptions``: plotting style and information for validphys, only utilized if
                           validphys is also installed.

The CommonMetaData defines how the CommonData file is to be loaded,
by modifying the CommonMetaData using one of the loaded Variants one can change the resulting
:py:class:`nnpdf_data.coredata.CommonData` object.
"""

import dataclasses
from functools import cache, cached_property
import logging
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
from validobj import ValidationError, parse_input
from validobj.custom import Parser

from .coredata import KIN_NAMES, CommonData
from .process_options import ValidProcess
from .utils import parse_yaml_inp, quick_yaml_load
from .validphys_compatibility import new_to_legacy_map, path_commondata

try:
    from validphys.plotoptions.plottingoptions import PlottingOptions, labeler_functions

    VP_AVAILABLE = True
except ModuleNotFoundError:
    # if validphys is not available, the __old__ plotting options from validphys
    # which we only still have because the world is a dark and horrible place
    # won't be loaded. Instead, the following file is loaded.
    from .validphys_compatibility import PlottingOptions, labeler_functions

    VP_AVAILABLE = False


# JCM:
# Some notes for developers
# The usage of `frozen` in the definitions of the dataclass is not strictly necessary
# however, changing the metadata can have side effects in many parts on validphys.
# By freezing the overall class (and leaving only specific attributes unfrozen) we have a more
# granular control. Please, use setter to modify frozen class instead of removing frozen.

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
    "EWK_RAP_ASY": ("$\\eta/y$", "$M^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
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
    "SHP_ASY": ("$\\eta$", "$p_T (GeV)$", "$\\sqrt{s} (GeV)$"),
    "JET_POL": ("$\\eta$", "$p_T^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "DIJET_POL": ("$\\m_{1,2} (GeV)", "$\\eta_1$", "$\\eta_2$"),
    "DY_Z_Y": ("$y_Z$", "$\\M^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "SINGLETOP": ("$y$", "$m_t^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "DY_MLL": ("$M_{ll} (GeV)$", "$M_{ll}^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
    "DY_W_ETA": ("$abs_\\eta$", "$M_W^2 (GeV^2)$", "$\\sqrt{s} (GeV)$"),
}

PROCESS_DESCRIPTION_LABEL = {
    "EWJ_JRAP": "Jet Rapidity Distribution",
    "EWK_RAP": "Drell-Yan Rapidity Distribution",
    "EWJ_RAP": "Jet Rapidity Distribution",
    "HQP_PTQ": "Heavy Quarks Production Single Quark Transverse Momentum Distribution",
    "JET": "Jets Rapidity Distribution",
    "HIG_RAP": "Higgs Rapidity Distribution",
    "HQP_YQ": "Heavy Quarks Production Single Quark Rapidity Distribution",
    "EWJ_JPT": "Jet Transverse Momentum Distribution",
    "DIS": "Deep Inelastic Scattering",
    "HQP_PTQQ": "Heavy Quarks Production Transverse Momentum Distribution",
    "EWK_PT": "Drell-Yan Transverse Momentum Distribution",
    "EWJ_PT": "Jet Transverse Momentum Distribution",
    "PHT": "Photon Production",
    "HQP_MQQ": "Heavy Quarks Production Mass Distribution",
    "EWK_PTRAP": "Drell-Yan Transverse Momentum Distribution",
    "HQP_YQQ": "Heavy Quarks Production Rapidity Distribution",
    "INC": "Heavy Quarks Total Cross Section",
    "EWJ_MLL": "Jet Mass Distribution",
    "EWK_MLL": "Drell-Yan Mass Distribution",
    "DIJET": "Dijets Invariant Mass and Rapidity Distribution",
    "DYP": "Fixed-Target Drell-Yan",
    "JET_POL": "Inclusive Jet longitudinal double-spin asymmetry",
    "DIJET_POL": "Dijets longitudinal double-spin asymmetry",
    "SHP_ASY": "double spin asymmetry in single hadron production",
    "DY_MLL": "Drell-Yan Mass Distribution of Lepton Pairs",
    "DY_W_ETA": "Drell-Yan W boson rapidity distribution",
}


def _get_ported_kinlabel(process_type):
    """Get the kinematic label for ported datasets
    In principle there is a one to one correspondance between the process label in the kinematic and
    ``KINLABEL_LATEX``, however, there were some special cases that need to be taken into account
    """
    process_type = str(process_type)
    if process_type in KINLABEL_LATEX:
        return KINLABEL_LATEX[process_type]
    # special case in which the process in DIS- or DYP-like
    if process_type[:3] in ("DIS", "DYP"):
        return _get_ported_kinlabel(process_type[:3])
    if len(process_type.split("_")) > 1:
        return _get_process_description(process_type.rsplit("_", 1)[0])
    raise KeyError(f"Label {process_type} not recognized in KINLABEL_LATEX")


def _get_process_description(process_type):
    """Get the process description string for a given process type
    Similarly to kinlabel, some special cases are taken into account.
    """
    try:
        return process_type.description
    except AttributeError:
        # This process needs to be updated
        pass

    if process_type in PROCESS_DESCRIPTION_LABEL:
        return PROCESS_DESCRIPTION_LABEL[process_type]
    # If not, is this a DYP- or DIS-like dataset?
    if process_type[:3] in ("DIS", "DYP"):
        return _get_process_description(process_type[:3])
    # Remove pieces of "_" until it is found
    if len(process_type.split("_")) > 1:
        return _get_process_description(process_type.rsplit("_", 1)[0])
    raise KeyError(f"Label {process_type} not found in PROCESS_DESCRIPTION_LABEL")


@Parser
def ValidPath(path_str: str) -> Path:
    """Parse strings into paths"""
    return Path(path_str)


### Theory metadata
@Parser
def ValidOperation(op_str: Optional[str]) -> str:
    """Ensures that the operation defined in the commondata file is implemented in validphys"""
    if op_str is None:
        op_str = "NONE"
    ret = op_str.upper()

    # TODO: move accepted operations to this module so that the convolution receives an operation to apply
    # instead of an operation to understand
    try:
        from validphys.convolution import OP

        if ret not in OP:
            raise ValidationError(f"The operation '{op_str}' is not implemented in validphys")
    except ModuleNotFoundError:
        # Don't perform any checks if VP is not available
        pass

    return str(ret)


@dataclasses.dataclass(frozen=True)
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

    This class is inmutable, what is read from the commondata metadata should be considered final

    Example
    -------
    >>> from nnpdf_data.commondataparser import TheoryMeta
    ... from validobj import parse_input
    ... from ruamel.yaml import YAML
    ... theory_raw = '''
    ... FK_tables:
    ...   - - fk1
    ...   - - fk2
    ...     - fk3
    ... operation: ratio
    ... '''
    ... theory = YAML(typ='safe').load(theory_raw)
    ... parse_input(theory, TheoryMeta)
    TheoryMeta(FK_tables=[['fk1'], ['fk2', 'fk3']], operation='RATIO', shifts = None, conversion_factor=1.0, comment=None, normalization=None))
    """

    FK_tables: list[tuple]
    operation: ValidOperation = "NULL"
    conversion_factor: float = 1.0
    shifts: Optional[dict] = None
    normalization: Optional[dict] = None
    comment: Optional[str] = None

    def fktables_to_paths(self, grids_folder):
        """Given a source for pineappl grids, constructs the lists of fktables
        to be loaded"""
        ret = []
        for operand in self.FK_tables:
            ret.append([grids_folder / f"{m}.{EXT}" for m in operand])
        return ret

    @classmethod
    def parser(cls, yaml_file):
        """The yaml databases in the server use "operands" as key instead of "FK_tables" """
        if not yaml_file.exists():
            raise FileNotFoundError(yaml_file)
        meta = quick_yaml_load(yaml_file)
        # Make sure the operations are upper-cased for compound-compatibility
        meta["operation"] = "NULL" if meta["operation"] is None else meta["operation"].upper()
        if "operands" in meta:
            meta["FK_tables"] = meta.pop("operands")
        return parse_input(meta, cls)

    def __hash__(self):
        """Include in the hash all pieces of information available for functions using a cache.
        This includes also any possible comment in the off-chance it could have a meaning."""
        to_be_hashed = [self.operation, self.conversion_factor]
        to_be_hashed.append(tuple([tuple(i) for i in self.FK_tables]))
        if self.shifts is not None:
            to_be_hashed.append(tuple(self.shifts.keys()))
            to_be_hashed.append(tuple(self.shifts.values()))
        if self.normalization is not None:
            to_be_hashed.append(tuple(self.normalization.keys()))
            to_be_hashed.append(tuple(self.normalization.values()))
        if self.comment is not None:
            to_be_hashed.append(self.comment)
        return hash(tuple(to_be_hashed))


## Theory end


@dataclasses.dataclass(frozen=True)
class Variant:
    """The new commondata format allow the usage of variants
    A variant can overwrite a number of keys, as defined by this dataclass:
        data_uncertainties
        theory
        data_central

    This class may overwrite *some* other keys for the benefit of reproducibility
    of old NNPDF fits, but the usage of these features is undocumented and discouraged.
    """

    data_uncertainties: Optional[list[ValidPath]] = None
    theory: Optional[TheoryMeta] = None
    data_central: Optional[ValidPath] = None
    # Undocumented feature for *_DW_* only, where the nuclear uncertainties were included
    # as part of the experimental uncertianties instead of as a separate type
    experiment: Optional[str] = None


ValidVariants = dict[str, Variant]


### Kinematic data
@dataclasses.dataclass(frozen=True)
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


@dataclasses.dataclass(frozen=True)
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
    variables: dict[str, ValidVariable]

    def get_label(self, var):
        """For the given variable, return the label as label (unit)
        If the label is an "extra" return the last one
        """
        if var.startswith("extra_"):
            return list(self.variables.values())[-1]
        return self.variables[var].full_label()

    def apply_label(self, var, value):
        """For a given value for a given variable, return the labels
        as label = value (unit)
        If the variable is not included in the list of variables, returns None
        as the variable could've been transformed by a kinematic transformation
        """
        if var not in self.variables:
            return None
        return self.variables[var].apply_label(value)


### kinematics end


### Observable and dataset definitions
@dataclasses.dataclass(frozen=True, eq=True)
class ObservableMetaData:
    observable_name: str
    observable: dict
    ndata: int
    # Plotting options
    plotting: PlottingOptions
    process_type: ValidProcess
    kinematic_coverage: list[str]

    # Data itself
    kinematics: ValidKinematics
    data_uncertainties: list[ValidPath]

    # The central data is optional _only_ for
    # positivity datasets, and will be checked as soon as the class is instantiated
    data_central: Optional[ValidPath] = None

    # Optional data
    theory: Optional[TheoryMeta] = None
    tables: Optional[list] = dataclasses.field(default_factory=list)
    npoints: Optional[list] = dataclasses.field(default_factory=list)
    variants: Optional[ValidVariants] = dataclasses.field(default_factory=dict)
    applied_variant: Optional[str] = None
    ported_from: Optional[str] = None

    # Derived quantities:
    # Note that an observable without a parent will fail in many different ways
    _parent: Optional[Any] = None

    def __post_init__(self):
        """
        Small modifications for better compatibility with the rest of validphys
        """
        # Since vp will rely on the kinematics being 3 variables,
        # fill the extra with whatever can be found in the kinematics dictionary
        # otherwise just fill with extra_x
        if len(self.kinematic_coverage) < 3:
            unused = list(set(self.kinematics.variables) - set(self.kinematic_coverage))
            diff_to_3 = 3 - len(self.kinematic_coverage)
            if unused:
                nkincov = self.kinematic_coverage + unused[diff_to_3:]
            else:
                nkincov = self.kinematic_coverage + [f"extra_{i}" for i in range(diff_to_3)]
            object.__setattr__(self, 'kinematic_coverage', nkincov)

    def __hash__(self):
        """ObservableMetaData is defined by:
        - the setname
        - the variant used
        - the data
        """
        return hash((self.name, self.applied_variant, self.data_central))

    def check(self):
        """Various checks to apply manually to the observable before it is used anywhere
        These are not part of the __post_init__ call since they can only happen after the metadata
        has been read, the observable selected and (likely) variants applied.
        """
        # Check whether the data central or the uncertainties are empty for a non-positivity/integrability set
        if not self.is_nnpdf_special:
            if self.data_central is None:
                raise ValidationError(f"Missing `data_central` field for {self.name}")

            if not self.data_uncertainties:
                ermsg = f"""Missing `data_uncertainties` for {self.name}.
                    Select one of the variants: {list(self.variants.keys())}"""
                raise ValidationError(ermsg)

        # Check that plotting.plot_x is being filled
        if self.plotting.plot_x is None:
            ermsg = f"No variable selected as x-axis in the plot for {self.name}. Please add `plotting::plot_x`."
            if self.plotting.x is not None:
                ermsg += "Please replace `plotting::x` with `plotting::plot_x`."
            raise ValidationError(ermsg)

        # Ensure that all variables in the kinematic coverage exist
        for var in self.kinematic_coverage:
            if var not in self.kinematics.variables and not var.startswith("extra_"):
                raise ValidationError(
                    f"Variable {var} is in `kinematic_coverage` but not included in `kinematics`"
                    f", nor it is a variable of type `extra_` for {self.name} dataset."
                )

        if len(self.kinematic_coverage) > 3:
            raise ValidationError(
                "Only a maximum of 3 variables can be used for `kinematic_coverage`"
            )

    def apply_variant(self, variant_name):
        """Return a new instance of this class with the variant applied

        This class also defines how the variant is applied to the commondata
        """
        try:
            variant = self.variants[variant_name]
        except KeyError as e:
            raise ValueError(f"The requested variant does not exist {variant_name}") from e

        variant_replacement = {}
        if variant.data_uncertainties is not None:
            variant_replacement["data_uncertainties"] = variant.data_uncertainties
        if variant.theory is not None:
            variant_replacement["theory"] = variant.theory
        if variant.data_central is not None:
            variant_replacement["data_central"] = variant.data_central

        # This section should only be used for the purposes of reproducibility
        # of legacy data, no new data should use these

        if variant.experiment is not None:
            new_nnpdf_metadata = dict(self._parent.nnpdf_metadata.items())
            new_nnpdf_metadata["experiment"] = variant.experiment
            setmetadata_copy = dataclasses.replace(self._parent, nnpdf_metadata=new_nnpdf_metadata)
            variant_replacement["_parent"] = setmetadata_copy
            variant_replacement["plotting"] = dataclasses.replace(self.plotting)

        return dataclasses.replace(self, applied_variant=variant_name, **variant_replacement)

    @property
    def is_positivity(self):
        return self.setname.startswith("NNPDF_POS")

    @property
    def is_integrability(self):
        return self.setname.startswith("NNPDF_INTEG")

    @property
    def is_nnpdf_special(self):
        """Is this an NNPDF special dataset used for e.g., Lagrange multipliers or QED fits"""
        return self.setname.startswith("NNPDF")

    @property
    def path_data_central(self):
        return self._parent.folder / self.data_central

    def load_data_central(self):
        """Loads the data for this commondata returns a dataframe

        Returns
        -------
        pd.DataFrame
            a dataframe containing the data
        """
        if self.is_nnpdf_special:
            data = np.zeros(self.ndata)
        else:
            datayaml = quick_yaml_load(self.path_data_central)
            data = datayaml["data_central"]

        if len(data) != self.ndata:
            raise ValueError(
                f"The number of bins in {self.path_data_central} does not match ndata={self.ndata}"
            )

        data_df = pd.DataFrame(data, index=range(1, self.ndata + 1), columns=["data"])
        data_df.index.name = _INDEX_NAME
        return data_df

    @property
    def paths_uncertainties(self):
        return [self._parent.folder / i for i in self.data_uncertainties]

    def load_uncertainties(self):
        """Returns a dataframe with all appropiate uncertainties

        Returns
        -------
        pd.DataFrame
            a dataframe containing the uncertainties
        """
        if self.is_nnpdf_special:
            return pd.DataFrame([{}] * self.ndata, index=range(1, self.ndata + 1))

        all_df = []
        for ufile in self.paths_uncertainties:
            uncyaml = quick_yaml_load(ufile)
            mindex = pd.MultiIndex.from_tuples(
                [(k, v["treatment"], v["type"]) for k, v in uncyaml["definitions"].items()],
                names=["name", "treatment", "type"],
            )
            bin_list = pd.DataFrame(uncyaml["bins"]).values.astype(float)
            if len(bin_list) != self.ndata:
                raise ValueError(f"The number of bins in {ufile} does not match ndata={self.ndata}")

            # I'm guessing there will be a better way of doing this than calling  dataframe twice for the same thing?
            final_df = pd.DataFrame(bin_list, columns=mindex, index=range(1, self.ndata + 1))
            final_df.index.name = _INDEX_NAME
            all_df.append(final_df)
        return pd.concat(all_df, axis=1)

    @property
    def path_kinematics(self):
        return self._parent.folder / self.kinematics.file

    def load_kinematics(self, fill_to_three=True, drop_minmax=True):
        """Returns a dataframe with the kinematic information

        Parameters
        ----------
        fill_to_three: bool
            ensure that there are always three columns (repeat the last one) in the kinematics

        drop_minmax: bool
            Drop the min and max value, necessary for legacy comparisons

        Returns
        -------
        pd.DataFrame
            a dataframe containing the kinematics
        """
        kinematics_file = self.path_kinematics
        kinyaml = quick_yaml_load(kinematics_file)

        kin_dict = {}
        for bin_index, dbin in enumerate(kinyaml["bins"], start=1):
            for d in dbin.values():
                if d["mid"] is None:
                    d["mid"] = 0.5 * (d["max"] + d["min"])

                if drop_minmax:
                    d["min"] = None
                    d["max"] = None
                else:
                    # If we are not dropping it, ensure that it has something!
                    d["min"] = d["min"] if d.get("min") is not None else d["mid"]
                    d["max"] = d["max"] if d.get("max") is not None else d["mid"]

            # The old commondata always had 3 kinematic variables and the code sometimes
            # relies on this fact
            # Add a fake one at the end repeating the last one
            if fill_to_three and (ncol := len(dbin)) < 3:
                for i in range(3 - ncol):
                    dbin[f"extra_{i}"] = d

            kin_dict[bin_index] = pd.DataFrame(dbin).stack()

        if len(kin_dict) != self.ndata:
            raise ValueError(
                f"The number of bins in {kinematics_file} does not match ndata={self.ndata}"
            )

        return pd.concat(kin_dict, axis=1, names=[_INDEX_NAME]).swaplevel(0, 1).T

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
        return self._parent.cm_energy

    @property
    def name(self):
        return f"{self.setname}_{self.observable_name}"

    @property
    def is_ported_dataset(self):
        """Return True if this is an automatically ported dataset that has not been updated"""
        if self.ported_from is None:
            return False

        # If not using a legacy variant, we consider it ported if the kin variables are still k1,k2,k3
        return {"k1", "k2", "k3"} == set(self.kinematic_coverage)

    @property
    def kinlabels(self):
        """Return the kinematic labels in the same order as they are set
        in ``kinematic_coverage`` (which in turns follow the key kinematic_coverage)
        If this is a ported dataset, rely on the process type using the legacy labels
        """
        if self.is_ported_dataset:
            return _get_ported_kinlabel(self.process_type)
        return [self.kinematics.get_label(i) for i in self.kinematic_coverage]

    def digest_plotting_variable(self, variable):
        """Digest plotting variables in the ``line_by`` or ``figure_by`` fields
        and return the appropiate ``kX`` or other label such that the plotting functions
        of validphys can understand it.

        These might be variables included as part of the kinematics or extra labels
        defined in the plotting dictionary.
        """
        if not VP_AVAILABLE:
            raise ModuleNotFoundError(
                "validphys, from the full nnpdf package, needs to be installed to use this functionality"
            )

        # If it is part of the coverage, just return the relevant KN
        if variable in self.kinematic_coverage:
            fig_idx = self.kinematic_coverage.index(variable)
            return f"k{fig_idx + 1}"

        # If it is not in the coverage, it might be a _known_ extra label
        if self.plotting.extra_labels is not None and variable in self.plotting.extra_labels:
            # In that case return it raw
            return variable

        # Or, it might be a variable that VP knows how to deal with automagically
        if variable in labeler_functions:
            return variable

        raise ValueError(f"Don't know what to do with plotting variable {variable} for {self.name}")

    def _plotting_options_set(self):
        """Set and return the PlottingOptions metadata

        Fill in missing information that can be learnt from the other variables (xlabel/ylabel)
        or that is shared by the whole dataset.
        """
        if self.plotting.already_digested:
            return self.plotting

        if self.plotting.nnpdf31_process is None:
            self.plotting.nnpdf31_process = self.nnpdf_metadata["nnpdf31_process"]

        if self.plotting.experiment is None:
            self.plotting.experiment = self.nnpdf_metadata["experiment"]

        if self.plotting.process_description is None:
            self.plotting.process_description = _get_process_description(self.process_type)

        ## Swap variables by the k_idx
        # Internally validphys takes the x/y to be "k1" "k2" or "k3"
        # Therefore, for the time being, swap the actual keys by k1/k2/k3
        try:
            x_idx = self.kinematic_coverage.index(self.plotting.plot_x)
            self.plotting.x = f"k{x_idx + 1}"

            if self.plotting.x_label is None and not self.is_ported_dataset:
                self.plotting.x_label = self.kinematics.get_label(self.plotting.plot_x)

        except ValueError:
            # it is possible that the x value is an "extra", if that's the case continue
            self.plotting.x = self.plotting.plot_x
            self.plotting.x_label = None

        # Swap the `figure_by` and `line_by` variables by k1/k2/k3
        # unless this is something coming from the "extra labels"
        if self.plotting.figure_by is not None:
            new_fig_by = []
            for var in self.plotting.figure_by:
                new_fig_by.append(self.digest_plotting_variable(var))
            self.plotting.figure_by = new_fig_by

        if self.plotting.line_by is not None:
            new_line_by = []
            for var in self.plotting.line_by:
                new_line_by.append(self.digest_plotting_variable(var))
            self.plotting.line_by = new_line_by

        # And do it also within the normalize dictionary
        if self.plotting.normalize is not None:
            # Copy the normalize dictionary and update the figure and line by
            tmp = dict(self.plotting.normalize)
            tmp["figure_by"] = []
            tmp["line_by"] = []
            for var in self.plotting.normalize.get("figure_by", []):
                tmp["figure_by"].append(self.digest_plotting_variable(var))
            for var in self.plotting.normalize.get("line_by", []):
                tmp["line_by"].append(self.digest_plotting_variable(var))
            self.plotting.normalize = tmp

        self.plotting.already_digested = True
        return self.plotting

    @cached_property
    def plotting_options(self):
        try:
            return self._plotting_options_set()
        except Exception as e:
            # There are many chances for failure here
            log.error(f"Failure for: {self.name}")
            raise e


@dataclasses.dataclass(frozen=True)
class ValidReference:
    """Holds literature information for the dataset"""

    url: str
    version: Optional[int] = None
    journal: Optional[str] = None
    tables: list[int] = dataclasses.field(default_factory=list)


@dataclasses.dataclass(frozen=True)
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
        return path_commondata / self.setname

    @property
    def cm_energy(self):
        """Return the center of mass energy as GeV if it can be understood from the name
        otherwise return None"""
        energy_string = self.setname.split("_")[2]
        if energy_string == "NOTFIXED":
            return None
        if energy_string.endswith("GEV"):
            factor = 1.0
        elif energy_string.endswith("TEV"):
            factor = 1000
        else:
            return None
        return float(energy_string[:-3].replace("P", ".")) * factor

    @cached_property
    def allowed_datasets(self):
        """Return the implemented datasets as a list <setname>_<observable>"""
        return [f"{self.setname}_{i.observable_name}" for i in self.implemented_observables]

    @cached_property
    def allowed_observables(self):
        """
        Returns the implemented observables as a {observable_name.upper(): observable} dictionary
        """
        return {o.observable_name.upper(): o for o in self.implemented_observables}

    def select_observable(self, obs_name_raw):
        """Check whether the observable is implemented and return said observable"""
        obs_name = obs_name_raw.upper()
        try:
            observable = self.allowed_observables[obs_name]
        except KeyError:
            raise ValueError(
                f"The selected observable {obs_name_raw} does not exist in {self.setname}"
            )

        # Now burn the _parent key into the observable and apply checks
        object.__setattr__(observable, "_parent", self)
        return observable


@cache
def parse_set_metadata(metadata_file):
    """Read the metadata file"""
    return parse_yaml_inp(metadata_file, SetMetaData)


@cache
def parse_new_metadata(metadata_file, observable_name, variant=None):
    """Given a metadata file in the new format and the specific observable to be read
    load and parse the metadata and select the observable. If any variants are selected, apply them.

    The triplet (metadata_file, observable_name, variant) define unequivocally the information
    to be parsed from the commondata library
    """
    set_metadata = parse_set_metadata(metadata_file)

    # Select one observable from the entire metadata
    metadata = set_metadata.select_observable(observable_name)

    # And apply variant if given
    if variant is not None:
        metadata = metadata.apply_variant(variant)

    return metadata


def load_commondata(metadata):
    """

    TODO: update this docstring since now the load_commondata_new takes the information from
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
    # Before loading, apply the checks
    metadata.check()

    # Now parse the data
    data_df = metadata.load_data_central()
    # the uncertainties
    uncertainties_df = metadata.load_uncertainties()
    # and the kinematics
    kin_df = metadata.load_kinematics()

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
    systypes = {"treatment": [], "name": []}
    for col in uncertainties_df.columns:
        if col[0].startswith("stat"):
            new_columns.append("stat")
        else:
            # if it is syst add the ADD/MULT information
            new_columns.append(col[1])
            systypes["treatment"].append(col[1])
            systypes["name"].append(col[2])

    uncertainties_df.columns = new_columns

    commondata_table = pd.concat([kin_df, data_df, uncertainties_df], axis=1)
    systype_table = pd.DataFrame(systypes, index=range(1, nsys + 1))
    systype_table.index.name = "sys_index"

    # TODO: Legacy compatibility
    # 1. Add a stat column if it doesn't exist
    # 2. Transform multiplicative uncertainties into % as it was done in the older version

    if "stat" not in commondata_table:
        commondata_table["stat"] = 0.0

    if "MULT" in commondata_table:
        commondata_table["MULT"] = commondata_table["MULT"].multiply(
            100 / commondata_table["data"], axis="index"
        )

    # The old -> new map is not bijective, as different old dataset can refer to the same new one
    # therefore "legacy_names", when filled, will be a list. With None otherwise.
    legacy_names = new_to_legacy_map(metadata.name, metadata.applied_variant)

    return CommonData(
        setname=metadata.name,
        ndata=metadata.ndata,
        commondataproc=procname,
        nkin=3,
        nsys=nsys,
        commondata_table=commondata_table,
        systype_table=systype_table,
        legacy_names=legacy_names,
        kin_variables=metadata.kinematic_coverage,
    )
