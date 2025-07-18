"""
Filters for NNPDF fits
"""

from collections.abc import Mapping
import dataclasses
import functools
from importlib.resources import read_text
import logging
import re
from typing import Union

import numpy as np

from nnpdf_data.coredata import generate_path_filtered_data
from nnpdf_data.process_options import PROCESSES
from reportengine.checks import check, make_check
import validphys.cuts
from validphys.utils import yaml_safe

log = logging.getLogger(__name__)

KIN_LABEL = {
    "DIS": ("x", "Q2", "y"),
    "DYP": ("y", "M2", "sqrts"),
    "JET": ("eta", "pT2", "sqrts"),
    "DIJET": ("eta", "m_12", "sqrts"),
    "PHT": ("eta_gamma", "E_{T,gamma)2", "sqrts"),
    "INC": ("0", "mu2", "sqrts"),
    "SHP_ASY": ("eta", "pT", "sqrts"),
    "EWK_RAP": ("etay", "M2", "sqrts"),
    "EWK_RAP_ASY": ("etay", "M2", "sqrts"),
    "EWK_PT": ("p_T", "M2", "sqrts"),
    "EWK_PTRAP": ("etay", "pT2", "sqrts"),
    "EWK_MLL": ("M_ll", "M_ll2", "sqrts"),
    "EWJ_RAP": ("etay", "M2", "sqrts"),
    "EWJ_PT": ("p_T", "M2", "sqrt(s)"),
    "EWJ_PTRAP": ("etay", "pT2", "sqrts"),
    "EWJ_JRAP": ("etay", "M2", "sqrts"),
    "EWJ_JPT": ("p_T", "M2", "sqrts"),
    "EWJ_MLL": ("M_ll", "M_ll2", "sqrts"),
    "HQP_YQQ": ("yQQ", "mu2", "sqrts"),
    "HQP_MQQ": ("MQQ", "mu2", "sqrts"),
    "HQP_PTQQ": ("p_TQQ", "mu2", "sqrts"),
    "HQP_YQ": ("yQ", "mu2", "sqrts"),
    "HQP_PTQ": ("p_TQ", "mu2", "sqrts"),
    "HIG_RAP": ("y", "M_H2", "sqrts"),
    "SIA": ("z", "Q2", "y"),
    "JET_POL": ("eta", "pT2", "sqrts"),
    "DIJET_POL": ("m_12", "eta_1", "eta_2"),
}


def _get_kinlabel_process_type(process_type):
    """Get KIN_LABEL from the dictionary above according
    to the process type
    This requires some extra digestion for DIS
    """
    if isinstance(process_type, str):
        process_type = PROCESSES.get(process_type.upper(), process_type.upper())
    if hasattr(process_type, "accepted_variables"):
        return process_type.accepted_variables
    process_type = str(process_type)
    if process_type[:3] == "DIS":
        return KIN_LABEL["DIS"]
    if process_type[:3] == "DYP":
        return KIN_LABEL["DYP"]
    return KIN_LABEL[process_type]


# TODO: in the new commondata instead of having this, let's always use the same
# variables
def _variable_understanding(variables_raw, process_vars):
    """Given a set of variable, check whether it might be a variation of existing
    variables for a process type"""
    variables = [i for i in variables_raw]

    def substitute(pr_v, cd_x):
        if pr_v in process_vars and cd_x in variables:
            variables[variables.index(cd_x)] = pr_v

    substitute("eta", "y")
    substitute("eta", "eta")
    substitute("etay", "eta")
    substitute("etay", "y")
    substitute("yQQ", "y_ttBar")
    substitute("yQ", "y_t")

    return variables


class RuleProcessingError(Exception):
    """Exception raised when we couldn't process a rule."""


class BadPerturbativeOrder(ValueError):
    """Exception raised when the perturbative order string is not
    recognized."""


class MissingRuleAttribute(RuleProcessingError, AttributeError):
    """Exception raised when a rule is missing required attributes."""


class FatalRuleError(Exception):
    """Exception raised when a rule application failed at runtime."""


@dataclasses.dataclass(frozen=True)
class FilterDefaults:
    """
    Dataclass carrying default values for filters (cuts) taking into
    account the values of ``q2min``, ``w2min`` and ``maxTau``.
    """

    q2min: float = None
    w2min: float = None
    maxTau: float = None

    def to_dict(self):
        return dataclasses.asdict(self)


@dataclasses.dataclass(frozen=True)
class FilterRule:
    """
    Dataclass which carries the filter rule information.

    """

    dataset: str = None
    process_type: str = None
    rule: str = None
    reason: str = dataclasses.field(
        default=None, hash=False, compare=False
    )  # Not relevant for hashing
    local_variables: Mapping[str, Union[str, float]] = dataclasses.field(
        default=None, hash=False
    )  # Avoid hash issues with caching
    PTO: str = None
    FNS: str = None
    IC: str = None

    def to_dict(self):
        rule_dict = dataclasses.asdict(self)
        filtered_dict = {k: v for k, v in rule_dict.items() if v is not None}
        return filtered_dict


@dataclasses.dataclass(frozen=True)
class AddedFilterRule(FilterRule):
    """
    Dataclass which carries extra filter rule that is added to the
    default rule.
    """

    pass


def default_filter_settings_input():
    """Return a FilterDefaults dataclass with the default hardcoded filter settings.
    These are defined in ``defaults.yaml`` in the ``validphys.cuts`` module.
    """
    return FilterDefaults(**yaml_safe.load(read_text(validphys.cuts, "defaults.yaml")))


def default_filter_rules_input():
    """
    Return a tuple of FilterRule objects.
    These are defined in ``filters.yaml`` in the ``validphys.cuts`` module.
    Similarly to `parse_added_filter_rules`, this function checks if the rules
    are unique, i.d. if there are no multiple rules for the same dataset of
    process with the same rule (`reason` and `local_variables` are not hashed).
    """
    list_rules = yaml_safe.load(read_text(validphys.cuts, "filters.yaml"))
    unique_rules = set(FilterRule(**rule) for rule in list_rules)
    if len(unique_rules) != len(list_rules):
        raise RuleProcessingError(
            "Detected repeated filter rules. Please, make sure that "
            " rules are not repeated in `filters.yaml`."
        )
    return tuple(unique_rules)


def check_nonnegative(var: str):
    """Ensure that `var` is positive"""

    @make_check
    def run_check(ns, **kwargs):
        val = ns[var]
        check(val >= 0, f"'{var}' must be positive or equal zero, but it is {val!r}.")

    return run_check


def make_dataset_dir(path):
    """Creates directory at path location."""
    if path.exists():
        log.warning(f"Dataset output folder exists: {path} Overwriting contents")
    else:
        path.mkdir(exist_ok=True)


def export_mask(path, mask):
    """Dump mask to file"""
    np.savetxt(path, mask, fmt='%d')


def filter_closure_data_by_experiment(
    filter_path, experiments_data, fakepdf, fakenoise, filterseed, data_index, sep_mult
):
    """
    Applies :py:func:`_filter_closure_data` on each experiment in the closure test.

    This function just peforms a ``for`` loop over ``experiments``, the reason
    we don't use ``reportengine.collect`` is that it can permute the order
    in which closure data is generate, which means that the pseudodata is
    not reproducible.

    """

    res = []
    for exp in experiments_data:
        experiment_index = data_index[data_index.isin([exp.name], level=0)]
        res.append(
            _filter_closure_data(
                filter_path,
                exp,
                fakepdf,
                fakenoise,
                filterseed,
                experiment_index,
                sep_mult,
                inconsistent_data_settings=None,
            )
        )

    return res


def filter_inconsistent_closure_data_by_experiment(
    filter_path,
    experiments_data,
    fakepdf,
    fakenoise,
    filterseed,
    data_index,
    sep_mult,
    inconsistent_data_settings,
):
    """
    Like :py:func:`filter_closure_data_by_experiment` except for inconsistent closure tests.
    """

    res = []
    for exp in experiments_data:
        experiment_index = data_index[data_index.isin([exp.name], level=0)]
        res.append(
            _filter_closure_data(
                filter_path,
                exp,
                fakepdf,
                fakenoise,
                filterseed,
                experiment_index,
                sep_mult,
                inconsistent_data_settings,
            )
        )

    return res


def filter_real_data(filter_path, data):
    """Filter real data, cutting any points which do not pass the filter rules."""
    log.info('Filtering real data.')
    return _filter_real_data(filter_path, data)


def filter(filter_data):
    """Summarise filters applied to all datasets"""
    total_data, total_cut_data = np.atleast_2d(filter_data).sum(axis=0)
    log.info(f'Summary: {total_cut_data}/{total_data} datapoints passed kinematic cuts.')


def _write_ds_cut_data(path, dataset):
    make_dataset_dir(path)
    all_dsndata = dataset.commondata.ndata
    datamask = dataset.cuts.load()
    if datamask is None:
        filtered_dsndata = all_dsndata
        log.info("All {all_ndata} points  in in {dataset.name} passed kinematic cuts.")
    else:
        filtered_dsndata = len(datamask)
        log.info(
            f"{len(datamask)}/{all_dsndata} datapoints " f"in {dataset.name} passed kinematic cuts."
        )
    # save to disk
    if datamask is not None:
        export_mask(path / f'FKMASK_{dataset.name}.dat', datamask)
    return all_dsndata, filtered_dsndata


def _filter_real_data(filter_path, data):
    """Filter real experimental data."""

    total_data_points = 0
    total_cut_data_points = 0
    for dataset in data.datasets:
        path = filter_path / dataset.name
        nfull, ncut = _write_ds_cut_data(path, dataset)
        total_data_points += nfull
        total_cut_data_points += ncut
        dataset.load_commondata().export(path)
    return total_data_points, total_cut_data_points


def _filter_closure_data(
    filter_path,
    data,
    fakepdf,
    fakenoise,
    filterseed,
    data_index,
    sep_mult,
    inconsistent_data_settings=None,
):
    """
    This function is accessed within a closure test only, that is, the fakedata
    namespace has to be True (If fakedata = False, the _filter_real_data function
    will be used to write the commondata files).

    The function writes commondata and systypes files within the
    name_closure_test/filter folder.
    If fakenoise is True, Level 1 type data is written to the filter folder, otherwise
    Level 0 data is written.

    Level 1 data is generated from the Level 0 data by adding noise sampled from
    the experimental covariance matrix using the validphys.pseudodata.make_replica
    function.

    Parameters
    ----------

    filter_path : str
                  path to filter folder

    data : validphys.core.DataGroupSpec

    fakepdf : validphys.core.PDF

    fakenoise : bool
                if fakenoise perform level1 shift of central data values

    filterseed : int
                 random seed used for the generation of
                 random noise added to Level 0 data


    data_index : pandas.MultiIndex

    inconsistent_data_settings: dict, default is None
        Settings for the inconsistent closure test. If None, the closure test is consistent.

    Returns
    -------
    tuple
         total data points and points passing the cuts

    """
    total_data_points = 0
    total_cut_data_points = 0

    # circular import generated @ core.py
    from validphys.pseudodata import level0_commondata_wc, make_level1_data

    closure_data = level0_commondata_wc(data, fakepdf)

    # Keep track of the original commondata, since it is what will be used to export
    # the data afterwards
    all_raw_commondata = {}

    for dataset in data.datasets:
        # == print number of points passing cuts, make dataset directory and write FKMASK  ==#
        path = filter_path / dataset.name
        nfull, ncut = _write_ds_cut_data(path, dataset)
        total_data_points += nfull
        total_cut_data_points += ncut
        all_raw_commondata[dataset.name] = dataset.commondata.load()

    if fakenoise:
        # ======= Level 1 closure test =======#

        closure_data = make_level1_data(data, closure_data, filterseed, data_index, sep_mult)

        if inconsistent_data_settings is not None:
            # avoid circular import
            from validphys.closuretest.inconsistent_closuretest.inconsistent_ct import (
                InconsistentCommonData,
            )

            # Convert the commondata to InconsistentCommonData
            closure_data = [
                InconsistentCommonData(
                    setname=cd.setname,
                    ndata=cd.ndata,
                    commondataproc=cd.commondataproc,
                    nkin=cd.nkin,
                    nsys=cd.nsys,
                    commondata_table=cd.commondata_table,
                    systype_table=cd.systype_table,
                )
                for cd in closure_data
            ]

            # Process the commondata
            closure_data = [
                cd.process_commondata(**inconsistent_data_settings) for cd in closure_data
            ]

        log.info("Writing Level1 data")

    else:
        log.info("Writing Level0 data")

    # ====== write commondata and systype files ======#
    for cd in closure_data:
        # Write the full dataset, not only the points that pass the filter
        data_path, unc_path = generate_path_filtered_data(filter_path.parent, cd.setname)
        data_path.parent.mkdir(exist_ok=True, parents=True)

        raw_cd = all_raw_commondata[cd.setname]

        data_range = np.arange(1, 1 + raw_cd.ndata)

        # Now put the closure data into the raw original commondata
        new_cv = cd.central_values.reindex(data_range, fill_value=0.0).values
        output_cd = raw_cd.with_central_value(new_cv)

        # And export it to file
        with open(data_path, "w", encoding="utf-8") as f:
            output_cd.export_data(f)

        if inconsistent_data_settings:
            if cd.setname in inconsistent_data_settings["inconsistent_datasets"]:
                # convert output_cd to InconsistentCommonData  (which has systematic_errors as a property and it's own export_uncertainties method)
                # this is done for the inconsistent datasets only as the systematics of the other datasets are not modified

                output_cd = InconsistentCommonData(
                    setname=output_cd.setname,
                    ndata=output_cd.ndata,
                    commondataproc=output_cd.commondataproc,
                    nkin=output_cd.nkin,
                    nsys=output_cd.nsys,
                    commondata_table=output_cd.commondata_table,
                    systype_table=output_cd.systype_table,
                )

                # put the inconsistent closure systematics into the raw systematics
                new_sys = cd.systematic_errors.reindex(data_range, fill_value=0.0)
                output_cd.systematic_errors = new_sys

        # export it to file
        with open(unc_path, "w", encoding="utf-8") as f:
            output_cd.export_uncertainties(f)

    return total_data_points, total_cut_data_points


def check_t0pdfset(t0pdfset):
    """T0 pdf check"""
    t0pdfset.load()
    log.info(f'{t0pdfset} T0 checked.')


def check_luxset(luxset):
    """Lux pdf check"""
    luxset.load()
    log.info(f'{luxset} Lux pdf checked.')


def check_unpolarized_bc(unpolarized_bc):
    """Check that unpolarized PDF bound can be loaded normally."""
    unpolarized_bc.load()
    log.info(f"Unpolarized boundary PDF '{unpolarized_bc}' pdf checked.")


def check_additional_errors(additional_errors):
    """Lux additional errors pdf check"""
    additional_errors.load()
    log.info(f'{additional_errors} Lux additional errors pdf checked.')


def check_positivity(posdatasets):
    """Verify positive datasets are ready for the fit."""
    log.info('Verifying positivity tables:')
    for pos in posdatasets:
        pos.load_commondata()
        log.info(
            f'{pos.name} checked, {len(pos.cuts.load())}/{pos.commondata.ndata} datapoints passed kinematic cuts.'
        )


def check_integrability(integdatasets):
    """Verify positive datasets are ready for the fit."""
    log.info('Verifying integrability tables:')
    for integ in integdatasets:
        integ.load_commondata()
        log.info(f'{integ.name} checked.')


class PerturbativeOrder:
    """Class that conveniently handles
    perturbative order declarations for use
    within the Rule class filter.


    Parameters
    ----------
    string: str
        A string in the format of NNLO or equivalently N2LO.
        This can be followed by one of ! + - or none.

        The syntax allows for rules to be executed only if the perturbative
        order is within a given range. The following enumerates all 4 cases
        as an example:

        NNLO+ only execute the following rule if the pto is 2 or greater
        NNLO- only execute the following rule if the pto is strictly less than 2
        NNLO! only execute the following rule if the pto is strictly not 2
        NNLO only execute the following rule if the pto is exactly 2

        Any unrecognized string will raise a BadPerturbativeOrder exception.

    Example
    -------
    >>> from validphys.filters import PerturbativeOrder
    >>> pto = PerturbativeOrder("NNLO+")
    >>> pto.numeric_pto
    2
    >>> 1 in pto
    False
    >>> 2 in pto
    True
    >>> 3 in pto
    True
    """

    def __init__(self, string):
        self.string = string.upper()
        self.parse()

    def parse(self):
        # Change an input like NNNLO or N3LO
        # to a numerical value for the pto.
        # In this example, we assign
        # self.numeric_pto to be 3.
        exp = re.compile(
            r"(N(?P<nnumber>\d+)|(?P<multiplens>N*))LO(?P<operator>[\+\-\!])?"
        ).fullmatch(self.string)
        if not exp:
            raise BadPerturbativeOrder(
                f"String {self.string!r} does not represent a valid perturbative order specfication."
            )
        if exp.group("multiplens") is None:
            self.numeric_pto = int(exp.group("nnumber"))
        else:
            self.numeric_pto = len(exp.group("multiplens"))

        self.operator = exp.group("operator")

    def __contains__(self, i):
        if self.operator == "!":
            return i != self.numeric_pto
        elif self.operator == "+":
            return i >= self.numeric_pto
        elif self.operator == "-":
            return i < self.numeric_pto
        else:
            return i == self.numeric_pto


class Rule:
    """Rule object to be used to generate cuts mask.

    A rule object is created for each rule in ./cuts/filters.yaml

    Old commondata relied on the order of the kinematical variables
    to be the same as specified in the `KIN_LABEL` dictionary set in this module.
    The new commondata specification instead defines explicitly the name of the
    variables in the metadata.
    Therefore, when using a new-format commondata, the KIN_LABEL dictionary
    will not be used and the variables defined in it will be used instead.

    Parameters
    ----------
    initial_data: dict
        A dictionary containing all the information regarding the rule.
        This contains the name of the dataset the rule to applies to
        and/or the process type the rule applies to. Additionally, the
        rule itself is defined, alongside the reason the rule is used.
        Finally, the user can optionally define their own custom local
        variables.

        By default these are defined in cuts/filters.yaml
    defaults: dict
        A dictionary containing default values to be used globally in
        all rules.

        By default these are defined in cuts/defaults.yaml
    theory_parameters:
        Dict containing pairs of (theory_parameter, value)
    loader: validphys.loader.Loader, optional
        A loader instance used to retrieve the datasets.
    """

    numpy_functions = {"sqrt": np.sqrt, "log": np.log, "fabs": np.fabs}

    def __init__(
        self, initial_data: FilterRule, *, defaults: dict, theory_parameters: dict, loader=None
    ):
        self.dataset = None
        self.process_type = None
        self._local_variables_code = {}

        # For compatibility with legacy code that passed a dictionary
        if isinstance(initial_data, FilterRule):
            initial_data = initial_data.to_dict()
        elif isinstance(initial_data, Mapping):
            initial_data = dict(initial_data)
        else:
            raise RuleProcessingError(
                "Expecting initial_data to be an instance of a FilterRule dataclass."
            )

        for key in initial_data:
            setattr(self, key, initial_data[key])

        if not hasattr(self, "rule"):
            raise MissingRuleAttribute("No rule defined.")

        if self.dataset is None and self.process_type is None:
            raise MissingRuleAttribute("Please define either a process type or dataset.")

        if self.process_type is None:
            from validphys.loader import Loader, LoaderError

            if loader is None:
                loader = Loader()

            try:
                cd = loader.check_commondata(self.dataset)
            except LoaderError as e:
                raise RuleProcessingError(f"Could not find dataset {self.dataset}") from e

            if cd.metadata.is_ported_dataset:
                self.variables = _get_kinlabel_process_type(cd.process_type)
            else:
                self.variables = cd.metadata.kinematic_coverage
        else:
            self.variables = _get_kinlabel_process_type(self.process_type)
            # TODO: for now this will be a string within this class

        if hasattr(self, "local_variables"):
            if not isinstance(self.local_variables, Mapping):
                raise RuleProcessingError(
                    f"Expecting local_variables to be a Mapping, not {type(self.local_variables)}."
                )
        else:
            self.local_variables = {}

        if hasattr(self, "PTO"):
            if not isinstance(self.PTO, str):
                raise RuleProcessingError(f"Expecting PTO to be a string, not {type(self.PTO)}.")
            try:
                self.PTO = PerturbativeOrder(self.PTO)
            except BadPerturbativeOrder as e:
                raise RuleProcessingError(e) from e

        self.rule_string = self.rule
        self.defaults = defaults
        self.theory_params = theory_parameters
        ns = {
            *self.numpy_functions,
            *self.defaults.to_dict().keys(),
            *self.variables,
            "idat",
            "central_value",
        }
        for k, v in self.local_variables.items():
            try:
                self._local_variables_code[k] = lcode = compile(
                    str(v), f"local variable {k}", "eval"
                )
            except Exception as e:
                raise RuleProcessingError(
                    f"Could not process local variable {k!r} ({v!r}): {e}"
                ) from e
            for name in lcode.co_names:
                if name not in ns:
                    raise RuleProcessingError(
                        f"Could not process local variable {k!r}: Unknown name {name!r}"
                    )
            ns.add(k)

        try:
            self.rule = compile(self.rule, "rule", "eval")
        except Exception as e:
            raise RuleProcessingError(f"Could not process rule {self.rule_string!r}: {e}") from e
        for name in self.rule.co_names:
            if name not in ns:
                raise RuleProcessingError(
                    f"Could not process rule {self.rule_string!r}: Unknown name {name!r}"
                )

        # Before returning, set the process type as a string for the rest of the filter
        if self.process_type is not None:
            self.process_type = str(self.process_type)

    @property
    def _properties(self):
        """Attributes of the Rule class that are defining. Two
        Rules with identical ``_properties`` are considered equal.
        """
        return (
            self.rule_string,
            self.dataset,
            self.process_type,
            self.theory_params['ID'],
            tuple(self.local_variables.items()),
        )

    def __eq__(self, other):
        return self._properties == other._properties

    def __hash__(self):
        return hash(self._properties)

    def __call__(self, dataset, idat):
        central_value = dataset.get_cv()[idat]
        process_name = dataset.commondataproc

        # We return None if the rule doesn't apply. This
        # is different to the case where the rule does apply,
        # but the point was cut out by the rule.
        if (
            not (
                dataset.setname == self.dataset
                or
                # for old rules, we might be using an old name of this dataset
                (dataset.legacy_names is not None and self.dataset in dataset.legacy_names)
            )
            and process_name != self.process_type
            and self.process_type != "DIS_ALL"
        ):
            return None

        # Handle the generalised DIS cut
        if self.process_type == "DIS_ALL" and not process_name.startswith("DIS"):
            return None

        ns = self._make_point_namespace(dataset, idat)
        for k, v in self.theory_params.items():
            if k == "PTO" and hasattr(self, "PTO"):
                if v not in self.PTO:
                    return None
            elif hasattr(self, k) and (getattr(self, k) != v):
                return None

        # Will return True if datapoint passes through the filter
        try:
            return eval(
                self.rule,
                self.numpy_functions,
                {**{"idat": idat, "central_value": central_value}, **self.defaults.to_dict(), **ns},
            )
        except Exception as e:  # pragma: no cover
            raise FatalRuleError(f"Error when applying rule {self.rule_string!r}: {e}") from e

    def __repr__(self):  # pragma: no cover
        return self.rule_string

    def _make_kinematics_dict(self, dataset, idat) -> dict:
        """Fill in a dictionary with the kinematics for each point"""
        # TODO
        # When applying a "process-type" rule the variables are as given
        # at the top of the module. However, for new commondata is important
        # that the variables come in the right order
        # This "understanding" should not be necessary and the process-variable
        # mapping in this module should only serve to check which variables are allowed
        kinematics = dataset.kinematics.values[idat]
        if dataset.legacy or "k1" in dataset.kin_variables:
            # For ported dataset, the naming is k1/k2/k3 and
            # thus it needs to rely on the process
            return dict(zip(self.variables, kinematics))

        # Use the order of the commondata and the sintax of KIN_LABEL
        new_vars = _variable_understanding(dataset.kin_variables, self.variables)
        return dict(zip(new_vars, kinematics))

    def _make_point_namespace(self, dataset, idat) -> dict:
        """Return a dictionary with kinematics and local
        variables evaluated for each point"""
        ns = self._make_kinematics_dict(dataset, idat)

        for key, value in self._local_variables_code.items():
            ns[key] = eval(value, {**self.numpy_functions, **ns})
        return ns


@functools.lru_cache
def get_cuts_for_dataset(commondata, rules) -> list:
    """Function to generate a list containing the index
    of all experimental points that passed kinematic
    cut rules stored in ./cuts/filters.yaml


    Parameters
    ----------
    commondata: nnpdf_data.coredata.CommonData
    rules: List[Rule]
        A list of Rule objects specifying the filters.

    Returns
    -------
    mask: list
        List object containing index of all passed experimental
        values

    Example
    -------
    >>> from validphys.filters import (get_cuts_for_dataset, Rule,
    ...     default_filter_settings, default_filter_rules_input)
    >>> from validphys.loader import Loader
    >>> l = Loader()
    >>> cd = l.check_commondata("NMC")
    >>> theory = l.check_theoryID(53)
    >>> filter_defaults = default_filter_settings()
    >>> params = theory.get_description()
    >>> rule_list = [Rule(initial_data=i, defaults=filter_defaults, theory_parameters=params)
    ...     for i in default_filter_rules_input()]
    >>> get_cuts_for_dataset(cd, rules=rule_list)
    """
    dataset = commondata.load()

    mask = []
    for idat in range(dataset.ndata):
        broken = False
        for rule in rules:
            rule_result = rule(dataset, idat)
            if rule_result is not None and not rule_result:
                broken = True
                break

        if not broken:
            mask.append(idat)

    return mask
