"""
Filters for NNPDF fits
"""

import logging
import re
from collections.abc import Mapping
from importlib.resources import read_text

import numpy as np

from NNPDF import CommonData
from reportengine.checks import make_argcheck, check, check_positive, make_check
from reportengine.compat import yaml
import validphys.cuts

log = logging.getLogger(__name__)

class RuleProcessingError(Exception):
    """Exception raised when we couldn't process a rule."""


class BadPerturbativeOrder(ValueError):
    """Exception raised when the perturbative order string is not
    recognized."""

class MissingRuleAttribute(RuleProcessingError, AttributeError):
    """Exception raised when a rule is missing required attributes."""


class FatalRuleError(Exception):
    """Exception raised when a rule application failed at runtime."""


def default_filter_settings_input():
    """Return a dictionary with the default hardcoded filter settings.
    These are defined in ``defaults.yaml`` in the ``validphys.cuts`` module.
    """
    return yaml.safe_load(read_text(validphys.cuts, "defaults.yaml"))


def default_filter_rules_input():
    """Return a dictionary with the input settings.
    These are defined in ``filters.yaml`` in the ``validphys.cuts`` module.
    """
    return yaml.safe_load(read_text(validphys.cuts, "filters.yaml"))


@make_argcheck
def check_rngalgo(rngalgo: int):
    """Check rngalgo content"""
    check(0 <= rngalgo < 17,
          "Invalid rngalgo. Must be int between [0, 16].")


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

@check_rngalgo
@check_nonnegative('filterseed')
@check_nonnegative('seed')
def prepare_nnpdf_rng(filterseed:int, rngalgo:int, seed:int):
    """Initialise the internal NNPDF RNG, specified by ``rngalgo`` which must
    be an integer between 0 and 16, seeded with ``filterseed``.
    The RNG can then be subsequently used to i.e generate pseudodata.
    """
    try:
        from NNPDF import RandomGenerator
    except ImportError as e:
        logging.error("Generating closure data needs a valid installation of libNNPDF")
        raise e

    log.warning("Importing libNNPDF")
    log.info("Initialising RNG")
    RandomGenerator.InitRNG(rngalgo, seed)
    RandomGenerator.GetRNG().SetSeed(filterseed)

@check_positive('errorsize')
def filter_closure_data(filter_path, data, t0pdfset, fakenoise, errorsize, prepare_nnpdf_rng):
    """Filter closure data. In addition to cutting data points, the data is
    generated from an underlying ``t0pdfset``, applying a shift to the data
    if ``fakenoise`` is ``True``, which emulates the experimental central values
    being shifted away from the underlying law.

    """
    log.info('Filtering closure-test data.')
    return _filter_closure_data(
        filter_path, data, t0pdfset, fakenoise, errorsize)


@check_positive("errorsize")
def filter_closure_data_by_experiment(
    filter_path, experiments_data, t0pdfset, fakenoise, errorsize, prepare_nnpdf_rng,
):
    """
    Like :py:func:`filter_closure_data` except filters data by experiment.

    This function just peforms a ``for`` loop over ``experiments``, the reason
    we don't use ``reportengine.collect`` is that it can permute the order
    in which closure data is generate, which means that the pseudodata is
    not reproducible.

    """
    return [
        _filter_closure_data(filter_path, exp, t0pdfset, fakenoise, errorsize)
        for exp in experiments_data
    ]


def filter_real_data(filter_path, data):
    """Filter real data, cutting any points which do not pass the filter rules.
    """
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
        log.info(f"{len(datamask)}/{all_dsndata} datapoints "
                 f"in {dataset.name} passed kinematic cuts.")
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
        dataset.load().Export(str(path))
    return total_data_points, total_cut_data_points


def _filter_closure_data(filter_path, data, fakepdfset, fakenoise, errorsize):
    """Filter closure test data."""
    total_data_points = 0
    total_cut_data_points = 0
    fakeset = fakepdfset.legacy_load()
    # Load data, don't cache result
    loaded_data = data.load.__wrapped__(data)
    # generate level 1 shift if fakenoise
    loaded_data.MakeClosure(fakeset, fakenoise)
    for j, dataset in enumerate(data.datasets):
        path = filter_path / dataset.name
        nfull, ncut = _write_ds_cut_data(path, dataset)
        total_data_points += nfull
        total_cut_data_points += ncut
        loaded_ds = loaded_data.GetSet(j)
        if errorsize != 1.0:
            loaded_ds.RescaleErrors(errorsize)
        loaded_ds.Export(str(path))
    return total_data_points, total_cut_data_points


def check_t0pdfset(t0pdfset):
    """T0 pdf check"""
    t0pdfset.load()
    log.info(f'{t0pdfset} T0 checked.')


def check_positivity(posdatasets):
    """Verify positive datasets are ready for the fit."""
    log.info('Verifying positivity tables:')
    for pos in posdatasets:
        pos.load()
        log.info(f'{pos.name} checked.')

def check_integrability(integdatasets):
    """Verify positive datasets are ready for the fit."""
    log.info('Verifying integrability tables:')
    for integ in integdatasets:
        integ.load()
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
        self,
        initial_data: dict,
        *,
        defaults: dict,
        theory_parameters: dict,
        loader=None,
    ):
        self.dataset = None
        self.process_type = None
        self._local_variables_code = {}
        for key in initial_data:
            setattr(self, key, initial_data[key])

        if not hasattr(self, "rule"):
            raise MissingRuleAttribute("No rule defined.")

        if self.dataset is None and self.process_type is None:
            raise MissingRuleAttribute(
                "Please define either a process type or dataset."
            )

        if self.process_type is None:
            from validphys.loader import Loader, LoaderError

            if loader is None:
                loader = Loader()
            try:
                cd = loader.check_commondata(self.dataset)
            except LoaderError as e:
                raise RuleProcessingError(
                    f"Could not find dataset {self.dataset}"
                ) from e
            if cd.process_type[:3] == "DIS":
                self.variables = CommonData.kinLabel["DIS"]
            else:
                self.variables = CommonData.kinLabel[cd.process_type]
        else:
            if self.process_type[:3] == "DIS":
                self.variables = CommonData.kinLabel["DIS"]
            else:
                self.variables = CommonData.kinLabel[self.process_type]

        if hasattr(self, "local_variables"):
            if not isinstance(self.local_variables, Mapping):
                raise RuleProcessingError(
                    f"Expecting local_variables to be a Mapping, not {type(self.local_variables)}."
                )
        else:
            self.local_variables = {}

        if hasattr(self, "PTO"):
            if not isinstance(self.PTO, str):
                raise RuleProcessingError(
                    f"Expecting PTO to be a string, not {type(self.PTO)}."
                )
            try:
                self.PTO = PerturbativeOrder(self.PTO)
            except BadPerturbativeOrder as e:
                raise RuleProcessingError(e) from e

        self.rule_string = self.rule
        self.defaults = defaults
        self.theory_params = theory_parameters
        ns = {
            *self.numpy_functions,
            *self.defaults,
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
            raise RuleProcessingError(
                f"Could not process rule {self.rule_string!r}: {e}"
            ) from e
        for name in self.rule.co_names:
            if name not in ns:
                raise RuleProcessingError(
                    f"Could not process rule {self.rule_string!r}: Unknown name {name!r}"
                )

    @property
    def _properties(self):
        """Attributes of the Rule class that are defining. Two
        Rules with identical ``_properties`` are considered equal.
        """
        return (self.rule_string, self.dataset, self.process_type, self.theory_params['ID'])

    def __eq__(self, other):
        return self._properties == other._properties

    def __hash__(self):
        return hash(self._properties)

    def __call__(self, dataset, idat):
        central_value = dataset.GetData(idat)
        # We return None if the rule doesn't apply. This
        # is different to the case where the rule does apply,
        # but the point was cut out by the rule.
        if (
            dataset.GetSetName() != self.dataset
            and dataset.GetProc(idat) != self.process_type
            and self.process_type != "DIS_ALL"
        ):
            return None

        # Handle the generalised DIS cut
        if self.process_type == "DIS_ALL" and dataset.GetProc(idat)[:3] != "DIS":
            return None

        ns = self._make_point_namespace(dataset, idat)
        for k, v in self.theory_params.items():
            if k == "PTO" and hasattr(self, "PTO"):
                if v not in self.PTO:
                    return None
            elif hasattr(self, k) and (
                getattr(self, k) != v
            ):
                return None

        # Will return True if datapoint passes through the filter
        try:
            return eval(
                self.rule,
                self.numpy_functions,
                {
                    **{"idat": idat, "central_value": central_value},
                    **self.defaults,
                    **ns,
                },
            )
        except Exception as e: # pragma: no cover
            raise FatalRuleError(
                f"Error when applying rule {self.rule_string!r}: {e}"
            ) from e

    def __repr__(self): # pragma: no cover
        return self.rule_string

    def _make_kinematics_dict(self, dataset, idat) -> dict:
        """Fill in a dictionary with the kinematics for each point"""
        kinematics = [dataset.GetKinematics(idat, j) for j in range(3)]
        return dict(zip(self.variables, kinematics))

    def _make_point_namespace(self, dataset, idat) -> dict:
        """Return a dictionary with kinematics and local
        variables evaluated for each point"""
        ns = self._make_kinematics_dict(dataset, idat)

        for key, value in self._local_variables_code.items():
            ns[key] = eval(value, {**self.numpy_functions, **ns})
        return ns

def get_cuts_for_dataset(commondata, rules) -> list:
    """Function to generate a list containing the index
    of all experimental points that passed kinematic
    cut rules stored in ./cuts/filters.yaml


    Parameters
    ----------
    commondata: NNPDF CommonData spec
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
    for idat in range(dataset.GetNData()):
        broken = False
        for rule in rules:
            rule_result = rule(dataset, idat)
            if rule_result is not None and not rule_result:
                broken = True
                break

        if not broken:
            mask.append(idat)

    return mask
