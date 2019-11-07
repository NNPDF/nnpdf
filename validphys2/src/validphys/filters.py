# -*- coding: utf-8 -*-
"""
Filters for NNPDF fits
"""

import logging
import numbers
import re
from collections.abc import Mapping

import numpy as np

from NNPDF import CommonData, RandomGenerator
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


@make_argcheck
def check_combocuts(combocuts: str):
    """Check combocuts content"""
    check(combocuts == 'NNPDF31',
          "Invalid combocut. Must be NNPDF31 (or implement it yourself).")


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
        log.warning(f"Dataset output folder exists: {path} Overwritting contents")
    else:
        path.mkdir(exist_ok=True)


def export_mask(path, mask):
    """Dump mask to file"""
    np.savetxt(path, mask, fmt='%d')


@check_combocuts
@check_rngalgo
@check_positive('errorsize')
@check_nonnegative('filterseed')
@check_nonnegative('seed')
def filter(experiments, theoryid, filter_path,
           fakedata: bool,
           filterseed:int, rngalgo:int, seed:int, fakenoise:bool,
           errorsize:numbers.Real, combocuts, t0pdfset,
           rules, defaults):
    """Apply filters to all datasets"""
    if not fakedata:
        log.info('Filtering real data.')
        total_data, total_cut_data = _filter_real_data(filter_path, experiments)
    else:
        log.info('Filtering closure-test data.')
        RandomGenerator.InitRNG(rngalgo, seed)
        RandomGenerator.GetRNG().SetSeed(filterseed)
        total_data, total_cut_data = _filter_closure_data(filter_path, experiments,
                                                          t0pdfset, fakenoise, errorsize)
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


def _filter_real_data(filter_path, experiments):
    """Filter real experimental data."""
    total_data_points = 0
    total_cut_data_points = 0
    for experiment in experiments:
        for dataset in experiment.datasets:
            path = filter_path / dataset.name
            nfull, ncut = _write_ds_cut_data(path, dataset)
            total_data_points += nfull
            total_cut_data_points += ncut
            dataset.load().Export(str(path))
    return total_data_points, total_cut_data_points


def _filter_closure_data(filter_path, experiments, fakepdfset, fakenoise, errorsize):
    """Filter closure test data."""
    total_data_points = 0
    total_cut_data_points = 0
    fakeset = fakepdfset.load()
    # Load experiments
    for experiment in experiments:
        #Don't want to save this in any cache since we are mutating it
        loaded_exp = experiment.load.__wrapped__(experiment)
        loaded_exp.MakeClosure(fakeset, fakenoise)
        for j, dataset in enumerate(experiment.datasets):
            path = filter_path / dataset.name
            nfull, ncut = _write_ds_cut_data(path, dataset)
            total_data_points += nfull
            total_cut_data_points += ncut
            loaded_ds = loaded_exp.GetSet(j)
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
    vfns: int
        Variable Flavour Number Scheme defined by the theory
    ic: int
        Defined by the theory
    pto: int
        Perturbative order defined by the theory
    loader: validphys.loader.Loader, optional
        A loader instance used to retrieve the datasets.
    """

    numpy_functions = {"sqrt": np.sqrt, "log": np.log, "fabs": np.fabs}

    def __init__(
        self,
        *,
        initial_data: dict,
        defaults: dict,
        theory_parameters: tuple,
        loader=None,
    ):
        self.dataset = None
        self.process_type = None
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
        try:
            self.rule = compile(self.rule, "rule", "eval")
        except Exception as e:
            raise RuleProcessingError(
                f"Could not process rule {self.rule_string!r}: {e}"
            ) from e
        ns = {
            *self.numpy_functions,
            *self.defaults,
            *self.local_variables,
            *self.variables,
            "idat",
            "central_value",
        }
        for name in self.rule.co_names:
            if not name in ns:
                raise RuleProcessingError(
                    f"Could not process rule {self.rule_string!r}: Unknown name {name!r}"
                )


    def __call__(self, dataset, idat):
        central_value = dataset.GetData(idat)
        # We return None if the rule doesn't apply. This
        # is different to the case where the rule does apply,
        # but the point was cut out by the rule.
        if (dataset.GetSetName() != self.dataset and
            dataset.GetProc(idat) != self.process_type and
            self.process_type != "DIS_ALL"):
            return None

        # Handle the generalised DIS cut
        if self.process_type == "DIS_ALL" and dataset.GetProc(idat)[:3] != "DIS":
            return None

        self._set_kinematics_dict(dataset, idat)
        self._set_local_variables_dict(dataset, idat)

        for parameter in self.theory_params:
            if parameter[0] == "PTO" and hasattr(self, "PTO"):
                if parameter[1] not in self.PTO:
                    return None
            elif hasattr(self, parameter[0]) and (getattr(self, parameter[0]) != parameter[1]):
                return None

        # Will return True if datapoint passes through the filter
        try:
            return eval(
                self.rule,
                self.numpy_functions,
                {
                    **{"idat": idat, "central_value": central_value},
                    **self.defaults,
                    **self.kinematics_dict,
                    **self.local_variables_dict,
                },
            )
        except Exception as e:
            raise FatalRuleError(f"Error when applyin rule {self.rule_string!r}: {e}") from e

    def __repr__(self):
        return self.rule_string

    def _set_kinematics_dict(self, dataset, idat) -> dict:
        kinematics = [dataset.GetKinematics(idat, j) for j in range(3)]
        self.kinematics_dict = dict(zip(self.variables, kinematics))

    def _set_local_variables_dict(self, dataset, idat) -> dict:
        local_variables = {}
        if not hasattr(self, "kinematics_dict"):
            self._set_kinematics_dict(dataset, idat)

        if hasattr(self, "local_variables"):
            for key, value in self.local_variables.items():
                local_variables[key] = eval(
                    str(value),
                    {**self.numpy_functions, **self.kinematics_dict, **local_variables},
                )
        self.local_variables_dict = local_variables

def get_cuts_for_dataset(commondata, rules, defaults) -> list:
    """Function to generate a list containing the index
    of all experimental points that passed kinematic
    cut rules stored in ./cuts/filters.yaml


    Parameters
    ----------
    commondata: NNPDF CommonData spec
    theoryid: NNPDF theoryID object

    Returns
    -------
    mask: list
        List object containing index of all passed experimental
        values

    Example
    -------
    >>> from importlib.resources import read_text
    >>> from reportengine.compat import yaml
    >>> from validphys.filters import get_cuts_for_dataset, Rule
    >>> from validphys.loader import Loader
    >>> import validphys.cuts
    >>> l = Loader()
    >>> cd = l.check_commondata("NMC")
    >>> filter_defaults = yaml.safe_load(read_text(validphys.cuts, "defaults.yaml")) 
    >>> filter_rules = yaml.safe_load(read_text(validphys.cuts, "filters.yaml"))
    >>> theory = l.check_theoryID(53)
    >>> params = tuple(theory.get_description().items())
    >>> rule_list = [Rule(initial_data=i, defaults=filter_defaults, theory_parameters=params) for i in filter_rules]
    >>> get_cuts_for_dataset(cd, rules=rule_list, defaults=filter_defaults)
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
