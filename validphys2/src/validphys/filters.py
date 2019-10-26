# -*- coding: utf-8 -*-
"""
Filters for NNPDF fits
"""

import functools
import logging
import numbers
import numpy as np
import re

from importlib.resources import read_binary

from NNPDF import DataSet, RandomGenerator, CommonData
from reportengine.checks import make_argcheck, check, check_positive, make_check
from reportengine.compat import yaml
import validphys.cuts

log = logging.getLogger(__name__)


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
           errorsize:numbers.Real, combocuts, t0pdfset):
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
                r'(N(?P<nnumber>\d+)|(?P<multiplens>N*))LO(?P<operator>[\+\-\!])?'
                ).fullmatch(self.string)
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
    """

    def __init__(self, *, initial_data: dict, defaults: dict, theory_parameters: tuple):
        self.dataset = None
        self.process_type = None
        for key in initial_data:
            setattr(self, key, initial_data[key])

        if not hasattr(self, "rule"):
            raise AttributeError("No rule defined.")

        if self.dataset is None and self.process_type is None:
            raise AttributeError("Please define either a process type or dataset.")

        if self.process_type is None:
            from validphys.loader import Loader

            l = Loader()
            cd = l.check_commondata(self.dataset)
            if cd.process_type[:3] == "DIS":
                self.variables = CommonData.kinLabel["DIS"]
            else:
                self.variables = CommonData.kinLabel[cd.process_type]
        else:
            if self.process_type[:3] == "DIS":
                self.variables = CommonData.kinLabel["DIS"]
            else:
                self.variables = CommonData.kinLabel[self.process_type]

        if hasattr(self, "PTO"):
            self.PTO = PerturbativeOrder(self.PTO)

        self.rule_string = self.rule
        self.rule = compile(self.rule, "rule", "eval")
        self.defaults = defaults
        self.theory_params = theory_parameters

    def __call__(self, dataset, idat):
        central_value = dataset.GetData(idat)
        self._set_kinematics_dict(dataset, idat)
        # TODO: check why division by zero happens
        try:
            self._set_local_variables_dict(dataset, idat)
        except ZeroDivisionError:
            pass

        # Handle the generalised DIS cut
        if self.process_type == "DIS_ALL" and dataset.GetProc(idat)[:3] == "DIS":
            return eval(
                self.rule,
                numpy_functions,
                {
                    **locals(),
                    **self.defaults,
                    **self.kinematics_dict,
                    **self.local_variables_dict,
                },
            )

        # We return None if the rule doesn't apply. This
        # is different to the case where the rule does apply,
        # but the point was cut out by the rule.
        if (dataset.GetSetName() != self.dataset and
            dataset.GetProc(idat) != self.process_type):
            return None

        for parameter in self.theory_params:
            if parameter[0] == "PTO" and hasattr(self, "PTO"):
                if parameter[1] not in self.PTO:
                    return None
            elif hasattr(self, parameter[0]) and (getattr(self, parameter[0]) != parameter[1]):
                return None

        # Will return True if datapoint passes through the filter
        return eval(
            self.rule,
            numpy_functions,
            {
                **locals(),
                **self.defaults,
                **self.kinematics_dict,
                **self.local_variables_dict,
            },
        )

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
                exec(key + "=" + str(value), {**numpy_functions, **self.kinematics_dict}, local_variables)
        self.local_variables_dict = local_variables

rules = yaml.safe_load(read_binary(validphys.cuts, "filters.yaml"))
defaults = yaml.safe_load(read_binary(validphys.cuts, "defaults.yaml"))

numpy_functions = {"sqrt": np.sqrt, "log": np.log, "fabs": np.fabs, "np": np}

@functools.lru_cache()
def get_rule(index, theory_parameters):
    return Rule(
        initial_data=rules[index],
        defaults=defaults,
        theory_parameters=theory_parameters,
    )

@functools.lru_cache()
def get_theory_parameters(theoryid):
    return tuple(theoryid.get_description().items())

def get_cuts_for_dataset(commondata, theoryid) -> list:
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
    >>> from validphys.loader import Loader
    >>> l = Loader()
    >>> cd = l.check_commondata("NMC")
    >>> theory = l.check_theoryID(53)
    >>> get_cuts_for_dataset(cd, theory)
    """
    dataset = commondata.load()

    theoryid_params = get_theory_parameters(theoryid)

    mask = []
    for idat in range(dataset.GetNData()):
        broken = False
        for i in range(len(rules)):
            rule = get_rule(i, theoryid_params)
            try:
                rule_result = rule(dataset, idat)
                if rule_result is not None and rule_result == False:
                    broken = True
                    break
            except Exception as e:
                print(e)
                return rule

        if not broken:
            mask.append(idat)

    return mask
