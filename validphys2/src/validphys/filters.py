# -*- coding: utf-8 -*-
"""
Filters for NNPDF fits
"""

import logging
import numbers
import numpy as np
import pathlib

from NNPDF import DataSet, RandomGenerator, CommonData
from reportengine.checks import make_argcheck, check, check_positive, make_check
from reportengine.compat import yaml

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

class Rule:
    """Rule object to be used to generate cuts mask.

    A rule object is created for each rule in ./cuts/filters.yaml
    """

    def __init__(self, *, initial_data: dict, defaults: dict, vfns, ic, pto):
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

        self.defaults = defaults
        self.pto = pto
        self.vfns = vfns
        self.ic = ic

    def __call__(self, dataset, idat):
        pto = self.pto
        vfns = self.vfns
        ic = self.ic
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
                globals(),
                {
                    **locals(),
                    **self.defaults,
                    **self.kinematics_dict,
                    **self.local_variables_dict,
                },
            )

        if (dataset.GetSetName() != self.dataset and
            dataset.GetProc(idat) != self.process_type):
            return

        if hasattr(self, "VFNS") and self.VFNS != vfns:
            return

        # Will return True if datapoint passes through the filter
        return eval(
            self.rule,
            globals(),
            {
                **locals(),
                **self.defaults,
                **self.kinematics_dict,
                **self.local_variables_dict,
            },
        )


    def _set_kinematics_dict(self, dataset, idat) -> dict:
        kinematics = [dataset.GetKinematics(idat, j) for j in range(3)]
        self.kinematics_dict = dict(zip(self.variables, kinematics))

    def _set_local_variables_dict(self, dataset, idat) -> dict:
        local_variables = {}
        if not hasattr(self, "kinematics_dict"):
            self._set_kinematics_dict(dataset, idat)

        if hasattr(self, "local_variables"):
            for key, value in self.local_variables.items():
                exec(key + "=" + str(value), {**globals(), **self.kinematics_dict}, local_variables)
        self.local_variables_dict = local_variables

path = pathlib.Path(__file__).resolve().parent
with open(path/"cuts/filters.yaml", "r") as rules_stream,\
     open(path/"cuts/defaults.yaml", "r") as defaults_stream:
    try:
        rules = yaml.safe_load(rules_stream)
        defaults = yaml.safe_load(defaults_stream)
    except yaml.YAMLError as exception:
        print(exception)

def get_cuts_for_dataset(commondata, theoryid) -> list:
    """Function to generate a list containing the index
    of all experimental points that passed kinematic
    cut rules stored in ./cuts/filters.yaml


    Parameters
    ----------
    commondata: NNPDF commondata object
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
    >>> ds = l.check_dataset("NMC", theoryid=53, cuts="nocuts")
    >>> theory = l.check_theoryID(53)
    >>> get_cuts_for_dataset(ds, theory)
    """
    dataset = commondata.load()

    theoryid_dict = theoryid.get_description()
    pto, ic, vfns = map(theoryid_dict.get, ["PTO", "IC", "VFNS"])

    rule_list = [Rule(initial_data=i, defaults=defaults,
                      pto=pto, vfns=vfns, ic=ic) for i in rules]

    mask = []
    for idat in range(dataset.GetNData()):
        broken = False
        for rule in rule_list:
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
