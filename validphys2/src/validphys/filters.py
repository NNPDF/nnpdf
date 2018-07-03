# -*- coding: utf-8 -*-
"""
Filters for NNPDF fits
"""

import logging
import numbers
import numpy as np
from NNPDF import DataSet, RandomGenerator
from reportengine.checks import make_argcheck, check, check_positive, make_check, CheckError

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


@make_argcheck
def check_rancutmethod(rancutmethod: int):
    """Check rancutmethod content
    If rancutmethod is 0: no random cut is applied to CT data.
    If rancutmethod is 1: cuts are pure random.
    If rancutmethod is 2: evently spread points (non-random).
    If rancutmethod is 3: random w/ exact 50:50 split.
    """
    check(0 <= rancutmethod <= 3,
          "Invalid rancutmethod. Must be int between [0, 3].")


def check_positive_int(var):
    """Ensure that `var` is positive"""
    @make_check
    def check(ns, **kwargs):
        val = ns[var]
        if not val>=0:
            raise CheckError(f"'{var}' must be positive or equal zero, but it is {val!r}.")
    return check


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
@check_rancutmethod
@check_positive('q2min')
@check_positive('w2min')
@check_positive('errorsize')
@check_positive_int('filterseed')
@check_positive_int('seed')
def filter(experiments, theoryid, filter_path,
           q2min:numbers.Real, w2min:numbers.Real, fakepdfset,
           fakedata, filterseed:int, rngalgo:int, seed:int,
           fakenoise:bool, rancutmethod:int, rancutprob:numbers.Real,
           errorsize:numbers.Real, rancuttrnval:bool, combocuts, t0set=None):
    """Apply filters to all datasets"""
    if not fakedata:
        log.info('Filtering real data.')
        total_data, total_cut_data = _filter_real_data(filter_path, experiments, theoryid, q2min, w2min)
    else:
        log.info('Filtering closure-test data.')
        RandomGenerator.InitRNG(rngalgo, seed)
        RandomGenerator.GetRNG().SetSeed(filterseed)
        total_data, total_cut_data = _filter_closure_data(filter_path, experiments, theoryid, q2min, w2min,
                                                          fakepdfset, fakenoise, rancutmethod, rancutprob,
                                                          rancuttrnval, errorsize)
    log.info(f'Summary: {total_cut_data}/{total_data} datapoints passed kinematic cuts.')


def _filter_real_data(filter_path, experiments, theoryid, q2min, w2min):
    """Filter real experimental data."""
    total_data_points = 0
    total_cut_data_points = 0
    for experiment in experiments:
        for dataset in experiment.datasets:
            path = filter_path / dataset.name
            make_dataset_dir(path)
            ds = dataset.load()
            total_data_points += ds.GetNData()
            # build data mask
            datamask = []
            for idat in range(ds.GetNData()):
                if pass_kincuts(ds, idat, theoryid, q2min, w2min):
                    datamask.append(idat)
            log.info(f'{len(datamask)}/{ds.GetNData()} datapoints in {dataset.name} passed kinematic cuts.')
            total_cut_data_points += len(datamask)
            # save to disk
            if len(datamask) != ds.GetNData():
                export_mask(path / ('FKMASK_%s.dat' % dataset.name), datamask)
                ds = DataSet(ds, datamask)
            ds.Export(str(path))
    return total_data_points, total_cut_data_points


def _filter_closure_data(filter_path, experiments, theoryid, q2min, w2min,
                         fakepdfset, fakenoise, rancutmethod, rancutprob,
                         rancuttrnval, errorsize):
    """Filter closure test data."""
    total_data_points = 0
    total_cut_data_points = 0
    fakeset = fakepdfset.load()
    # Load experiments
    for experiment in experiments:
        uncut_exp = experiment.load()
        uncut_exp.MakeClosure(fakeset, fakenoise)
        for j, dataset in enumerate(experiment.datasets):
            path = filter_path / dataset.name
            make_dataset_dir(path)
            ds = uncut_exp.GetSet(j)
            total_data_points += ds.GetNData()
            # build data mask
            datamask = []
            for idat in range(ds.GetNData()):
                if pass_kincuts(ds, idat, theoryid, q2min, w2min):
                    datamask.append(idat)
            if 1 <= rancutmethod <= 3:
                random_cut(datamask, rancutmethod, rancutprob, rancuttrnval)
            log.info(f'{len(datamask)}/{ds.GetNData()} datapoints in {dataset.name} passed kinematic cuts.')
            total_cut_data_points += len(datamask)
            # save to disk
            if len(datamask) != ds.GetNData():
                export_mask(path / ('FKMASK_%s.dat' % dataset.name), datamask)
                ds = DataSet(ds, datamask)
            if errorsize != 1.0:
                ds.RescaleErrors(errorsize)
            ds.Export(str(path))
    return total_data_points, total_cut_data_points


def check_positivity(posdatasets):
    """Verify positive datasets are ready for the fit."""
    log.info('Verifying positivity tables:')
    for pos in posdatasets:
        pos.load()
        log.info(f'{pos.name} checked.')


def pass_kincuts(dataset, idat, theoryid, q2min, w2min):
    """Applies cuts as in C++ for NNPDF3.1 combo cuts.
    This function replicas the c++ code but should be upgraded as
    discussed several times.
    """
    pto = theoryid.get_description().get('PTO')
    vfns = theoryid.get_description().get('FNS')
    ic = theoryid.get_description().get('IC')

    check(0 <= pto <= 2, "Invalid PTO. Must be 0, 1 or 2.")
    check(vfns in ('FONLL-A', 'FONLL-B', 'FONLL-C'), "Invalid FNS. Must be FONLL-A/B/C.")
    check(ic in (True, False), "Invalid IC. Must be True or False.")

    if dataset.GetSetName() == 'ATLAS1JET11':
        # allow only first rapidity bin of ATLAS1JET11
        return dataset.GetKinematics(idat, 0) < 0.3

    if dataset.GetSetName() in ('LHCBWZMU8TEV', 'LHCBWZMU7TEV'):
        if pto == 2:
            return dataset.GetKinematics(idat, 0) >= 2.25

    if dataset.GetSetName() in ('D0WMASY', 'D0WEASY'):
        if pto == 2:
            return dataset.GetData(idat) >= 0.03

    if dataset.GetSetName() == 'ATLASZPT7TEV':
        pt = np.sqrt(dataset.GetKinematics(idat, 1))
        if pt < 30 or pt > 500:
            return False
        return True

    if dataset.GetSetName() == 'ATLASZPT8TEVMDIST':
        return dataset.GetKinematics(idat, 0) >= 30

    if dataset.GetSetName() == 'ATLASZPT8TEVYDIST':
        pt = np.sqrt(dataset.GetKinematics(idat, 1))
        if pt < 30 or pt > 150:
            return False
        return True

    if dataset.GetSetName() == 'CMSZDIFF12':
        pt = np.sqrt(dataset.GetKinematics(idat, 1))
        y = dataset.GetKinematics(idat, 0)
        if pt < 30 or pt > 170 or y > 1.6:
            return False
        return True

    if dataset.GetSetName() == 'ATLASWPT31PB':
        return dataset.GetKinematics(idat, 0) > 30

    if dataset.GetProc(idat)[0:3] in ('EWK', 'DYP'):
        # building rapidity and pT or Mll
        y = dataset.GetKinematics(idat, 0)
        pTmv = np.sqrt(dataset.GetKinematics(idat, 1))

        # generalized cuts
        maxCMSDY2Dy = 2.2
        maxCMSDY2Dminv = 200.0
        minCMSDY2Dminv = 30.0
        maxTau = 0.080
        maxY = 0.663

        if dataset.GetSetName() in ('CMSDY2D11', 'CMSDY2D12'):
            if pto == 0 or pto == 1:
                if pTmv > maxCMSDY2Dminv or pTmv < minCMSDY2Dminv or y > maxCMSDY2Dy:
                    return False
            if pto == 2:
                if pTmv > maxCMSDY2Dminv or y > maxCMSDY2Dy:
                    return False
            return True

        if dataset.GetSetName() in ('ATLASZHIGHMASS49FB', 'LHCBLOWMASS37PB'):
            if pTmv > maxCMSDY2Dminv:
                return False
            return True

        if dataset.GetSetName() == 'ATLASLOMASSDY11':
            if pto == 0 or pto == 1:
                if idat < 6:
                    return False
            return True

        if dataset.GetSetName() == 'ATLASLOMASSDY11EXT':
            if pto == 0 or pto == 1:
                if idat < 2:
                    return False
            return True

        # new cuts for the fixed target DY
        if dataset.GetSetName() in ('DYE886P', 'DYE605'):
            rapidity = dataset.GetKinematics(idat, 0)
            invM2 = dataset.GetKinematics(idat, 1)
            sqrts = dataset.GetKinematics(idat, 2)
            tau = invM2 / sqrts**2
            ymax = -0.5 * np.log(tau)

            if tau > maxTau or np.fabs(rapidity / ymax) > maxY:
                return False
            return True

    # DIS cuts
    if dataset.GetProc(idat)[0:3] == 'DIS':
        # load kinematics
        x = dataset.GetKinematics(idat, 0)
        Q2 = dataset.GetKinematics(idat, 1)
        W2 = Q2 * (1 - x) / x

        # basic cuts
        if W2 <= w2min or Q2 <= q2min:
            return False

        if dataset.GetSetName() in ('EMCF2P', 'EMCF2D'):
            return x > 0.1

        # additional F2C cuts in case of FONLL-A
        if dataset.GetProc(idat) == 'DIS_NCP_CH' and vfns == 'FONLL-A':
            Q2cut1_f2c = 4
            Q2cut2_f2c = 10
            xcut_f2c = 1e-3

            if Q2 <= Q2cut1_f2c: # cut if Q2 <= 4
                return False

            if Q2 <= Q2cut2_f2c and x <= xcut_f2c: # cut if Q2 <= 10 and x <= 10 ^ -3
                return False

        # additional F2C cut in case of FONLL-C + IC
        if dataset.GetProc(idat) == 'DIS_NCP_CH' and vfns == 'FONLL-C' and ic:
            Q2cut1_f2c = 8
            if Q2 <= Q2cut1_f2c:
                return False
    return True


def random_cut(datamask, rancutmethod, rancutprob, rancuttrnval):
    """Apply random cuts for closure test data"""
    log.info(f'Applying random cuts to data using method {rancutmethod}')
    log.info(f'Cutting to {rancutprob}%')
    valdatamask = []
    if rancutmethod == 1: # Option 1: pure random
        i = 0
        while i < len(datamask):
            if RandomGenerator.GetRNG().GetRandomUniform() > datamask and len(datamask) > 2:
                valdatamask.append(datamask[i])
                datamask.pop(i)
                i -= 1
            i += 1
    elif rancutmethod == 2: # Option 2: evently spread points (non-random)
        counter = 0.0
        i = 0
        while i < len(datamask):
            counter += (1.0 - rancutprob)
            if counter >= 1 and len(datamask) > 2:
                valdatamask.append(datamask[i])
                datamask.pop(i)
                i -= 1
                counter -= 1
            i += 1
    elif rancutmethod == 3: # Option 3: random w/ exact 50:50 split
        Ndatremove = int(min(len(datamask) * (1.0 - rancutprob), len(datamask) - 2))
        for i in range(Ndatremove):
            position = RandomGenerator.GetRNG().GetRandomUniform(len(datamask))
            valdatamask.append(datamask[position])
            datamask.pop(position)
    if rancuttrnval:
        datamask = valdatamask

