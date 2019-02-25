# -*- coding: utf-8 -*-
"""
Filters for NNPDF fits
"""

import logging
import numbers
import numpy as np

from NNPDF import DataSet, RandomGenerator
from reportengine.checks import make_argcheck, check, check_positive, make_check

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


def get_cuts_for_dataset(commondata, theoryid, q2min, w2min):
    """Return cut mask for dataset"""
    datamask = []
    ds = commondata.load()
    for idat in range(ds.GetNData()):
        if pass_kincuts(ds, idat, theoryid, q2min, w2min):
            datamask.append(idat)
    return datamask


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


def pass_kincuts(dataset, idat, theoryid, q2min, w2min):
    """Applies cuts as in C++ for NNPDF3.1 combo cuts.
    This function replicas the c++ code but should be upgraded as
    discussed several times.
    """
    pto = theoryid.get_description().get('PTO')
    vfns = theoryid.get_description().get('FNS')
    ic = theoryid.get_description().get('IC')

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
