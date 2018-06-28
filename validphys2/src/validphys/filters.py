# -*- coding: utf-8 -*-
"""
Filters for NNPDF fits
"""

import logging
import numpy as np
from reportengine.checks import make_argcheck, check
from NNPDF import DataSet

log = logging.getLogger(__name__)


@make_argcheck
def check_combocuts(combocuts: str):
    """Check combocuts content"""
    check(combocuts in 'NNPDF31',
          "Invalid combocut. Must be NNPDF31 (or implement it yourself).")


def make_dataset_dir(path):
    """check if results folder exists"""
    if path.exists():
        log.warning(f"Dataset output folder exists: {path} Overwritting contents")
    else:
        path.mkdir(exist_ok=True)


def export_mask(path, mask):
    """Dump mask to file"""
    np.savetxt(path, mask, fmt='%d')


@check_combocuts
def filter(experiments, theoryid, filter_path, q2min, w2min, combocuts, t0set=None):
    """Apply filters to all datasets"""
    pto = theoryid.get_description().get('PTO')
    vfns = theoryid.get_description().get('FNS')
    ic = theoryid.get_description().get('IC')

    total_data_points = 0
    total_cut_data_points = 0
    # Load experiments
    for experiment in experiments:
        for dataset in experiment.datasets:
            path = filter_path / dataset.name
            make_dataset_dir(path)
            ds = dataset.load()
            total_data_points += ds.GetNData()

            # build data mask
            datamask = []
            for idat in range(ds.GetNData()):
                if pass_kincuts(ds, idat, pto, q2min, w2min, vfns, ic):
                    datamask.append(idat)

            log.info('%d/%d datapoints in %s passed kinematic cuts.' %
                    (len(datamask), ds.GetNData(), dataset.name))
            total_cut_data_points += len(datamask)

            # save to disk
            if len(datamask) != ds.GetNData():
                export_mask(path / ('FKMASK_%s.dat' % dataset.name), datamask)
                ## in principle one should do:
                # dataset.cuts = datamask
                # ds = dataset.load()
                ## however the lru cache ignores the self.cuts change
                ## thus this is the quick solution:
                ds = DataSet(ds, datamask)

            ds.Export(str(path))

    log.info('TOTAL: %d/%d datapoints passed kinematic cuts.' % (total_cut_data_points, total_data_points))


def pass_kincuts(set, idat, pto, q2min, w2min, vfns, ic):
    """Applies cuts as in C++ for NNPDF3.1 combo cuts.
    This function replicas the c++ code but should be upgraded as
    discussed several times.
    """
    if set.GetSetName() in 'ATLAS1JET11':
        # allow only first rapidity bin of ATLAS1JET11
        return set.GetKinematics(idat, 0) < 0.3

    if set.GetSetName() in ('LHCBWZMU8TEV', 'LHCBWZMU7TEV'):
        if pto == 2:
            return set.GetKinematics(idat, 0) >= 2.25

    if set.GetSetName() in ('D0WMASY', 'D0WEASY'):
        if pto == 2:
            return set.GetData(idat) >= 0.03

    if set.GetSetName() in 'ATLASZPT7TEV':
        pt = np.sqrt(set.GetKinematics(idat, 1))
        if pt < 30 or pt > 500:
            return False
        return True

    if set.GetSetName() in 'ATLASZPT8TEVMDIST':
        return set.GetKinematics(idat, 0) >= 30

    if set.GetSetName() in 'ATLASZPT8TEVYDIST':
        pt = np.sqrt(set.GetKinematics(idat, 1))
        if pt < 30 or pt > 150:
            return False
        return True

    if set.GetSetName() in 'CMSZDIFF12':
        pt = np.sqrt(set.GetKinematics(idat, 1))
        y = set.GetKinematics(idat, 0)
        if pt < 30 or pt > 170 or y > 1.6:
            return False
        return True

    if set.GetSetName() in 'ATLASWPT31PB':
        return set.GetKinematics(idat, 0) > 30

    if set.GetProc(idat)[0:3] in ('EWK', 'DYP'):
        # building rapidity and pT or Mll
        y = set.GetKinematics(idat, 0)
        pTmv = np.sqrt(set.GetKinematics(idat, 1))

        # generalized cuts
        maxCMSDY2Dy = 2.2
        maxCMSDY2Dminv = 200.0
        minCMSDY2Dminv = 30.0
        maxTau = 0.080
        maxY = 0.663

        if set.GetSetName() in ('CMSDY2D11', 'CMSDY2D12'):
            if pto == 0 or pto == 1:
                if pTmv > maxCMSDY2Dminv or pTmv < minCMSDY2Dminv or y > maxCMSDY2Dy:
                    return False
            if pto == 2:
                if pTmv > maxCMSDY2Dminv or y > maxCMSDY2Dy:
                    return False
            return True

        if set.GetSetName() in ('ATLASZHIGHMASS49FB', 'LHCBLOWMASS37PB'):
            if pTmv > maxCMSDY2Dminv:
                return False
            return True

        if set.GetSetName() in 'ATLASLOMASSDY11':
            if pto == 0 or pto == 1:
                if idat < 6:
                    return False
            return True

        if set.GetSetName() in 'ATLASLOMASSDY11EXT':
            if pto == 0 or pto == 1:
                if idat < 2:
                    return False
            return True

        # new cuts for the fixed target DY
        if set.GetSetName() in ('DYE886P', 'DYE605'):
            rapidity = set.GetKinematics(idat, 0)
            invM2 = set.GetKinematics(idat, 1)
            sqrts = set.GetKinematics(idat, 2)
            tau = invM2 / sqrts**2
            ymax = -0.5 * np.log(tau)

            if tau > maxTau or np.fabs(rapidity / ymax) > maxY:
                return False
            return True

    # DIS cuts
    if set.GetProc(idat)[0:3] in 'DIS':
        # load kinematics
        x = set.GetKinematics(idat, 0)
        Q2 = set.GetKinematics(idat, 1)
        W2 = Q2 * (1 - x) / x

        # basic cuts
        if W2 <= w2min or Q2 <= q2min:
            return False

        if set.GetSetName() in ('EMCF2P', 'EMCF2D'):
            return x > 0.1

        # additional F2C cuts in case of FONLL-A
        if set.GetProc(idat) in 'DIS_NCP_CH' and vfns in 'FONLL-A':
            Q2cut1_f2c = 4
            Q2cut2_f2c = 10
            xcut_f2c = 1e-3

            if Q2 <= Q2cut1_f2c: # cut if Q2 <= 4
                return False

            if Q2 <= Q2cut2_f2c and x <= xcut_f2c: # cut if Q2 <= 10 and x <= 10 ^ -3
                return False

        # additional F2C cut in case of FONLL-C + IC
        if set.GetProc(idat) in 'DIS_NCP_CH' and vfns in 'FONLL-C' and ic:
            Q2cut1_f2c = 8
            if Q2 <= Q2cut1_f2c:
                return False

    return True