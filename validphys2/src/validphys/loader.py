# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:40:38 2016

@author: Zahari Kassabov

Resolve paths to useful objects, and query the existence of different resources
within the specified paths.
"""
import pathlib
import functools
import logging

import numpy as np

from NNPDF import CommonData, FKTable

from validphys.core import CommonDataSpec, FitSpec, TheoryIDSpec

log = logging.getLogger(__name__)

class DataNotFoundError(FileNotFoundError): pass

class SysNotFoundError(FileNotFoundError): pass

class Loader():
    """Load various resources from the NNPDF data path."""

    def __init__(self, datapath, resultspath):
        datapath = pathlib.Path(datapath)
        resultspath = pathlib.Path(resultspath)
        self.datapath = datapath
        self.resultspath = resultspath

    @property
    @functools.lru_cache()
    def available_theories(self):
        """Return a string token for each of the available theories"""
        theory_token  = 'theory_'
        return {folder.name[len(theory_token):]
                for folder in self.datapath.glob(theory_token+'*')}

    @property
    @functools.lru_cache()
    def available_datasets(self):

        data_str = "DATA_"
        return {file.stem[len(data_str):] for
                file in self.commondata_folder.glob(data_str+"*.dat")}

    @property
    def commondata_folder(self):
        return self.datapath / 'commondata'

    def check_commondata(self, setname, sysnum):
        datafile = self.commondata_folder / ('DATA_' + setname + '.dat')
        if not datafile.exists():
            raise DataNotFoundError(("Could not find Commondata set: '%s'. "
                  "File '%s' does not exist.")
                 % (setname, datafile))
        sysfile = (self.commondata_folder / 'systypes' /
                   ('SYSTYPE_%s_%d.dat' % (setname, sysnum)))

        plotfiles = list(self.commondata_folder.glob('PLOTTING_' + setname +
                                                         '*' + '.y[a]ml'))

        if not sysfile.exists():
            raise SysNotFoundError(("Could not find systype %d for "
                 "dataset '%s'. File %s does not exist.") % (sysnum, setname,
                  sysfile))

        return CommonDataSpec(datafile, sysfile, plotfiles)

    @functools.lru_cache()
    def check_theoryID(self, theoryID):
        theoryID = str(theoryID)
        theopath = self.datapath / ('theory_%s' % theoryID)
        if not theopath.exists():
            raise FileNotFoundError(("Could not find theory %s. "
                  "Folder '%s' not found") % (theoryID, theopath) )
        return TheoryIDSpec(theoryID, theopath)

    def get_commondata(self, setname, sysnum, plotfiles=None):
        """Get a Commondata from the set name and number.
           The plotfiles argument is accepted to keep symmetry with
           the commondataSpec,
           returned by check_commondata, but it doesn't do anything."""
        datafile, sysfile, *_ = self.check_commondata(setname, sysnum)
        return CommonData.ReadFile(str(datafile), str(sysfile))

    @functools.lru_cache()
    def check_fktable(self, theoryID, setname):
        _, theopath = self.check_theoryID(theoryID)
        fkpath = theopath/ 'fastkernel' / ('FK_%s.dat' % setname)
        if not fkpath.exists():
          raise FileNotFoundError(("Could not find FKTable for set '%s'. "
          "File '%s' not found") % (setname, fkpath) )
        return fkpath

    def get_fktable(self, theoryID, setname):

        fkpath = self.check_fktable(theoryID, setname)
        return FKTable(str(fkpath), [])

    def check_cfactor(self, theoryID, setname, cfactors):
        theopath = self.check_theoryID(theoryID)
        cf = []
        for cfactor in cfactors:
            cfactorpath = (theopath / 'cfactor' /
                           'CF_{cfactor}_{setname}.dat'.format(**locals()))
            if not cfactorpath.exists():
                msg = ("Could not find cfactor '{cfactor}' for set {setname} "
                       "in theory {theoryID}. File {cfactorpath} does not "
                       "exist.").format(**locals())
                raise FileNotFoundError(msg)
            cf.append(cfactorpath)

        return cf

    def check_fit(self, fitname):
        resultspath = self.resultspath
        p = resultspath / fitname
        if p.is_dir():
            return FitSpec(fitname, p)
        if not p.exists():
            msg = ("Could not find fit '{fitname}' in '{resultspath}'. "
                   "Folder '{p}' not found").format(**locals())
            raise FileNotFoundError(msg)
        msg = ("Could not load fit '{fitname}' from '{resultspath}. "
                   "'{p}' must be a folder").format(**locals())
        raise IOError(msg)

    @property
    def available_fits(self):
        return [p.name for p in self.resultspath.iterdir() if p.is_dir()]

    def get_cuts(self, setname, fit):
        fitname, fitpath = fit
        p = (fitpath/'filter')/setname/('FKMASK_' + setname+ '.dat')
        if not p.parent.exists():
            raise FileNotFoundError("Bad filter configuration. "
            "Could not find: %s" % p.parent)
        if not p.exists():
            return None
        cuts = np.loadtxt(str(p), dtype=int)
        log.debug("Loading cuts for %s" % setname)
        return cuts
