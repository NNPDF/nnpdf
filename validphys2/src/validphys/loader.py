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
import re

import numpy as np
import yaml

from NNPDF import CommonData

from validphys.core import (CommonDataSpec, FitSpec, TheoryIDSpec, FKTableSpec,
                            PositivitySetSpec, DataSetSpec)

log = logging.getLogger(__name__)

class LoadFailedError(FileNotFoundError): pass

class DataNotFoundError(LoadFailedError): pass

class SysNotFoundError(LoadFailedError): pass

class FKTableNotFound(LoadFailedError): pass

class CfactorNotFound(LoadFailedError): pass

class CompoundNotFound(LoadFailedError): pass

class TheoryNotFound(LoadFailedError): pass

class FitNotFound(LoadFailedError): pass

class CutsNotFound(LoadFailedError): pass


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

    def check_commondata(self, setname, sysnum=None):
        datafile = self.commondata_folder / ('DATA_' + setname + '.dat')
        if not datafile.exists():
            raise DataNotFoundError(("Could not find Commondata set: '%s'. "
                  "File '%s' does not exist.")
                 % (setname, datafile))
        if sysnum is None:
            sysnum = 'DEFAULT'
        sysfile = (self.commondata_folder / 'systypes' /
                   ('SYSTYPE_%s_%s.dat' % (setname, sysnum)))

        #Good luck debugging this!
        plotfiles = tuple(p for p in self.commondata_folder.iterdir()
                        if re.match('PLOTTING_%s(_.*)?\.ya?ml'%setname, p.name))

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
            raise TheoryNotFound(("Could not find theory %s. "
                  "Folder '%s' not found") % (theoryID, theopath) )
        return TheoryIDSpec(theoryID, theopath)

    def get_commondata(self, setname, sysnum, plotfiles=None):
        """Get a Commondata from the set name and number.
           The plotfiles argument is accepted to keep symmetry with
           the commondataSpec,
           returned by check_commondata, but it doesn't do anything."""
        datafile, sysfile, *_ = self.check_commondata(setname, sysnum)
        return CommonData.ReadFile(str(datafile), str(sysfile))

    #   @functools.lru_cache()
    def check_fktable(self, theoryID, setname, cfac):
        _, theopath = self.check_theoryID(theoryID)
        fkpath = theopath/ 'fastkernel' / ('FK_%s.dat' % setname)
        if not fkpath.exists():
          raise FKTableNotFound(("Could not find FKTable for set '%s'. "
          "File '%s' not found") % (setname, fkpath) )

        cfactors = self.check_cfactor(theoryID, setname, cfac)
        return FKTableSpec(fkpath, cfactors)

    def check_compound(self, theoryID, setname, cfac):
        thid, theopath = self.check_theoryID(theoryID)
        compound_spec_path = theopath / 'compound' / ('FK_%s-COMPOUND.dat' % setname)
        try:
            with compound_spec_path.open() as f:
                #Drop first line with comment
                next(f)
                txt = f.read()
        except FileNotFoundError as e:
            msg = ("Could not find COMPOUND set '%s' for theory %d: %s" %
                   (setname, int(thid), e))
            raise CompoundNotFound(msg)
        #This is a little bit stupid, but is the least amount of thinking...
        yaml_format = 'FK:\n' + re.sub('FK:', ' - ', txt)
        data = yaml.load(yaml_format)
        #we have to split out 'FK_' the extension to get a name consistent
        #with everything else
        tables = [self.check_fktable(theoryID, name[3:-4], cfac)
                  for name in data['FK']]
        op = data['OP']
        return tuple(tables), op


    def get_fktable(self, theoryID, setname, cfac):

        fkspec= self.check_fktable(theoryID, setname, cfac)
        return fkspec.load()

    def check_cfactor(self, theoryID, setname, cfactors):
        _, theopath = self.check_theoryID(theoryID)
        cf = []
        for cfactor in cfactors:
            cfactorpath = (theopath / 'cfactor' /
                           'CF_{cfactor}_{setname}.dat'.format(**locals()))
            if not cfactorpath.exists():
                msg = ("Could not find cfactor '{cfactor}' for FKTable {setname} "
                       "in theory {theoryID}. File {cfactorpath} does not "
                       "exist.").format(**locals())
                raise CfactorNotFound(msg)
            cf.append(cfactorpath)

        return tuple(cf)

    def check_posset(self, theiryID, setname, postlambda):
        cd = self.check_commondata(setname, 0)
        fk = self.check_fktable(theiryID, setname, [])
        th =  self.check_theoryID(theiryID)
        return PositivitySetSpec(cd, fk, postlambda, th)

    def get_posset(self, theiryID, setname, postlambda):
        return self.check_posset(theiryID, setname, postlambda).load()

    def check_fit(self, fitname):
        resultspath = self.resultspath
        p = resultspath / fitname
        if p.is_dir():
            return FitSpec(fitname, p)
        if not p.exists():
            msg = ("Could not find fit '{fitname}' in '{resultspath}'. "
                   "Folder '{p}' not found").format(**locals())
            raise FitNotFound(msg)
        msg = ("Could not load fit '{fitname}' from '{resultspath}. "
                   "'{p}' must be a folder").format(**locals())
        raise FitNotFound(msg)

    def check_dataset(self, *, name, sysnum=None,
                     theoryid, cfac=(),
                     use_cuts, fit=None):

        if not isinstance(theoryid, TheoryIDSpec):
            theoryid = self.check_theoryID(theoryid)

        theoryno, _ = theoryid

        op = None

        commondata = self.check_commondata(name, sysnum)

        try:
            fkspec, op = self.check_compound(theoryno, name, cfac)
        except CompoundNotFound:
            fkspec = self.check_fktable(theoryno, name, cfac)

        if use_cuts:
            cuts = self.get_cuts(name, fit)
        else:
            cuts = None

        return DataSetSpec(name=name, commondata=commondata,
                           fkspecs=fkspec, thspec=theoryid, cuts=cuts, op=op)

    @property
    def available_fits(self):
        try:
            return [p.name for p in self.resultspath.iterdir() if p.is_dir()]
        except OSError:
            return []

    def get_cuts(self, setname, fit):
        fitname, fitpath = fit
        p = (fitpath/'filter')/setname/('FKMASK_' + setname+ '.dat')
        if not p.parent.exists():
            raise CutsNotFound("Bad filter configuration. "
            "Could not find: %s" % p.parent)
        if not p.exists():
            return None
        cuts = np.loadtxt(str(p), dtype=int)
        log.debug("Loading cuts for %s" % setname)
        return cuts
