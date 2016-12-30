# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:40:38 2016

@author: Zahari Kassabov

Resolve paths to useful objects, and query the existence of different resources
within the specified paths.
"""
import sys
import pathlib
import functools
import logging
import re
import tempfile
import shutil
import os.path as osp

import numpy as np
import yaml
import requests

from NNPDF import CommonData

from validphys.core import (CommonDataSpec, FitSpec, TheoryIDSpec, FKTableSpec,
                            PositivitySetSpec, DataSetSpec, PDF)
from validphys import lhaindex

log = logging.getLogger(__name__)

class LoaderError(Exception): pass

class LoadFailedError(FileNotFoundError, LoaderError): pass

class DataNotFoundError(LoadFailedError): pass

class SysNotFoundError(LoadFailedError): pass

class FKTableNotFound(LoadFailedError): pass

class CfactorNotFound(LoadFailedError): pass

class CompoundNotFound(LoadFailedError): pass

class TheoryNotFound(LoadFailedError): pass

class FitNotFound(LoadFailedError): pass

class CutsNotFound(LoadFailedError): pass

class PDFNotFound(LoadFailedError): pass

class RemoteLoaderError(LoaderError): pass

class LoaderBase:

    def __init__(self, datapath, resultspath):
        datapath = pathlib.Path(datapath)
        resultspath = pathlib.Path(resultspath)
        self.datapath = datapath
        self.resultspath = resultspath


#TODO: Deprecate get methods?
class Loader(LoaderBase):
    """Load various resources from the NNPDF data path."""

    @property
    def available_fits(self):
        try:
            return [p.name for p in self.resultspath.iterdir() if p.is_dir()]
        except OSError:
            return []

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
    @functools.lru_cache()
    def available_pdfs(self):
        return lhaindex.expand_local_names('*')

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
            raise SysNotFoundError(("Could not find systype %s for "
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

        commondata = self.check_commondata(name, sysnum)

        try:
            fkspec, op = self.check_compound(theoryno, name, cfac)
        except CompoundNotFound:
            fkspec = self.check_fktable(theoryno, name, cfac)
            op = None

        if use_cuts:
            cuts = self.get_cuts(name, fit)
        else:
            cuts = None

        return DataSetSpec(name=name, commondata=commondata,
                           fkspecs=fkspec, thspec=theoryid, cuts=cuts, op=op)

    def check_pdf(self, name):
        if lhaindex.isinstalled(name):
            return PDF(name)
        raise PDFNotFound(name)

    def get_pdf(self, name):
        return self.check_pdf(name).load()


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



#http://stackoverflow.com/a/15645088/1007990
def _download_file(url, stream):
    response = requests.get(url, stream=True)
    total_length = response.headers.get('content-length')

    if total_length is None or not log.isEnabledFor(logging.INFO):
        stream.write(response.content)
    else:
        dl = 0
        total_length = int(total_length)
        for data in response.iter_content(chunk_size=4096):
            dl += len(data)
            stream.write(data)
            done = int(50 * dl / total_length)
            sys.stdout.write("\r[%s%s] (%d%%)" % ('=' * done, ' ' * (50-done),
                             done*2))
            sys.stdout.flush()
        sys.stdout.write('\n')

def download_file(url, stream_or_path):
    #TODO: change to os.fspath
    if isinstance(stream_or_path, (str,bytes)):
        log.info("Downloading %s to %s." , url,stream_or_path)
        with open(stream_or_path, 'wb') as f:
            return _download_file(url, f)
    else:
        log.info("Downloading %s." , url,)
        return _download_file(url, stream_or_path)


def download_and_extract(url, local_path):
    name = url.split('/')[-1]
    with tempfile.NamedTemporaryFile(delete=False, suffix=name) as t:
        log.info("Saving data to %s" , t.name)
        download_file(url, t)
    log.info("Extracting archive to %s" , local_path)
    shutil.unpack_archive(t.name, extract_dir=str(local_path))



#TODO: Make this async someday
class RemoteLoader(LoaderBase):

    #TODO: Move these to nnprofile
    fit_url = 'http://pcteserver.mi.infn.it/~nnpdf/fits/'
    fit_index = 'fitdata.json'

    theory_url = 'http://pcteserver.mi.infn.it/~apfelcomb/commondatatheory/'
    theory_index = 'theorydata.json'

    lhapdf_url = 'http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/'

    @property
    @functools.lru_cache()
    def remote_theories(self):
        index_url = self.theory_url + self.theory_index
        token = 'theory_'
        try:
            resp = requests.get(index_url)
        except Exception as e:
            raise RemoteLoaderError("Failed to fetch remote fits index %s: %s" % (index_url,e)) from e

        try:
            info = resp.json()['files']
        except Exception as e:
            raise RemoteLoaderError("Malformed index %s. Expecting json with a key 'files': %s" % (index_url, e)) from e

        return {th[len(token):].split('.')[0] : self.theory_url+th for th in info}

    @property
    @functools.lru_cache()
    def remote_fits(self):
        index_url = self.fit_url+self.fit_index
        try:
            resp = requests.get(index_url)
        except Exception as e:
            raise RemoteLoaderError("Failed to fetch remote fits index %s: %s" % (index_url,e)) from e

        try:
            info = resp.json()['files']
        except Exception as e:
            raise RemoteLoaderError("Malformed index %s. Expecting json with a key 'files': %s" % (index_url, e)) from e

        return {fit.split('.')[0] : self.fit_url+fit for fit in info}

    @property
    def downloadable_fits(self):
        return list(self.remote_fits)

    @property
    def downloadable_theories(self):
        return list(self.remote_theories)

    @property
    def lhapdf_pdfs(self):
        return lhaindex.expand_index_names('*')

    @property
    def downloadable_pdfs(self):
        return set((*self.lhapdf_pdfs, *self.downloadable_fits))


    def download_fit(self, fitname):
        if not fitname in self.remote_fits:
            raise FitNotFound("Could not find fit '{}' in remote index {}".format(fitname, self.fit_index))

        #TODO: Change the crazy paths in fitmanager. Why depend on results??
        download_and_extract(self.remote_fits[fitname], self.resultspath.parent)

        fitpath = self.resultspath / fitname
        if lhaindex.isinstalled(fitname):
            log.warn("The PDF corresponding to the downloaded fit '%s' "
             "exists in the LHAPDF path."
             " Will be erased and replaced with the new one.", fitname)
            p = pathlib.Path(lhaindex.finddir(fitname))
            if p.is_symlink():
                p.unlink()
            else:
                shutil.rmtree(str(p))
        else:
            p = osp.join(lhaindex.get_lha_datapath(), fitname)
        gridpath = fitpath / 'nnfit' / fitname
        shutil.copytree(str(gridpath), str(p))


    def download_pdf(self, name):
        #It would be good to use the LHAPDF command line, except that it does
        #stupid things like returning 0 exit status when it fails to download
        if name in self.lhapdf_pdfs:
            url = self.lhapdf_url + name + '.tar.gz'
            download_and_extract(url, lhaindex.get_lha_datapath())
        elif name in self.downloadable_fits:
            self.download_fit(name)
        else:
            raise PDFNotFound("PDF %s is neither an uploaded fit nor a LHAPDF set.")

    def download_theoryID(self, thid):
        thid = str(thid)
        remote = self.remote_theories
        if thid not in remote:
            raise TheoryNotFound("Theory %s not available." % thid)
        download_and_extract(remote[thid], self.datapath)


class FallbackLoader(Loader, RemoteLoader):
    """A loader that first tries to find resources locally
    (calling Loader.check_*) and if it fails, it tries to download them
    (calling RemoteLoader.download_*)."""

    def make_checker(self, resource):
        #We are intercepting the check_
        orig = super().__getattribute__('check_' + resource)
        download = getattr(self, 'download_' + resource)

        @functools.wraps(orig)
        def f(*args, **kwargs):
            try:
                return orig(*args, **kwargs)
            except LoadFailedError as e:
                saved_exception = e
                log.info("Could not find a resource (%s): %s. "
                "Attempting to download it.", resource, saved_exception)
                try:
                    download(*args, **kwargs)
                except RemoteLoaderError as e:
                    log.error("Failed to download resource: %s",  e)
                    raise e
                except LoadFailedError as e:
                    log.error("Resource not in the remote repository: %s", e)
                    raise saved_exception
                except requests.RequestException as e:
                    log.error("There was a problem with the connection: %s", e)
                    raise saved_exception from e

                except Exception as e:
                    #Simply raise these for now so we can find and fix them
                    raise e
                else:
                    return orig(*args, **kwargs)
        return f



    def __getattribute__(self, attr):
        token = 'check_'
        if attr.startswith(token):
            resname = attr[len(token):]
            if hasattr(RemoteLoader, 'download_' + resname):
                return super().__getattribute__('make_checker')(resname)
        return super().__getattribute__(attr)
