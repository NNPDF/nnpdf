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
import os
import os.path as osp
import urllib.parse as urls
import mimetypes

import yaml
import requests

from reportengine import filefinder
from validphys.core import (CommonDataSpec, FitSpec, TheoryIDSpec, FKTableSpec,
                            PositivitySetSpec, DataSetSpec, PDF, Cuts,
                            peek_commondata_metadata)
from validphys import lhaindex
import NNPDF as nnpath


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

    def __init__(self, datapath=None, resultspath=None):
        if datapath is None:
            datapath = nnpath.get_data_path()

        if resultspath is None:
            resultspath = nnpath.get_results_path()
        datapath = pathlib.Path(datapath)
        resultspath = pathlib.Path(resultspath)

        if not datapath.exists():
            raise LoaderError(f"The data path {datapath} does not exist.")

        if not resultspath.exists():
             raise LoaderError(f"The results path {resultspath} does not exist.")

        self.datapath = datapath
        self.resultspath = resultspath
        self._nnprofile = None

    @property
    def nnprofile(self):
        if self._nnprofile is None:
            profile_path = nnpath.get_profile_path()
            if not osp.exists(profile_path):
                raise LoaderError(f"Could not find the profile path at "
                                  "{profile_path}. Check your libnnpdf configuration")
            with open(profile_path) as f:
                try:
                    self._nnprofile = yaml.safe_load(f)
                except yaml.YAMLError as e:
                    raise LoaderError(f"Could not parse profile file "
                                      f"{profile_path}: {e}") from e
        return self._nnprofile
    def _vp_cache(self):
        """Return the vp-cache path, and create it if it doesn't exist"""
        vpcache = pathlib.Path(self.nnprofile['validphys_cache_path'])
        if not vpcache.exists():
            try:
                log.info(f"Creating validphys cache directory: {vpcache}")
                vpcache.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                raise LoaderError("Could not create the cache directory "
                              f"at {vpcache}") from e
        return vpcache


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

        if not sysfile.exists():
            raise SysNotFoundError(("Could not find systype %s for "
                 "dataset '%s'. File %s does not exist.") % (sysnum, setname,
                  sysfile))

        plotfiles = []


        metadata = peek_commondata_metadata(datafile)
        process_plotting_root = self.commondata_folder/f'PLOTTINGTYPE_{metadata.process_type}'
        type_plotting = (process_plotting_root.with_suffix('.yml'),
                         process_plotting_root.with_suffix('.yaml'),)

        data_plotting_root = self.commondata_folder/f'PLOTTING_{setname}'

        data_plotting = (data_plotting_root.with_suffix('.yml'),
                         data_plotting_root.with_suffix('.yaml'),
                        )
        #TODO: What do we do when both .yml and .yaml exist?
        for tp in (type_plotting, data_plotting):
            for p in tp:
                if p.exists():
                    plotfiles.append(p)

        return CommonDataSpec(datafile, sysfile, plotfiles, name=setname, metadata=metadata)

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
        cd = self.check_commondata(setname, sysnum, plotfiles)
        return cd.load()

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
        cd = self.check_commondata(setname, 'DEFAULT')
        fk = self.check_fktable(theiryID, setname, [])
        th =  self.check_theoryID(theiryID)
        return PositivitySetSpec(setname, cd, fk, postlambda, th)

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
            cuts = self.check_cuts(name, fit)
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

    def check_cuts(self, setname, fit):
        if fit is None:
            raise TypeError("Must specify a fit to use the cuts.")
        if not isinstance(fit, FitSpec):
            fit = self.check_fit(fit)
        fitname, fitpath = fit
        p = (fitpath/'filter')/setname/('FKMASK_' + setname+ '.dat')
        if not p.parent.exists():
            raise CutsNotFound("Bad filter configuration. "
            "Could not find: %s" % p.parent)
        if not p.exists():
            return None
        return Cuts(setname, p)


    def get_cuts(self, setname, fit):
        cuts = self.check_cuts(setname, fit)
        if cuts:
            return cuts.load()
        return None

    def check_vp_output_file(self, filename, extra_paths=('.',)):
        """Find a file in the vp-cache folder, or (with higher priority) in
        the ``extra_paths``."""
        try:
            vpcache = self._vp_cache()
        except KeyError as e:
            log.warn("Entry validphys_cache_path expected but not found "
                     "in the nnprofile.")
        else:
            extra_paths = (*extra_paths, vpcache)

        finder = filefinder.FallbackFinder(extra_paths)
        try:
            path, name = finder.find(filename)
        except FileNotFoundError as e:
            raise LoadFailedError(f"Could not find '{filename}'") from e
        except filefinder.FinderError as e:
            raise LoaderError(e) from e
        return path/name


#http://stackoverflow.com/a/15645088/1007990
def _download_and_show(response, stream):

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

def download_file(url, stream_or_path, make_parents=False):
    """Download a file and show a progress bar if the INFO log level is
    enabled. If ``make_parents`` is ``True`` ``stream_or_path``
    is path-like, all the parent folders will
    be created."""
    #There is a bug in CERN's
    #Apache that incorrectly sets the Content-Encodig header to gzip, even
    #though it doesn't compress two times.
    # See: http://mail-archives.apache.org/mod_mbox/httpd-dev/200207.mbox/%3C3D2D4E76.4010502@talex.com.pl%3E
    # and e.g. https://bugzilla.mozilla.org/show_bug.cgi?id=610679#c30
    #If it looks like the url is already encoded, we do not request
    #it to be compressed
    headers = {}
    if mimetypes.guess_type(url)[1] is not None:
        headers['Accept-Encoding'] = None

    response = requests.get(url, stream=True, headers=headers)

    response.raise_for_status()

    if isinstance(stream_or_path, (str, bytes, os.PathLike)):
        log.info("Downloading %s to %s." , url, stream_or_path)
        if make_parents:
            pathlib.Path(stream_or_path).parent.mkdir(exist_ok=True, parents=True)

        with open(stream_or_path, 'wb') as f:
            return _download_and_show(response, f)
    else:
        log.info("Downloading %s." , url,)
        return _download_and_show(response, stream_or_path)


def download_and_extract(url, local_path):
    """Download a compressed archive and then extract it to the given path"""
    name = url.split('/')[-1]
    with tempfile.NamedTemporaryFile(delete=False, suffix=name) as t:
        log.debug("Saving data to %s" , t.name)
        download_file(url, t)
    log.info("Extracting archive to %s" , local_path)
    shutil.unpack_archive(t.name, extract_dir=str(local_path))




def _key_or_loader_error(f):
    @functools.wraps(f)
    def f_(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except KeyError as e:
            log.error(f"nnprofile is configured "
                      f"improperly: Key {e} is missing! "
                      f"Fix it at {nnpath.get_profile_path()}")
            raise LoaderError("Cannot attempt download because "
                              "nnprofile is configured improperly: "
                              f"Missing key '{e}'") from e
    return f_


#TODO: Make this async someday
class RemoteLoader(LoaderBase):

    @property
    @_key_or_loader_error
    def fit_urls(self):
        return self.nnprofile['fit_urls']

    @property
    @_key_or_loader_error
    def fit_index(self):
        return self.nnprofile['fit_index']

    @property
    @_key_or_loader_error
    def theory_urls(self):
        return self.nnprofile['theory_urls']

    @property
    @_key_or_loader_error
    def theory_index(self):
        return self.nnprofile['theory_index']

    @property
    @_key_or_loader_error
    def nnpdf_pdfs_urls(self):
        return self.nnprofile['nnpdf_pdfs_urls']

    @property
    @_key_or_loader_error
    def nnpdf_pdfs_index(self):
        return self.nnprofile['nnpdf_pdfs_index']

    @property
    @_key_or_loader_error
    def lhapdf_urls(self):
        return self.nnprofile['lhapdf_urls']


    def _remote_files_from_url(self, url, index, thing='files'):
        index_url = url + index
        try:
            resp = requests.get(index_url)
            resp.raise_for_status()
        except Exception as e:
            raise RemoteLoaderError("Failed to fetch remote %s index %s: %s" % (thing, index_url,e)) from e

        try:
            info = resp.json()['files']
        except Exception as e:
            raise RemoteLoaderError("Malformed index %s. Expecting json with a key 'files': %s" % (index_url, e)) from e

        return {file.split('.')[0] : url+file for file in info}

    def remote_files(self, urls, index, thing='files'):
        d = {}
        for url in urls:
            d.update(self._remote_files_from_url(url, index, thing))
        return d

    @property
    @functools.lru_cache()
    def remote_fits(self):
        return self.remote_files(self.fit_urls, self.fit_index, thing="fits")

    @property
    @functools.lru_cache()
    def remote_theories(self):
        token = 'theory_'
        rt = self.remote_files(self.theory_urls, self.theory_index, thing="theories")
        return {k[len(token):]: v for k,v in rt.items()}

    @property
    @functools.lru_cache()
    def remote_nnpdf_pdfs(self):
        return self.remote_files(self.nnpdf_pdfs_urls, self.nnpdf_pdfs_index,
                                 thing="PDFs")

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
    def nnpdf_pdfs(self):
        return list(self.remote_nnpdf_pdfs)

    @property
    def downloadable_pdfs(self):
        return set((*self.lhapdf_pdfs, *self.downloadable_fits, *self.nnpdf_pdfs))


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
        #Check if the pdf is an existing fit first
        try:
            #We don't want to download the fit here
            fit = Loader.check_fit(self, name)
        except FitNotFound:
            pass
        else:
            if (fit.path/'nnfit').exists():
                p = osp.join(lhaindex.get_lha_datapath(), fit.name)
                log.info("Found existing fit with the same name as the "
                "requested PDF (%s). Copying the grid to the LHAPDF path (%s).",
                name, p)

                gridpath = fit.path / 'nnfit' / fit.name
                shutil.copytree(str(gridpath), str(p))
                return

        #It would be good to use the LHAPDF command line, except that it does
        #stupid things like returning 0 exit status when it fails to download
        if name in self.lhapdf_pdfs:
            url = self.lhapdf_url + name + '.tar.gz'
            download_and_extract(url, lhaindex.get_lha_datapath())
        elif name in self.downloadable_fits:
            self.download_fit(name)
        elif name in self.remote_nnpdf_pdfs:
            download_and_extract(self.remote_nnpdf_pdfs[name], lhaindex.get_lha_datapath())
        else:
            raise PDFNotFound("PDF '%s' is neither an uploaded fit nor a LHAPDF set." % name)

    def download_theoryID(self, thid):
        thid = str(thid)
        remote = self.remote_theories
        if thid not in remote:
            raise TheoryNotFound("Theory %s not available." % thid)
        download_and_extract(remote[thid], self.datapath)

    def download_vp_output_file(self, filename, **kwargs):
        try:
            root_url = self.nnprofile['reports_root_url']
        except KeyError as e:
            raise LoadFailedError('Key report_root_url not found in nnprofile')
        try:
            url = root_url  + filename
        except Exception as e:
            raise LoadFailedError(e) from e
        try:
            filename = pathlib.Path(filename)

            download_file(url, self._vp_cache()/filename, make_parents=True)
        except requests.HTTPError as e:
            if e.response.status_code == requests.codes.not_found:
                raise RemoteLoaderError(f"Ressource {filename} could not "
                                        f"be found on the validphys "
                                        f"server {url}") from e
            elif e.response.status_code == requests.codes.unauthorized:
                log.error("Could not access the validphys reports page "
                          "because the authentification is not provided. "
                          "Please, update your ~/.netrc file to contain the "
                          "following:\n\n"
                          f"machine {urls.urlsplit(root_url).netloc}\n"
                          f"    login nnpdf\n"
                          f"    password <PASSWORD>\n"
                 )
            raise

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
                log.info("Could not find a resource "
                    f"({resource}): {saved_exception}. "
                    f"Attempting to download it.")
                try:
                    download(*args, **kwargs)
                except RemoteLoaderError as e:
                    log.error(f"Failed to download resource: {e}")
                    raise e
                except LoadFailedError as e:
                    log.error(f"Resource not in the remote repository: {e}")
                    raise saved_exception
                except requests.RequestException as e:
                    log.error(f"There was a problem with the connection: {e}")
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
