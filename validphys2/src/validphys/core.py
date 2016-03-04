# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 10:21:36 2016

@author: Zahari Kassabov
"""
import sys
import pathlib
import functools
import logging
import argparse
import itertools
from collections import OrderedDict
#We can't use this for type checking yet
#import typing

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats

from reportengine.configparser import ConfigError
from reportengine import configparser
#from reportengine.broadcast import broadcast

from NNPDF import CommonData, FKTable, ThPredictions, LHAPDFSet
from NNPDF.fkset import FKSet
from NNPDF.dataset import DataSet

from validphys import lhaindex

log = logging.getLogger(__name__)

class PDFDoesNotExist(Exception): pass

class _PDFSETS():
    """Convenient way to access installed PDFS, by e.g. tab completing
       in ipython."""
    def __getattr__(self, attr):
        if lhaindex.isinstalled(attr):
            return PDF(attr)
        raise AttributeError()

    def __dir__(self):
        return lhaindex.expand_local_names('*')


PDFSETS = _PDFSETS()

class PDF:

    def __init__(self, name):
        self.name = name

    def __getattr__(self, attr):
        try:
            return lhaindex.parse_info(self.name)[attr]
        except KeyError:
            raise AttributeError("'%r' has no attribute '%s'" % (type(self),
                                                                 attr))
        except IOError:
            raise PDFDoesNotExist(self.name)

    @property
    def stats_class(self):
        """Return the stats calculator for this error type"""
        error = self.ErrorType
        klass = STAT_TYPES[error]
        if hasattr(self, 'ErrorConfLevel'):
            klass = functools.partial(klass, rescale_factor=self.rescale_factor)
        return klass

    @property
    def infopath(self):
        return lhaindex.infofilename(self.name)

    @property
    def rescale_factor(self):
        if hasattr(self, "ErrorConfLevel"):
            if self.ErrorType == 'replicas':
                raise ValueError("Attribute at %s 'ErrorConfLevel' doesn't "
                "make sense with 'replicas' error type" % self.infopath)
            val = scipy.stats.norm.isf((1 - 0.01*self.ErrorConfLevel)/2)
            if np.isnan(val):
                raise ValueError("Invalid 'ErrorConfLevel' of PDF %s: %s" %
                                 (self, val))
            return val
        else:
            return 1

    def load(self):
        return LHAPDFSet(self.name, self.nnpdf_error)



    @property
    def nnpdf_error(self):
        """Return the NNPDF error tag, used to build the `LHAPDFSet` objeect"""
        error = self.ErrorType
        if error == "replicas":
            return LHAPDFSet.ER_MC
        if error == "hessian":
            if hasattr(self, 'ErrorConfLevel'):
                cl = self.ErrorConfLevel
                if cl == 90:
                    return LHAPDFSet.ER_EIG90
                elif cl == 68:
                    return LHAPDFSet.ER_EIG
                else:
                    raise NotImplementedError("No hessian errors with confidence"
                                              " interval %s" % (cl,) )
            else:
                return LHAPDFSet.ER_EIG

        if error == "symmhessian":
            if hasattr(self, 'ErrorConfLevel'):
                cl = self.ErrorConfLevel
                if cl == 68:
                    return LHAPDFSet.ER_SYMEIG
                else:
                    raise NotImplementedError("No symmetric hessian errors "
                                              "with confidence"
                                              " interval %s" % (cl,) )
            else:
                return LHAPDFSet.ER_EIG

        raise NotImplementedError("Error type for %s: '%s' is not implemented" %
                                  (self.name, error))

class DataSetSpec:

    def __init__(self, *, name, commondata, cfac, fkpath, thspec, cuts):
        self.name = name
        self.commondata = commondata
        self.cfac = cfac
        self.fkpath = fkpath
        self.thspec = thspec
        self.cuts = cuts


class Stats:

    def __init__(self, data):
        """`data `should be N_pdf*N_bins"""
        self.data = np.atleast_2d(data)

    def central_value(self):
        raise NotImplementedError()

    def std(self):
        raise NotImplementedError()

    def sample_values(self, size):
        raise NotImplementedError()

    def std_interval(self, nsigma):
        raise NotImplementedError()

    def errorbar68(self):
        raise NotImplementedError()

    #TODO...
    ...

class MCStats(Stats):
    """Result obtained from a Monte Carlo sample"""

    def central_value(self):
        return np.mean(self.data, axis=0)

    def std_error(self):
        return np.std(self.data, axis=0)

    def sample_values(self, size):
        return np.random.choice(self, size=size)


class SymmHessianStats(Stats):
    """Compute stats in the 'assymetric' hessian format: The first index (0)
    is the
    central value. The rest of the indexes are results for each eigenvector.
    A 'rescale_factor is allowed in case the eigenvector confidence interval
    is not 68%'."""
    def __init__(self, data, rescale_factor=1):
        super().__init__(data)
        self.rescale_factor = rescale_factor

    def central_value(self):
        return self.data[0]

    def std_error(self):
        data = self.data
        diffsq = (data[0] - data[1:])**2
        return np.sqrt(diffsq.sum(axis=0))/self.rescale_factor

class HessianStats(SymmHessianStats):
    """Compute stats in the 'assymetric' hessian format: The first index (0)
    is the
    central value. The odd indexes are the
    results for lower eigenvectors and the
    even are the upper eigenvectors.A 'rescale_factor is allowed in
    case the eigenvector confidence interval
    is not 68%'."""
    def std_error(self):
        data = self.data
        diffsq = (data[1::2] - data[2::2])**2
        return np.sqrt(diffsq.sum(axis=0))/self.rescale_factor


STAT_TYPES = dict(
                    symmhessian = SymmHessianStats,
                    hessian = HessianStats,
                    replicas   = MCStats,
                   )


class Environment:
    """Container for information to be filled at run time"""
    def __init__(self,*, data_path, results_path ,output_path, this_folder):
        self.deta_path = pathlib.Path(data_path)
        self.results_path = pathlib.Path(results_path)

        self.output_path = pathlib.Path(output_path)
        self.this_folder = pathlib.Path(this_folder)

        self.loader = Loader(data_path, resultspath=self.results_path)



class Config(configparser.Config):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    @property
    def loader(self):
        return self.environment.loader

    def parse_pdf(self, name:str):
        if lhaindex.isinstalled(name):
            pdf = PDF(name)
            try:
                pdf.nnpdf_error
            except NotImplementedError as e:
                raise ConfigError(str(e))
            return pdf
        raise ConfigError("Bad PDF: {} not installed".format(name), name,
                          lhaindex.expand_local_names('*'))

    def parse_theoryid(self, theoryID: (str, int)):
        try:
            return (str(theoryID), self.loader.check_theoryID(theoryID))
        except FileNotFoundError as e:
            raise ConfigError(str(e), theoryID,
                              self.loader.available_theories,
                              display_alternatives='all')

    def parse_use_cuts(self, use_cuts:bool, *, fit=None):
        if use_cuts and not fit:
            raise ConfigError("Setting 'use_cuts' true requires "
            "specifying a fit, e.g. 'fit' on which filter "
            "has been executed, e.g.\nfit:NNPDF30_nlo_as_0118")
        return use_cuts



    #TODO: load fit config from here
    def parse_fit(self, fit:str):
        try:
            return fit, self.loader.check_fit(fit)
        except Exception as e:
            raise ConfigError(str(e), self.loader.available_fits)


    @configparser.element_of('datasets')
    def parse_dataset(self, dataset:dict, * ,theoryid, use_cuts, fit=None):
        """We check data related things here, and theory related things when
        making the FKSet"""
        theoryno, theopath = theoryid
        try:
            name, sysnum = dataset['dataset'], dataset['sys']
        except KeyError:
            raise ConfigError("'dataset' must be a mapping with "
                              "'name' and 'sysnum'")

        try:
            commondata = self.loader.check_commondata(name, sysnum)
        except DataNotFoundError as e:
            raise ConfigError(str(e), name, self.loader.available_datasets)
        except SysNotFoundError as e:
            raise ConfigError(str(e))

        try:
            if 'cfac' in dataset:
                cfac = dataset['cfac']
                cfac = self.loader.check_cfactor(theoryno, name, cfac)
            else:
                cfac = []

            fkpath = self.loader.check_fktable(theoryno, name)

        except FileNotFoundError as e:
            raise ConfigError(e)


        if use_cuts:
            cuts = self.loader.get_cuts(name, fit)
        else:
            cuts = None

        return DataSetSpec(name=name, commondata=commondata, cfac=cfac,
                           fkpath=fkpath, thspec=theoryid, cuts=cuts)


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
        sysfile = (self.datapath / 'commondata' / 'systypes' /
                   ('SYSTYPE_%s_%d.dat' % (setname, sysnum)))

        if not sysfile.exists():
            raise SysNotFoundError(("Could not find systype %d for "
                 "dataset '%s'. File %s does not exist.") % (sysnum, setname,
                  sysfile))

        return datafile, sysfile

    @functools.lru_cache()
    def check_theoryID(self, theoryID):
        theoryID = str(theoryID)
        theopath = self.datapath / ('theory_%s' % theoryID)
        if not theopath.exists():
            raise FileNotFoundError(("Could not find theory %s. "
                  "Folder '%s' not found") % (theoryID, theopath) )
        return theopath

    def get_commondata(self, setname, sysnum):
        """Get a Commondata from the set name and number"""
        datafile, sysfile = self.check_commondata(setname, sysnum)
        return CommonData.ReadFile(str(datafile), str(sysfile))

    @functools.lru_cache()
    def check_fktable(self, theoryID, setname):
        theopath = self.check_theoryID(theoryID)
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
            return p
        if not p.exists():
            msg = ("Could not find fit '{firname}' in '{resultspath}. "
                   "Folder '{p}' not found").format(**locals())
            raise FileNotFoundError(msg)
        msg = ("Could not load fit '{firname}' from '{resultspath}. "
                   "'{p}' must be a folder").format(**locals())
        raise IOError(msg)

    @property
    def available_fits(self):
        return [p for p in self.resultspath.iterdir() if p.is_dir()]

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

def get_error_type(pdfname):
    info = lhaindex.parse_info(pdfname)
    error = info["ErrorType"]
    if error == "replicas":
        return LHAPDFSet.ER_MC
    if error == "hessian":
        if 'ErrorConfLevel' in info:
            cl = info['ErrorConfLevel']
            if cl == 90:
                return LHAPDFSet.ER_EIG90
            elif cl == 68:
                return LHAPDFSet.ER_EIG
            else:
                raise NotImplementedError("No hessian errors with confidence"
                                          " interval %s" % (cl,) )
        else:
            return LHAPDFSet.ER_EIG

    if error == "symmhessian":
        if 'ErrorConfLevel' in info:
            cl = info['ErrorConfLevel']
            if cl == 68:
                return LHAPDFSet.ER_SYMEIG
            else:
                raise NotImplementedError("No symmetric hessian errors "
                                          "with confidence"
                                          " interval %s" % (cl,) )
        else:
            return LHAPDFSet.ER_EIG

    raise NotImplementedError("Error type for %s: '%s' is not implemented" %
                              (pdfname, error))


#==============================================================================
#
# def load_fkset(dataset, theoryID):
#
#     return FKSet(FKSet.parseOperator("NULL"), [fktable])
#
# def get_dataset(commondata, fktable):
#     fkset = fktable_to_fkset(fktable)
#     return DataSet(commondata, fkset)
#
#
#
# def make_results(loader, setname, theoryid, pdfname):
#     """High level functions that produces results by parsing the config"""
#     fktable = loader.get_fktable(theoryid, setname)
#==============================================================================

def results(dataset:DataSetSpec, pdf:PDF):
    cdpath,syspth = dataset.commondata
    cd = CommonData.ReadFile(str(cdpath), str(syspth))
    thlabel, thpath = dataset.thspec

    fktable = FKTable(str(dataset.fkpath), [str(factor) for factor in dataset.cfac])
    #IMPORTANT: We need to tell the python garbage collector to NOT free the
    #memory owned by the FKTable on garbage collection.
    #TODO: Do this automatically
    fktable.thisown = 0
    fkset = FKSet(FKSet.parseOperator("NULL"), [fktable])

    data = DataSet(cd, fkset)

    if dataset.cuts is not None:
        #ugly need to convert from numpy.int64 to int, so we can pass
        #it happily to the vector to the SWIG wrapper.
        #Do not do this (or find how to enable in SWIG):
        #data = DataSet(data, list(dataset.cuts))
        intmask = [int(ele) for ele in dataset.cuts]
        data = DataSet(data, intmask)

    nnpdf_pdf = pdf.load()
    th_predictions = ThPredictions(nnpdf_pdf, data)

    stats = pdf.stats_class


    return DataResult(data), ThPredictionsResult(thlabel, th_predictions,
                                                 stats)

def chi2_data(results):
    data_result, th_result = results
    diffs = th_result._rawdata.T - data_result.central_value
    #chi²_i = diff_ij @ invcov_jk @ diff_ki
    result =  np.einsum('ij, jk, ik -> i',
                     diffs, data_result.invcovmat, diffs)/len(data_result)

    return th_result.stats_class(result[:, np.newaxis])

def chi2_stats(chi2_data):
    return OrderedDict([
            (r'$\left< \chi^2 \right>$'  , chi2_data.central_value().mean() ),
            (r'std($\chi^2$)'            , chi2_data.std_error().mean()     ),
           ])


def plot_chi2dist(results, setlabel, chi2_data, chi2_stats, pdf):
    fig, ax = plt.subplots()
    label = pdf.name
    if not isinstance(chi2_data, MCStats):
        ax.set_axis_bgcolor("#ffcccc")
        log.warn("Chi² distribution plots have a different meaning for non MC sets.")
        label += " (%s!)" % pdf.ErrorType
    label += '\n'+ ', '.join(str(k)+(' %.2f' % v) for (k,v) in chi2_stats.items())
    ax.set_title("$\chi^2$ distribution for %s" % setlabel)
    ax.hist(chi2_data.data, label=label, zorder=10000)
    ax.legend()
    return fig


def plot_results(results, setlabel, normalize_to = None):

    cvs = [r.central_value for r in results]
    errors = [r.std_error for r in results]

    y_label = "Observable value"

    if normalize_to is not None:
        y_label = "Ratio to %s" % normalize_to.label
        norm_cv = normalize_to.central_value
        cvs = [cv/norm_cv for cv in cvs]
        errors = [e/norm_cv for e in errors]

    l = len(results)
    if l < 5:
        delta = iter(np.linspace(-0.05*l, 0.05*l, l))
    else:
        delta = iter(np.linspace(-0.25, 0.25, l))

    figure, ax = plt.subplots()

    ax.set_title(setlabel)
    ax.set_ylabel(y_label)

    for result, cv, error in zip(results, cvs, errors):
        x = np.arange(1, len(result) + 1) + next(delta)
        ax.errorbar(x, cv, yerr=error,
                         linestyle='none',
                         label=result.label, elinewidth = 2,
                         capsize=10)

    ax.set_xticks(range(1,len(result) + 1), len(result) // 15 + 1)

    ax.grid(axis='x')

    #TODO: Abtract this out
    dfs = ax.get_yticks() - 1
    l = len(dfs) // 2  + 1 - ((len(dfs) // 2) % 2)
    mdiff = np.max(dfs)
    ax.set_yticks(np.linspace(-mdiff, mdiff, l) + 1)



    ax.legend()

    return figure


class Result:
    def __init__(self, dataobj):
        self.dataobj = dataobj

    @property
    @functools.lru_cache()
    def std_error(self):
        return np.sqrt(np.diag(self.covmat))

    @property
    @functools.lru_cache()
    def central_value(self):
        return self.dataobj.get_cv()

    def __len__(self):
        return len(self.dataobj)

    def __getattr__(self, attr):
        return getattr(self.dataobj, attr)


class DataResult(Result):

    @property
    def label(self):
        return "CommonData"

    @property
    @functools.lru_cache()
    def covmat(self):
        return self.dataobj.get_covmat()

    @property
    @functools.lru_cache()
    def invcovmat(self):
        return self.dataobj.get_invcovmat()


class ThPredictionsResult(Result):

    def __init__(self, thlabel, dataobj, stats_class):
        self.thlabel = thlabel
        self.stats_class = stats_class
        super().__init__(dataobj)

    @property
    @functools.lru_cache()
    def std_error(self):
        return self.dataobj.get_error()

    @property
    @functools.lru_cache()
    def _rawdata(self):
        return self.dataobj.get_data()

    @property
    def data(self):
        return self.stats_class(self._rawdata)

    @property
    def label(self):
        return "<Theory %s>@%s" % (self.thlabel, self.dataobj.GetPDFName())


def main():
    #TODO: Oberhaul this to use reportengine properly

    parser = argparse.ArgumentParser(
             description = "Validphys developer preview",
             )

    parser.add_argument('config_yml',
                        help = "path to the configuration file")

    parser.add_argument('-o','--output', help="output folder where to "
                                         "store resulting plots and tables",
                        default='output')

    loglevel = parser.add_mutually_exclusive_group()

    loglevel.add_argument('-q','--quiet', help="supress INFO messages and C output",
                        action='store_true')

    loglevel.add_argument('-d', '--debug', help = "show debug info",
                          action='store_true')

    parser.add_argument('-p','--datapath', help="path where the NNPDF "
                        "data is located",
                        default='../nnpdfcpp/data')

    parser.add_argument('--resultspath', help="path where the fit results "
                          "are located. Calculated from 'datapath' by default",
                         )

    args = parser.parse_args()

    if args.quiet:
        level = logging.WARN
    elif args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(format='%(levelname)s: %(message)s', level=level)

    if not args.resultspath:
        args.resultspath = pathlib.Path(args.datapath).parent / 'nnpdfbuild' / 'results'

    environment = Environment(data_path = args.datapath,
                              results_path = args.resultspath,
                              output_path = args.output,
                              this_folder = pathlib.Path(__file__).parent
                              )

    try:
        with open(args.config_yml) as f:
            try:
                c = Config.from_yaml(f, environment=environment)
            except ConfigError as e:
                print("Bad configuration encountered:")
                print(e)
                print(e.alternatives_text())
                sys.exit(1)
    except OSError as e:
        print("Could not open configuration file: %s" % e)
        sys.exit(1)

    environment.output_path.mkdir(exist_ok = True)
    plt.style.use(str(environment.this_folder / 'small.mplstyle'))

    #TODO; Use resourcebuilder instead of this horrible code
    n = itertools.count()
    for dataset in c['datasets']:
        results_ = (dataresult, theoresult)= results(dataset, c['pdf'])


        fig = plot_results(results_, setlabel=dataset.name,
                       normalize_to=dataresult)

        fig.savefig(str(environment.output_path / ("result_%2d.pdf"%next(n))),
                bbox_inches='tight')

        chi2dt = chi2_data(results_)
        chi2st = chi2_stats(chi2dt)
        fig = plot_chi2dist(results_, dataset.name,
                            chi2dt, chi2st, c['pdf'])
        fig.savefig(str(environment.output_path / ("chi2_%2d.pdf"%next(n))),
                bbox_inches='tight')
