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
#We can't use this for type checking yet
#import typing

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from reportengine.configparser import ConfigError
from reportengine import configparser
#from reportengine.broadcast import broadcast

from NNPDF import CommonData, FKTable, ThPredictions, LHAPDFSet
from NNPDF.fkset import FKSet
from NNPDF.dataset import DataSet

import lhaindex


class Environment:
    """Container for information to be filled at run time"""
    def __init__(self,*, data_path, output_path, this_folder):
        self.deta_path = pathlib.Path(data_path)
        self.output_path = pathlib.Path(output_path)
        self.this_folder = pathlib.Path(this_folder)
        self.loader = Loader(data_path)



class Config(configparser.Config):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    @property
    def loader(self):
        return self.environment.loader

    def check_pdf(self, name:str):
        if lhaindex.isinstalled(name):
            try:
                get_error_type(name)
            except NotImplementedError:
                raise ConfigError(str(e))
            return name
        raise ConfigError("Bad PDF: {} not installed".format(name), name,
                          lhaindex.expand_local_names('*'))

    def check_theoryid(self, theoryID: (str, int)):
        try:
            return (str(theoryID), self.loader.check_theoryID(theoryID))
        except FileNotFoundError as e:
            raise ConfigError(str(e), theoryID,
                              self.loader.available_theories,
                              display_alternatives='all')

    @configparser.element_of('datasets')
    def check_dataset(self, dataset:dict, * ,theoryid):
        """We check data related things here, and theory related things when
        making the FKSet"""
        theoryid, theopath = theoryid
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
                cfac = self.loader.check_cfactor(theoryid, name, cfac)
            else:
                cfac = []
    
            fkpath = self.loader.check_fktable(theoryid, name)
        
        except FileNotFoundError as e:
            raise ConfigError(e)
        
        return {'setname': name, 'commondata': commondata, 'cfac': cfac, 
                'fk': fkpath}


class DataNotFoundError(FileNotFoundError): pass

class SysNotFoundError(FileNotFoundError): pass

class Loader():
    """Load various resources from the NNPDF dapa path."""

    def __init__(self, datapath):
        if isinstance(datapath, str):
            datapath = pathlib.Path(datapath)
        self.datapath = datapath

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



def load_pdf(name):
    """Return a LHAPDFSet with the correct errro type"""
    err_type = get_error_type(name)
    return LHAPDFSet(name, err_type)


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

def results(setname, commondata, cfac, fk ,theoryid, pdf):
    cdpath,syspth = commondata
    cd = CommonData.ReadFile(str(cdpath), str(syspth))
    thlabel, thpath = theoryid
    
    fktable = FKTable(str(fk), [str(factor) for factor in cfac])
    #IMPORTANT: We need to tell the python garbage collector to NOT free the
    #memory owned by the FKTable on garbage collection.
    #TODO: Do this automatically
    fktable.thisown = 0
    fkset = FKSet(FKSet.parseOperator("NULL"), [fktable])
    data = DataSet(cd, fkset)
    pdf = load_pdf(pdf)
    th_predictions = ThPredictions(pdf, data)

    return DataResult(data), ThPredictionsResult(thlabel, th_predictions)

def chi2_data(results):
    data_result, th_result = results
    diffs = th_result.rawdata.T - data_result.central_value
    #chi²_i = diff_ij @ invcov_jk @ diff_ki
    return np.einsum('ij, jk, ik -> i',
                     diffs, data_result.invcovmat, diffs)/len(data_result)

#TODO; Implement for hessian
def chi2_stats(chi2_data):
    return {
            r'$\left \chi^2 \right>$': chi2_data.mean(),
            r'std(\chi²)'            : chi2_data.std(),
           }


def plot_chi2dist(results, setlabel, chi2_data, chi2_stats):
    fig, ax = plt.subplots()
    ax.set_title("$\chi^2$ distribution for %s" % setlabel)
    ax.hist(chi2_data)
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



class DataResult:

    def __init__(self, dataobj):
        self.dataobj = dataobj

    @property
    def label(self):
        return "CommonData"

    @property
    @functools.lru_cache()
    def central_value(self):
        return self.dataobj.get_cv()

    @property
    @functools.lru_cache()
    def std_error(self):
        return np.sqrt(np.diag(self.dataobj.get_covmat()))

    def __len__(self):
        return len(self.dataobj)

    def __getattr__(self, attr):
        return getattr(self.dataobj, attr)


class ThPredictionsResult(DataResult):

    def __init__(self, thlabel, dataobj):
        self.thlabel = thlabel
        super().__init__(dataobj)

    @property
    @functools.lru_cache()
    def std_error(self):
        return self.dataobj.get_error()

    @property
    @functools.lru_cache()
    def rawdata(self):
        return self.dataobj.get_data()

    def label(self):
        return "<Theory %s>: %s" % (self.thlabel, self.dataobj.GetPDFName())


if __name__ == '__main__':
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

    loglevel.add_argument('-q','--quiet', help="Supress INFO messages",
                        action='store_true')

    loglevel.add_argument('-d', '--debug', help = "Show debug info",
                          action='store_true')

    parser.add_argument('-p','--datapath', help="Path where the NNPDF "
                        "data is located",
                        default='../nnpdfcpp/data')

    args = parser.parse_args()

    if args.quiet:
        level = logging.WARN
    elif args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(format='%(levelname)s: %(message)s', level=level)



    environment = Environment(data_path = args.datapath,
                              output_path = args.output,
                              this_folder = pathlib.Path(__file__).parent
                              )
    with open(args.config_yml) as f:
        try:
            c = Config.from_yaml(f, environment=environment)
        except ConfigError as e:
            print("Bad configuration encountered:")
            print(e)
            print(e.alternatives_text())
            sys.exit(1)
    
    environment.output_path.mkdir(exist_ok = True)
    plt.style.use(str(environment.this_folder / 'small.mplstyle'))

    #TODO; Use resourcebuilder instead of this horrible code
    n = itertools.count()
    for dataset in c['datasets']:
        inp = {**dataset, 'theoryid': c['theoryid'], 'pdf': c['pdf']}
        results = (dataresult, theoresult)= results(**inp)


        fig = plot_results(results, setlabel=dataset['setname'],
                       normalize_to=None)

        fig.savefig(str(environment.output_path / ("result_%2d.pdf"%next(n))),
                bbox_inches='tight')

        chi2dt = chi2_data(results)
        chi2st = chi2_stats(chi2dt)
        fig = plot_chi2dist(results, dataset['setname'],
                            chi2dt, chi2st)
        fig.savefig(str(environment.output_path / ("chi2_%2d.pdf"%next(n))),
                bbox_inches='tight')
