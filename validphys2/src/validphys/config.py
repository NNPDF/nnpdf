# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:43:10 2016

@author: Zahari Kassabov
"""
import logging
import pathlib
import functools
import inspect

from  reportengine import configparser
from reportengine.environment import Environment
from reportengine.configparser import ConfigError, element_of, _parse_func
from reportengine.helputils import get_parser_type
from reportengine import report

from validphys.core import ExperimentSpec, DataSetInput
from validphys.loader import (Loader, LoadFailedError ,DataNotFoundError,
                              PDFNotFound, FallbackLoader)

log = logging.getLogger(__name__)



class Environment(Environment):
    """Container for information to be filled at run time"""

    def __init__(self,*, datapath, resultspath, this_folder, net=True,
                 upload=False, **kwargs):
        self.deta_path = pathlib.Path(datapath)
        self.results_path = pathlib.Path(resultspath)
        self.this_folder = pathlib.Path(this_folder)

        if net:
            loader_class = FallbackLoader
        else:
            loader_class = Loader

        self.loader = loader_class(self.deta_path, resultspath=self.results_path)

        self.upload = upload
        super().__init__(**kwargs)


def _id_with_label(f):
    f = _parse_func(f)
    def parse_func(self, item, **kwargs):
        if not isinstance(item, dict):
            return f(self, item, **kwargs)
        keydiff =  {'id', 'label'} - item.keys()

        if  keydiff and 'id' in keydiff:
            raise ConfigError("'%s' must be a single id, or a mapping "
                              "with keys 'id', 'label'"%(item,))
        id = item['id']
        val = f(self, id, **kwargs)
        if 'label' in item:
            val.label = str(item['label'])
        return val

    currsig = inspect.signature(parse_func)
    origsig = inspect.signature(f)
    parse_func = functools.wraps(f)(parse_func)

    params = [*list(currsig.parameters.values())[:2],
              *list(origsig.parameters.values())[2:]]

    parse_func.__signature__ = inspect.Signature(
                 parameters=params)

    labeldoc =  ((" Either just an id %s, or a mapping "
                          "with 'id' and 'label'.") %
                          (get_parser_type(f),))
    if parse_func.__doc__ is None:
        parse_func.__doc__ = labeldoc
    else:
        parse_func.__doc__ += labeldoc

    return parse_func

class Config(report.Config):


    @property
    def loader(self):
        return self.environment.loader

    @element_of('pdfs')
    @_id_with_label
    def parse_pdf(self, name:str):
        """A PDF set installed in LHAPDF."""
        try:
            pdf = self.loader.check_pdf(name)
        except PDFNotFound as e:
            raise ConfigError("Bad PDF: {} not installed".format(name), name,
                          self.loader.available_pdfs) from e

        #Check that we know how to compute errors
        try:
            pdf.nnpdf_error
        except NotImplementedError as e:
            raise ConfigError(str(e))
        return pdf


    @element_of('theoryids')
    @_id_with_label
    def parse_theoryid(self, theoryID: (str, int)):
        """A number corresponding to the database theory ID where the
        corresponding theory folder is installed in te data directory."""
        try:
            return self.loader.check_theoryID(theoryID)
        except FileNotFoundError as e:
            raise ConfigError(str(e), theoryID,
                              self.loader.available_theories,
                              display_alternatives='all')

    def parse_use_cuts(self, use_cuts:bool, *, fit=None):
        """Whether to use the filtered points in the fit, or the whole
        data in the dataset."""
        if use_cuts and not fit:
            raise ConfigError("Setting 'use_cuts' true requires "
            "specifying a fit on which filter "
            "has been executed, e.g.\nfit : NNPDF30_nlo_as_0118")
        return use_cuts

    #TODO: load fit config from here
    @element_of('fits')
    @_id_with_label
    def parse_fit(self, fit:str):
        """A fit in the results folder, containing at least a valid filter result."""
        try:
            return self.loader.check_fit(fit)
        except LoadFailedError as e:
            raise ConfigError(str(e), fit ,self.loader.available_fits)

    def produce_fitcontext(self, fit):
        """Set PDF, theory ID and experiments from the fit config"""


        _, pdf         = self.parse_from_('fit', 'pdf', write=False)
        _, theory      = self.parse_from_('fit', 'theory', write=False)
        _, experiments = self.parse_from_('fit', 'experiments', write=False)

        #TODO: parse we need multilevel from to do theoruid nicely
        thid = theory['theoryid']

        return {'pdf': pdf, 'theoryid':thid, 'experiments': experiments}



    def parse_dataset_input(self, dataset):
        try:
            name = dataset['dataset']
        except KeyError:
            raise ConfigError("'dataset' must be a mapping with "
                              "'name' and 'sysnum'")


        sysnum = dataset.get('sys')
        cfac = dataset.get('cfac', tuple())
        return DataSetInput(name=name, sys=sysnum, cfac=cfac)


    def produce_dataset(self, *, dataset_input ,theoryid, use_cuts, fit=None,
                      check_plotting:bool=False):
        """Dataset specification from the theory and CommonData.
           Use the cuts from the fit, if provided. If check_plotting is set to
           True, attempt to lod and check the PLOTTING files
           (note this may cause a noticeable slowdown in general)."""
        name = dataset_input.name
        sysnum = dataset_input.sys
        cfac = dataset_input.cfac



        try:
            ds =  self.loader.check_dataset(name=name, sysnum=sysnum,
                                             theoryid=theoryid, cfac=cfac,
                                             use_cuts=use_cuts, fit=fit)
        except DataNotFoundError as e:
            raise ConfigError(str(e), name, self.loader.available_datasets)

        except LoadFailedError as e:
            raise ConfigError(e)

        if check_plotting:
            from validphys.plotoptions import get_infos
            #normalize=True should check for more stuff
            get_infos(ds, normalize=True)
            if not ds.commondata.plotfiles:
                log.warn("Plotting files not found for: %s" % (ds,))
        return ds


    @configparser.element_of('experiments')
    def parse_experiment(self, experiment:dict, *, theoryid, use_cuts,
                         fit=None, check_plotting:bool=False):
        """A set of datasets where correlated systematics are taken
           into account. It is a mapping where the keys are the experiment
           name 'experiment' and a list of datasets."""
        try:
            name, datasets = experiment['experiment'], experiment['datasets']
        except KeyError as e:
            raise ConfigError("'experiment' must be a mapping with "
                              "'name' and 'datasets', but %s is missing" % e)

        dsinputs = [self.parse_dataset_input(ds) for ds in datasets]
        #autogenerated func, from elemet_of
        datasets = [self.produce_dataset(dataset_input=dsinp, theoryid=theoryid,
                                       use_cuts=use_cuts, fit=fit,
                                       check_plotting=check_plotting)
                                       for dsinp in dsinputs]

        return ExperimentSpec(name=name, datasets=datasets, dsinputs=dsinputs)

    #TODO: Worth it to do some black magic to not pass params explicitly?
    #Note that `parse_experiments` doesn't exist yet.
    def parse_reweighting_experiments(self, experiments, *, theoryid,
                                      use_cuts, fit=None):
        """A list of experiments to be used for reweighting."""
        return self.parse_experiments(experiments,
                                     theoryid=theoryid,
                                     use_cuts=use_cuts, fit=fit)
    def parse_t0pdfset(self, name):
        """PDF set used to generate the t0 covmat."""
        return self.parse_pdf(name)

    def parse_use_t0(self, do_use_t0:bool, t0pdfset=None):
        """Whether to use the t0 PDF set to generate covariance matrices."""
        if do_use_t0 and not t0pdfset:
            raise ConfigError("Setting use_t0 requires specifying a valid t0pdfset")

        return do_use_t0

    #TODO: Find a good name for this
    def produce_t0set(use_t0=False, t0pdfset=None):
        """Return the t0set if use_t0 is True and None otherwise."""
        if use_t0:
            return t0pdfset
        else:
            return None

    @element_of('posdatasets')
    def parse_posdataset(self, posset:dict, * ,theoryid):
        """An observable used as positivity constrain in the fit.
        It is a mapping containing 'dataset' and 'poslambda'."""
        bad_msg = ("posset must be a mapping with a name ('dataset') and "
                   "a float multiplier(poslambda)")

        theoryno, theopath = theoryid
        try:
            name = posset['dataset']
            poslambda = float(posset['poslambda'])
        except KeyError as e:
            raise ConfigError(bad_msg, e.args[0], posset.keys()) from e
        except ValueError as e:
            raise ConfigError(bad_msg) from e

        try:
            return self.loader.check_posset(theoryno, name, poslambda)
        except FileNotFoundError as e:
            raise ConfigError(e) from e
