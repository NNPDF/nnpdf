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

from validphys import lhaindex
from validphys.core import PDF, DataSetSpec
from validphys.loader import Loader, DataNotFoundError, SysNotFoundError

log = logging.getLogger(__name__)



class Environment(Environment):
    """Container for information to be filled at run time"""
    def __init__(self,*, data_path, results_path, this_folder, **kwargs):
        self.deta_path = pathlib.Path(data_path)
        self.results_path = pathlib.Path(results_path)
        self.this_folder = pathlib.Path(this_folder)

        self.loader = Loader(data_path, resultspath=self.results_path)
        super().__init__(**kwargs)


def _id_with_label(f):
    f = _parse_func(f)
    def parse_func(self, item, **kwargs):
        if not isinstance(item, dict):
            return f(self, item, **kwargs)
        keydiff =  {'id', 'label'} - item.keys()

        if  keydiff and 'id' in keydiff:
            raise ConfigError("'%s' must be a single id, or a mapping with keys 'id', 'label'"%(item,))
        id = item['id']
        val = f(self, id, **kwargs)
        if 'label' in item:
            val.label = str(item['label'])
        return val

    currsig = inspect.signature(parse_func)
    origsig = inspect.signature(f)
    parse_func = functools.wraps(f)(parse_func)

    params = [*list(currsig.parameters.values())[:2], *list(origsig.parameters.values())[2:]]

    parse_func.__signature__ = inspect.Signature(
                 parameters=params)


    return parse_func

class Config(configparser.Config):


    @property
    def loader(self):
        return self.environment.loader

    @element_of('pdfs')
    @_id_with_label
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

    @element_of('theoryids')
    @_id_with_label
    def parse_theoryid(self, theoryID: (str, int)):
        try:
            return self.loader.check_theoryID(theoryID)
        except FileNotFoundError as e:
            raise ConfigError(str(e), theoryID,
                              self.loader.available_theories,
                              display_alternatives='all')

    def parse_use_cuts(self, use_cuts:bool, *, fit=None):
        if use_cuts and not fit:
            raise ConfigError("Setting 'use_cuts' true requires "
            "specifying a fit on which filter "
            "has been executed, e.g.\nfit : NNPDF30_nlo_as_0118")
        return use_cuts

    #TODO: load fit config from here
    def parse_fit(self, fit:str):
        try:
            return self.loader.check_fit(fit)
        except Exception as e:
            #TODO: maybe drop pathlib because it's too annoying to use

            raise ConfigError(str(e), fit ,self.loader.available_fits)


    @configparser.element_of('datasets')
    def parse_dataset(self, dataset:dict, * ,theoryid, use_cuts, fit=None):
        """Load a dataset specification from the corrsponding theory.
        Use the cuts from the fit, if provided."""
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