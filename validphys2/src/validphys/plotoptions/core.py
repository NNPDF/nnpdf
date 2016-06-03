# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 19:27:44 2016

@author: Zahari Kassabov
"""
import logging

import numpy as np
import pandas as pd


from reportengine.configparser import Config, ConfigError, named_element_of
from reportengine.utils import get_functions, ChainMap
from NNPDF import CommonData

from validphys.plotoptions.utils import apply_to_all_columns
from validphys.plotoptions import labelers, kintransforms, resulttransforms

log = logging.getLogger(__name__)

default_labels = ('idat', 'k1', 'k2', 'k3')

labeler_functions = get_functions(labelers)
transform_functions = get_functions(kintransforms)
result_functions = get_functions(resulttransforms)


kinlabels_latex = CommonData.kinLabel_latex.asdict()

def get_plot_kinlabels(commondata):
    """Return the LaTex kinematic labels for a given Commondata"""
    #Since there is no 1:1 correspondence between latex keys and GetProc,
    #we match the first key such that the proc label starts with it.

    l = commondata.GetProc(0)
    try:
        key = next(k for k in kinlabels_latex if l.startswith(k))
    except StopIteration:
        raise ValueError("Could not find a set of kinematic "
                         "variables matching  the process %s Check the "
                         "labels defined in commondata.cc. " % (l))
    return kinlabels_latex[key]

def get_info(commondata, file=None, cuts=None):
    try:
        return PlotInfo.from_commondata(commondata, file=file, cuts=cuts)
    except ConfigError:
        log.error("Problem processing file %s" % getattr(file, 'name', file))
        raise

class PlotInfo:
    def __init__(self, kinlabels, x=None ,extra_labels=None, func_labels=None,
                 figure_by=None, line_by=None, kinematics_override=None,
                 result_transform=None, y_label=None, x_label=None,
                 x_scale=None, y_scale=None):
        self.kinlabels = kinlabels
        if x is None:
            x = 'idat'
        self.x = x
        self.extra_labels = extra_labels
        self.func_labels = func_labels
        self.figure_by = figure_by
        self.line_by = line_by
        self.kinematics_override = kinematics_override
        self.result_transform = result_transform
        self._x_label = x_label
        self.y_label = y_label
        self.x_scale = x_scale
        self.y_scale = y_scale

    def name_to_label(self, name):
        if name in labeler_functions:
            func = labeler_functions[name]
            return getattr(func, 'label', name)
        try:
            ix = ('k1', 'k2', 'k3').index(name)
        except ValueError:
            return name
        return self.kinlabels[ix]

    @property
    def xlabel(self):
        if self._x_label is not None:
            return self._x_label
        return self.name_to_label(self.x)



    def group_label(self, same_vals, groupby):
        if not groupby:
            return ''
        pieces = []
        for column, val in zip(groupby, same_vals):
            label = self.name_to_label(column)
            pieces.append('%s = %s' % (label, val))
        return '%s' % ' '.join(pieces)



    @classmethod
    def from_commondata(cls, commondata, file=None, cuts=None):

        if file:
            config = PlotConfigParser.from_yaml(file, commondata, cuts=cuts)
            plot_params = config.process_all_params()
        else:
            plot_params = {}

        if 'kinematics_override' in plot_params:
            kinlabels = plot_params['kinematics_override'].new_labels
        else:
            kinlabels = get_plot_kinlabels(commondata)
        return cls(kinlabels=kinlabels, **plot_params)


class PlotConfigParser(Config):

    def __init__(self, input_params ,commondata, cuts=None, **kwargs):
        self.commondata = commondata
        self.cuts = cuts
        super().__init__(input_params, **kwargs)

    def process_all_params(self, input_params=None):
        self._output_params = ChainMap({'func_labels':{}})
        return super().process_all_params(input_params=input_params,
                                              ns=self._output_params)


    @named_element_of('extra_labels')
    def parse_label(self, elems:list):
        if self.cuts is not None:
            elems = [elems[c] for c in self.cuts]
        if len(elems) != len(self.commondata):
            raise ConfigError("The number of elements in %s (%d) must be the same as "
                              "the number of points in the CommonData (%d)" %
                              (elems, len(elems), len(self.commondata)))
        return elems

    def resolve_name(self, val, extra_labels):
        if extra_labels is None:
            all_labels = list(default_labels)
        else:
            all_labels = list(extra_labels.keys()) + list(default_labels)
        if val in all_labels:
            return val
        if val in labeler_functions:
            self._output_params['func_labels'][val] = labeler_functions[val]
            return val

        raise ConfigError("Unknown label %s" % val, val, all_labels +
                          list(labeler_functions),
                              display_alternatives='all')

    def parse_x(self, x:str, extra_labels=None):
        return self.resolve_name(x, extra_labels)

    def parse_figure_by(self, gb:list, extra_labels=None):
        return [self.resolve_name(val, extra_labels) for val in gb]

    def parse_line_by(self, lb:list, extra_labels=None):
        return self.parse_figure_by(lb, extra_labels)

    def parse_kinematics_override(self, tr:str):
        if not tr in transform_functions:
            raise ConfigError("Unknown transform function '%s'" % tr, tr,
                              transform_functions)
        return transform_functions[tr]

    def parse_result_transform(self, tr:str):
        if not tr in result_functions:
            raise ConfigError("Unknown transform function '%s'" % tr, tr,
                              result_functions)
        return result_functions[tr]

    def parse_y_label(self, label:str):
        return label

    def parse_x_label(self, label:str):
        return label

    def _parse_scale(self, scale:str):
        if not (scale == 'linear' or scale=='log'):
            raise ConfigError("Scale must be 'linear' or 'log'")
        return scale

    def parse_x_scale(self, scale:str):
        return self._parse_scale(scale)

    def parse_y_scale(self, scale:str):
        return self._parse_scale(scale)


def kitable(commondata, info):
    def _expand(f,x):
        return f(**x)
    table = pd.DataFrame(commondata.get_kintable(), columns=default_labels[1:])
    table.index.name = default_labels[0]
    if info.kinematics_override:
        transform = apply_to_all_columns(table, info.kinematics_override)
        table = pd.DataFrame(np.array(transform).T, columns=table.columns)

    #TODO: This is a little bit ugly. We want to call the functions
    #with all the
    #extra labels
    if info.extra_labels:
        vals = tuple(info.extra_labels.items())
    else:
        vals = ()

    if info.func_labels:
        funcs = tuple(info.func_labels.items())
    else:
        funcs = ()

    for label, value in vals:
        table[label] = value

    nreal_labels = len(table.columns)

    for label, func in funcs:
        #Pass only the "real" labels and not the derived functions
        table[label] = apply_to_all_columns(table.iloc[:,:nreal_labels], func)

    return table

def transform_result(cv, error, kintable, info):
    if not info.result_transform:
        return cv, error
    f = info.result_transform

    df = pd.DataFrame({'cv':cv, 'error':error})
    newcv, newerror = apply_to_all_columns(pd.concat([df,kintable], axis=1),f)

    return np.array(newcv), np.array(newerror)
