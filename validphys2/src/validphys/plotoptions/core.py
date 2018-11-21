# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 19:27:44 2016

@author: Zahari Kassabov
"""
import logging

import numpy as np
import pandas as pd
import numbers

from reportengine.floatformatting import format_number
from reportengine.configparser import Config, ConfigError, named_element_of
from reportengine.utils import get_functions, ChainMap

from validphys.core import CommonDataSpec, DataSetSpec, Cuts
from validphys.plotoptions.utils import apply_to_all_columns, get_subclasses
from validphys.plotoptions import labelers, kintransforms, resulttransforms

log = logging.getLogger(__name__)

default_labels = ('idat', 'k1', 'k2', 'k3')

labeler_functions = get_functions(labelers)
transform_functions = get_subclasses(kintransforms, kintransforms.Kintransform)
result_functions = get_functions(resulttransforms)


def get_info(data, *, normalize=False, cuts=None, use_plotfiles=True):
    """Retrieve and process the plotting information for the input data (which could
    be a DatasetSpec or a CommonDataSpec).

    If ``use_plotfiles`` is ``True`` (the default), the PLOTTING files will be
    used to retrieve the infromation. Otherwise the default configuration
    (which depends of the process type) will be used.

    If cuts is None, the cuts of the dataset will be used, but no cuts for
    commondata.

    If cuts is False, no cuts will be used.

    If cuts is an instance of Cuts, it will be used.

    If normalize is True, the specialization for ratio plots will be used to
    generate the PlotInfo objects.
    """
    if cuts is None:
        if isinstance(data, DataSetSpec):
            cuts = data.cuts.load() if data.cuts else None
    elif isinstance(cuts, Cuts):
        cuts = cuts.load()
    elif not cuts:
        cuts = None

    if isinstance(data, DataSetSpec):
        data = data.commondata
    if not isinstance(data, CommonDataSpec):
        raise TypeError("Unrecognized data type: %s" % type(data) )

    info = PlotInfo.from_commondata(data, cuts=cuts, normalize=normalize)
    return info



class PlotInfo:
    def __init__(self, kinlabels,dataset_label,* ,x=None ,extra_labels=None, func_labels=None,
                 figure_by=None, line_by=None, kinematics_override=None,
                 result_transform=None, y_label=None, x_label=None,
                 x_scale=None, y_scale=None, process_description='-', **kwargs):
        self.kinlabels = kinlabels
        if x is None:
            x = 'idat'
        self.x = x
        self.extra_labels = extra_labels
        self.func_labels = func_labels
        self.figure_by = figure_by
        self.line_by = line_by
        if kinematics_override is None:
            raise ValueError(f'A kinematics_override must be set for {dataset_label}')
        self.kinematics_override = kinematics_override
        self.result_transform = result_transform
        self._x_label = x_label
        self.y_label = y_label
        self.x_scale = x_scale
        self.y_scale = y_scale
        self.dataset_label = dataset_label
        self.process_description = process_description

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

    def get_xcol(self, table):
        """Return a numpy array with the x column or the index as appropriate"""
        if self.x == 'idat':
            return np.array(table.index)
        else:
            return np.asarray(table[self.x])



    def group_label(self, same_vals, groupby):
        if not groupby:
            return ''
        if len(same_vals) == 1 and isinstance(same_vals[0], str):
            return f'({same_vals[0]})'
        pieces = []
        for column, val in zip(groupby, same_vals):
            label = self.name_to_label(column)
            if isinstance(val, numbers.Real):
                val = format_number(val)
            pieces.append('%s = %s' % (label, val))
        return '%s' % ' '.join(pieces)



    @classmethod
    def from_commondata(cls, commondata, cuts=None, normalize=False):

        #The only reason to call the parser once per config file is to
        #give better error messages and stricter checks
        plot_params = ChainMap()
        if commondata.plotfiles:
            for file in commondata.plotfiles:
                with open(file) as f:
                    config = PlotConfigParser.from_yaml(f, commondata, cuts=cuts)
                try:
                    config_params = config.process_all_params()
                except ConfigError:
                    log.error(f"Error in plotting file: {file}")
                    raise

                plot_params = plot_params.new_child(config_params)
            if normalize and 'normalize' in plot_params:
                #We might need to use reportengine.namespaces.resolve here
                plot_params = plot_params.new_child(config_params['normalize'])
            if not 'dataset_label' in plot_params:
                log.warning(f"'dataset_label' key not found in {file}")
                plot_params['dataset_label'] = commondata.name

        else:
            plot_params = {'dataset_label':commondata.name}

        kinlabels = commondata.plot_kinlabels
        if 'kinematics_override' in plot_params:
            kinlabels = plot_params['kinematics_override'].new_labels(*kinlabels)

        return cls(kinlabels=kinlabels, **plot_params)


class PlotConfigParser(Config):

    allowed_keys = {'normalize':dict}


    #TODO: Remove any reference to cuts from here entirely
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
        ndata = self.commondata.ndata
        if len(elems) != ndata:
            raise ConfigError("The number of elements in %s (%d) must be the same as "
                              "the number of points in the CommonData (%d)" %
                              (elems, len(elems), (ndata)))
        if self.cuts is not None:
            elems = [elems[c] for c in self.cuts]
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

    def parse_process_description(self, desc:str):
        return desc

    def parse_dataset_label(self, lb:str):
        return lb

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
        return transform_functions[tr]()

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
        if not (scale == 'linear' or scale=='log' or scale == 'symlog'):
            raise ConfigError("Scale must be 'linear', 'log' or 'symlog'")
        return scale

    def parse_x_scale(self, scale:str):
        return self._parse_scale(scale)

    def parse_y_scale(self, scale:str):
        return self._parse_scale(scale)

    def parse_theory_reference(self, ref:str):
        return ref

    def parse_data_reference(self, ref:str):
        return ref


def kitable(commondata, info):
    if isinstance(commondata, (DataSetSpec, CommonDataSpec)):
        commondata = commondata.load()
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

def get_xq2map(kintable, info):
    """Return a tuple of (x,Q²) from the kinematic values defined in kitable
    (usually obtained by calling ``kitable``) using machinery specified in
    ``info``"""
    return apply_to_all_columns(kintable, info.kinematics_override.xq2map)
