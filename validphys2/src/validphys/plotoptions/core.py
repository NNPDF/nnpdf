# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 19:27:44 2016

@author: Zahari Kassabov
"""
import dataclasses
import enum
import logging
import numbers
import typing

import numpy as np
import pandas as pd

from reportengine.floatformatting import format_number
from reportengine.utils import ChainMap
from validphys.core import CommonDataSpec, DataSetSpec
from validphys.coredata import CommonData
from validphys.plotoptions.plottingoptions import PlottingOptions, default_labels, labeler_functions
from validphys.plotoptions.utils import apply_to_all_columns
from validphys.utils import parse_yaml_inp

log = logging.getLogger(__name__)


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
    elif hasattr(cuts, 'load'):
        cuts = cuts.load()

    if cuts is not None and not len(cuts):
        raise NotImplementedError("No point passes the cuts. Cannot retieve info")

    if isinstance(data, DataSetSpec):
        data = data.commondata
    if not isinstance(data, CommonDataSpec):
        raise TypeError("Unrecognized data type: %s" % type(data))

    info = PlotInfo.from_commondata(data, cuts=cuts, normalize=normalize)
    return info


class PlotInfo:
    def __init__(
        self,
        kinlabels,
        dataset_label,
        *,
        experiment,
        x=None,
        extra_labels=None,
        func_labels=None,
        figure_by=None,
        line_by=None,
        kinematics_override=None,
        result_transform=None,
        y_label=None,
        x_label=None,
        x_scale=None,
        y_scale=None,
        ds_metadata=None,
        process_description='-',
        nnpdf31_process,
        **kwargs,
    ):
        self.kinlabels = kinlabels
        self.experiment = experiment
        self.nnpdf31_process = nnpdf31_process
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
        self.ds_metadata = ds_metadata  # For new commondata format

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
            if self.ds_metadata is not None and column in ('k1', 'k2', 'k3'):
                # new-style commondata, we can have a nicer label!
                ix = ('k1', 'k2', 'k3').index(column)
                var_key = self.ds_metadata.kinematic_coverage[ix]
                pieces.append(self.ds_metadata.kinematics.apply_label(var_key, val))
            else:
                label = self.name_to_label(column)
                if isinstance(val, numbers.Real):
                    val = format_number(val)
                pieces.append('%s = %s' % (label, val))

        return '%s' % ' '.join(pieces)

    @classmethod
    def from_commondata(cls, commondata, cuts=None, normalize=False):
        plot_params = ChainMap()
        kinlabels = commondata.plot_kinlabels

        if commondata.legacy:
            if commondata.plotfiles:
                for file in commondata.plotfiles:
                    pf = parse_yaml_inp(file, PlottingFile)
                    config_params = dataclasses.asdict(pf, dict_factory=dict_factory)
                    plot_params = plot_params.new_child(config_params)
                if normalize and 'normalize' in plot_params:
                    plot_params = plot_params.new_child(config_params['normalize'])
                if 'dataset_label' not in plot_params:
                    log.warning(f"'dataset_label' key not found in {file}")
                    plot_params['dataset_label'] = commondata.name

            else:
                plot_params = {'dataset_label': commondata.name}

        else:
            pcd = commondata.metadata.plotting_options
            config_params = dataclasses.asdict(pcd, dict_factory=dict_factory)
            plot_params = plot_params.new_child(config_params)
            # Add a reference to the metadata to the plot_params so that it is stored in PlotInfo
            plot_params["ds_metadata"] = commondata.metadata
            # If normalize, we need to update some of the parameters
            if normalize and pcd.normalize is not None:
                plot_params = plot_params.new_child(pcd.normalize)

        kinlabels = plot_params['kinematics_override'].new_labels(*kinlabels)

        if "extra_labels" in plot_params and cuts is not None:
            cut_extra_labels = {
                k: [v[i] for i in cuts] for k, v in plot_params["extra_labels"].items()
            }
            plot_params["extra_labels"] = cut_extra_labels

        return cls(kinlabels=kinlabels, **plot_params)


def dict_factory(key_value_pairs):
    """A dictionary factory to be used in conjunction with dataclasses.asdict
    to remove nested None values and convert enums to their name.

    https://stackoverflow.com/questions/59481989/dict-from-nested-dataclasses
    """
    new_kv_pairs = []
    for key, value in key_value_pairs:
        if isinstance(value, enum.Enum):
            new_kv_pairs.append((key, value.name))
        elif value is not None:
            new_kv_pairs.append((key, value))

    return dict(new_kv_pairs)


class KinLabel(enum.Enum):
    k1 = enum.auto()
    k2 = enum.auto()
    k3 = enum.auto()


@dataclasses.dataclass
class PlottingFile(PlottingOptions):
    normalize: typing.Optional[PlottingOptions] = None


def kitable(data, info, *, cuts=None):
    """Obtain a DataFrame with the kinematics for each data point

    Parameters
    ----------
    data: (DataSetSpec, CommonDataSpec, Dataset, CommonData)
        A data object to extract the kinematics from.
    info: PlotInfo
        The description of the transformations to apply to the kinematics.
        See :py:func:`get_info`
    cuts: Cuts or None, default=None
        An object to load the cuts from. It is an error to set this if ``data``
        is a dataset. If `data` is a CommonData, these **must** be the same as
        those passed to :py:func:`get_info`.

    Returns
    -------
    table: pd.DataFrame
       A DataFrame containing the kinematics for all points after cuts.
    """
    if isinstance(data, (DataSetSpec)) and cuts is not None:
        raise TypeError("Cuts must be None when a dataset is given")

    if isinstance(data, DataSetSpec):
        data = data.load_commondata()
    elif isinstance(data, CommonDataSpec):
        data = data.load()

    table = pd.DataFrame(data.get_kintable(), columns=default_labels[1:])
    if isinstance(data, CommonData) and cuts is not None:
        table = table.loc[cuts.load()]
    table.index.name = default_labels[0]
    if info.kinematics_override:
        transform = apply_to_all_columns(table, info.kinematics_override)
        table = pd.DataFrame(np.array(transform).T, columns=table.columns, index=table.index)

    # TODO: This is a little bit ugly. We want to call the functions
    # with all the
    # extra labels
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
        # Pass only the "real" labels and not the derived functions
        table[label] = apply_to_all_columns(table.iloc[:, :nreal_labels], func)

    return table


def transform_result(cv, error, kintable, info):
    if not info.result_transform:
        return cv, error
    f = info.result_transform

    df = pd.DataFrame({'cv': cv, 'error': error})
    newcv, newerror = apply_to_all_columns(pd.concat([df, kintable], axis=1), f)

    return np.array(newcv), np.array(newerror)


def get_xq2map(kintable, info):
    """Return a tuple of (x,Q²) from the kinematic values defined in kitable
    (usually obtained by calling ``kitable``) using machinery specified in
    ``info``"""
    return apply_to_all_columns(kintable, info.kinematics_override.xq2map)
