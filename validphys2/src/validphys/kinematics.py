"""
Provides information on the kinematics involved in the data.

Uses the PLOTTING file specification.
"""
from collections import namedtuple
import logging

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.checks import check_positive
from reportengine.table import table

from validphys.core import CutsPolicy
from validphys.plotoptions import core as plotoptions_core

from validphys import plotutils
from validphys import plotoptions
from validphys.core import CutsPolicy
from validphys.closuretest import fits_normed_dataset_central_delta
from validphys.closuretest import dataset_xi
from validphys.closuretest import dataset_replica_and_central_diff

from reportengine.figure import figuregen


log = logging.getLogger(__name__)


@check_positive('titlelevel')
def describe_kinematics(commondata, titlelevel: int = 1):
    """Output a markdown text describing the stored metadata for a given
    commondata.

    titlelevel can be used to control the header level of the title.
    """
    import inspect

    cd = commondata
    info = plotoptions_core.get_info(cd)
    proc = cd.load_commondata().commondataproc
    src = inspect.getsource(info.kinematics_override.xq2map)
    titlespec = '#' * titlelevel
    return f"""
{titlespec} {cd}

{info.dataset_label}

Stored data:

 - Process type: **{proc}** ({info.process_description})

 - variables:
     * k1: {info.kinlabels[0]}
     * k2: {info.kinlabels[1]}
     * k3: {info.kinlabels[2]}



Map:

```python
{src}
```

"""


describe_kinematics.highlight = 'markdown'


nfittedlabel = '$N_{fitted}$'
ndatalabel = '$N_{data}$'


def kinlimits(commondata, cuts, use_cuts, use_kinoverride: bool = True):
    """Return a mapping containing the number of fitted and used datapoints, as
    well as the label, minimum and maximum value for each of the three
    kinematics. If ``use_kinoverride`` is set to False, the PLOTTING files will
    be ignored and the kinematics will be interpred based on the process type
    only. If use_cuts is 'CutsPolicy.NOCUTS', the information on the total
    number of points will be displayed, instead of the fitted ones."""
    info = plotoptions_core.get_info(commondata, cuts=None, use_plotfiles=use_kinoverride)

    kintable = plotoptions_core.kitable(commondata, info)
    ndata = len(kintable)
    if cuts:
        kintable = kintable.loc[cuts.load()]
        nfitted = len(kintable)
    elif use_cuts is not CutsPolicy.NOCUTS:
        nfitted = len(kintable)
    else:
        nfitted = '-'

    d = {'dataset': commondata, ndatalabel: ndata, nfittedlabel: nfitted}
    for i, key in enumerate(['k1', 'k2', 'k3']):
        kmin = kintable[key].min()
        kmax = kintable[key].max()
        label = info.kinlabels[i]
        d[key] = label
        d[key + ' min'] = kmin
        d[key + ' max'] = kmax
    return d


all_kinlimits = collect(kinlimits, ('dataset_inputs',))

def xq2_dataset_map(commondata, cuts,internal_multiclosure_dataset_loader,
                        _internal_max_reps=None,
                        _internal_min_reps=20):
    """
    Load in a dictionary all the specs of a dataset meaning:
    - ds name
    - ds coords
    - standard deviation (in multiclosure)
    - mean (in multiclosure again)
    - (x,Q^2) coords
    """
    xq2_map_obj = xq2map_with_cuts(commondata, cuts)
    coords = xq2_map_obj[2]
    central_deltas = fits_normed_dataset_central_delta(internal_multiclosure_dataset_loader)
    std_devs = np.std(central_deltas, axis = 0)
    means = np.mean(central_deltas, axis = 0)
    xi = dataset_xi(dataset_replica_and_central_diff(internal_multiclosure_dataset_loader,False))

    # for case of DY observables we have 2 (x,Q) for each experimental point
    if coords[0].shape[0] != std_devs.shape[0]:
        std_devs = np.concatenate((std_devs,std_devs))
        means = np.concatenate((means,means))
        xi = np.concatenate((xi,xi))
    return {'x_coords':coords[0], 'Q_coords':coords[1],'std_devs':std_devs,'name':commondata.name,'process':commondata.process_type, 'means': means, 'xi': xi}

@table
def all_kinlimits_table(all_kinlimits, use_kinoverride: bool = True):
    """Return a table with the kinematic limits for the datasets given as input
    in dataset_inputs. If the PLOTTING overrides are not used, the information on
    sqrt(k2) will be displayed."""

    table = pd.DataFrame(
        all_kinlimits,
        columns=[
            'dataset',
            '$N_{data}$',
            '$N_{fitted}$',
            'k1',
            'k1 min',
            'k1 max',
            'k2',
            'k2 min',
            'k2 max',
            'k3',
            'k3 min',
            'k3 max',
        ],
    )

    # We really want to see the square root of the scale
    if not use_kinoverride:
        table['k2'] = 'sqrt(' + table['k2'] + ')'
        table['k2 min'] = np.sqrt(table['k2 min'])
        table['k2 max'] = np.sqrt(table['k2 max'])
        # renaming the columns is overly complicated
        cols = list(table.columns)
        cols[6:9] = ['sqrt(k2)', 'sqrt(k2) min', 'sqrt(k2) max']
        table.columns = cols

    return table


@table
def all_commondata_grouping(all_commondata, metadata_group):
    """Return a table with the grouping specified
    by `metadata_group` key for each dataset for all available commondata.
    """
    records = []
    for cd in all_commondata:
        records.append(
            {
                'dataset': str(cd),
                metadata_group: getattr(plotoptions_core.get_info(cd), metadata_group),
            }
        )
    df = pd.DataFrame.from_records(records, index='dataset')
    # sort first by grouping alphabetically and then dataset name
    return df.sort_values([metadata_group, 'dataset'])


def total_fitted_points(all_kinlimits_table) -> int:
    """Print the total number of fitted points in a given set of data"""
    tb = all_kinlimits_table
    return int(tb[nfittedlabel].sum())


XQ2Map = namedtuple('XQ2Map', ('experiment', 'commondata', 'fitted', 'masked', "group"))


def xq2map_with_cuts(commondata, cuts, group_name=None):
    """Return two (x,Q²) tuples: one for the fitted data and one for the
    cut data. If `display_cuts` is false or all data passes the cuts, the second
    tuple will be empty."""
    info = plotoptions_core.get_info(commondata)
    kintable = plotoptions_core.kitable(commondata, info)
    if cuts:
        mask = cuts.load()
        boolmask = np.zeros(len(kintable), dtype=bool)
        boolmask[mask] = True
        fitted_kintable = kintable.loc[boolmask]
        masked_kitable = kintable.loc[~boolmask]
        xq2fitted = plotoptions_core.get_xq2map(fitted_kintable, info)
        xq2masked = plotoptions_core.get_xq2map(masked_kitable, info)
    else:
        xq2fitted = plotoptions_core.get_xq2map(kintable, info)
        xq2masked = (np.array([]), np.array([]))

    return XQ2Map(info.experiment, commondata, xq2fitted, xq2masked, group_name)


dataset_inputs_by_groups_xq2map = collect(
    xq2map_with_cuts, ('group_dataset_inputs_by_metadata', 'data_input')
)
xq2_data_map = collect("xq2_dataset_map",("data",))

@figuregen
def xq2_data_prcs_maps(xq2_data_map):
    import matplotlib.pyplot as plt
    import matplotlib.colors
    from matplotlib.colors import ListedColormap
    keys = ["std_devs","xi"]
    for elem in xq2_data_map:
        for k in keys:
            fig, ax = plotutils.subplots()
            im = ax.scatter(elem['x_coords'],elem['Q_coords'],
                                c=(np.asarray(elem[k])), 
                                cmap='viridis',
                                label = elem['name'])
            fig.colorbar(im,label=k)
            ax.set_xscale('log')  # Set x-axis to log scale
            ax.set_yscale('log')  # Set y-axis to log scale
            ax.set_xlabel('x')
            ax.set_ylabel('Q2')
            ax.set_title(elem["name"]+" "+k)
            yield fig


def kinematics_table_notable(commondata, cuts, show_extra_labels: bool = False):
    """
    Table containing the kinematics of a commondata object,
    indexed by their datapoint id. The kinematics will be tranfsormed as per the
    PLOTTING file of the dataset or process type, and the column headers will
    be the labels of the variables defined in the metadata.

    If ``show_extra_labels`` is ``True`` then extra label defined in the
    PLOTTING files will be displayed. Otherwise only the original three
    kinematics will be shown.
    """
    info = plotoptions_core.get_info(commondata, cuts=cuts)
    res = plotoptions_core.kitable(commondata, info, cuts=cuts)
    res.columns = [*info.kinlabels, *res.columns[3:]]
    if not show_extra_labels:
        res = res.iloc[:, :3]
    return res


@table
def kinematics_table(kinematics_table_notable):
    """Same as kinematics_table_notable but writing the table to file"""
    return kinematics_table_notable
