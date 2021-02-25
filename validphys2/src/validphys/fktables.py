"""
Actions for investigating properties of fktables. Module exposes the utils for
loading fktables found in :py:mod:`validphys.fkparser` as providers and contains
actions which produce plots/tables for use in reports.

"""
import numpy as np
from matplotlib import transforms
import matplotlib.collections as mcollections
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd

from reportengine import collect
from reportengine.figure import figuregen

from validphys.convolution import (
    central_fk_differential_predictions,
    load_fktable
)
from validphys.checks import (
    check_pdf_normalize_to,
    check_scale,
    check_pdfs_noband,
    make_argcheck,
    CheckError
)
from validphys.pdfplots import BandPDFPlotter
from validphys.plotoptions import get_info


def check_op_is_null(dataset):
    if dataset.op != "NULL":
        raise CheckError("only datasets with no operation are supported right now.")

def check_dataset_is_dis(dataset):
    for loaded_fk in map(load_fktable, dataset.fkspecs):
        if loaded_fk.hadronic:
            raise CheckError("Only DIS datasets are supported right now.")


def normalised_averaged_differential_prediction(dataset, pdf, kinematics_table_notable):
    """Return the xgrid of a ``dataset`` s fktable."""
    info = get_info(dataset)
    cuts = dataset.cuts.load() if dataset.cuts is not None else None
    dfs = []
    k1list = info.get_xcol(kinematics_table_notable)
    for fkspec in dataset.fkspecs:
        fk_df = central_fk_differential_predictions(
            load_fktable(fkspec).with_cuts(cuts), pdf
        )
        for datapoint, datapoint_df in fk_df.groupby(level=0):
            k1 = k1list[datapoint]
            new_index = pd.MultiIndex.from_product(
                [[k1], datapoint_df.index.get_level_values(1)])
            out = datapoint_df.set_index(new_index).fillna(0)
            # divide by total contribution to prediction (sum over flavour and x)
            out /= out.to_numpy().sum()
            # Just get absolute contribution
            out = out.abs()
            dfs.append(out)
    return dfs


class BandPDFWithFKXPlotter(BandPDFPlotter):
    """Overload the __call__ method of BandPDFPlotter to add rug plot of
    fk xgrid to each figure. The y position and height of the rugplot is
    based on the y-limits set when plotting the PDF bands.

    """
    def __init__(self, mean_diff_pred, plotinfo, *args, **kwargs):
        self.mean_diff_pred = mean_diff_pred
        self.plotinfo = plotinfo
        super().__init__(*args, **kwargs)

    def __call__(self,):
        for fig, partonname in super().__call__():
            ax = fig.gca()
            title = (
                ax.get_title() +
                f" with partial predictions from {self.plotinfo.dataset_label}"
            )
            ax.set_title(title)
            total_df = self.mean_diff_pred
            kinvals = [df.index.get_level_values(0).unique() for df in total_df]
            vmin, vmax = min(kinvals), max(kinvals)
            if self.plotinfo.x_scale == 'log':
                norm = mcolors.LogNorm(vmin, vmax)
            else:
                norm = mcolors.Normalize(vmin, vmax)
            sm = plt.cm.ScalarMappable(
                cmap=plt.cm.viridis,
                norm=norm
            )
            colors = sm.to_rgba(kinvals)
            for df, color in zip(total_df, colors):
                df = df.droplevel(0, axis=0).T
                # might need to get name from basis here.
                segment_data_series = df.get(partonname)
                if segment_data_series is not None:
                    xgrid = segment_data_series.index.values
                    segments = np.c_[
                        xgrid,
                        np.zeros_like(xgrid),
                        xgrid,
                        segment_data_series.values,
                    ].reshape(-1, 2, 2)
                    rugs = mcollections.LineCollection(
                        segments,
                        # Make the x coordinate refer to the data but the y (the height
                        # relative to the plot height.
                        # https://matplotlib.org/tutorials/advanced/transforms_tutorial.html?highlight=transforms%20blended_transform_factory#blended-transformations
                        transform=transforms.blended_transform_factory(ax.transData, ax.transAxes),
                        colors=color,
                    )
                    ax.add_collection(rugs)
            plt.colorbar(sm, ax=ax, label=self.plotinfo.xlabel)
            yield fig, partonname

fk_xplotting_grids = collect(
    "xplotting_grid", ('pdfs', "fk_basis_flavour"))

@figuregen
@check_pdf_normalize_to
@check_pdfs_noband
@check_scale("xscale", allow_none=True)
def plot_pdfs_fktable_xgrids(
    normalised_averaged_differential_prediction,
    dataset,
    pdfs,
    fk_xplotting_grids,
    xscale: (str, type(None)) = None,
    normalize_to: (int, str, type(None)) = None,
    ymin=None,
    ymax=None,
    pdfs_noband: (list, type(None)) = None,
    show_mc_errors: bool = True,
):
    """Like :py:func:`plot_pdfs` with additional rug plot for x grid points
    in ``data`` FK tables. Restricted to the evolution basis,
    and plots all flavours used in fk tables.

    """
    plotinfo = get_info(dataset)
    yield from BandPDFWithFKXPlotter(
        normalised_averaged_differential_prediction,
        plotinfo,
        pdfs,
        fk_xplotting_grids,
        xscale,
        normalize_to,
        ymin,
        ymax,
        pdfs_noband=pdfs_noband,
        show_mc_errors=show_mc_errors,
    )
