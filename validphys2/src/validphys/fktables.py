"""
Actions for investigating properties of fktables. Module exposes the utils for
loading fktables found in :py:mod:`validphys.fkparser` as providers and contains
actions which produce plots/tables for use in reports.

"""
import numpy as np
from matplotlib import transforms
import matplotlib.collections as mcollections
from matplotlib.lines import Line2D
import pandas as pd

from reportengine import collect
from reportengine.figure import figuregen

from validphys.convolution import central_differential_predictions
from validphys.checks import (
    check_pdf_normalize_to,
    check_scale,
    check_pdfs_noband,
    make_argcheck,
    CheckError
)
from validphys.fkparser import load_fktable
from validphys.pdfplots import BandPDFPlotter
import validphys.plotutils as plotutils

def check_op_is_null(dataset):
    if dataset.op != "NULL":
        raise CheckError("only datasets with no operation are supported right now.")

def check_dataset_is_dis(dataset):
    for loaded_fk in map(load_fktable, dataset.fkspecs):
        if loaded_fk.hadronic:
            raise CheckError("Only DIS datasets are supported right now.")


# TODO: do work so these checks aren't required.
@make_argcheck(check_op_is_null)
@make_argcheck(check_dataset_is_dis)
def normalised_averaged_differential_prediction(dataset, pdf):
    """Return the xgrid of a ``dataset`` s fktable."""
    df = central_differential_predictions(dataset, pdf).fillna(0)
    unstacked = df.unstack()
    normed = unstacked / unstacked.to_numpy().sum(axis=1, keepdims=True)
    # mean of absolute values across data points and then restack
    norm_abs_mean_vals = normed.abs().to_numpy().mean(axis=0, keepdims=True)
    # dummy index to stack but then instantly drop.
    return pd.DataFrame(
        norm_abs_mean_vals, columns=normed.columns, index=["dummy"]).stack().droplevel(0)


class BandPDFWithFKXPlotter(BandPDFPlotter):
    """Overload the __call__ method of BandPDFPlotter to add rug plot of
    fk xgrid to each figure. The y position and height of the rugplot is
    based on the y-limits set when plotting the PDF bands.

    """
    def __init__(self, mean_diff_pred, *args, **kwargs):
        self.mean_diff_pred = mean_diff_pred
        super().__init__(*args, **kwargs)

    def setup_flavour(self, flstate):
        # make proxy for adding rug plot to legend
        proxy = Line2D([0], [0], color="k", visible="False")
        flstate.handles=[proxy]
        flstate.labels=["FKTable xgrid points"]
        flstate.hatchit=plotutils.hatch_iter()

    def __call__(self,):
        for fig, partonname in super().__call__():
            ax = fig.gca()

            df = self.mean_diff_pred.T
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
                    colors="k",
                )
                ax.add_collection(rugs)
            yield fig, partonname

fk_xplotting_grids = collect(
    "xplotting_grid", ('pdfs', "fk_basis_flavour"))

@figuregen
@check_pdf_normalize_to
@check_pdfs_noband
@check_scale("xscale", allow_none=True)
def plot_pdfs_fktable_xgrids(
    normalised_averaged_differential_prediction,
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
    yield from BandPDFWithFKXPlotter(
        normalised_averaged_differential_prediction,
        pdfs,
        fk_xplotting_grids,
        xscale,
        normalize_to,
        ymin,
        ymax,
        pdfs_noband=pdfs_noband,
        show_mc_errors=show_mc_errors,
    )
