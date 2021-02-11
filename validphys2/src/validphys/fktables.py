"""
Actions for investigating properties of fktables. Module exposes the utils for
loading fktables found in :py:mod:`validphys.fkparser` as providers and contains
actions which produce plots/tables for use in reports.

"""
import numpy as np
from matplotlib import transforms
import matplotlib.collections as mcollections
from matplotlib.lines import Line2D

from reportengine import collect
from reportengine.figure import figuregen

from validphys.checks import check_pdf_normalize_to, check_scale, check_pdfs_noband
from validphys.fkparser import load_fktable
from validphys.pdfplots import BandPDFPlotter

def dataset_input_fktables(dataset):
    """Loads the fktables associated with a dataset."""
    cuts = dataset.cuts.load() if dataset.cuts is not None else None
    all_tables = [
        load_fktable(fk).with_cuts(cuts) for fk in dataset.fkspecs
    ]
    return all_tables

def dataset_input_fktable_xgrid(dataset_input_fktables):
    """Return the xgrid of a ``dataset`` s fktable."""
    all_xgrids = [fktab.xgrid for fktab in dataset_input_fktables]
    return np.concatenate(all_xgrids)

dataset_inputs_fktable_xgrid = collect(
    "dataset_input_fktable_xgrid", ("data_input",))

def total_data_fktable_xgrid(dataset_inputs_fktable_xgrid):
    """Concatenate the fktable xgrids from each dataset into a single array."""
    return np.concatenate(dataset_inputs_fktable_xgrid)

class BandPDFWithFKXPlotter(BandPDFPlotter):
    """Overload the __call__ method of BandPDFPlotter to add rug plot of
    fk xgrid to each figure. The y position and height of the rugplot is
    based on the y-limits set when plotting the PDF bands.

    """
    def __init__(self, data_fktable_xgrid: np.array, *args, **kwargs):
        self.data_fktable_xgrid = data_fktable_xgrid
        super().__init__(*args, **kwargs)

    def __call__(self,):
        for fig, partonname in super().__call__():
            ax = fig.gca()
            segment_data = self.data_fktable_xgrid
            segments = np.c_[
                segment_data,
                np.zeros_like(segment_data),
                segment_data,
                np.full_like(segment_data, 0.1),
            ].reshape(-1, 2, 2)
            rugs = mcollections.LineCollection(
                segments,
                # Make the x coordinate refer to the data but the y (the height
                # relative to the plot height.
                # https://matplotlib.org/tutorials/advanced/transforms_tutorial.html?highlight=transforms%20blended_transform_factory#blended-transformations
                transform=transforms.blended_transform_factory(ax.transData, ax.transAxes),
                colors="k",
            )
            # make proxy for legend
            proxy = Line2D([0], [0], color="k")
            ax.add_collection(rugs)
            # add label for rug plot
            handles, labels = ax.get_legend_handles_labels()
            handles.append(proxy)
            labels.append("FKTable xgrid points")
            ax.legend(handles, labels)
            yield fig, partonname

fk_xplotting_grids = collect(
    "xplotting_grid", ('pdfs', "fk_basis_flavour"))

@figuregen
@check_pdf_normalize_to
@check_pdfs_noband
@check_scale("xscale", allow_none=True)
def plot_pdfs_fktable_xgrids(
    total_data_fktable_xgrid,
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
        total_data_fktable_xgrid,
        pdfs,
        fk_xplotting_grids,
        xscale,
        normalize_to,
        ymin,
        ymax,
        pdfs_noband=pdfs_noband,
        show_mc_errors=show_mc_errors,
    )
