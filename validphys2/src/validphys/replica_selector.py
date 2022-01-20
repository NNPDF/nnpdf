"""
replica_selector.py

Tools for filtering replica sets based on criteria on the replicas.
"""
import logging
import numbers
import shutil
import warnings
import re

import numpy as np
import pandas as pd

from reportengine.checks import check_positive, make_argcheck, check
from reportengine.compat import yaml
from reportengine.table import table
from reportengine.figure import figuregen

from validphys.pdfbases import flavour
from validphys.pdfoutput import pdfset
from validphys.lhio import new_pdf_from_indexes
from validphys.checks import check_pdf_is_montecarlo, check_scale
from validphys.core import PDF
from validphys.pdfplots import ReplicaPDFPlotter
from validphys.renametools import rename_pdf
from validphys.utils import tempfile_cleaner

log = logging.getLogger(__name__)

def _fixup_new_replica(alphas_pdf: PDF, new_replica_file):
    """Helper function that takes in a
    :py:class:`validphys.core.PDF` object as well as
    the path to the central replica corresponding to the
    PDF and handles the writing of the alphas values
    to the header file.
    """
    AlphaS_MZ = alphas_pdf.AlphaS_MZ
    AlphaS_Vals = alphas_pdf.AlphaS_Vals
    with open(new_replica_file, 'rb') as in_stream:
        data = in_stream.read()
    with open(new_replica_file, 'wb') as out_stream:
        # Add the AlphaS_MZ and AlphaS_Vals keys
        out_stream.write(f"AlphaS_MZ: {AlphaS_MZ}\nAlphaS_Vals: {AlphaS_Vals}\n".encode())
        out_stream.write(data)

@make_argcheck
def _check_target_name(target_name):
    """Make sure this specifies a name and not some kid of path"""
    if target_name is None:
        return
    check(
        re.fullmatch(r'[\w]+', target_name),
        "`target_name` must contain alphnumeric characters and underscores only",
    )

@_check_target_name
def alpha_s_bundle_pdf(pdf, pdfs, output_path, target_name: (str, type(None)) = None):
    """Action that bundles PDFs for distributing to the LHAPDF
    format. The baseline pdf is declared as the ``pdf`` key
    and the PDFs from which the replica 0s are to be added is
    declared as the ``pdfs`` list.

    The bundled PDF set is stored inside the ``output`` directory.

    Parameters
    ----------
    pdf: :py:class:`validphys.core.PDF`
        The baseline PDF to which the new replicas will be added
    pdfs: list of :py:class:`validphys.core.PDF`
        The list of PDFs from which replica0 will be appended
    target_name: str or None
        Optional argument specifying the name of the output PDF.
        If ``None``, then the name of the original pdf is used
        but with ``_pdfas`` appended
    """
    base_pdf_path = pdf.infopath.parent
    nrep = len(pdf)

    target_name = target_name or pdf.name + '_pdfas'
    target_path = output_path / target_name

    alphas_paths = [i.infopath.parent for i in pdfs]
    alphas_replica0s = [path / f'{p}_0000.dat' for path, p in zip(alphas_paths, pdfs)]
    new_nrep = nrep + len(alphas_replica0s)
    alphas_values = [str(p.AlphaS_MZ) for p in pdfs]

    if target_path.exists():
        log.warning(f"{target_path} already exists. Deleting contents.")
        shutil.rmtree(target_path)

    # We create a temporary directory to handle the manipulations inside.
    # We move the files to the new directory at the end.
    with tempfile_cleaner(
        root=output_path, exit_func=shutil.rmtree, exc=KeyboardInterrupt
    ) as tempdir:
        # Copy the base pdf into the temporary directory
        temp_pdf = shutil.copytree(base_pdf_path, tempdir / pdf.name)

        # Copy the alphas PDF replica0s into the new PDF
        for i, (alphas_pdf, rep) in enumerate(zip(pdfs, alphas_replica0s)):
            to = temp_pdf / f'{pdf.name}_{str(i + nrep).zfill(4)}.dat'
            shutil.copy(rep, to)
            _fixup_new_replica(alphas_pdf, to)

        #  Fixup the info file
        info_file = (temp_pdf / temp_pdf.name).with_suffix('.info')

        with open(info_file, 'r') as stream:
            yaml_obj = yaml.YAML()
            info_yaml = yaml_obj.load(stream)
        info_yaml['NumMembers'] = new_nrep
        info_yaml['error_type'] += '+as'
        extra_desc = '; '.join(
            f"mem={i} => alphas(MZ)={val}"
            for val, i in zip(alphas_values, range(nrep, new_nrep))
        )
        info_yaml['SetDesc'] += f"; {extra_desc}"
        with open(info_file, 'w') as stream:
            yaml_obj.dump(info_yaml, stream)

        # Rename the base pdf to the final name
        rename_pdf(temp_pdf, pdf.name, target_name)
        # This is the pdf path after the above renaming
        # i.e new_pdf.exists() == True
        new_pdf = temp_pdf.with_name(target_name)
        # Move the final pdf outside the temporary directory
        new_pdf = new_pdf.rename(target_path)
    log.info(f"alpha_s bundle written at {new_pdf}")
    return target_name



@check_positive('Q')
def gluon_values(pdf, Q, xgrid):
    """Return the x*gluon values of the PDF at Q and for each element in xgrid."""
    scale, x = xgrid
    grid = flavour.grid_values(pdf, ['g'], x, Q)
    #Remove Q and flavour axes
    return grid.reshape((grid.shape[0], grid.shape[2]))


@check_positive('range_percent')
def unconstrained_region_index(
        gluon_values, range_percent: numbers.Real = 68) -> int:
    """Return the first index where the mean of the gluon values is smaller than the size
    of the centered percentile range of size ``range_percent``.
    If no value satisfies the condition, return -1"""
    s = (100 - range_percent)/2
    p = np.percentile(gluon_values, q=[s, 50, 100 - s], axis=0)
    good = p[1] > p[2] - p[0]
    where = np.argwhere(good)
    if not where.size:
        log.warning(
            f"Could not any unconstrained region at "
            f"percentile {range_percent}")
        return -1
    first_good_index = where[0]
    return int(first_good_index)


def unconstrained_limit(xgrid, unconstrained_region_index):
    return xgrid[1][unconstrained_region_index]


@check_positive('nsigma')
def not_outlier_mask(
        xplotting_grid, unconstrained_region_index, nsigma: numbers.Real = 4):
    """Return a boolean mask with the replicas that are never outside of the given
    percentile in the constrained region"""
    lim = unconstrained_region_index
    gv = xplotting_grid.grid_values[:, :, lim:]
    delta = nsigma * np.std(gv, axis=0)
    mean = np.mean(gv, axis=0)
    return ((gv >= mean - delta) & (gv <= mean+delta)).all(axis=2).all(axis=1)


def growing_gluon_mask(gluon_values, unconstrained_region_index):
    """Return a boolean mask marking the replicas where the gluon
    is monothonic in the unconstrained region"""
    lim = unconstrained_region_index
    return (np.diff(gluon_values[:, :lim], axis=1) < 0).all(axis=1)


def _get_exponent_along_xgrid(values, xgrid):
    """Return the small-x exponent computed at each point in x, defined as:
    exp = log(x*f(x))/log(x)
    """
    #PDF can be negative and we just return nan for the corresponding entries
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.log(values)/np.log(xgrid)


def _filter_exponents(values, xgrid, max_allowed_diff):
    """Filter out the values that differ from the best fit of
    x^alpha at any point. This is determined by checking whether the
    log(values)/log(x) differs from the least squares fit to a constant .
    This is the implementation of
    filter_exponents_mask"""
    if len(values.shape) == 2:
        values = values[:, np.newaxis, :]
    exps = _get_exponent_along_xgrid(values, xgrid)
    reshaped = exps.reshape(
        values.shape[0] * values.shape[1], values.shape[2]).T
    #Throw away negative stuff
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        #polyfit throws an MKL error if the inputs are NaNs. While the outputs
        #are also nans as expected, and while probably an overkill,
        #I couldn't find in the documentation that
        #this must be the case, so just set bogus values manually.
        nanmask = np.isnan(reshaped).any(axis=0)
        #bogus values.
        reshaped[:, nanmask] = -1000
        fitted_exps = (np.polyfit(xgrid, reshaped, 0))
        fitted_exps[:, nanmask] = np.nan
        fitted_exps = fitted_exps.reshape(values.shape[0], values.shape[1])
        mask = (
            np.abs(exps - fitted_exps[..., np.newaxis]).max(axis=2) <
            max_allowed_diff).all(axis=1)
    return mask


@check_positive('max_allowed_diff')
def frozen_exponents_mask(
        xplotting_grid,
        unconstrained_region_index,
        max_allowed_diff: numbers.Real = 0.3):
    """Return a boolean mask with the replicas that can be approximated by
    x^ in the uncinstrained region."""
    lim = unconstrained_region_index
    x = xplotting_grid.xgrid[:lim]
    gv = xplotting_grid.grid_values[..., :lim]
    return _filter_exponents(gv, x, max_allowed_diff=max_allowed_diff)


def mcpdf_total_mask(
        not_outlier_mask, growing_gluon_mask, frozen_exponents_mask):
    """Mask of the replicas that fullfill all conditions."""
    return np.all(
        [frozen_exponents_mask, growing_gluon_mask, not_outlier_mask], axis=0)


#TODO: I don't like passing around these things everywhere.
#Find a way to group them.
@table
def mcpdf_stats_table(
        not_outlier_mask, growing_gluon_mask, frozen_exponents_mask,
        mcpdf_total_mask):
    """Show the number of replicas that satisfy each condition, the number
    of replicas that satisfy all the conditions and the starting number of
    replicas."""
    d = {
        "Total replicas": len(not_outlier_mask),
        "Pass frozen exponents": np.sum(frozen_exponents_mask),
        "Pass growing gluon": np.sum(growing_gluon_mask),
        "Pass not outlier": np.sum(not_outlier_mask),
        "Pass all": np.sum(mcpdf_total_mask),
    }
    #This meeses with the index order
    #return pd.Series(d).to_frame("# Replicas")
    s = pd.Series(d, index=d.keys())
    s.name = "# Replicas"
    return s


def discarded_mcpdf_replicas(
        not_outlier_mask, growing_gluon_mask, frozen_exponents_mask,
        mcpdf_total_mask):
    """Return a dictionary of masks where ``Selected`` corresponds to the
    replicas passing all the criteria, and is complementary with the union of
    the rest of the values
    of the entries, which correpond to the replicas that failed each of the
    criteria."""
    return {
        "Discarded: outlier": ~not_outlier_mask,
        "Discarded: exponents not frozen": ~frozen_exponents_mask,
        "Discarded: gluon not monothonic": ~growing_gluon_mask,
        "Selected": mcpdf_total_mask,
    }


@pdfset
@check_pdf_is_montecarlo
def make_montecarlo_pdf(
        pdf,
        output_path,
        set_name: str,
        mcpdf_total_mask,
        installgrid: bool = True,
):
    """Make a "Monte Carlo" PDF set by filtering the replicas from a given
    initial set."""
    indexes = np.arange(1, len(mcpdf_total_mask)+1)[mcpdf_total_mask]
    new_pdf_from_indexes(pdf, indexes, set_name=set_name, folder=output_path)


class FilteredReplicaPlotter(ReplicaPDFPlotter):
    def __init__(self, filter_map, *args, **kwargs):
        self.filter_map = filter_map
        super().__init__(*args, **kwargs)

    def draw(self, pdf, grid, flstate):
        limits = super().draw(pdf, grid, flstate)
        ax = flstate.ax
        for label, filt in self.filter_map.get(pdf, {}).items():
            #TODO: Ideally this woud not modify the color of the subsequent
            #PDFs, and instead would do something different
            next_prop = next(ax._get_lines.prop_cycler)
            color = next_prop['color']
            gv = grid.grid_values[filt, flstate.flindex, :]
            ax.plot(
                grid.xgrid,
                gv.T,
                alpha=1,
                linewidth=0.7,
                color=color,
                zorder=1)
            stats = pdf.stats_class(gv)
            ax.plot(
                grid.xgrid,
                stats.central_value(),
                color=color,
                linewidth=2,
                label=f'{label}')

        return limits


@figuregen
@check_scale('xscale', allow_none=True)
def plot_filtered_replicas(
        pdf,
        xplotting_grid,
        discarded_mcpdf_replicas,
        xscale: (str, type(None)) = None,
        normalized: bool = False,
):
    """
    A replica plot where each replica is marked according to
    ``discarded_mcppdf_replicas``. That is, each replica is marked either as
    selected or as failing one of the criteria.
    """
    if normalized:
        normalize_to = 0
    else:
        normalize_to = None
    yield from FilteredReplicaPlotter(
        filter_map={pdf: discarded_mcpdf_replicas},
        pdfs=[pdf],
        xplotting_grids=[xplotting_grid],
        xscale=xscale,
        normalize_to=normalize_to)
