# -*- coding: utf-8 -*-
"""
Tools for computing and plotting effective exponents.
"""
from __future__ import generator_stop

import logging
import numbers
import random
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.compat import yaml
from reportengine.figure import figuregen
from reportengine.floatformatting import format_number, significant_digits
from reportengine.table import table

from validphys.checks import check_positive, check_pdf_normalize_to, make_argcheck, check_xlimits
from validphys.core import PDF, FitSpec
from validphys.pdfbases import check_basis, Basis
from validphys.pdfplots import BandPDFPlotter, PDFPlotter

import validphys.pdfgrids as pdfgrids

log = logging.getLogger(__name__)

INTERNAL_LINESTYLE = ['-.', ':']
INTERNAL_COLOR = plt.rcParams['axes.prop_cycle'].by_key()["color"]


@check_positive('Q')
@make_argcheck(check_basis)
@check_xlimits
def alpha_eff(pdf: PDF, *,
              xmin: numbers.Real = 1e-6,
              xmax: numbers.Real = 1e-3,
              npoints: int = 200,
              Q: numbers.Real = 1.65,
              basis: (str, Basis),
              flavours: (list, tuple, type(None)) = None):
    """Return a list of xplotting_grids containing the value of the effective
    exponent alpha at the specified values of x and flavour.
    alpha is relevant at small x, hence the linear scale.

    basis: Is one of the bases defined in pdfbases.py. This includes 'flavour'
    and 'evolution'.

    flavours: A set of elements from the basis.
    If None, the defaults for that basis will be selected.

    Q: The PDF scale in GeV.
    """
    #Loading the filter map of the fit/PDF
    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']

    if npoints == 2:
        xGrid = np.array([xmin, xmax])
    else:
        xGrid = pdfgrids.xgrid(xmin, xmax, 'log', npoints)

    pdfGrid = pdfgrids.xplotting_grid(
        pdf, Q, xgrid=xGrid, basis=basis, flavours=flavours)
    pdfGrid_values = pdfGrid.grid_values
    # NOTE: without this I get "setting an array element with a sequence"
    xGrid = pdfGrid.xgrid
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        alphaGrid_values = -np.log(abs(pdfGrid_values/xGrid))/np.log(xGrid)
        alphaGrid_values[alphaGrid_values == - np.inf] = np.nan  # when PDF_i =0
    alphaGrid = pdfGrid._replace(grid_values=alphaGrid_values)
    return alphaGrid

@check_positive('Q')
@make_argcheck(check_basis)
@check_xlimits
def beta_eff(pdf, *,
             xmin: numbers.Real = 0.6,
             xmax: numbers.Real = 0.9,
             npoints: int = 200,
             Q: numbers.Real = 1.65,
             basis: (str, Basis),
             flavours: (list, tuple, type(None)) = None):
    """Return a list of xplotting_grids containing the value of the effective
    exponent beta at the specified values of x and flavour.
    beta is relevant at large x, hence the linear scale.

    basis: Is one of the bases defined in pdfbases.py. This includes 'flavour'
    and 'evolution'.

    flavours: A set of elements from the basis.
    If None, the defaults for that basis will be selected.

    Q: The PDF scale in GeV.
    """
    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']

    if npoints == 2:
        xGrid = np.array([xmin, xmax])
    else:
        xGrid = pdfgrids.xgrid(xmin, xmax, 'linear', npoints)


    pdfGrid = pdfgrids.xplotting_grid(
        pdf, Q, xgrid=xGrid, basis=basis, flavours=flavours)
    pdfGrid_values = pdfGrid.grid_values
    # NOTE: without this I get "setting an array element with a sequence"
    xGrid = pdfGrid.xgrid
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        betaGrid_values = np.log(abs(pdfGrid_values/xGrid))/np.log(1-xGrid)
        betaGrid_values[betaGrid_values == -np.inf] = np.nan  # when PDF_i =0
    betaGrid = pdfGrid._replace(grid_values=betaGrid_values)

    return betaGrid  # .grid_values

class PreprocessingPlotter(PDFPlotter):
    """ Class inherenting from BandPDFPlotter, changing title and ylabel to reflect the effective
    exponent being plotted.
    """

    def __init__(self, exponent, *args,  **kwargs):
        self.exponent = exponent
        super().__init__(*args, **kwargs)

    def get_title(self, parton_name):
        return fr"$\{self.exponent}_e$ for ${parton_name}$ at {format_number(self.Q, 3)} Gev"

    def get_ylabel(self, parton_name):
        if self.normalize_to is not None:
            return "Ratio to {}".format(self.normalize_pdf.label)
        else:
            return fr"$\{self.exponent}_e$ for ${parton_name}$"

def get_alpha_lines(effective_exponents_table_internal):
    """Given an effective_exponents_table_internal returns the rows with bounds
    of the alpha effective exponent for all flavours, used to plot horizontal
    lines on the alpha effective exponent plots.

    """
    return effective_exponents_table_internal.iloc[0::2, :]

def get_beta_lines(effective_exponents_table_internal):
    """Same as `get_alpha_lines` but for beta"""
    return effective_exponents_table_internal.iloc[1::2, :]

pdfs_alpha_lines = collect('get_alpha_lines', ("pdfs",))
pdfs_beta_lines = collect('get_beta_lines', ("pdfs",))

fits_alpha_lines = collect('get_alpha_lines', ('fits', 'fitpdf'))
fits_beta_lines = collect('get_beta_lines', ('fits', 'fitpdf'))

class ExponentBandPlotter(BandPDFPlotter, PreprocessingPlotter):
    def __init__(self, hlines, exponent, *args,  **kwargs):
        super().__init__(exponent, *args,  **kwargs)
        self.hlines = hlines

    def draw(self, pdf, grid, flstate):
        """Overload :py:meth:`BandPDFPlotter.draw` to plot bands of the
        effective exponent calculated from the replicas and horizontal lines
        for the effective exponents of the previous/next fits, if possible.

        ``flstate`` is an element of the flavours for the first pdf specified in
        pdfs. If this flavour doesn't exist in the current pdf's fitbasis or
        the set of flavours for which the preprocessing exponents exist for the
        current pdf no horizontal lines are plotted.

        """
        pdf_index = self.pdfs.index(pdf)
        hlines = self.hlines[pdf_index]
        # get the correct index label - don't assume table ordering. Basis must
        # be same for all fits so assuming flavour exists in table is valid.
        table_fl_index = f"${grid.basis.elementlabel(flstate.fl)}$"

        errdown, errup = super().draw(pdf, grid, flstate)
        col_label = hlines.columns.get_level_values(0).unique()
        # need to have plotted bands before getting x limit
        xmin, xmax = flstate.ax.get_xlim()

        for i, label in enumerate(col_label):
            # wrap color index since number of pdfs could in theory exceed
            # number of colors
            handle = flstate.ax.hlines(
                hlines.loc[table_fl_index, label].values,
                xmin=xmin,
                xmax=xmax,
                linestyle=INTERNAL_LINESTYLE[i],
                color=INTERNAL_COLOR[pdf_index % len(INTERNAL_COLOR)]
            )
            flstate.handles.append(handle)
            flstate.labels.append(label)
        # need to return xgrid shaped object but with hlines taken into account to get plots nice
        hline_positions = hlines.loc[table_fl_index, :].values.flatten()
        new_errdown = min(
            [*errdown, *hline_positions,])
        new_errup = max(
            [*errup, *hline_positions,])
        return new_errdown*np.ones_like(errdown), new_errup*np.ones_like(errup)


alpha_eff_pdfs = collect('alpha_eff', ('pdfs',))

@figuregen
@check_pdf_normalize_to
def plot_alpha_eff_internal(
        pdfs, alpha_eff_pdfs, pdfs_alpha_lines,
        normalize_to: (int, str, type(None)) = None,
        ybottom=None, ytop=None):
    """Plot the central value and the uncertainty of a list of effective
    exponents as a function of x for a given value of Q. If normalize_to
    is given, plot the ratios to the corresponding alpha effective.
    Otherwise, plot absolute values.
    See the help for ``xplotting_grid`` for information on how to set basis,
    flavours and x ranges. Yields one figure per PDF flavour.

    normalize_to:  Either the name of one of the alpha effective or its
    corresponding index in the list, starting from one, or None to plot
    absolute values.
    """
    yield from ExponentBandPlotter(
        pdfs_alpha_lines, 'alpha', pdfs, alpha_eff_pdfs, 'log', normalize_to, ybottom, ytop)

alpha_eff_fits = collect('alpha_eff', ('fits', 'fitpdf',))

@figuregen
def plot_alpha_eff(
        fits_pdf, alpha_eff_fits, fits_alpha_lines,
        normalize_to: (int, str, type(None)) = None,
        ybottom=None, ytop=None):
    """Plot the central value and the uncertainty of a list of effective
    exponents as a function of x for a given value of Q. If normalize_to
    is given, plot the ratios to the corresponding alpha effective.
    Otherwise, plot absolute values.
    See the help for ``xplotting_grid`` for information on how to set basis,
    flavours and x ranges. Yields one figure per PDF flavour.

    normalize_to:  Either the name of one of the alpha effective or its
    corresponding index in the list, starting from one, or None to plot
    absolute values.

    xscale: One of the matplotlib allowed scales. If undefined, it will be
    set based on the scale in xgrid, which should be used instead.
    """
    return plot_alpha_eff_internal(
        fits_pdf, alpha_eff_fits, fits_alpha_lines, normalize_to, ybottom, ytop)

beta_eff_pdfs = collect('beta_eff', ('pdfs',))

@figuregen
@check_pdf_normalize_to
def plot_beta_eff_internal(
        pdfs, beta_eff_pdfs, pdfs_beta_lines,
        normalize_to: (int, str, type(None)) = None,
        ybottom=None, ytop=None):
    """ Same as plot_alpha_eff_internal but for beta effective exponent """
    yield from ExponentBandPlotter(
        pdfs_beta_lines, 'beta', pdfs, beta_eff_pdfs, 'linear', normalize_to, ybottom, ytop)

beta_eff_fits = collect('beta_eff', ('fits', 'fitpdf',))

@figuregen
def plot_beta_eff(
        fits_pdf, beta_eff_fits, fits_beta_lines,
        normalize_to: (int, str, type(None)) = None,
        ybottom=None, ytop=None):
    """ Same as plot_alpha_eff but for beta effective exponents """
    return plot_beta_eff_internal(
        fits_pdf, beta_eff_fits, fits_beta_lines, normalize_to, ybottom, ytop)


def previous_effective_exponents(basis:str, fit: (FitSpec, type(None)) = None):
    """If provided with a fit, check that the `basis` is the basis which was fitted
    if so then return the previous effective exponents read from the fit runcard.
    """
    if fit is None:
        return None
    else:
        fitting = fit.as_input()["fitting"]
        if fitting["fitbasis"] == basis:
            return fitting["basis"]
        else:
            return None

@table
def previous_effective_exponents_table(fit: FitSpec):
    """Given a fit, reads the previous exponents from the fit runcard"""
    fitting = fit.as_input()["fitting"]
    checked = check_basis(
        fitting["fitbasis"],
        [runcard_fl['fl'] for runcard_fl in fitting["basis"]])
    basis = checked["basis"]
    flavours = checked["flavours"]
    prev_a_bounds = [runcard_fl['smallx'] for runcard_fl in fitting["basis"]]
    prev_b_bounds = [runcard_fl['largex'] for runcard_fl in fitting["basis"]]
    # make single list alternating alpha and beta bounds
    data = [vals for pair in zip(prev_a_bounds, prev_b_bounds) for vals in pair]
    flavours_label = [f"${basis.elementlabel(fl)}$" for fl in flavours]
    ind = pd.MultiIndex.from_product([flavours_label, [r"$\alpha$", r"$\beta$"]])
    columns = pd.MultiIndex.from_product([[f"prev ({fit.label})"], ["Min", "Max"]])
    return pd.DataFrame(data, index=ind, columns=columns)

@table
@make_argcheck(check_basis)
def next_effective_exponents_table(
    pdf: PDF,
    *,
    x1_alpha: numbers.Real = 1e-6,
    x2_alpha: numbers.Real = 1e-3,
    x1_beta: numbers.Real = 0.65,
    x2_beta: numbers.Real = 0.95,
    basis:(str, Basis),
    flavours: (list, tuple, type(None)) = None,
):
    """Given a PDF, calculate the next effective exponents

    By default `x1_alpha = 1e-6`, `x2_alpha = 1e-3`, `x1_beta = 0.65`, and
    `x2_beta = 0.95`, but different values can be specified in the runcard. The
    values control where the bounds of alpha and beta are evaluated:

    alpha_min:
        singlet/gluon: the 2x68% c.l. lower value evaluated at x=`x1_alpha`
        others  : min(2x68% c.l. lower value evaluated at x=`x1_alpha` and x=`x2_alpha`)

    alpha_max:
        singlet/gluon: min(2 and the 2x68% c.l. upper value evaluated at x=`x1_alpha`)
        others    : min(2 and max(2x68% c.l. upper value evaluated at x=`x1_alpha`
                    and x=`x2_alpha`))

    beta_min:
        max(0 and min(2x68% c.l. lower value evaluated at x=`x1_beta` and x=`x2_beta`))
    beta_max:
        max(2x68% c.l. upper value evaluated at x=`x1_beta` and x=`x2_beta`)

    """
    Qmin = pdf.q_min

    alpha_effs = alpha_eff(
        pdf, xmin=x1_alpha, xmax=x2_alpha, npoints=2, Q=Qmin, basis=basis, flavours=flavours)
    beta_effs = beta_eff(
        pdf, xmin=x1_beta, xmax=x2_beta, npoints=2, Q=Qmin, basis=basis, flavours=flavours)

    eff_exp_data = []

    alphagrid = alpha_effs.grid_values
    betagrid = beta_effs.grid_values

    alphastats = pdf.stats_class(alphagrid)
    betastats = pdf.stats_class(betagrid)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)

        alpha_cv = np.nanmean(alphagrid, axis=0)
        beta_cv = np.nanmean(betagrid, axis=0)
        #tuple of low and high values repectively
        alpha68 = alphastats.errorbar68()
        beta68 = betastats.errorbar68()

    alpha_sigup = alpha68[1] - alpha_cv
    beta_sigup = beta68[1] - beta_cv
    alpha_sigdown = -alpha68[0] + alpha_cv
    beta_sigdown = -beta68[0] + beta_cv
    flavours_label = []
    for (j, fl) in enumerate(flavours):
        # the gluon/singlet case
        if fl in (r"\Sigma", "g"):
            new_alpha_bounds = [
                alpha_cv[j, 0] - 2*alpha_sigdown[j, 0],
                min(2, alpha_cv[j, 0] + 2*alpha_sigup[j, 0])]
        else:
            new_alpha_bounds = [
                min(alpha_cv[j, :] - 2*alpha_sigdown[j, :]),
                min(2, max(alpha_cv[j, :] + 2*alpha_sigup[j, :]))]

        new_beta_bounds = [
            max(0, min(beta_cv[j, :] - 2*beta_sigdown[j, :])),
            max(beta_cv[j, :] + 2*beta_sigup[j, :])]

        eff_exp_data.extend((new_alpha_bounds, new_beta_bounds))
        flavours_label.append(f"${basis.elementlabel(fl)}$")
    ind = pd.MultiIndex.from_product([flavours_label, [r"$\alpha$", r"$\beta$"]])
    eff_exp_columns = pd.MultiIndex.from_product([[f"next ({pdf.label})"], ["Min", "Max"]])
    df = pd.DataFrame(eff_exp_data, index=ind,
                      columns=eff_exp_columns)
    return df

@table
def effective_exponents_table_internal(
    next_effective_exponents_table,
    *,
    fit=None,
    basis,
):
    """Returns a table which concatenates previous_effective_exponents_table
    and next_effective_exponents_table if both tables contain effective exponents
    in the same basis.

    If the previous exponents are in a different basis, or no fit was given to
    read the previous exponents from, then only the next exponents table is
    returned, for plotting purposes.

    """
    if fit is not None and fit.as_input()["fitting"]["fitbasis"] == basis:
        # have to call action here in case fit is None
        previous_table = previous_effective_exponents_table(fit)
        df = pd.concat((previous_table, next_effective_exponents_table), axis=1)
    else:
        df = next_effective_exponents_table
    return df


effective_exponents_table = collect(
    'effective_exponents_table_internal', ('fitpdfandbasis',))
fmt = lambda a: float(significant_digits(a, 4))

next_fit_eff_exps_table = collect("next_effective_exponents_table", ("fitpdfandbasis",))


def iterate_preprocessing_yaml(
    fit, next_fit_eff_exps_table, _flmap_np_clip_arg=None):
    """Using py:func:`next_effective_exponents_table` update the preprocessing
    exponents of the input ``fit``. This is part of the usual pipeline referred
    to as "iterating a fit", for more information see: :ref:`run-iterated-fit`.
    A fully iterated runcard can be obtained from the action
    :py:func:`iterated_runcard_yaml`.

    This action can be used in a report but should be wrapped in a code block
    to be formatted correctly, for example:

    ```yaml
    {@iterate_preprocessing_yaml@}
    ```

    Alternatively, using the API, the yaml dump returned by this function can
    be written to a file e.g

    >>> from validphys.api import API
    >>> yaml_output = API.iterate_preprocessing_yaml(fit=<fit name>)
    >>> with open("output.yml", "w+") as f:
    ...     f.write(yaml_output)

    Parameters
    ----------
    fit: validphys.core.FitSpec
        Whose preprocessing range will be iterated, the output runcard will be
        the same as the one used to run this fit, except with new preprocessing
        range.
    next_fit_eff_exps_table: pd.DataFrame
        Table outputted by :py:func:`next_fit_eff_exps_table` containing
        the next preprocessing ranges.
    _flmap_np_clip_arg: dict
        Internal argument used by ``vp-nextfitruncard``. Dictionary containing
        a mapping like
        ``{<flavour>: {<largex/smallx>: {a_min: <min value>, a_max: <max value>}}}``.
        If a flavour is present in ``_flmap_np_clip_arg`` then the preprocessing
        ranges will be passed through ``np.clip`` with the arguments supplied
        in the mapping.

    """
    (df_effexps,) = next_fit_eff_exps_table
    # Use round trip loader rather than safe_load in fit.as_input()
    with open(fit.path / "filter.yml", "r") as f:
        filtermap = yaml.load(f, yaml.RoundTripLoader)
    previous_exponents = filtermap["fitting"]["basis"]
    basis = filtermap["fitting"]["fitbasis"]
    checked = check_basis(basis, None)
    basis = checked["basis"]

    # use order defined in runcard.
    runcard_flavours = [
        f"{basis.elementlabel(ref_fl['fl'])}" for ref_fl in previous_exponents]
    for i, fl in enumerate(runcard_flavours):
        alphas = df_effexps.loc[(f"${fl}$", r"$\alpha$")].values
        betas = df_effexps.loc[(f"${fl}$", r"$\beta$")].values
        flmap_key = previous_exponents[i]["fl"]
        if _flmap_np_clip_arg is not None and _flmap_np_clip_arg.get(flmap_key) is not None:
            smallx_args = _flmap_np_clip_arg[flmap_key].get("smallx")
            largex_args = _flmap_np_clip_arg[flmap_key].get("largex")
            if smallx_args is not None:
                alphas = np.clip(alphas, **smallx_args)
            if largex_args is not None:
                betas = np.clip(betas, **largex_args)
        previous_exponents[i]["smallx"] = [
            fmt(alpha) for alpha in alphas
        ]
        previous_exponents[i]["largex"] = [
            fmt(beta) for beta in betas
        ]
    return yaml.dump(filtermap, Dumper=yaml.RoundTripDumper)


def update_runcard_description_yaml(
    iterate_preprocessing_yaml, _updated_description=None):
    """Take the runcard with iterated preprocessing and update the description
    if ``_updated_description`` is provided. As with
    :py:func:`iterate_preprocessing_yaml` the result can be used in a report
    but should be wrapped in a code block to be formatted correctly,
    for example:

    ```yaml
    {@update_runcard_description_yaml@}
    ```

    """
    filtermap = yaml.load(iterate_preprocessing_yaml, yaml.RoundTripLoader)

    # update description if necessary
    if _updated_description is not None:
        filtermap["description"] = _updated_description

    return yaml.dump(filtermap, Dumper=yaml.RoundTripDumper)


def iterated_runcard_yaml(
    fit, update_runcard_description_yaml):
    """
    Takes the runcard with preprocessing iterated and description updated then

    - Updates the t0 pdf set to be ``fit``
    - Modifies the random seeds (to random unsigned long ints)

    This should facilitate running a new fit with identical input settings
    as the specified ``fit`` with the t0, seeds and preprocessing iterated. For
    more information see: :ref:`run-iterated-fit`

    This action can be used in a report but should be wrapped in a code block
    to be formatted correctly, for example:

    ```yaml
    {@iterated_runcard_yaml@}
    ```

    alternatively, using the API, the yaml dump returned by this function can
    be written to a file e.g

    >>> from validphys.api import API
    >>> yaml_output = API.iterated_runcard_yaml(
    ...     fit=<fit name>,
    ...     _updated_description="My iterated fit"
    ... )
    >>> with open("output.yml", "w+") as f:
    ...     f.write(yaml_output)

    """
    filtermap = yaml.load(update_runcard_description_yaml, yaml.RoundTripLoader)
    # iterate t0
    filtermap["datacuts"]["t0pdfset"] = fit.name

    # Update seeds with valid pseudorandom unsigned long int
    # Check if seeds exist especially since extra seeds needed in n3fit vs nnfit
    # Start with seeds in "fitting" section of runcard
    fitting_seeds = ["seed", "trvlseed", "nnseed", "mcseed"]
    fitting_data = filtermap.get("fitting")
    maxint = np.iinfo('int32').max
    for seed in fitting_seeds:
        if seed in filtermap:
            filtermap[seed] = random.randrange(0, maxint)
        elif fitting_data is not None and seed in fitting_data:
            #BCH
            # For older runcards the seeds are inside the `fitting` namespace
            fitting_data[seed] = random.randrange(0, maxint)

    # Next "closuretest" section of runcard
    if "closuretest" in filtermap:
        closuretest_data = filtermap["closuretest"]
        if "filterseed" in closuretest_data:
            closuretest_data["filterseed"] = random.randrange(0, maxint)

    return yaml.dump(filtermap, Dumper=yaml.RoundTripDumper)
