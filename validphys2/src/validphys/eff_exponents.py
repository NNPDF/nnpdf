# -*- coding: utf-8 -*-
"""
Tools for computing and plotting effective exponents.
"""
from __future__ import generator_stop

import logging
import warnings
import numbers
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import NNPDF as nnpath

from reportengine.figure import figuregen
from reportengine.table import table
from reportengine.floatformatting import format_number
from reportengine.compat import yaml

from validphys.checks import check_positive, check_pdf_normalize_to
from validphys.pdfplots import BandPDFPlotter, PDFPlotter, FlavourState
from validphys.pdfbases import check_basis
from validphys.core import PDF
from validphys import plotutils


import validphys.pdfgrids as pdfgrids

log = logging.getLogger(__name__)


@check_positive('Q')
@pdfgrids._check_limits
def alpha_eff(pdfs,
              xmin: numbers.Real = 1e-6,
              xmax: numbers.Real = 1e-3,
              npoints: int = 200,
              Q: numbers.Real = 1.65):
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
    any_pdf = pdfs[0]
    pdfpath = nnpath.get_results_path()+any_pdf.name
    filtermap = yaml.safe_load(open(pdfpath+'/filter.yml'))

    #Extracting initial scale from the theory info
    infomap = yaml.safe_load(open(pdfpath+'/nnfit/'+any_pdf.name+'.info'))
    Q = infomap['QMin']
    basis = filtermap['fitting']['fitbasis']+'FitBasis'

    flavours = None
    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']

    alphaGrids = []
    xGrid = pdfgrids.xgrid(xmin, xmax, 'log', npoints)

    for pdf in pdfs:
        pdfGrid = pdfgrids.xplotting_grid(
            pdf, Q, xgrid=xGrid, basis=basis, flavours=flavours)
        pdfGrid_values = pdfGrid.grid_values
        # NOTE: without this I get "setting an array element with a sequence"
        xGrid = pdfGrid.xgrid
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            alphaGrid_values = -np.log(abs(pdfGrid_values/xGrid))/np.log(xGrid)
            alphaGrid_values[alphaGrid_values == -
                             np.inf] = np.nan  # when PDF_i =0
        alphaGrid = pdfGrid._replace(grid_values=alphaGrid_values)
        alphaGrids.append(alphaGrid)

    return alphaGrids  # .grid_values


@check_positive('Q')
@pdfgrids._check_limits
def beta_eff(pdfs,
             xmin: numbers.Real = 0.6,
             xmax: numbers.Real = 0.9,
             npoints: int = 200,
             Q: numbers.Real = 1.65):
    """Return a list of xplotting_grids containing the value of the effective
    exponent beta at the specified values of x and flavour.
    beta is relevant at large x, hence the linear scale.

    basis: Is one of the bases defined in pdfbases.py. This includes 'flavour'
    and 'evolution'.

    flavours: A set of elements from the basis.
    If None, the defaults for that basis will be selected.

    Q: The PDF scale in GeV.
    """
    #Loading the filter map of the fit/PDF
    any_pdf = pdfs[0]
    pdfpath = nnpath.get_results_path()+any_pdf.name
    filtermap = yaml.safe_load(open(pdfpath+'/filter.yml'))

    #Extracting initial scale from the theory info
    infomap = yaml.safe_load(open(pdfpath+'/nnfit/'+any_pdf.name+'.info'))
    Q = infomap['QMin']
    basis = filtermap['fitting']['fitbasis']+'FitBasis'

    flavours = None
    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']

    betaGrids = []
    xGrid = pdfgrids.xgrid(xmin, xmax, 'linear', npoints)

    for pdf in pdfs:
        pdfGrid = pdfgrids.xplotting_grid(
            pdf, Q, xgrid=xGrid, basis=basis, flavours=flavours)
        pdfGrid_values = pdfGrid.grid_values
        # NOTE: without this I get "setting an array element with a sequence"
        xGrid = pdfGrid.xgrid
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            betaGrid_values = np.log(abs(pdfGrid_values/xGrid))/np.log(1-xGrid)
            betaGrid_values[betaGrid_values == -
                            np.inf] = np.nan  # when PDF_i =0
        betaGrid = pdfGrid._replace(grid_values=betaGrid_values)
        betaGrids.append(betaGrid)

    return betaGrids  # .grid_values


class PreprocessingPlotter(PDFPlotter):
    """ Class inherenting from BandPDFPlotter, has the same functionality
    but with overloaded title and ylabel to take into account the effective
    exponents names.
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

    def __call__(self,):
        if not self.xplotting_grids:
            return

        basis = self.firstgrid.basis

        for flindex, fl in enumerate(self.firstgrid.flavours):
            fig, ax = plt.subplots()
            parton_name = basis.elementlabel(fl)
            flstate = FlavourState(flindex=flindex, fl=fl, fig=fig, ax=ax,
                                   parton_name=parton_name)

            self.setup_flavour(flstate)
            ax.set_title(self.get_title(parton_name))

            all_vals = []
            dfs = []
            for pdf, grid in zip(self.pdfs, self.xplotting_grids):
                limits = self.draw(pdf, grid, flstate)
                if limits is not None:
                    all_vals.append(np.atleast_2d(limits))
                #here, by default x1_alpha=1e-6,x2_alpha=1e-3,x1_beta=0.65,x2_beta=0.95
                dfs.append(effective_exponents_table(pdf))

            #Note these two lines do not conmute!
            ax.set_xscale(self.xscale)
            plotutils.frame_center(
                ax, self.firstgrid.xgrid, np.concatenate(all_vals))
            if (self.ymin is not None):
                ax.set_ylim(ymin=self.ymin)
            if (self.ymax is not None):
                ax.set_ylim(ymax=self.ymax)

            ax.set_xlabel('$x$')
            ax.set_xlim(self.firstgrid.xgrid[0])

            ax.set_ylabel(self.get_ylabel(parton_name))

            ax.set_axisbelow(True)

            for df in dfs:
                if self.exponent == 'alpha':
                    #prev
                    ax.axhline(df.iat[2*flindex, 1],
                               linestyle='--', label='prev')
                    #self.labels.append('prev')
                    ax.axhline(df.iat[2*flindex, 2], linestyle='--')
                    #next
                    ax.axhline(df.iat[2*flindex, 3], label='next')
                    ax.axhline(df.iat[2*flindex, 4])

                elif self.exponent == 'beta':
                    #prev
                    ax.axhline(df.iat[2*flindex+1, 1],
                               linestyle='--', label='prev')
                    ax.axhline(df.iat[2*flindex+1, 2], linestyle='--')
                    #next
                    ax.axhline(df.iat[2*flindex+1, 3], label='next')
                    ax.axhline(df.iat[2*flindex+1, 4])

            self.legend(flstate)

            yield fig, parton_name


class ExponentBandPlotter(BandPDFPlotter, PreprocessingPlotter):
    pass


@figuregen
@check_pdf_normalize_to
def plot_alphaEff(pdfs,
                  alpha_eff,
                  normalize_to: (int, str, type(None)) = None,
                  ybottom=None,
                  ytop=None):
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
    yield from ExponentBandPlotter('alpha', pdfs, alpha_eff, 'log', normalize_to, ybottom, ytop)


@figuregen
@check_pdf_normalize_to
def plot_betaEff(pdfs,
                 beta_eff,
                 normalize_to: (int, str, type(None)) = None,
                 ybottom=None,
                 ytop=None):
    """ Same as plot_alphaEff but for beta effective exponent """
    yield from ExponentBandPlotter('beta', pdfs, beta_eff, 'linear', normalize_to, ybottom, ytop)


@table
def effective_exponents_table(pdf: PDF,
                              x1_alpha: numbers.Real = 1e-6,
                              x2_alpha: numbers.Real = 1e-3,
                              x1_beta: numbers.Real = 0.65,
                              x2_beta: numbers.Real = 0.95):
    """Returns a table with the effective exponents for the next fit"""

    #Reading from the filter
    pdfpath = nnpath.get_results_path()+pdf.name
    filtermap = yaml.safe_load(open(pdfpath+'/filter.yml'))
    infomap = yaml.safe_load(open(pdfpath+'/nnfit/'+pdf.name+'.info'))
    Q = infomap['QMin']
    basis = filtermap['fitting']['fitbasis']+'FitBasis'

    flavours = None
    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']

    pdfs = [pdf]
    alphamin_grids = alpha_eff(pdfs,
                               xmin=x1_alpha,
                               xmax=x1_alpha,
                               npoints=1,
                               Q=Q)
    alphamax_grids = alpha_eff(pdfs,
                               xmin=x2_alpha,
                               xmax=x2_alpha,
                               npoints=1,
                               Q=Q)
    betamin_grids = beta_eff(pdfs,
                             xmin=x1_beta,
                             xmax=x1_beta,
                             npoints=1, Q=Q)
    betamax_grids = beta_eff(pdfs,
                             xmin=x2_beta,
                             xmax=x2_beta,
                             npoints=1, Q=Q)

    eff_exp_data = []

    alphamin_grid_values = alphamin_grids[0].grid_values
    alphamax_grid_values = alphamax_grids[0].grid_values
    betamin_grid_values = betamin_grids[0].grid_values
    betamax_grid_values = betamax_grids[0].grid_values

    alphamin_stats = pdf.stats_class(alphamin_grid_values)
    alphamax_stats = pdf.stats_class(alphamax_grid_values)
    betamin_stats = pdf.stats_class(betamin_grid_values)
    betamax_stats = pdf.stats_class(betamax_grid_values)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)

        alphamin_cv = np.nanmean(alphamin_grid_values, axis=0)
        alphamax_cv = np.nanmean(alphamax_grid_values, axis=0)
        betamin_cv = np.nanmean(betamin_grid_values, axis=0)
        betamax_cv = np.nanmean(betamax_grid_values, axis=0)

        alphamin_err68down, alphamin_err68up = alphamin_stats.errorbar68()
        alphamax_err68down, alphamax_err68up = alphamax_stats.errorbar68()
        betamin_err68down, betamin_err68up = betamin_stats.errorbar68()
        betamax_err68down, betamax_err68up = betamax_stats.errorbar68()

        alphamin_sigup = alphamin_err68up - alphamin_cv
        alphamax_sigup = alphamax_err68up - alphamax_cv
        betamin_sigup = betamin_err68up - betamin_cv
        betamax_sigup = betamax_err68up - betamax_cv

        alphamin_sigdown = -alphamin_err68down + alphamin_cv
        alphamax_sigdown = -alphamax_err68down + alphamax_cv
        betamin_sigdown = -betamin_err68down + betamin_cv
        betamax_sigdown = -betamax_err68down + betamax_cv

    flavours_label = []

    for (j, fl) in enumerate(flavours):

        prev_amin_bound = None
        prev_amax_bound = None
        prev_bmin_bound = None
        prev_bmax_bound = None
        YAMLlabel = ""
        YAMLflaliases = {r'\Sigma': 'sng', 'V': 'v', 'T3': 't3',
                         'V3': 'v3', 'T8': 't8', 'V8': 'v8', 'gluon': 'g', r'c^+': 'cp'}

        #Matching the flavour between vp and the runcards
        for k, ref_fl in enumerate(filtermap['fitting']['basis']):
            if basis.elementlabel(fl) in YAMLflaliases.keys():
                YAMLlabel = YAMLflaliases[basis.elementlabel(fl)]
            else:
                YAMLlabel = basis.elementlabel(fl)

            if ref_fl['fl'] == YAMLlabel:
                prev_amin_bound = filtermap['fitting']['basis'][k]['smallx'][0]
                prev_amax_bound = filtermap['fitting']['basis'][k]['smallx'][1]
                prev_bmin_bound = filtermap['fitting']['basis'][k]['largex'][0]
                prev_bmax_bound = filtermap['fitting']['basis'][k]['largex'][1]

        # the gluon/singlet case
        #Defining the bounds
        # alpha_min = singlet/gluon: the 2x68% c.l. lower value evaluated at x=1e-6.
        #                  others  : min(2x68% c.l. lower value evaluated at x=1e-6 and x=1e-3)
        # alpha_max = singlet/gluon: min(2 and the 2x68% c.l. upper value evaluated at x=1e-6)
        #                others    : min(2 and max(2x68% c.l. upper value evaluated at x=1e-6 and x=1e-3))
        if fl == r'\Sigma' or fl == "g":
            new_min_bound = alphamin_cv[j][0]-2*alphamin_sigdown[j][0]
            new_max_bound = round(
                min(2, alphamin_cv[j][0]+2*alphamin_sigup[j][0]), 3)
            alpha_line = [r"$\alpha$", prev_amin_bound,
                          prev_amax_bound, new_min_bound, new_max_bound]
        else:
            new_min_bound = round(
                min(alphamin_cv[j][0]-2*alphamin_sigdown[j][0],
                    alphamax_cv[j][0]-2*alphamax_sigdown[j][0]), 3)
            new_max_bound = round(
                min(2, max(alphamin_cv[j][0]+2*alphamin_sigup[j][0],
                           alphamax_cv[j][0]+2*alphamax_sigup[j][0])), 3)

            alpha_line = [r"$\alpha$", prev_amin_bound,
                          prev_amax_bound, new_min_bound, new_max_bound]

        # beta_min  =  max(0 and min(2x68% c.l. lower value evaluated at x=0.65 and x=0.95))
        # beta_max  =  max(2x68% c.l. upper value evaluated at x=0.65 and x=0.95)
        new_min_bound = round(
            max(0,
                min(betamin_cv[j][0]-2*betamin_sigdown[j][0],
                    betamax_cv[j][0]-2*betamax_sigdown[j][0])), 3)

        new_max_bound = round(
            max(betamin_cv[j][0]+2*betamin_sigup[j][0],
                betamax_cv[j][0]+2*betamax_sigup[j][0]), 3)

        beta_line = [r"$\beta$", prev_bmin_bound,
                     prev_bmax_bound, new_min_bound, new_max_bound]

        eff_exp_data.append(alpha_line)
        eff_exp_data.append(beta_line)
        #flavours_label.append(f'${basis.elementlabel(fl)}$')
        flavours_label.append(f'${fl}$')
        flavours_label.append("")

    eff_exp_columns = [r"$\alpha/\beta$", "prev Min",
                       "prev Max", "next Min", "next Max"]
    df = pd.DataFrame(eff_exp_data, index=flavours_label,
                      columns=eff_exp_columns)
    df.name = pdf.name
    return df


def next_effective_exponents_yaml(pdf: PDF,
                                  x1_alpha: numbers.Real = 1e-6,
                                  x2_alpha: numbers.Real = 1e-3,
                                  x1_beta: numbers.Real = 0.65,
                                  x2_beta: numbers.Real = 0.95):
    """-Returns a table in yaml format called NextEffExps.yaml in output
       -Prints the yaml table in the report"""

    df_effexps = effective_exponents_table(pdf,
                                           x1_alpha,
                                           x2_alpha,
                                           x1_beta,
                                           x2_beta)
    #Reading from the filter
    pdfpath = nnpath.get_results_path()+pdf.name
    filtermap = yaml.safe_load(open(pdfpath+'/filter.yml'))
    basis = filtermap['fitting']['fitbasis']+'FitBasis'

    flavours = None
    checked = check_basis(basis, flavours)
    basis = checked['basis']
    flavours = checked['flavours']

    YAMLflaliases = {r'\Sigma': 'sng', 'V': 'v', 'T3': 't3',
                     'V3': 'v3', 'T8': 't8', 'V8': 'v8', 'g': 'g', r'c^+': 'cp'}

    import io
    s = io.StringIO()
    inputs = []

    with open("output/NextEffExps.yaml", 'w') as out:
        out.write("basis:\n")

        for (j, fl) in enumerate(flavours):

            YAMLlabel = " "
            if basis.elementlabel(fl) in YAMLflaliases.keys():
                YAMLlabel = YAMLflaliases[basis.elementlabel(fl)]
            else:
                YAMLlabel = basis.elementlabel(fl)

            for k, ref_fl in enumerate(filtermap['fitting']['basis']):
                if basis.elementlabel(fl) in YAMLflaliases.keys():
                    YAMLlabel = YAMLflaliases[basis.elementlabel(fl)]
                else:
                    YAMLlabel = basis.elementlabel(fl)

                if ref_fl['fl'] == YAMLlabel:
                    pos = filtermap['fitting']['basis'][k]['pos']
                    mutsize = "["+str(filtermap['fitting']
                                      ['basis'][k]['mutsize'][0])+"]"
                    mutprob = "["+str(filtermap['fitting']
                                      ['basis'][k]['mutprob'][0])+"]"
                    smallx = "["+str(df_effexps.iat[2*j, 3])+", " + \
                        str(df_effexps.iat[2*j, 4])+"]"
                    largex = "["+str(df_effexps.iat[2*j+1, 3]) + \
                        ", "+str(df_effexps.iat[2*j+1, 4])+"]"

                    out.write(" - {")
                    out.write("fl: "+YAMLlabel)
                    out.write(", pos: "+pos)
                    out.write(", mutsize: "+mutsize)
                    out.write(", mutprob: "+mutprob)
                    out.write(", smallx: "+smallx)
                    out.write(", largex: "+largex)
                    out.write("")
                    out.write("}\n")

                    d = {'fl': YAMLlabel, 'pos': pos, 'mutsize': mutsize,
                         'mutprob': mutprob, 'smallx': smallx, 'largex': largex}
                    inputs.append(d)

    yaml.dump(inputs, s)

    return s.getvalue()
