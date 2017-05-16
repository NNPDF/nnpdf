# -*- coding: utf-8 -*-
"""
Figures for visualizing results
"""
from __future__ import generator_stop

import logging
import functools
import warnings
import abc
from types import SimpleNamespace
from collections import defaultdict, Sequence
import copy
import numbers

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors, ticker as mticker
import scipy.stats as stats

from reportengine.figure import figure, figuregen
from reportengine.checks import make_check, CheckError, make_argcheck

from validphys.core import MCStats
from validphys.results import chi2_stat_labels
from validphys.plotoptions import get_info, kitable, transform_result
from validphys.checks import check_scale
from validphys import plotutils
from validphys.utils import sane_groupby_iter, split_ranges

from validphys.gridvalues import LUMI_CHANNELS

log = logging.getLogger(__name__)

@figure
def plot_chi2dist(results, dataset, abs_chi2_data, chi2_stats, pdf):
    """Plot the distribution of chi²s of the members of the pdfset."""
    setlabel = dataset.name
    fig, ax = plt.subplots()
    label = pdf.name
    alldata, central, npoints = abs_chi2_data
    if not isinstance(alldata, MCStats):
        ax.set_facecolor("#ffcccc")
        log.warn("Chi² distribution plots have a different meaning for non MC sets.")
        label += " (%s!)" % pdf.ErrorType
    label += '\n'+ '\n'.join(str(chi2_stat_labels[k])+(' %.2f' % v) for (k,v) in chi2_stats.items())
    ax.set_title("$\chi^2$ distribution for %s" % setlabel)

    ax.hist(alldata.data, label=label, zorder=100)
    l = ax.legend()
    l.set_zorder(1000)
    return fig

#TODO: This should be simplified if at all possible. For now some more examples
#are needed for a spec to emerge.
@make_check
def _check_normalize_to(ns, **kwargs):
    """Transforn normalize_to into an index."""

    msg = ("normalize_to should be either 'data', a pdf id or an index of the "
           "result (0 for the data, and i for the ith pdf)")

    val = ns.get('normalize_to', None)
    if val is None:
        return

    if 'pdf' in ns:
        names = ['data', ns['pdf'].name]
    else:
        names = ['data', *(pdf.name for pdf in ns['pdfs'])]
    if isinstance(val, int):
        if not val < len(names):
            raise CheckError(msg)
        return

    if isinstance(val, str):
        try:
            val = names.index(val)
        except ValueError:
            raise CheckError(msg, val, alternatives=names)
        ns['normalize_to'] = val
        return

    raise RuntimeError("Should not be here")


@_check_normalize_to
@figuregen
def plot_fancy(one_or_more_results, dataset,
               normalize_to:(int,str,type(None)) = None):
    """
    Read the PLOTTING configuration for the dataset and generate the
    corrspondig data theory plot.

    The input results are assumed to be such that the first one is the data,
    and the subsequent ones are the predictions for the PDFfs. See
    ``one_or_more_results``. The labelling of the predictions can be
    influenced by setting ``label`` attribute of theories and pdfs.

    normalize_to: should be either 'data', a pdf id or an index of the
    result (0 for the data, and i for the ith pdf). None means plotting
    absolute values.

    See docs/plotting_format.md for details on the format of the PLOTTING
    files.
    """


    results = one_or_more_results

    info = get_info(dataset, normalize=(normalize_to is not None))

    table = kitable(dataset, info)
    nkinlabels = len(table.columns)

    if normalize_to is not None:
        norm_result = results[normalize_to]
        #We modify the table, so we pass only the label columns
        norm_cv, _ = transform_result(norm_result.central_value,
                                   norm_result.std_error,
                                   table.iloc[:,:nkinlabels], info)

    for i,result in enumerate(results):
        #We modify the table, so we pass only the label columns
        cv, err = transform_result(result.central_value, result.std_error,
                                   table.iloc[:,:nkinlabels], info)
        #By doing tuple keys we avoid all possible name collisions
        if normalize_to is None:
            table[('cv', i)] = cv
            table[('err', i)] = err
        else:
            table[('cv', i)] = cv/norm_cv
            table[('err', i)] = err/norm_cv

    figby = sane_groupby_iter(table, info.figure_by)


    for samefig_vals, fig_data in figby:
        #For some reason matplotlib doesn't set the axis right
        min_vals = []
        max_vals = []
        fig, ax = plt.subplots()
        plotutils.setup_ax(ax)
        ax.set_title("%s %s"%(info.dataset_label,
                     info.group_label(samefig_vals, info.figure_by)))

        lineby = sane_groupby_iter(fig_data, info.line_by)

        first = True


        for (sameline_vals, line_data) in lineby:
            ax.set_prop_cycle(None)
            if first:
                labels = True
            else:
                labels = False
            first = False

            offset_iter = plotutils.offset_xcentered(len(results), ax)

            x = info.get_xcol(line_data)

            try:
                x = np.asanyarray(x, np.float)
            except ValueError:
                xticklabels = x
                npoints = len(x)
                x = np.arange(npoints)
                ax.set_xticks(x)
                ax.set_xticklabels(xticklabels)
                #TODO: Remove this when mpl stops doing the idiotic thing
                #(in v2?)
                ax.set_xlim(-npoints/20, npoints - 1+ npoints/20)


            #Use black for the first iteration (data),
            #and follow the cycle for
            #the rest.
            color = '#262626'
            for i, res in enumerate(results):

                if labels:
                    label = res.label
                else:
                    label = None
                cv = line_data[('cv', i)].as_matrix()
                err = line_data[('err', i)].as_matrix()



                max_vals.append(np.nanmax(cv+err))
                min_vals.append(np.nanmin(cv-err))


                ax.errorbar(x, cv, yerr=err,
                     linestyle='--',
                     lw=0.25,
                     label= label,
                     #elinewidth = 2,
                     capsize=2,
                     marker = 's',
                     markeredgewidth=0.25,
                     c=color,
                     zorder=1000,
                     transform=next(offset_iter))

                color = 'C'+str(i)

            glabel = info.group_label(sameline_vals, info.line_by)

            #Use some anchor that is not in y=1 for ratio plots
            if normalize_to is not None:
                next_after_normalize = (normalize_to + 1) % len(results)
                annotate_point = x[-1], line_data[('cv', next_after_normalize)].as_matrix()[-1]
            else:
                annotate_point = x[-1], line_data[('cv', 0)].as_matrix()[-1]
            ax.annotate(glabel, annotate_point, xytext=(15 ,-10),
                             size='xx-small',
                             textcoords='offset points', zorder=10000)



        if info.x_scale:
            ax.set_xscale(info.x_scale)

        if info.y_scale:
            ax.set_yscale(info.y_scale)

        if normalize_to is None:
            if info.y_label:
                ax.set_ylabel(info.y_label)
        else:
            ax.set_ylabel("Ratio to %s" % norm_result.label)


        ax.legend().set_zorder(100000)
        ax.set_xlabel(info.xlabel)
        fig.tight_layout()


        yield fig

def _scatter_marked(ax, x, y, marked_dict, *args, **kwargs):
    kwargs['s'] = kwargs.get('s', 30) + 10
    x = np.array(x, copy=False)
    y = np.array(y, copy=False)
    for label, indexes in marked_dict.items():
        ax.scatter(x[indexes],y[indexes], *args, **kwargs, label=label,
                   facecolors='none', linewidth=0.5, edgecolor='red')
        kwargs['s'] += 10


#I need to use the informations contained in experiments_chi2_table

@figure
def experiments_chi2_plot(experiments, experiments_chi2):
    """Return a plot with the chi² of all the experiments"""
    exchi2 = []
    xvalues = []
    xticks = []
    for experiment, expres in zip(experiments, experiments_chi2):
        exchi2.append(expres.central_result/expres.ndata)
        xvalues.append(experiments_chi2.index(expres))
        xticks.append(experiment.name) 
    fig, ax = plt.subplots()
    width = 0.5
    plt.xticks(xvalues, xticks,rotation=80)
    ax.bar(xvalues, exchi2, width)
    return fig

@figure
def datasets_chi2_plot(experiments, experiments_chi2,each_dataset_chi2):
    """Return a plot with the chi² of all the datasets"""
    ds = iter(each_dataset_chi2)
    dschi2 = []
    xvalues = []
    xticks = []
    for experiment, expres in zip(experiments, experiments_chi2):
        for dataset, dsres in zip(experiment, ds):
            stats = dsres.central_result
            dschi2.append(dsres.central_result/dsres.ndata)
            xticks.append(dataset.name)
            xvalues.append(xticks.index(dataset.name))
    fig, ax = plt.subplots()
    width = 0.5
    plt.xticks(xvalues, xticks,rotation=80)
    ax.bar(xvalues, dschi2, width)
    return fig

@figure
def plot_training_validation(fit, replica_data, replica_filters=None):
    """Scatter plot with the training and validation chi² for each replica
    in the fit. The mean is also displayed"""
    training, valid = zip(*((dt.training, dt.validation) for dt in replica_data))
    fig, ax = plt.subplots()
    ax.plot(training, valid, marker='o', linestyle='none', markersize=5, zorder=100)
    if replica_filters:
        _scatter_marked(ax, training,valid, replica_filters, zorder=90)
        ax.legend().set_zorder(10000)

    ax.set_title(getattr(fit, 'label', fit.name))

    ax.set_xlabel(r'$\chi^2/N_{dat}$ train')
    ax.set_ylabel(r'$\chi^2/N_{dat}$ valid')

    ax.plot(np.mean(training), np.mean(valid),
         marker='s', color='red', markersize=7, zorder=1000)

    return fig


@figure
def plot_trainvaliddist(fit, replica_data):
    """KDEs for the trainning and validation distributions for
    each replica in the fit."""
    training, valid = zip(*((dt.training, dt.validation) for dt in replica_data))
    fig, ax = plt.subplots()

    kde_train = stats.gaussian_kde(training, bw_method='silverman')
    kde_valid = stats.gaussian_kde(valid, bw_method='silverman')
    mean = (np.array(training) + np.array(valid))*0.5
    kde_mean = stats.gaussian_kde(mean, bw_method='silverman')

    x = np.linspace(np.min([training,valid]),
                    np.max([training, valid]), 150)
    ax.plot(x, kde_train(x), label="Training")
    ax.plot(x, kde_valid(x), label="Validation")
    ax.plot(x, kde_mean(x), label="Mean")

    ax.set_xlabel(r"$\chi^2/N_{dat}$")
    ax.set_title("KDE of the fit distributions for %s" % getattr(fit,
                                                                 'label', fit.name))

    ax.set_ylim(0, None)
    ax.legend()
    return fig

@figure
def plot_covmat_eigs(experiment):
    """Plot the eigenvalues of the covariance matrix for a given experiment."""
    eigs = la.eigvalsh(experiment.load().get_covmat())
    fig,ax = plt.subplots()
    x = np.arange(1,len(eigs) + 1)
    ax.plot(x, eigs, 'o', markersize=10)
    ax.set_yscale('log')
    ax.yaxis.grid(False)
    plt.title("Covmat eigenvalues for %s" % experiment.name)
    plt.xlabel("# Eigenvector")
    return fig

@figure
def plot_corrmat_eigs(experiment):
    """Plot the eigenvalues of the correlation matrix for a given experiment."""
    covmat = experiment.load().get_covmat()
    stds = np.sqrt(np.diag(covmat))
    corrmat = covmat/np.outer(stds,stds)
    eigs = la.eigvalsh(corrmat)
    fig,ax = plt.subplots()
    ax.plot(eigs, 'o')
    ax.set_yscale('log')
    return fig


#The indexing to one instead of zero is so that we can be consistent with
#how plot_fancy works, so normalize_to: 1 would normalize to the first pdf
#for both.
@make_argcheck
def _check_pdf_normalize_to(pdfs, normalize_to):
    """Transforn normalize_to into an index."""

    msg = ("normalize_to should be, a pdf id or an index of the "
           "pdf (starting from one)")

    if normalize_to is None:
        return

    names = [pdf.name for pdf in pdfs]
    if isinstance(normalize_to, int):
        normalize_to -= 1
        if not normalize_to < len(names) or normalize_to<0:
            raise CheckError(msg)
        return {'normalize_to': normalize_to}

    if isinstance(normalize_to, str):
        try:
            normalize_to = names.index(normalize_to)
        except ValueError:
            raise CheckError(msg, normalize_to, alternatives=names)
        return {'normalize_to': normalize_to}


    raise RuntimeError("Should not be here")



def _scale_from_grid(grid):
    return 'linear' if grid.scale == 'linear' else 'log'


class FlavourState(SimpleNamespace):
    """This is the namespace for the pats specific for each flavour"""
    pass


class PDFPlotter(metaclass=abc.ABCMeta):
    """Stateful object breaks plotting grids by favour, as a function of x and
    for fixed Q.

    This class has a lot of state, but it should all be defined at
    initialization time. Things that change e.g. per flavour should be passed
    explicitly as arguments.
    """

    def __init__(self, pdfs, xplotting_grids, xscale, normalize_to):
        self.pdfs = pdfs
        self._xplotting_grids = xplotting_grids
        self._xscale = xscale
        self.normalize_to = normalize_to
        self.xplotting_grids = self.normalize()


    def setup_flavour(self, flstate):
        pass


    def normalize(self):
        normalize_to = self.normalize_to
        if normalize_to is not None:
            normalize_pdf = self.normalize_pdf
            normalize_grid = self._xplotting_grids[normalize_to]
            normvals = normalize_pdf.stats_class(
                            normalize_grid.grid_values).central_value()

            #Handle division by zero more quietly
            def fp_error(tp, flag):
                log.warn("Invalid values found computing normalization to %s: "
                 "Floating point error (%s).", normalize_pdf, tp)
                #Show warning only once
                np.seterr(all='ignore')

            newgrids = []
            with np.errstate(all='call'):
                np.seterrcall(fp_error)
                for grid in self._xplotting_grids:
                    newvalues = grid.grid_values/normvals
                    #newgrid is like the old grid but with updated values
                    newgrid = type(grid)(**{**grid._asdict(),
                                             'grid_values':newvalues})
                    newgrids.append(newgrid)

            return newgrids
        return self._xplotting_grids

    @property
    def normalize_pdf(self):
        if self.normalize_to is None:
            raise AttributeError("Need to set a normalize_to index.")
        return self.pdfs[self.normalize_to]

    def get_ylabel(self, parton_name):
        if self.normalize_to is not None:
            return "Ratio to {}".format(self.normalize_pdf.label)
        else:
            return '$x{}(x)$'.format(parton_name)

    def get_title(self, parton_name):
        return "$%s$ at %.1f GeV" % (parton_name, self.Q)

    @property
    def xscale(self):
        if self._xscale is None:
            return _scale_from_grid(self.firstgrid)
        return self._xscale

    @property
    def Q(self):
        return self.firstgrid.Q

    @property
    def firstgrid(self):
        if self.xplotting_grids:
            return self.xplotting_grids[0]
        raise AttributeError("Need at least one xgrid")


    @abc.abstractmethod
    def draw(self, pdf, grid, flstate):
        """Plot the desired function of the grid and return the array to be
        used for autoscaling"""
        pass

    def legend(self, flstate):
        return flstate.ax.legend()

    def __iter__(self):
        yield from self()


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
            for pdf, grid in zip(self.pdfs, self.xplotting_grids):
                all_vals.append(np.atleast_2d(self.draw(pdf, grid, flstate)))

            #Note these two lines do not conmute!
            ax.set_xscale(self.xscale)
            plotutils.frame_center(ax, self.firstgrid.xgrid, np.concatenate(all_vals))

            ax.set_xlabel('x')


            ax.set_ylabel(self.get_ylabel(parton_name))

            ax.set_axisbelow(True)

            self.legend(flstate)
            yield fig, parton_name



@functools.lru_cache()
def _warn_pdf_not_montecarlo(pdf):
    et = pdf.ErrorType
    if et != 'replicas':
        log.warn("Plotting members of a non-Monte Carlo PDF set:"
        " %s with error type '%s'.", pdf.name, et)

#Cant't add the lru_cache here because pdfs is not hashable at the moment
@make_argcheck
def _warn_any_pdf_not_montecarlo(pdfs):
    for pdf in pdfs:
        _warn_pdf_not_montecarlo(pdf)


class ReplicaPDFPlotter(PDFPlotter):
    def draw(self, pdf, grid, flstate):
        ax = flstate.ax
        next_prop = next(ax._get_lines.prop_cycler)
        color = next_prop['color']
        gv = grid.grid_values[:,flstate.flindex,:]
        ax.plot(grid.xgrid, gv.T, alpha=0.2, linewidth=0.5,
                color=color, zorder=1)
        stats = pdf.stats_class(gv)
        ax.plot(grid.xgrid, stats.central_value(), color=color,
                linewidth=2,
                label=pdf.label)
        return gv

@figuregen
@_check_pdf_normalize_to
@check_scale('xscale', allow_none=True)
@_warn_any_pdf_not_montecarlo
def plot_pdfreplicas(pdfs, xplotting_grids, xscale:(str,type(None))=None,
                      normalize_to:(int,str,type(None))=None):
    """Plot the replicas of the specifid PDFs. Otherise it works the same as
    plot_pdfs.

    - xscale sets the scale of the plot. E.g. 'linear' or 'log'. Default is
    deduced from the xplotting_grid, which in turn is 'log' by default.

    - normalize_to should be, a pdf id or an index of the pdf (starting from one).
    """
    yield from ReplicaPDFPlotter(pdfs=pdfs, xplotting_grids=xplotting_grids,
                                 xscale=xscale, normalize_to=normalize_to)


class UncertaintyPDFPlotter(PDFPlotter):

    def get_ylabel(self, parton_name):
        if self.normalize_to is not None:
            return r"$\sigma($%s$)$" % super().get_ylabel(parton_name)
        return r"$\sigma/\sigma_{ref}$"

    def draw(self, pdf, grid, flstate):
        ax = flstate.ax
        flindex = flstate.flindex
        gv = grid.grid_values[:,flindex,:]
        stats = pdf.stats_class(gv)

        res = stats.std_error()

        ax.plot(grid.xgrid, res, label=pdf.label)

        return res

@figuregen
@_check_pdf_normalize_to
@check_scale('xscale', allow_none=True)
def plot_pdf_uncertainties(pdfs, xplotting_grids, xscale:(str,type(None))=None,
                      normalize_to:(int,str,type(None))=None):
    """Plot the PDF standard deviations as a function of x.
    If normalize_to is set, the ratio to that
    PDF's central value is plotted. Otherwise it is the absolute values."""
    yield from UncertaintyPDFPlotter(pdfs, xplotting_grids, xscale, normalize_to)

class DistancePDFPlotter(PDFPlotter):

    def normalize(self):
        #if self.normalize_to is None:
            #throw an error and exit
        normalize_to = self.normalize_to
        if normalize_to is not None:
            normalize_pdf = self.normalize_pdf
            normalize_grid = self._xplotting_grids[normalize_to]
            normvals_central = normalize_pdf.stats_class(
                            normalize_grid.grid_values).central_value()
            normvals_sigma = normalize_pdf.stats_class(
                            normalize_grid.grid_values).std_error()

            #need sigma of normalize_pdf
            #need sigma of each pdfs
            #self.pdfs

            #distance = |central_i - central_norm  |/(sigma_i^2 + sigma_norm^2)^1/2

            #Handle division by zero more quietly
            def fp_error(tp, flag):
                log.warn("Invalid values found computing normalization to %s: "
                 "Floating point error (%s).", normalize_pdf, tp)
                #Show warning only once
                np.seterr(all='ignore')

            newgrids = []
            with np.errstate(all='call'):
                np.seterrcall(fp_error)
                for grid,pdf in zip(self._xplotting_grids,self.pdfs):
                    numerator = pow(pdf.stats_class(grid.grid_values).central_value()-normvals_central,2)
                    denominator = pow(pdf.stats_class(grid.grid_values).std_error(),2)+pow(normvals_sigma,2)
                    newvalues = np.sqrt(numerator/denominator)
                    #newgrid is like the old grid but with updated values
                    newgrid = type(grid)(**{**grid._asdict(),
                                             'grid_values':newvalues})
                    newgrids.append(newgrid)

            return newgrids
        return self._xplotting_grids

    def get_ylabel(self, parton_name):
        #if self.normalize_to is None:
            #throw an error and exit
        return "distance from {}".format(self.normalize_pdf.label)


    def draw(self, pdf, grid, flstate):
        ax = flstate.ax
        flindex = flstate.flindex
        gv = grid.grid_values[flindex,:]

        p,=ax.plot(grid.xgrid, gv, label=pdf.label)
        color=p.get_color()
        if(self.normalize_pdf.ErrorType=="replicas" and pdf.ErrorType=="replicas"):
            draw_line = np.sqrt((1./(len(self.normalize_pdf)-1)+1./(len(pdf)-1))/2)
            ax.axhline(draw_line,color=color,alpha=0.5,linestyle="--")


        return gv

@figuregen
@_check_pdf_normalize_to
@check_scale('xscale', allow_none=True)
def plot_pdfdistance(pdfs, xplotting_grids, xscale:(str,type(None))=None,
                      normalize_to:(int,str,type(None))=None):
    """Plot the PDF standard deviations as a function of x.
    If normalize_to is set, the ratio to that
    PDF's central value is plotted. Otherwise it is the absolute values."""
    yield from DistancePDFPlotter(pdfs, xplotting_grids, xscale, normalize_to)

class BandPDFPlotter(PDFPlotter):
    def setup_flavour(self, flstate):
        flstate.handles=[]
        flstate.labels=[]
        flstate.hatchit=plotutils.hatch_iter()

    def draw(self, pdf, grid, flstate):
        ax = flstate.ax
        flindex = flstate.flindex
        hatchit = flstate.hatchit
        labels = flstate.labels
        handles = flstate.handles
        stats = pdf.stats_class(grid.grid_values[:,flindex,:])
        pcycler = ax._get_lines.prop_cycler
        #This is ugly but can't think of anything better

        next_prop = next(pcycler)
        cv = stats.central_value()
        xgrid = grid.xgrid
        #Ignore spurious normalization warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            err68down, err68up = stats.errorbar68()

        color = next_prop['color']
        ax.plot(xgrid, cv, color=color)
        alpha = 0.5
        ax.fill_between(xgrid, err68up, err68down, color=color, alpha=alpha,
                        zorder=1)
        #http://stackoverflow.com/questions/5195466/matplotlib-does-not-display-hatching-when-rendering-to-pdf
        hatch = next(hatchit)
        ax.fill_between(xgrid, err68up, err68down, facecolor='None', alpha=alpha,
                        edgecolor=color,
                        hatch=hatch,
                        zorder=1)
        if isinstance(stats, MCStats):
            errorstdup, errorstddown = stats.errorbarstd()
            ax.plot(xgrid, errorstdup, linestyle='--', color=color)
            ax.plot(xgrid, errorstddown, linestyle='--', color=color)
            label  = "%s ($68\%%$ c.l.+$1\sigma$)" % pdf.label
            outer = True
        else:
            outer = False
            label = "%s ($68\%%$ c.l.)" % pdf.label
        handle = plotutils.HandlerSpec(color=color, alpha=alpha,
                                               hatch=hatch,
                                               outer=outer)
        handles.append(handle)
        labels.append(label)

        return [err68down, err68up]

    def legend(self, flstate):
        return flstate.ax.legend(flstate.handles, flstate.labels,
                                 handler_map={plotutils.HandlerSpec:
                                             plotutils.ComposedHandler()
                                             }
                                 )

@figuregen
@_check_pdf_normalize_to
@check_scale('xscale', allow_none=True)
def plot_pdfs(pdfs, xplotting_grids, xscale:(str,type(None))=None,
                      normalize_to:(int,str,type(None))=None):
    """Plot the central value and the uncertainty of a list of pdfs as a
    function of x for a given value of Q. If normalize_to is given, norma.
    See the help for ``xplotting_grid`` for information on how to set basis,
    flavours and x ranges. Yields one figure per PDF flavour.

    normalize_to:  Either the name of one of the PDFs or its corresponding
    index in the list, starting from one, or None to plot absolute values.

    xscale: One of the matplotlib allowed scales. If undefined, it will be
    set based on the scale in xgrid, which should be used instead.

    """
    yield from BandPDFPlotter(pdfs, xplotting_grids, xscale, normalize_to)

class FLavoursPlotter(BandPDFPlotter):

    def setup_flavour(self, flstate):
        flstate.handles= self.handles
        flstate.labels= self.doesnothing
        flstate.hatchit= self.hatchit

    def __call__(self,):
        if not self.xplotting_grids:
            return

        self.handles=[]
        self.doesnothing=[]
        self.labels = []
        self.hatchit=plotutils.hatch_iter()


        basis = self.firstgrid.basis
        fig, ax = plt.subplots()
        ax.set_xlabel('x')
        ax.set_xscale(self.xscale)
        ax.set_title(f'{self.pdfs[0]} Q={self.Q : .1f} GeV  ')

        all_vals = []
        for flindex, fl in enumerate(self.firstgrid.flavours):

            parton_name = basis.elementlabel(fl)
            self.labels.append(f'${parton_name}$')
            flstate = FlavourState(flindex=flindex, fl=fl, fig=fig, ax=ax,
                                    parton_name=parton_name)
            self.setup_flavour(flstate)



            for pdf, grid in zip(self.pdfs, self.xplotting_grids):
                all_vals.append(np.atleast_2d(self.draw(pdf, grid, flstate)))

        plotutils.frame_center(ax, self.firstgrid.xgrid, np.concatenate(all_vals))
        ax.set_axisbelow(True)
        ax.set_xlim(self.firstgrid.xgrid[0])
        flstate.labels = self.labels
        self.legend(flstate)
        return fig

@figure
@check_scale('xscale', allow_none=True)
def plot_flavours(pdf, xplotting_grid, xscale:(str,type(None))=None,
                      normalize_to:(int,str,type(None))=None):
    """Plot the absolute central value and the uncertainty of all the flavours
    of a pdf as a function of x for a given value of Q.

    xscale: One of the matplotlib allowed scales. If undefined, it will be
    set based on the scale in xgrid, which should be used instead.

    """
    obj = FLavoursPlotter([pdf], [xplotting_grid], xscale, normalize_to=None)
    return obj()


@figuregen
def plot_smpdf(pdf, dataset, obs_pdf_correlations, mark_threshold:float=0.9):
    """
    Plot the correlations between the change in the observable and the change
    in the PDF in (x,fl) space.

    mark_threshold is the proportion of the maximum absolute correlation
    that will be used to mark the corresponding area in x in the
    background of the plot. The maximum absolute values are used for
    the comparison."""
    info = get_info(dataset)

    table = kitable(dataset, info)
    figby = sane_groupby_iter(table, info.figure_by)

    basis = obs_pdf_correlations.basis

    fullgrid = obs_pdf_correlations.grid_values

    fls = obs_pdf_correlations.flavours
    x = obs_pdf_correlations.xgrid
    nf = len(fls)

    plotting_var = info.get_xcol(table)

    #TODO: vmin vmax should be global or by figure?
    vmin,vmax = min(plotting_var), max(plotting_var)
    if info.x_scale == 'log':
        norm = mcolors.LogNorm(vmin, vmax)
    else:
        norm = mcolors.Normalize(vmin, vmax)
    #http://stackoverflow.com/a/11558629/1007990
    sm = cm.ScalarMappable(cmap=cm.viridis, norm=norm)
    sm._A = []

    for same_vals, fb in figby:
        grid = fullgrid[ np.asarray(fb.index),...]


        #Use the maximum absolute correlation for plotting purposes
        absgrid = np.max(np.abs(grid), axis=0)
        mark_mask = absgrid > np.max(absgrid)*mark_threshold

        label = info.group_label(same_vals, info.figure_by)
        #TODO: PY36ScalarMappable
        #TODO Improve title?
        title = "%s %s\n[%s]" % (info.dataset_label, '(%s)'%label if label else '' ,pdf.label)

        #Start plotting
        w,h = plt.rcParams["figure.figsize"]
        h*=2.5
        fig,axes = plt.subplots(nrows=nf ,sharex=True, figsize=(w,h), sharey=True)
        fig.suptitle(title)
        colors = sm.to_rgba(info.get_xcol(fb))
        for flindex, (ax, fl) in enumerate(zip(axes, fls)):
            for i,color in enumerate(colors):
                ax.plot(x, grid[i,flindex,:].T, color=color)


            flmask = mark_mask[flindex,:]
            ranges = split_ranges(x, flmask, filter_falses=True)
            for r in ranges:
                ax.axvspan(r[0], r[-1], color='#eeeeff')

            ax.set_ylabel("$%s$"%basis.elementlabel(fl))
            ax.set_xscale(_scale_from_grid(obs_pdf_correlations))
            ax.set_ylim(-1,1)
            ax.set_xlim(x[0], x[-1])
        ax.set_xlabel('$x$')
        #fig.subplots_adjust(hspace=0)

        fig.colorbar(sm, ax=axes.ravel().tolist(), label=info.xlabel,
                     aspect=100)
        #TODO: Fix title for this
        #fig.tight_layout()
        yield fig

@figure
def plot_obscorrs(corrpair_datasets, obs_obs_correlations, pdf):
    """NOTE: EXPERIMENTAL. Plot the correlation matrix between a pair of datasets."""
    fig, ax = plt.subplots()

    ds1, ds2 = corrpair_datasets
    in1,in2 = get_info(ds1), get_info(ds2)

    im = ax.imshow(obs_obs_correlations, cmap=cm.Spectral_r, vmin=-1, vmax=1)

    ax.set_ylabel(str(ds1))
    ax.set_xlabel(str(ds2))
    fig.colorbar(im, [ax])
    return fig

@figure
def plot_positivity(pdfs, positivity_predictions_for_pdfs, posdataset):
    """Plot the value of a positivity observable on a symlog scale as a
    function of the data point index."""
    fig,ax = plt.subplots()
    ax.axhline(0, color='red')
    offsets = plotutils.offset_xcentered(len(pdfs), ax)
    minscale = np.inf
    for i,(pdf, pred) in enumerate(zip(pdfs, positivity_predictions_for_pdfs)):
        cv = pred.central_value
        ax.errorbar(np.arange(len(cv)), cv, yerr=pred.std_error,
                    linestyle='--',
                    marker='s',
                    label=pdf.label, lw=0.5,transform=next(offsets))
        minscale = min(minscale, np.abs(np.min(cv)))
    ax.legend()
    ax.set_title(str(posdataset))
    ax.set_xlabel('idat')
    ax.set_ylabel('Observable Value')
    ax.set_yscale('symlog', linthreshy=minscale)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    return fig


@make_argcheck
def _check_display_cuts_requires_use_cuts(display_cuts, use_cuts):
    if display_cuts and not use_cuts:
        raise CheckError("The display_cuts option requires setting use_cuts to True")

@make_argcheck
def _check_marker_by(marker_by):
    markers = ('process type', 'experiment', 'dataset')
    if marker_by not in markers:
        raise CheckError("Unknown marker_by value", marker_by, markers)

#TODO: Right now this is hackish Could we turn it into a permanent interface?
@make_argcheck
def _check_highlights(experiments, highlight_datasets):
    if highlight_datasets:
        values = frozenset(highlight_datasets)
        names_set = {ds.name for experiment in experiments for ds in experiment }
        diff = values - names_set
        if diff:
            raise CheckError(f"The following highlight elements are "
                             "not dataset names: {diff}")
        return {'highlight_datasets': values}


@make_argcheck
def _check_aspect(aspect):
    aspects = ('landscape', 'portrait', 'square')
    if aspect not in aspects:
        raise CheckError(f"Unknown aspect {aspect}", aspect, aspects)

@figure
@_check_display_cuts_requires_use_cuts
@_check_marker_by
@_check_highlights
@_check_aspect
def plot_xq2(experiments_xq2map, use_cuts ,display_cuts:bool=True,
                 marker_by:str='process type', highlight_label:str='highlight',
                 highlight_datasets:(Sequence,type(None))=None,
                 aspect:str='landscape'):
    """Plot the (x,Q²) coverage based of the data based on some LO
    approximations. These are governed by the relevant kintransform.

    The representation of the filtered data depends on the `display_cuts` and
    `use_cuts` options:

     - If `use_cuts` is False, all the data will be plotted (and setting
    `display_cuts` to True is an error).

     - If `use_cuts` is True and `display_cuts` is False, the masked points
    will be ignored.

     - If `use_cuts` is True and `display_cuts` is True, the filtered points
    will be displaed and marked.

    The points are grouped according to the `marker_by` option. The possible
    values are: "process type", "experiment" or "dataset".
    """

    w,h = plt.rcParams["figure.figsize"]
    rescaling_factor = 1.6
    w *= rescaling_factor
    h *= rescaling_factor
    if aspect=='landscape':
        figsize = w, h
    elif aspect=='portrait':
        figsize = h, w
    elif aspect=='square':
        figsize = h, h
    else:
        raise ValueError(f"Unknown aspect {aspect}")
    fig, ax = plt.subplots(figsize=figsize)

    filteredx = []
    filteredq2 = []

    x = defaultdict(list)
    q2 = defaultdict(list)

    xh = defaultdict(list)
    q2h = defaultdict(list)

    if not highlight_datasets:
        highlight_datasets = set()

    def next_options():
        #Get the colors
        prop_settings = plt.rcParams['axes.prop_cycle']
        #Apparently calling the object gives us an infinite cycler
        settings_cycler = prop_settings()
        #So far, I don't understand how this is done with mpl "cycler"
        #objects, or wether  I like it. So far this is godd enough
        for  markeropts, settings in zip(plotutils.marker_iter_plot(), settings_cycler):
            #Override last with first
            options = {
                'linestyle': 'none',
                **markeropts,
                **settings,
            }
            yield options

    next_opts = next_options()
    key_options = {}

    for experiment, commondata, fitted, masked in experiments_xq2map:
        info = get_info(commondata)
        if marker_by == 'process type':
            key = info.process_description
        elif marker_by == 'experiment':
            key = str(experiment)
        elif marker_by == 'dataset':
            key = info.dataset_label
        else:
            raise ValueError('Unknown marker_by value')

        #TODO: This is an ugly check. Is there a way to do it with .setdefault
        # or defaultdict?
        if key not in key_options:
            key_options[key] = next(next_opts)

        if commondata.name in highlight_datasets:
            xdict = xh
            q2dict = q2h
        else:
            xdict = x
            q2dict = q2

        xdict[key].append(fitted[0])
        q2dict[key].append(fitted[1])
        if display_cuts:
            xdict[key].append(masked[0])
            q2dict[key].append(masked[1])
            filteredx.append(masked[0])
            filteredq2.append(masked[1])

    for key in key_options:
        if key in x:
            coords = np.concatenate(x[key]), np.concatenate(q2[key])
        else:
            #This is to get the label key
            coords = [], []
        ax.plot(*coords,
            label=key,
            markeredgewidth=1,
            markeredgecolor=None,
            **key_options[key],
        )

    #Iterate again so highlights are printed on top.
    for key in xh:
        ax.plot(np.concatenate(xh[key]), np.concatenate(q2h[key]),
            markeredgewidth=0.6,
            markeredgecolor="black",
            **key_options[key],
        )
    if xh:
        #Get legend key
        ax.plot([], [], marker='s', markeredgewidth=0.6, color='none',
            markersize=5,
            markeredgecolor="black", label= f'Black edge: {highlight_label}',
        )

    if display_cuts:
        ax.scatter(np.concatenate(filteredx), np.concatenate(filteredq2),
            marker='o',
            facecolors='none', edgecolor='red', s=40, lw=0.8, label="Cut"
        )

    ax.set_title("Kinematic coverage")
    ax.legend()
    ax.set_xlabel('x')
    ax.set_ylabel(r'$Q^2 (GeV^2)$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return fig


@figuregen
def plot_lumi1d(pdf, lumi_channel, lumigrid1d, sqrts:numbers.Real):
    """Plot PDF luminosities at a given center of mass energy.
    sqrts is the center of mass energy (GeV).
    """

    fig, ax = plt.subplots()
    mx = lumigrid1d.m
    gv = lumigrid1d.grid_values

    cv = gv.central_value()
    err = gv.std_error()

    ax.fill_between(mx, 1-err/cv, 1+err/cv, alpha=0.5)
    ax.plot(mx, cv/cv, label='%s' % pdf.label)
    ax.legend(loc='best')
    ax.set_xlabel('$M_{X}$ (GeV)')
    ax.set_xscale('log')
    ax.grid(False)
    ax.set_title("$%s$ luminosity\n%s - "
                 "$\\sqrt{s}=%.1f$ GeV" % (LUMI_CHANNELS[lumi_channel],
                                           pdf.label, sqrts))

    yield fig


#TODO: Move these to utils somewhere? Find better implementations?
def _reflect_matrl(mat, odd=False):
    """Reflect a matrix with positive values in the first axis to have the
    same balues for the nwgative axis. The first value is not reflected.

    If ``odd`` is set, the negative part will be multiplied by -1.

    """
    mat = np.asarray(mat)
    res = np.empty(shape=(mat.shape[0]*2-1, *mat.shape[1:]),dtype=mat.dtype)
    neglen = mat.shape[0]-1
    fact = -1 if odd else 1
    res[:neglen,...] = fact*mat[:0:-1,...]
    res[neglen:,...] = mat
    return res

def _reflect_matud(mat, odd=False):
    """Reflect a matrix with positive values in the second axis to have the
    same balues for the nwgative axis. The first value is not reflected.

    If ``odd`` is set, the negative part will be multiplied by -1.

    """
    mat = np.asarray(mat)
    res = np.empty(shape=(mat.shape[0], mat.shape[1]*2-1, *mat.shape[2:]),
        dtype=mat.dtype)
    neglen = mat.shape[1]-1
    fact = -1 if odd else 1
    res[:,:neglen,...] = fact*mat[:,:0:-1,...]
    res[:,neglen:,...] = mat
    return res


@figure
def plot_lumi2d(pdf, lumi_channel, lumigrid2d, sqrts,
                display_negative:bool=True):
    """Plot the absolute luminosity on a grid of invariant mass and
    rapidity for a given center of mass energy `sqrts`.
    The color scale is logarithmic.
    If `display_negative` is True, mark the negative values.

    The luminosity is calculated for positive rapidity, and reflected for
    negative rapidity for display purposes.

    """


    cmap = copy.copy(cm.viridis_r)
    cmap.set_bad("white", alpha=0)
    fig, ax = plt.subplots()
    gv = lumigrid2d.grid_values
    mat = gv.central_value()

    fig, ax = plt.subplots()

    mat = _reflect_matud(mat)
    y = _reflect_matrl(lumigrid2d.y, odd=True)
    masked_weights = np.ma.masked_invalid(mat, copy=False)

    #TODO: SymLogNorm is really the right thing to do here, but I can't be
    #bothered to make it work. Mostly the ticks around zero are completely
    #broken and looks like it takes a lot of fidlling wirh the mpl internals
    #to fix it.

    with np.errstate(invalid='ignore'):
        positive_mask = masked_weights>0
    linlim = np.nanpercentile(masked_weights[positive_mask],90)/1e5

    #norm = mcolors.SymLogNorm(linlim, vmin=None)

    norm = mcolors.LogNorm(vmin=linlim)
    with np.errstate(invalid='ignore'):
        masked_weights[masked_weights<linlim] = linlim

    mesh = ax.pcolormesh(y, lumigrid2d.m, masked_weights, cmap=cmap,
        shading='gouraud',
        linewidth=0,
        edgecolor='None',
        rasterized=True,
        norm=norm,
    )
    #Annoying code because mpl does the defaults horribly
    #loc = mticker.SymmetricalLogLocator(base=10, linthresh=linlim,)
    #loc.numticks = 5

    #fig.colorbar(mesh, ticks=loc)

    if display_negative:
        cmap_neg =  mcolors.ListedColormap(['red', (0,0,0,0)])
        neg_norm = mcolors.BoundaryNorm([0,0.5,1],2)
        ax.pcolormesh(y, lumigrid2d.m, positive_mask, cmap=cmap_neg,
            shading='gouraud',
            linewidth=0,
            edgecolor='None',
            rasterized=True,
            norm=neg_norm,
        )

    fig.colorbar(mesh, extend='min', label="Differential luminosity ($GeV^{-1}$)")
    ax.set_ylabel('$M_{X}$ (GeV)')
    ax.set_xlabel('y')
    ax.set_yscale('log')
    ax.grid(False)

    ax.set_title("$%s$ luminosity\n%s - "
             "$\\sqrt{s}=%.1f$ GeV" % (LUMI_CHANNELS[lumi_channel],
                     pdf.label, sqrts))

    return fig

@figure
def plot_lumi2d_uncertainty(pdf, lumi_channel, lumigrid2d, sqrts:numbers.Real):
    """
    Plot 2D luminosity unciertainty plot at a given center of mass energy.
    Porting code from https://github.com/scarrazza/lumi2d.

    The luminosity is calculated for positive rapidity, and reflected for
    negative rapidity for display purposes.
    """
    norm = mcolors.LogNorm(vmin=1, vmax=50)
    cmap = copy.copy(cm.viridis_r)
    cmap.set_bad("white", alpha=0)

    grid = lumigrid2d
    channel = lumi_channel

    gv = grid.grid_values
    mat = gv.std_error()/np.abs(gv.central_value())*100

    fig, ax = plt.subplots()

    mat = _reflect_matud(mat)
    y = _reflect_matrl(grid.y, odd=True)


    masked_weights = np.ma.masked_invalid(mat, copy=False)

    mesh = ax.pcolormesh(y, grid.m, masked_weights, norm=norm, cmap=cmap,
                         shading='gouraud',
                         linewidth=0,
                         edgecolor='None',
                         rasterized=True)

    # some extra options
    extup = np.nanmax(masked_weights) > 50
    extdown = np.nanmin(masked_weights) < 1

    #TODO: Wrap this somewhere
    if extup:
        if extdown:
            extend = 'both'
        else:
            extend = 'max'
    elif extdown:
        extend = 'min'
    else:
        extend = None

    fig.colorbar(mesh, label="Relative uncertainty (%)",
        ticks=[1,5,10,25,50], format='%.0f', extend=extend)
    ax.set_yscale('log')
    ax.set_title("Relative uncertainty for $%s$-luminosity\n%s - "
                 "$\\sqrt{s}=%.1f$ GeV" % (LUMI_CHANNELS[channel],
                         pdf.label, sqrts))
    ax.set_ylabel('$M_{X}$ (GeV)')
    ax.set_xlabel('y')
    ax.grid(False)

    return fig
