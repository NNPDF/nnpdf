# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:19:52 2016

@author: Zahari Kassabov
"""
from __future__ import generator_stop

from collections import OrderedDict, namedtuple, Sequence
import itertools
import logging

import numpy as np
import scipy.linalg as la
import pandas as pd

from NNPDF import ThPredictions, CommonData
from NNPDF.experiments import Experiment
from reportengine.checks import require_one, remove_outer
from reportengine.table import table
from reportengine import collect

from validphys.checks import make_check, CheckError
from validphys.core import DataSetSpec, PDF, ExperimentSpec

log = logging.getLogger(__name__)



class Result: pass


#TODO: Eventually,only one of (NNPDFDataResult, StatsResult) should survive
class NNPDFDataResult(Result):
    """A result fills its values from a libnnpf data object"""
    def __init__(self, dataobj):
        self._central_value = dataobj.get_cv()

    @property
    def central_value(self):
        return self._central_value

    def __len__(self):
        return len(self.central_value)

class StatsResult(Result):
    def __init__(self, stats):
        self.stats = stats

    @property
    def central_value(self):
        return self.stats.central_value()

    @property
    def std_error(self):
        return self.stats.std_error()




class DataResult(NNPDFDataResult):

    def __init__(self, dataobj):
        super().__init__(dataobj)
        self._covmat = dataobj.get_covmat()
        self._sqrtcovmat = dataobj.get_sqrtcovmat()

    @property
    def label(self):
        return "CommonData"

    @property
    def std_error(self):
        return np.sqrt(np.diag(self.covmat))

    @property
    def covmat(self):
        return self._covmat

    @property
    def sqrtcovmat(self):
        """Lower part of the Cholesky decomposition"""
        return self._sqrtcovmat


class ThPredictionsResult(NNPDFDataResult):

    def __init__(self, dataobj, stats_class, label=None):
        self.stats_class = stats_class
        self.label = label
        self._std_error = dataobj.get_error()
        self._rawdata = dataobj.get_data()
        super().__init__(dataobj)

    @property
    def std_error(self):
        return self._std_error

    @staticmethod
    def make_label(pdf, dataset):
        """Deduce a reasonsble label for the result based on pdf and dataspec"""
        th = dataset.thspec
        if hasattr(pdf,'label'):
            if hasattr(th, 'label'):
                label = ' '.join((pdf.label, th.label))
            else:
                label = pdf.label
        elif hasattr(th, 'label'):
            label = th.label
        else:
            label = ('%s@<Theory %s>' % (pdf, th.id))
        return label


    @classmethod
    def from_convolution(cls, pdf, dataset, loaded_pdf=None, loaded_data=None):
        if loaded_pdf  is None:
            loaded_pdf = pdf.load()
        if loaded_data is None:
            loaded_data = dataset.load()
        th_predictions = ThPredictions(loaded_pdf, loaded_data)


        label = cls.make_label(pdf, dataset)


        return cls(th_predictions, pdf.stats_class, label)

class PositivityResult(StatsResult):
    @classmethod
    def from_convolution(cls, pdf, posset):
        loaded_pdf = pdf.load()
        loaded_pos = posset.load()
        data = loaded_pos.GetPredictions(loaded_pdf)
        stats = pdf.stats_class(data.T)
        return cls(stats)

    @property
    def rawdata(self):
        return self.stats.data



def experiments_index(experiments):
    """Return an empy dataframe with index
       per experiment per dataset per point"""
    records = []
    for exp_index, experiment in enumerate(experiments):
        loaded_exp = experiment.load()
        set_lens = [len(loaded_exp.GetSet(i)) for i in
                    range(len(experiment.datasets))]
        #TODO: This code is very ugly and slow...
        cum_sum = [sum(set_lens[:i+1]) for i in range(len(set_lens))]
        curr_ds_domain = iter(enumerate(cum_sum))
        index_offset = 0
        ds_id, curr_ds_len = next(curr_ds_domain)
        for index in range(cum_sum[-1]):
            if index >= curr_ds_len:
                index_offset = curr_ds_len
                ds_id, curr_ds_len = next(curr_ds_domain)
            dataset = experiment.datasets[ds_id]

            records.append(OrderedDict([
                                 ('experiment', str(experiment.name)),
                                 ('dataset', str(dataset.name)),
                                 ('id', index - index_offset),
                                  ]))

    columns = ['experiment', 'dataset', 'id']
    df = pd.DataFrame(records, columns=columns)
    df.set_index(columns, inplace=True)
    return df.index

#TODO: Find a good name for this
def t0set(use_t0=False, t0pdfset=None):
    """A conveninece provider to not have to drag the two
    related parameters everywhere"""
    if use_t0:
        return t0pdfset
    else:
        return None



@table
def experiment_result_table(experiments, pdf, experiments_index):
    """Generate a table containing the data central value, the central prediction,
    and the prediction for each PDF member."""

    result_records = []
    for exp_index, experiment in enumerate(experiments):
        loaded_exp = experiment.load()



        data_result = DataResult(loaded_exp)
        th_result = ThPredictionsResult.from_convolution(pdf, experiment,
                                                         loaded_data=loaded_exp)


        for index in range(len(data_result.central_value)):
            replicas = (('rep_%05d'%(i+1), th_result._rawdata[index,i]) for
                        i in range(th_result._rawdata.shape[1]))

            result_records.append(OrderedDict([
                                 ('data_central', data_result.central_value[index]),
                                 ('theory_central', th_result.central_value[index]),
                                  *replicas
                                 ]))

    if not result_records:
        log.warn("Empty records for experiment results")
        return pd.DataFrame()
    df =  pd.DataFrame(result_records, columns=result_records[0].keys(),
                       index=experiments_index)
    return df

@table
def experiments_covmat(experiments, experiments_index):
    """Export the covariance matrix for the experiments. It exports the full
    (symmetric) matrix, with the 3 first rows and columns being:

        - experiment name

        - dataset name

        - index of the point within the dataset.
    """
    data = np.zeros((len(experiments_index),len(experiments_index)))
    df = pd.DataFrame(data, index=experiments_index, columns=experiments_index)
    for experiment in experiments:
        name = experiment.name
        loaded_exp = experiment.load()
        mat = loaded_exp.get_covmat()
        df.loc[[name],[name]] = mat
    return df

@table
def experiments_invcovmat(experiments, experiments_index):
    """Compute and export the inverse covariance matrix.
    Note that this inverts the matrices with the LU method which is
    suboptimal."""
    data = np.zeros((len(experiments_index),len(experiments_index)))
    df = pd.DataFrame(data, index=experiments_index, columns=experiments_index)
    for experiment in experiments:
        name = experiment.name
        loaded_exp = experiment.load()
        #Improve this inversion if this method tuns out to be important
        mat = la.inv(loaded_exp.get_covmat())
        df.loc[[name],[name]] = mat
    return df

@table
def closure_pseudodata_replicas(experiments, pdf, nclosure:int,
                                experiments_index, nnoisy:int=0):
    """Generate closure pseudodata replicas from the given pdf.

    nclosure: Number of Level 1 pseudodata replicas.

    nnoisy:   Number of Level 2 replicas generated out of each pseudodata replica.

    The columns of the table are of the form (clos_0, noise_0_n0 ..., clos_1, ...)
    """

    #TODO: Do this somewhere else
    from NNPDF import randomgenerator
    randomgenerator.RandomGenerator.InitRNG(0,0)
    data = np.zeros((len(experiments_index), nclosure*(1+nnoisy)))

    cols = []
    for i in range(nclosure):
        cols += ['clos_%04d'%i, *['noise_%04d_%04d'%(i,j) for j in range(nnoisy)]]


    loaded_pdf = pdf.load()

    for exp in experiments:
        #Since we are going to modify the experiments, we copy them
        #(and work on the copies) to avoid all
        #sorts of weirdness with other providers. We don't want this to interact
        #with ExperimentSpec at all, because it could do funny things with the
        #cache when calling load(). We need to copy this yet again, for each
        # of the noisy replicas.
        closure_exp = Experiment(exp.load())

        #TODO: This is probably computed somewhere else... All this code is
        #very error prone.
        #The predictions are for the unmodified experiment.
        predictions = [ThPredictions(loaded_pdf, d.load()) for d in exp]


        exp_location = experiments_index.get_loc(closure_exp.GetExpName())

        index = itertools.count()
        for i in range(nclosure):
            #Generate predictions with experimental noise, a different for
            #each closure set.
            closure_exp.MakeClosure(predictions, True)
            data[exp_location, next(index)] = closure_exp.get_cv()
            for j in range(nnoisy):
                #If we don't copy, we generate noise on top of the noise,
                #which is not what we want.
                replica_exp = Experiment(closure_exp)
                replica_exp.MakeReplica()

                data[exp_location, next(index)] = replica_exp.get_cv()


    df = pd.DataFrame(data, index=experiments_index,
                      columns=cols)

    return df


def results(dataset:(DataSetSpec), pdf:PDF, t0set:(PDF, type(None))=None):
    """Tuple of data and theory results for a single pdf.
    The theory is specified as part of the dataset.
    An experiment is also allowed.
    (as a result of the C++ code layout)."""

    data = dataset.load()

    if t0set:
        #Copy data to avoid chaos
        data = type(data)(data)
        log.debug("Setting T0 predictions for %s" % dataset)
        data.SetT0(t0set.load())

    return DataResult(data), ThPredictionsResult.from_convolution(pdf, dataset,
                                                 loaded_data=data)

def experiment_results(experiment, pdf:PDF, t0set:(PDF, type(None))=None):
    return results(experiment, pdf, t0set)

#It's better to duplicate a few lines than to complicate the logic of
#``results`` to support this.
#TODO: The above comment doesn't make sense after adding T0
def pdf_results(dataset:(DataSetSpec,  ExperimentSpec), pdfs:Sequence, t0set:(PDF, type(None))):
    """Return a list of results, the first for the data and the rest for
    each of the PDFs."""

    data = dataset.load()

    if t0set:
        #Copy data to avoid chaos
        data = type(data)(data)
        log.debug("Setting T0 predictions for %s" % dataset)
        data.SetT0(t0set.load())

    th_results = []
    for pdf in pdfs:
        th_result = ThPredictionsResult.from_convolution(pdf, dataset,
                                                         loaded_data=data)
        th_results.append(th_result)


    return (DataResult(data), *th_results)

@require_one('pdfs', 'pdf')
@remove_outer('pdfs', 'pdf')
def one_or_more_results(dataset:(DataSetSpec, ExperimentSpec),
                        pdfs:(type(None), Sequence)=None,
                        pdf:(type(None), PDF)=None,
                        t0set:(PDF, type(None))=None):
    """Generate a list of results, where the first element is the data values,
    and the next is either the prediction for pdf or for each of the pdfs.
    Which of the two is selected intelligently depending on the namespace,
    when executing as an action."""
    if pdf:
        return results(dataset, pdf, t0set)
    else:
        return pdf_results(dataset, pdfs, t0set)
    raise ValueError("Either 'pdf' or 'pdfs' is required")


def _calc_chi2(sqrtcov, diffs):
    """Elementary function to compute the chi², given a Cholesky decomposed
    lower triangular part and a vector of differences"""
    #Note la.cho_solve doesn't really improve things here
    vec = la.solve_triangular(sqrtcov, diffs, lower=True)
    #This sums up the result for the chi² for any input shape.
    #Sum the squares over the first dimension and leave the others alone
    return np.einsum('i...,i...->...', vec,vec)

def _all_chi2(results):
    """Return the chi² for all elements in the result"""
    data_result, th_result = results
    diffs = th_result._rawdata - data_result.central_value[:,np.newaxis]
    return _calc_chi2(sqrtcov=data_result.sqrtcovmat, diffs=diffs)


Chi2Data = namedtuple('Chi2data', ('replica_result', 'central_result', 'ndata'))

def abs_chi2_data(results):
    """Return a tuple (member_chi², central_chi², numpoints)"""
    data_result, th_result = results

    chi2s = _all_chi2(results)

    central_diff = th_result.central_value - data_result.central_value
    central_result = _calc_chi2(data_result.sqrtcovmat, central_diff)

    return Chi2Data(th_result.stats_class(chi2s[:, np.newaxis]),
                    central_result, len(data_result))

def abs_chi2_data_experiment(experiment_results):
    return abs_chi2_data(experiment_results)

def _chs_per_replica(chs):
    th, _, l = chs
    return th.data.ravel()/l


@table
def experiments_chi2_table(experiments, pdf, experiments_chi2,
                           each_dataset_chi2):
    dschi2 = iter(each_dataset_chi2)
    records = []
    for experiment, expres in zip(experiments, experiments_chi2):
        stats = chi2_stats(expres)
        stats['experiment'] = experiment.name
        records.append(stats)
        for dataset, dsres in zip(experiment, dschi2):
            stats = chi2_stats(dsres)
            stats['experiment'] = dataset.name
            records.append(stats)
    return pd.DataFrame(records)

@table
def correlate_bad_experiments(experiments, replica_data, pdf):
    """Generate a table for each replica with entries
    ("Replica_mean_chi2", "Worst_dataset","Worst_dataset_chi2")."""
    datasets = [ds for exp in experiments for ds in exp.datasets]
    mchi2 = [0.5*(val.training + val.validation) for val in replica_data]

    chs = [_chs_per_replica(abs_chi2_data(results(ds, pdf))) for ds in datasets]
    worst_indexes = np.argmax(chs, axis=0)
    mchi2 = np.mean(chs, axis=0)
    print(worst_indexes)
    worst_values = np.max(chs, axis=0)
    worst_ds = [datasets[i].name for i in worst_indexes]
    v = np.array([mchi2, worst_ds, worst_values])
    print(v)
    df = pd.DataFrame(v.T,
                      index=np.arange(1, len(pdf)),
                      columns=["Replica_mean_chi2", "Worst_dataset",
                      "Worst_dataset_chi2"])
    df.sort_values(df.columns[0], inplace=True, ascending=False)
    return df

#TODO: Compute results
@table
def perreplica_chi2_table(experiments, pdf):
    """Chi² per point for each replica for each experiment.
    Also outputs the total chi² per replica."""

    chs = [abs_chi2_data(results(exp, pdf)) for exp in experiments]

    total_chis = np.zeros((len(experiments) + 1, len(pdf)))
    ls = []
    for i,ch in enumerate(chs, 1):
        th, central, l = ch
        total_chis[i]= [central, *th.error_members()]
        ls.append(l)

    #total_chis/=total_l
    total_chis[0] = np.sum(total_chis[1:,:], axis=0)
    total_chis[0]/= np.sum(ls)
    total_chis[1:,:]/= np.array(ls)[:, np.newaxis]

    return pd.DataFrame(total_chis.T, columns = ['Total', *[exp.name for exp in experiments]])

@make_check
def _assert_use_cuts_true(ns, **kwargs):
    if not ns['use_cuts']:
        raise CheckError("use_cuts needs to be True for this action.")

@_assert_use_cuts_true
@table
def closure_shifts(experiments_index, fit, use_cuts, experiments):
    """Save the differenve between the fitted data and the real commondata
    values.

    Actually shifts is what should be saved in the first place, rather than
    thi confusing fiddling with Commondata, but until we can implement this at
    the C++ level, we just dave it here.
    """
    name, fitpath = fit
    result = np.zeros(len(experiments_index))
    for experiment in experiments:
        for dataset in experiment:
            dspath = fitpath/'filter'/dataset.name
            cdpath = dspath/("DATA_" + dataset.name + ".dat")
            try:
                syspath = next( (dspath/'systypes').glob('*.dat'))
            except StopIteration as e:
                raise FileNotFoundError("No systype "
                "file found in filter folder %s" % (dspath/'systypes')) from e
            cd = CommonData.ReadFile(str(cdpath), str(syspath))
            loc = experiments_index.get_loc((experiment.name, dataset.name))
            result[loc] = cd.get_cv() - dataset.load().get_cv()
    return pd.DataFrame(result, index=experiments_index)




def positivity_predictions(pdf, positivityset):
    return PositivityResult.from_convolution(pdf, positivityset)

#TODO: Replace this with reportengine.collect
def possets_predictions(pdf, posdatasets):
    return [positivity_predictions(pdf, pos) for pos in posdatasets]

def count_negative_points(possets_predictions):
    return np.sum([(r.rawdata < 0).sum(axis=1) for r in possets_predictions], axis=0)




chi2_stat_labels = {
    'central_mean': r'$<\chi^2_{0}>_{data}$',
    'npoints': r'$N_{data}$',
    'perreplica_mean': r'$\left< \chi^2 \right>_{rep,data}$',
    'perreplica_std': r'$\left<std_{rep}(\chi^2)\right>_{data}$',
    'chi2_per_data': r'$\frac{\chi^2}{N_{data}}$'
}

def chi2_stats(abs_chi2_data):
    """Compute severa estimators from the chi²:

     - central_mean

     - npoints

     - perreplica_mean

     - perreplica_std

     - chi2_per_data
    """
    rep_data, central_result, npoints = abs_chi2_data
    m = central_result.mean()
    rep_mean = rep_data.central_value().mean()
    return OrderedDict([
            ('central_mean'        ,  m),
            ('npoints'             , npoints),
            ('chi2_per_data', m/npoints),
            ('perreplica_mean', rep_mean),
            ('perreplica_std',  rep_data.std_error().mean()),
           ])

experiments_results = collect(experiment_results, ('experiments',))
each_dataset_results = collect(results, ('experiments', 'experiment'))

experiments_chi2 = collect(abs_chi2_data_experiment, ('experiments',))
each_dataset_chi2 = collect(abs_chi2_data, ('experiments', 'experiment'))
