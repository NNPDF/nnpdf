# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:19:52 2016

@author: Zahari Kassabov
"""
from collections import OrderedDict
import itertools
import logging

import numpy as np
import pandas as pd

from NNPDF import ThPredictions
from NNPDF.experiments import Experiment
from reportengine.checks import require_one
from reportengine.table import table

from validphys.core import DataSetSpec, PDF

log = logging.getLogger(__name__)

#TODO: Is this abstraction any useful?
class Result:
    def __init__(self, dataobj):
        self._central_value = dataobj.get_cv()

    @property
    def std_error(self):
        return np.sqrt(np.diag(self.covmat))

    @property
    def central_value(self):
        return self._central_value

    def __len__(self):
        return len(self.central_value)


class DataResult(Result):

    def __init__(self, dataobj):
        super().__init__(dataobj)
        self._covmat = dataobj.get_covmat()
        self._invcovmat = dataobj.get_invcovmat()

    @property
    def label(self):
        return "CommonData"

    @property
    def covmat(self):
        return self._covmat

    @property
    def invcovmat(self):
        return self._invcovmat


class ThPredictionsResult(Result):

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
                label = ' '.join((pdf.label, dataset.label))
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
    """Export the inverse covariance matrix. See ``experiments_covmat``."""
    data = np.zeros((len(experiments_index),len(experiments_index)))
    df = pd.DataFrame(data, index=experiments_index, columns=experiments_index)
    for experiment in experiments:
        name = experiment.name
        loaded_exp = experiment.load()
        mat = loaded_exp.get_invcovmat()
        df.loc[[name],[name]] = mat
    return df

@table
def closure_pseudodata_replicas(experiments, pdf, nclosure:int,
                                experiments_index, nnoisy:int=0):
    """Generate closure pseudodata replicas from the given pdf.

    nclosure: Number of Level 1 pseudodata replicas.
    nnoisy:   Number of Level 2 replicas generated out of each pseudodata replica.

    The columns of the table are of the form (clos_0, noise_0_n0 ... ,clos_1, ...)
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

#TODO: don't do computations here
@table
def experiment_chi2_table(experiments, pdf):

    records = []
    for experiment in experiments:

        #TODO: This is probably computed somewhere else....
        r = data, theory = results(experiment, pdf)
        data = chi2_data(r)
        stats = chi2_stats(data)
        stats['experiment'] = experiment.name
        records.append(stats)

        for dataset in experiment:
            r = data, theory = results(dataset, pdf)
            data = chi2_data(r)
            stats = chi2_stats(data)
            stats['experiment'] = dataset.name
            records.append(stats)
    return pd.DataFrame(records)




def results(dataset:(DataSetSpec), pdf:PDF, t0set:(PDF, type(None))):
    """Tuple of data and theory results for a single pdf.
    The theory is specified as part of the dataset
    (as a result of the C++ code layout)."""

    data = dataset.load()

    if t0set:
        #Copy data to avoid chaos
        data = type(data)(data)
        t0_preds = ThPredictions(t0set.load(), data)
        log.debug("Setting T0 predictions for %s" % dataset)
        data.SetT0(t0_preds)

    return DataResult(data), ThPredictionsResult.from_convolution(pdf, dataset,
                                                 loaded_data=data)

#It's better to duplicate a few lines than to complicate the logic of
#``results`` to support this.
def pdf_results(dataset:DataSetSpec, pdfs:list):
    """Return a list of results, the first for the data and the rest for
    each of the PDFs."""

    data = dataset.load()

    if t0set:
        #Copy data to avoid chaos
        data = type(data)(data)
        t0_preds = ThPredictions(t0set.load(), data)
        log.debug("Setting T0 predictions for %s" % dataset)
        data.SetT0(t0_preds)

    th_results = []
    for pdf in pdfs:
        th_result = ThPredictionsResult.from_convolution(pdf, dataset,
                                                         loaded_data=data)
        th_results.append(th_result)


    return (DataResult(data), *th_results)

@require_one('pdfs', 'pdf')
def one_or_more_results(dataset:DataSetSpec, pdfs:list=None, pdf:PDF=None,
                        t0set:(PDF, type(None))=None):
    if pdf:
        return results(dataset, pdf, t0set)
    else:
        return pdf_results(dataset, pdfs, t0set)
    raise ValueError("Either 'pdf' or 'pdfs' is required")



def chi2_data(results):
    data_result, th_result = results
    diffs = th_result._rawdata.T - data_result.central_value
    #chi²_i = diff_ij @ invcov_jk @ diff_ki
    result =  np.einsum('ij, jk, ik -> i',
                     diffs, data_result.invcovmat, diffs)/len(data_result)

    print(result)

    central_diff = th_result.central_value - data_result.central_value

    central_result = (central_diff@data_result.invcovmat@central_diff)/len(data_result)


    return (th_result.stats_class(result[:, np.newaxis]), central_result)


chi2_stat_labels = {
    'central_mean': r'$<\chi^2_{0}>_{data}$',
    'central_std' : r'$std_{data}(chi^2_{0})$',
    'perreplica_mean': r'$\left< \chi^2 \right>_{rep,data}$',
    'perreplica_std': r'$\left<std_{rep}(\chi^2)\right>_{data}$',
}

def chi2_stats(chi2_data):
    rep_data, central_result = chi2_data

    return OrderedDict([
            ('central_mean'        ,  central_result.mean()),
            ('central_std'        ,  central_result.std()),
            ('perreplica_mean', rep_data.central_value().mean()),
            ('perreplica_std',  rep_data.std_error().mean()),
           ])
