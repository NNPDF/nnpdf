# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:19:52 2016

@author: Zahari Kassabov
"""
from collections import OrderedDict
import logging

import numpy as np
import pandas as pd

from NNPDF import ThPredictions
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
        covmat = loaded_exp.get_covmat()
        df.ix[name][name] = covmat
    return df

@table
def experiments_invcovmat(experiments, experiments_index):
    """Export the inverse covariance matrix. See ``experiments_covmat``."""
    data = np.zeros((len(experiments_index),len(experiments_index)))
    df = pd.DataFrame(data, index=experiments_index, columns=experiments_index)
    for experiment in experiments:
        name = experiment.name
        loaded_exp = experiment.load()
        covmat = loaded_exp.get_invcovmat()
        df.ix[name][name] = covmat
    return df




def results(dataset:DataSetSpec, pdf:PDF):
    """Tuple of data and theory results for a single pdf.
    The theory is specified as part of the dataset
    (as a result of the C++ code layout)."""

    data = dataset.load()

    return DataResult(data), ThPredictionsResult.from_convolution(pdf, dataset,
                                                 loaded_data=data)

#It's better to duplicate a few lines than to complicate the logic of
#``results`` to support this
def pdf_results(dataset:DataSetSpec, pdfs:list):
    """Return a list of results, the first for the data and the rest for
    each of the PDFs."""

    data = dataset.load()

    th_results = []
    for pdf in pdfs:
        th_result = ThPredictionsResult.from_convolution(pdf, dataset,
                                                         loaded_data=data)
        th_results.append(th_result)


    return (DataResult(data), *th_results)

@require_one('pdfs', 'pdf')
def one_or_more_results(dataset:DataSetSpec, pdfs:list=None, pdf:PDF=None):
    if pdf:
        return results(dataset, pdf)
    else:
        return pdf_results(dataset, pdfs)
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
