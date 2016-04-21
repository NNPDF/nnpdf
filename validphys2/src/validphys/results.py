# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:19:52 2016

@author: Zahari Kassabov
"""
from collections import OrderedDict
import functools

import numpy as np

from NNPDF import ThPredictions
from reportengine.checks import require_one

from validphys.core import DataSetSpec, PDF

class Result:
    def __init__(self, dataobj):
        self.dataobj = dataobj

    @property
    @functools.lru_cache()
    def std_error(self):
        return np.sqrt(np.diag(self.covmat))

    @property
    @functools.lru_cache()
    def central_value(self):
        return self.dataobj.get_cv()

    def __len__(self):
        return len(self.dataobj)

    def __getattr__(self, attr):
        return getattr(self.dataobj, attr)


class DataResult(Result):

    @property
    def label(self):
        return "CommonData"

    @property
    @functools.lru_cache()
    def covmat(self):
        return self.dataobj.get_covmat()

    @property
    @functools.lru_cache()
    def invcovmat(self):
        return self.dataobj.get_invcovmat()


class ThPredictionsResult(Result):

    def __init__(self, thlabel, dataobj, stats_class):
        self.thlabel = thlabel
        self.stats_class = stats_class
        super().__init__(dataobj)

    @property
    @functools.lru_cache()
    def std_error(self):
        return self.dataobj.get_error()

    @property
    @functools.lru_cache()
    def _rawdata(self):
        return self.dataobj.get_data()

    @property
    def data(self):
        return self.stats_class(self._rawdata)

    @property
    def label(self):
        return "<Theory %s>@%s" % (self.thlabel, self.dataobj.GetPDFName())

def results(dataset:DataSetSpec, pdf:PDF):
    """Tuple of data and theory results for a single pdf.
    The theory is specified as part of the dataset
    (as a result of the C++ code layout)."""

    nnpdf_pdf = pdf.load()
    data = dataset.load()
    th_predictions = ThPredictions(nnpdf_pdf, data)

    stats = pdf.stats_class

    thlabel, thpath = dataset.thspec
    return DataResult(data), ThPredictionsResult(thlabel, th_predictions,
                                                 stats)

#It's better to duplicate a few lines than to complicate the logic of
#``results`` to support this
def pdf_results(dataset:DataSetSpec, pdfs:list):
    """Return a list of results, the first for the data and the rest for
    each of the PDFs."""

    data = dataset.load()
    thlabel, thpath = dataset.thspec

    th_results = []
    for pdf in pdfs:
        nnpdf_pdf = pdf.load()
        th_predictions = ThPredictions(nnpdf_pdf, data)
        stats = pdf.stats_class
        th_result = ThPredictionsResult(thlabel, th_predictions,
                                                 stats)
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

    return th_result.stats_class(result[:, np.newaxis])

def chi2_stats(chi2_data):
    return OrderedDict([
            (r'$\left< \chi^2 \right>$'  , chi2_data.central_value().mean() ),
            (r'std($\chi^2$)'            , chi2_data.std_error().mean()     ),
           ])
