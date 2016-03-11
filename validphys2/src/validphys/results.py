# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:19:52 2016

@author: Zahari Kassabov
"""
from collections import OrderedDict
import functools

import numpy as np

from NNPDF import CommonData, FKTable, ThPredictions
from NNPDF.fkset import FKSet
from NNPDF.dataset import DataSet

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
    cdpath,syspth = dataset.commondata
    cd = CommonData.ReadFile(str(cdpath), str(syspth))
    thlabel, thpath = dataset.thspec

    fktable = FKTable(str(dataset.fkpath), [str(factor) for factor in dataset.cfac])
    #IMPORTANT: We need to tell the python garbage collector to NOT free the
    #memory owned by the FKTable on garbage collection.
    #TODO: Do this automatically
    fktable.thisown = 0
    fkset = FKSet(FKSet.parseOperator("NULL"), [fktable])

    data = DataSet(cd, fkset)

    if dataset.cuts is not None:
        #ugly need to convert from numpy.int64 to int, so we can pass
        #it happily to the vector to the SWIG wrapper.
        #Do not do this (or find how to enable in SWIG):
        #data = DataSet(data, list(dataset.cuts))
        intmask = [int(ele) for ele in dataset.cuts]
        data = DataSet(data, intmask)

    nnpdf_pdf = pdf.load()
    th_predictions = ThPredictions(nnpdf_pdf, data)

    stats = pdf.stats_class


    return DataResult(data), ThPredictionsResult(thlabel, th_predictions,
                                                 stats)

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
