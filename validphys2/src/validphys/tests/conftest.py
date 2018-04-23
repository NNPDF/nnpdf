#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 18:05:51 2018

@author: zah
"""
import pytest
from validphys.loader import FallbackLoader as Loader
from validphys.core import ExperimentSpec
from validphys import results

#Fortunately py.test works much like reportengine and providers are
#connected by argument names.
@pytest.fixture(scope='module')
def data():
    l = Loader()
    ds = l.check_dataset(name='NMC', theoryid=162, use_cuts=False)
    exp = ExperimentSpec('NMC Experiment', [ds])
    pdf = l.check_pdf("NNPDF31_nnlo_as_0118")
    exps = [exp]
    return pdf, exps

@pytest.fixture(scope='module')
def convolution_results(data):
    pdf, exps = data
    return [results.experiment_results(exp, pdf, pdf) for exp in exps]

@pytest.fixture
def dataset_t0_convolution_results(data):
    pdf, exps = data
    ds = exps[0].datasets[0]
    return results.results(ds, pdf, t0set=pdf)

@pytest.fixture(scope='module')
def chi2data(convolution_results):
    return [results.abs_chi2_data_experiment(r) for r in convolution_results]