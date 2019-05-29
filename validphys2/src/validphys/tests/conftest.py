"""
conftest.py

Pytest fixtures.
"""
import pathlib

import pytest
from hypothesis import settings

from validphys.loader import FallbackLoader as Loader
from validphys.core import ExperimentSpec
from validphys import results

#Adding this here to change the time of deadline from default (200ms) to 1000ms
settings.register_profile("extratime", deadline=1000)
settings.load_profile("extratime")

#Fortunately py.test works much like reportengine and providers are
#connected by argument names.
@pytest.fixture
def tmp(tmpdir):
    """A tempdir that is manipulated like pathlib Paths"""
    return pathlib.Path(tmpdir)

@pytest.fixture(scope='module')
def data():
    l = Loader()
    dataset_inputs = [{'name': 'NMC'},
                      {'name':'ATLASTTBARTOT', 'cfac':['QCD']},
                      {'name':'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sysnum':10}]
    ds = [l.check_dataset(**x, theoryid=162, cuts=None) for x in dataset_inputs]
    exps = [ExperimentSpec(x.name, [x]) for x in ds]
    pdf = l.check_pdf("NNPDF31_nnlo_as_0118")
    return pdf, exps

@pytest.fixture(scope='module')
def exps_covariance_matrices(data):
    """produces a list of covariance matrix outputs for each experiment"""
    _, exps = data
    covs = [results.experiment_covariance_matrix(exp, False, None) for exp in exps]
    return covs

@pytest.fixture(scope='module')
def t0_exps_covariance_matrices(data):
    """produces a list of covariance matrix outputs for each experiment"""
    pdf, exps = data
    covs = [results.experiment_covariance_matrix(exp, False, pdf) for exp in exps]
    return covs

def convolution_results_implement(data):
    pdf, exps = data
    #no theory covmat here
    covs = [results.experiment_covariance_matrix(exp, False, pdf) for exp in exps]
    return [results.experiment_results(exp, pdf, cov) for exp, cov in zip(exps, covs)]

@pytest.fixture(scope='module')
def theory_data():
    l = Loader()
    names = ['NMC', 'ATLASTTBARTOT']
    theoryids = [163, 180, 173]
    ds1 = [l.check_dataset(name=names[0], theoryid=x, cuts=None) for x in theoryids]
    ds2 = [l.check_dataset(name=names[1], theoryid=x, cuts=None) for x in theoryids]
    exp1 = [ExperimentSpec(x.name, [x]) for x in ds1]
    exp2 = [ExperimentSpec(x.name, [x]) for x in ds2]
    exps_by_theoryid = [exp1, exp2]
    exps_central_theory = [exp1[0], exp2[0]]
    pdf = l.check_pdf("NNPDF31_nnlo_as_0118")
    return pdf, exps_by_theoryid, exps_central_theory, theoryids

@pytest.fixture(scope='module')
def convolution_results(data):
    return convolution_results_implement(data)

@pytest.fixture
def dataset_t0_convolution_results(data):
    pdf, exps = data
    ds = [x.datasets[0] for x in exps]
    covs = [results.covariance_matrix(x, False, pdf) for x in ds]
    return [results.results(x, pdf, cov) for x, cov in zip(ds, covs)]

@pytest.fixture(scope='module')
def single_exp_data():
    l = Loader()
    dataset_inputs = [{'name': 'NMC'},
                      {'name':'ATLASTTBARTOT', 'cfac':['QCD']},
                      {'name':'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sysnum':10}]
    ds = [l.check_dataset(**x, theoryid=162, cuts=None) for x in dataset_inputs]
    exp = ExperimentSpec('pseudo experiment', ds)
    pdf = l.check_pdf("NNPDF31_nnlo_as_0118")
    return pdf, exp

@pytest.fixture(scope='module')
def dataset_convolution_results(single_exp_data):
    pdf, exp = single_exp_data
    covs = [results.covariance_matrix(ds, False, pdf) for ds in exp.datasets]
    return [results.results(ds, pdf, cov) for ds, cov in zip(exp.datasets, covs)]

@pytest.fixture(scope='module')
def dataset_chi2data(dataset_convolution_results):
    return [results.abs_chi2_data(r) for r in dataset_convolution_results]

def chi2data_implement(convolution_results):
    return [results.abs_chi2_data_experiment(r) for r in convolution_results]

@pytest.fixture(scope='module')
def chi2data(convolution_results):
    return chi2data_implement(convolution_results)

@pytest.fixture(scope='module')
def weighted_data():
    l = Loader()
    ds = l.check_dataset(name='NMC', theoryid=162, cuts=None)
    wds = l.check_dataset(name='NMC', theoryid=162, cuts=None, weight=100)
    exp = ExperimentSpec('NMC Experiment', [ds])
    wexp = ExperimentSpec('Weighted', [wds])
    pdf = l.check_pdf("NNPDF31_nnlo_as_0118")
    exps = [exp, wexp]
    return pdf, exps

@pytest.fixture(scope='module')
def convolution_results_with_weights(weighted_data):
    return convolution_results_implement(weighted_data)

@pytest.fixture(scope='module')
def weighted_chi2data(convolution_results_with_weights):
    return chi2data_implement(convolution_results_with_weights)
