"""
conftest.py

Pytest fixtures.
"""
import pathlib

import pytest

from validphys.loader import FallbackLoader as Loader
from validphys.core import ExperimentSpec
from validphys import results



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
def theory_data():
    l = Loader()
    names = ['NMC', 'ATLASTTBARTOT']
    theoryids = [52, 180, 173]
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
    pdf, exps = data
    return [results.experiment_results(exp, pdf, pdf) for exp in exps]

@pytest.fixture
def dataset_t0_convolution_results(data):
    pdf, exps = data
    ds = [x.datasets[0] for x in exps]
    return [results.results(x, pdf, t0set=pdf) for x in ds]

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
    return [results.results(ds, pdf, pdf) for ds in exp.datasets]

@pytest.fixture(scope='module')
def dataset_chi2data(dataset_convolution_results):
    return [results.abs_chi2_data(r) for r in dataset_convolution_results]

@pytest.fixture(scope='module')
def chi2data(convolution_results):
    return [results.abs_chi2_data_experiment(r) for r in convolution_results]

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
    return convolution_results(weighted_data)

@pytest.fixture(scope='module')
def weighted_chi2data(convolution_results_with_weights):
    return chi2data(convolution_results_with_weights)
