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
def data_config():
    experiment_list = [
        {
            'experiment': 'NMC',
            'datasets': [{'dataset': 'NMC'}]},
        {
            'experiment': 'ATLASTTBARTOT',
            'datasets': [{'dataset': 'ATLASTTBARTOT', 'cfac':['QCD']}]},
        {
            'experiment': 'CMSZDIFF12',
            'datasets': [{'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10}]}
        ]
    config_dict = dict(
        pdf="NNPDF31_nnlo_as_0118",
        use_cuts='nocuts',
        experiments=experiment_list,
        theoryid=162,
        use_t0=False,
        use_fitthcovmat=False
    )
    return config_dict

@pytest.fixture(scope='module')
def data_witht0_config():
    experiment_list = [
        {
            'experiment': 'NMC',
            'datasets': [{'dataset': 'NMC'}]},
        {
            'experiment': 'ATLASTTBARTOT',
            'datasets': [{'dataset': 'ATLASTTBARTOT', 'cfac':['QCD']}]},
        {
            'experiment': 'CMSZDIFF12',
            'datasets': [{'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10}]}
        ]
    config_dict = dict(
        pdf="NNPDF31_nnlo_as_0118",
        use_cuts='nocuts',
        experiments=experiment_list,
        theoryid=162,
        use_t0=True,
        t0pdfset="NNPDF31_nnlo_as_0118",
        use_fitthcovmat=False
    )
    return config_dict

@pytest.fixture(scope='module')
def data_singleexp_witht0_config():
    experiment_list = [
        {
            'experiment': 'pseudo experiment',
            'datasets': [
                {'dataset': 'NMC'},
                {'dataset': 'ATLASTTBARTOT', 'cfac':['QCD']},
                {'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10}]}]
    config_dict = dict(
        pdf="NNPDF31_nnlo_as_0118",
        use_cuts='nocuts',
        experiments=experiment_list,
        theoryid=162,
        use_t0=True,
        t0pdfset="NNPDF31_nnlo_as_0118",
        use_fitthcovmat=False
    )
    return config_dict

@pytest.fixture(scope='module')
def weighted_data_witht0_config():
    experiment_list = [
        {
            'experiment': 'NMC Experiment',
            'datasets': [{'dataset': 'NMC'}]},
        {
            'experiment': 'Weighted',
            'datasets': [{'dataset': 'NMC', 'weight': 100}]},
        ]
    config_dict = dict(
        pdf="NNPDF31_nnlo_as_0118",
        use_cuts='nocuts',
        experiments=experiment_list,
        theoryid=162,
        use_t0=True,
        t0pdfset="NNPDF31_nnlo_as_0118",
        use_fitthcovmat=False
    )
    return config_dict
