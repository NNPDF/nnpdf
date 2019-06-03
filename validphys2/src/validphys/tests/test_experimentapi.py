import pytest

from validphys.api import API
from validphys.core import ExperimentAPI

def test_load():
    experiments = [
        {
            'experiment': 'NMCexp',
            'datasets': [{'dataset': 'NMC'}]},
        {
            'experiment': 'ATLASxp',
            'datasets': [{'dataset': 'ATLASTTBARTOT', 'cfac':['QCD']}]},
        {
            'experiment': 'CMSexp',
            'datasets': [{'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10}]}
        ]
    ExperimentAPI(experiments_list=experiments)
