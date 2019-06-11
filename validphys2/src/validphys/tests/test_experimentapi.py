import pytest
from validphys.api import API


def test_load():
    """Loading datasets without libnnpdf from commondata"""
    data = [
        {'name': 'NMC'},
        {'name': 'ATLASTTBARTOT', 'cfac': ['QCD']},
        {'name': 'CMSZDIFF12', 'cfac': ['QCD', 'NRM'], 'sys': 10}
    ]
    for d in data:
        API.commondata(dataset_input=d).load()
