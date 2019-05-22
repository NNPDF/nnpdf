import matplotlib
#This is to fix a weird bug in LHAPDF
matplotlib.use('agg')

import pytest

from validphys.api import API

@pytest.mark.mpl_image_compare
def test_plotpdfs():
    pdfs = ['NNPDF31_nnlo_as_0118']
    Q = 10
    flavours = ['g']
    #plot_pdfs returns a generator with (figure, name_hint)
    return next(API.plot_pdfs(pdfs=pdfs, Q=Q, flavours=flavours))[0]

@pytest.mark.mpl_image_compare
def test_dataspecschi2():
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
    dataspecs = [
        {'pdf': 'NNPDF31_nnlo_as_0118', 'theoryid': 162, 'speclabel': 'no t0'},
        {'pdf': 'NNPDF31_nnlo_as_0118', 'theoryid': 162, 'use_t0': False, 'speclabel': 'with t0'}
    ]
    return API.plot_dataspecs_datasets_chi2(
        experiments=experiments, dataspecs=dataspecs, use_cuts='nocuts')
