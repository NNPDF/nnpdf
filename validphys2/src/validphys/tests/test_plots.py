import matplotlib
#This is to fix a weird bug in LHAPDF
matplotlib.use('agg')

import pytest

from validphys.api import API
from validphys.tests.conftest import PDF, THEORYID

@pytest.mark.linux
@pytest.mark.mpl_image_compare
def test_plotpdfs():
    pdfs = [PDF]
    Q = 10
    flavours = ['g']
    #plot_pdfs returns a generator with (figure, name_hint)
    return next(API.plot_pdfs(pdfs=pdfs, Q=Q, flavours=flavours))[0]

@pytest.mark.linux
@pytest.mark.mpl_image_compare
def test_dataspecschi2():
    dsinpts = [
        {'dataset': 'NMC'},
        {'dataset': 'ATLASTTBARTOT', 'cfac':['QCD']},
        {'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10}
    ]
    dataspecs = [
        {'pdf': PDF, 'theoryid': THEORYID, 'speclabel': 'no t0'},
        {'pdf': PDF, 'theoryid': THEORYID, 'use_t0': False, 'speclabel': 'with t0'}
    ]
    return API.plot_dataspecs_datasets_chi2(
        dataset_inputs=dsinpts,
        dataspecs=dataspecs,
        use_cuts='internal',
        metadata_group='experiment'
    )

@pytest.mark.linux
@pytest.mark.mpl_image_compare
def test_plot_xq2():
    theoryid = THEORYID
    use_cuts = "nocuts"
    display_cuts = False
    marker_by = "process type"
    metadata_group = "experiment"
    dataset_inputs = [
        {'dataset': 'NMC'},
        {'dataset': 'ATLASTTBARTOT', 'cfac':['QCD']},
        {'dataset': 'CMSZDIFF12', 'cfac':('QCD', 'NRM'), 'sys':10}    
        ]    
        
    return API.plot_xq2(
            theoryid=theoryid,
            use_cuts=use_cuts,
            dataset_inputs=dataset_inputs,
            display_cuts=display_cuts,
            marker_by=marker_by,
            metadata_group=metadata_group)
