import matplotlib

# This is to fix a weird bug in LHAPDF
matplotlib.use('agg')

import pytest

from validphys.api import API
from validphys.tests.conftest import DATA, PDF, SINGLE_DATASET, THEORYID

TOLERANCE_VALUE = 18


@pytest.mark.linux
@pytest.mark.mpl_image_compare()
def test_plotpdfs():
    pdfs = [PDF]
    Q = 10
    flavours = ['g']
    # plot_pdfs returns a generator with (figure, name_hint)
    return next(iter(API.plot_pdfs(pdfs=pdfs, Q=Q, flavours=flavours)))[0]


@pytest.mark.linux
@pytest.mark.mpl_image_compare()
def test_dataspecschi2():
    dsinpts = [
        {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy'},
        {'dataset': 'ATLAS_TTBAR_7TEV_TOT_X-SEC', 'variant': 'legacy'},
        {'dataset': 'CMS_Z0J_8TEV_PT-Y', 'cfac': ['NRM'], 'variant': 'legacy'},
    ]
    dataspecs = [
        {'pdf': PDF, 'theoryid': THEORYID, 'speclabel': 'no t0'},
        {'pdf': PDF, 'theoryid': THEORYID, 'use_t0': False, 'speclabel': 'with t0'},
    ]
    return API.plot_dataspecs_datasets_chi2(
        dataset_inputs=dsinpts,
        dataspecs=dataspecs,
        use_cuts='internal',
        metadata_group='experiment',
    )


@pytest.mark.linux
@pytest.mark.mpl_image_compare()
def test_plotfancy():
    fig = API.plot_fancy(dataset_input=DATA[2], theoryid=THEORYID, pdfs=[PDF], use_cuts='internal')[
        0
    ]
    fig.tight_layout()
    return fig


@pytest.mark.linux
@pytest.mark.mpl_image_compare()
def test_plot_smpdf(single_data_internal_cuts_config):
    return next(iter(API.plot_smpdf(**single_data_internal_cuts_config)))


@pytest.mark.linux
@pytest.mark.mpl_image_compare()
def test_plot_smpdf_categorical(single_data_categorical_internal_cuts_config):
    return next(iter(API.plot_smpdf(**single_data_categorical_internal_cuts_config)))


@pytest.mark.linux
@pytest.mark.mpl_image_compare()
def test_plot_obscorrs(single_data_internal_cuts_config):
    corrpair = [{"corrpair": (i["dataset"],)} for i in DATA[:2]]
    return API.plot_obscorrs(**single_data_internal_cuts_config, corrpair=corrpair)


@pytest.mark.linux
@pytest.mark.mpl_image_compare()
def test_plot_xq2():
    theoryid = THEORYID
    use_cuts = "nocuts"
    display_cuts = False
    marker_by = "process type"
    metadata_group = "experiment"
    dataset_inputs = [
        {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy'},
        {'dataset': 'ATLAS_TTBAR_7TEV_TOT_X-SEC', 'variant': 'legacy'},
        {'dataset': 'CMS_Z0J_8TEV_PT-Y', 'cfac': ['NRM']},
    ]

    return API.plot_xq2(
        theoryid=theoryid,
        use_cuts=use_cuts,
        dataset_inputs=dataset_inputs,
        display_cuts=display_cuts,
        marker_by=marker_by,
        metadata_group=metadata_group,
    )


@pytest.mark.linux
@pytest.mark.mpl_image_compare()
def test_plot_xq2_custom():
    theoryid = THEORYID
    use_cuts = "nocuts"
    display_cuts = False

    marker_by = "group"
    metadata_group = "custom_group"

    dataset_inputs = [
        {'dataset': 'NMC_NC_NOTFIXED_P_EM-SIGMARED', 'variant': 'legacy', 'custom_group': 'one'},
        {'dataset': 'ATLAS_TTBAR_7TEV_TOT_X-SEC', 'variant': 'legacy', 'custom_group': 'one'},
        {'dataset': 'CMS_Z0J_8TEV_PT-Y', 'cfac': ['NRM'], 'custom_group': 'two'},
    ]

    return API.plot_xq2(
        theoryid=theoryid,
        use_cuts=use_cuts,
        dataset_inputs=dataset_inputs,
        display_cuts=display_cuts,
        marker_by=marker_by,
        metadata_group=metadata_group,
    )
