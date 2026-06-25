import matplotlib

# This is to fix a weird bug in LHAPDF
matplotlib.use('agg')

from io import BytesIO

from PIL import Image, ImageFilter
import matplotlib.pyplot as plt
import pytest

from validphys.api import API
from validphys.tests.conftest import DATA, PDF, THEORYID

BLUR_RADIUS = 1.5


def _blurred_figure(fig, blur_radius=BLUR_RADIUS):
    """Return a blurred raster copy of a figure for less pixel-sensitive mpl comparisons."""
    buffer = BytesIO()
    try:
        fig.savefig(buffer, format='png')
    finally:
        plt.close(fig)

    buffer.seek(0)
    image = Image.open(buffer).convert('RGB').filter(ImageFilter.GaussianBlur(blur_radius))

    blurred_fig = plt.figure(figsize=(image.width / 100, image.height / 100), dpi=100)
    ax = blurred_fig.add_axes([0, 0, 1, 1])
    ax.imshow(image)
    ax.set_axis_off()
    return blurred_fig


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
    fig = API.plot_dataspecs_datasets_chi2(
        dataset_inputs=dsinpts,
        dataspecs=dataspecs,
        use_cuts='internal',
        metadata_group='experiment',
    )
    return _blurred_figure(fig)


@pytest.mark.linux
@pytest.mark.mpl_image_compare()
def test_plotfancy():
    fig = API.plot_fancy(
        dataset_input=DATA[2], theoryid=THEORYID, pdfs=[PDF], use_cuts='internal', with_shift=True
    )[0]
    fig.tight_layout()
    return _blurred_figure(fig)


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
