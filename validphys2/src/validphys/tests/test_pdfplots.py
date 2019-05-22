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
