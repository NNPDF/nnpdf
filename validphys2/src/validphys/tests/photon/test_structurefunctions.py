import pytest
from validphys.photon.structure_functions import F2LO
import numpy as np

def test_zero_pdfs():
    class fake_pdfs:
        def xfxQ(self, x, Q):
            res = {}
            for i in range(1, 6+1):
                res[i] = res[-i] = 0.
            return res

    pdfs = fake_pdfs()

    fake_theory = {
        "Qmc": 1.3,
        "Qmb": 5. ,
        "Qmt": 172.,
        "MaxNfPdf": 5,
    }

    f2lo = F2LO(pdfs, fake_theory)

    np.testing.assert_equal(f2lo.Qmt, np.inf)

    for x in np.geomspace(1e-4, 1., 10):
        for Q in np.geomspace(10, 1000000, 10):
            np.testing.assert_allclose(f2lo.FxQ(x, Q), 0.)

