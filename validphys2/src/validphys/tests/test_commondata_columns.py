"""
Test that the multiplicative and additive columns of commondata
are consistent.

"""
import numpy as np
import pytest

from validphys.commondataparser import load_commondata
from validphys.loader import Loader

def _cd_column_comparison(cd):
    """Given a py:func:`validphys.coredata.CommonData` object, test
    if the multiplicative and additive columns are consistent. In general
    the ``mult`` columns and the ``add`` columns should be related by:

        add = mult * central_value * 1e-2

    """
    add = cd.systematics_table.loc[:, ["ADD"]]
    mult = cd.systematics_table.loc[:, ["MULT"]]
    cv = cd.central_values
    np.testing.assert_allclose(
        cv.to_numpy()[:, np.newaxis] * mult.to_numpy() * 1e-2,
        add.to_numpy(),
        rtol=1e-05, # use the allclose values (less stringent tolerances)
        atol=1e-08,
    )

DS_NAMES = [
    "LHCBZEE2FB_40",
    "LHCBW36PB_40",
    "CDFZRAP_NEW",
    "CMSTTBARTOT7TEV",
    "CMSTTBARTOT8TEV",
    "CMSTTBARTOT13TEV",
    "CMS_2JET_3D_8TEV",
    "D0ZRAP_40",
    "ATLASPHT15",
]

L = Loader()

@pytest.mark.parametrize("dataset_name", DS_NAMES)
def test_dataset_commondata_columns_match(dataset_name):
    """Apply :py:func:`_cd_column_comparison` to list of
    datasets in ``DS_NAMES``. Checking if the columns in the
    commondata file are self consistent.

    """
    # Note: If the dataset should have additional contributions to the
    # covmat, like sys: 10 then DS_NAMES should be changed to be proper dataset
    # inputs and we should use API here. Use loader to keep test fast.
    cdspec = L.check_commondata(dataset_name)
    lcd = load_commondata(cdspec)
    _cd_column_comparison(lcd)
