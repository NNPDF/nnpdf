import numpy as np
import pandas as pd
import pytest

from validphys.api import API
from validphys.commondataparser import load_commondata
from validphys.covmats_utils import construct_covmat
from validphys.loader import FallbackLoader as Loader
from validphys.tests.conftest import FIT, THEORYID_NEW


def test_basic_commondata_loading():
    l = Loader()
    cd = l.check_commondata(setname="SLAC_NC_NOTFIXED_D_EM-F2", variant="legacy_dw")
    res = load_commondata(cd)
    # Test commondata loading
    assert res.ndata == 211
    assert isinstance(res.commondata_table, pd.DataFrame)
    # Test systype loading
    assert res.nsys == 103
    assert isinstance(res.systype_table, pd.DataFrame)
    rules = API.rules(
        **{
            "dataset_input": "SLAC_NC_NOTFIXED_D_DW_EM-F2",
            "theoryid": THEORYID_NEW,
            "use_cuts": "internal",
        }
    )
    # Test a dataset with no systematics
    emptysyscd = l.check_posset(
        theoryID=THEORYID_NEW, setname='NNPDF_POS_2P24GEV_XDQ', postlambda=1e-10, rules=rules
    )
    emptysysres = load_commondata(emptysyscd.commondata)
    assert emptysysres.nsys == 0
    assert emptysysres.systype_table.empty is True


def test_commondata_with_cuts():
    l = Loader()
    setname = "NMC_NC_NOTFIXED_P_EM-SIGMARED"

    cd = l.check_commondata(setname=setname, variant="legacy")
    loaded_cd = load_commondata(cd)

    fit_cuts = l.check_fit_cuts(fit=FIT, commondata=cd)
    internal_cuts = l.check_internal_cuts(cd, API.rules(theoryid=THEORYID_NEW, use_cuts="internal"))

    loaded_cd_fit_cuts = loaded_cd.with_cuts(fit_cuts)
    # We must do these - 1 subtractions due to the fact that cuts indexing
    # starts at 0 while commondata indexing starts at 1
    assert all(loaded_cd_fit_cuts.commondata_table.index - 1 == fit_cuts.load())
    assert all(loaded_cd_fit_cuts.additive_errors.index - 1 == fit_cuts.load())
    assert all(loaded_cd_fit_cuts.multiplicative_errors.index - 1 == fit_cuts.load())

    loaded_cd_internal_cuts = loaded_cd.with_cuts(internal_cuts)
    assert all(loaded_cd_internal_cuts.commondata_table.index - 1 == internal_cuts.load())

    loaded_cd_nocuts = loaded_cd.with_cuts(None)
    assert all(loaded_cd_nocuts.commondata_table.index == range(1, cd.ndata + 1))

    preloaded_fit_cuts = fit_cuts.load()
    loaded_cd_preloaded_cuts = loaded_cd.with_cuts(fit_cuts)
    assert all(loaded_cd_preloaded_cuts.commondata_table.index - 1 == preloaded_fit_cuts)

    assert all(loaded_cd.with_cuts([1, 2, 3]).commondata_table.index - 1 == [1, 2, 3])

    # Check that giving cuts for another dataset raises the correct ValueError exception
    cd_bad = l.check_commondata(setname="NMC_NC_NOTFIXED_EM-F2", variant="legacy")
    bad_cuts = l.check_fit_cuts(fit=FIT, commondata=cd_bad)
    with pytest.raises(ValueError):
        loaded_cd.with_cuts(bad_cuts)


def test_commondata_load_write_load(tmp):
    """Test that we can a commondata, write it down, and load it again"""
    l = Loader()

    # Select a dataset that we know mixes ADD and MULT (so that the ordering is checked)
    setname = "ATLAS_2JET_7TEV_R06_M12Y"
    # And a complicated variant
    variant = "legacy"

    # Get a reference to the commondata
    cd = l.check_commondata(setname=setname, variant=variant)

    # Load it up, save the covmat, and write it down
    original_data = cd.load()
    data_path, unc_path = original_data.export(tmp)

    # Now, reload it with the new data/unc paths
    new_data = cd.with_modified_data(data_path, uncertainties_file=unc_path).load()

    # central value!
    original_cv = original_data.central_values.to_numpy()
    new_cv = new_data.central_values.to_numpy()
    np.testing.assert_allclose(original_cv, new_cv)

    # stats!
    original_stats = original_data.stat_errors.to_numpy()
    new_stats = new_data.stat_errors.to_numpy()
    np.testing.assert_allclose(original_stats, new_stats)

    # Create fake data in order to check whether the covmats are truly the same
    # the fake data ensures that the MULT and ADD are treated in the same way in both
    # otherwise, since the data is saved wrt to the original central value, the test will always pass
    fake_unc = np.diag(construct_covmat(original_stats, original_data.systematic_errors()))
    fake_data = np.random.rand(len(original_stats)) * fake_unc + original_cv

    new_covmat = construct_covmat(new_stats, new_data.systematic_errors(fake_data))
    original_covmat = construct_covmat(original_stats, original_data.systematic_errors(fake_data))
    np.testing.assert_allclose(new_covmat, original_covmat)

def test_variant_nnpdf_metadata():
    """Tests the undocumented feature of a variant which updates the key `experiment`
          within the nnpdf_metadata
    """
    l = Loader()
    set_name = "SLAC_NC_NOTFIXED_D_EM-F2"

    for v1, v2 in [("legacy", "legacy_dw"), ("legacy_dw", "legacy")]:
        cd1 = l.check_commondata(setname=set_name, variant=v1)
        pcd1 = cd1.metadata.plotting_options
        cd2 = l.check_commondata(setname=set_name, variant=v2)
        pcd2 = cd2.metadata.plotting_options

        # ensure the nnpdf_metadata and the plotting are changed
        assert cd1.metadata.nnpdf_metadata["experiment"] != cd2.metadata.nnpdf_metadata["experiment"]
        assert pcd2.experiment != pcd1.experiment
        # but the real experiment is the same
        assert cd1.metadata.experiment == cd2.metadata.experiment
        # And check that the legacy names are _not_ the same
        assert cd1.legacy_names != cd2.legacy_names
