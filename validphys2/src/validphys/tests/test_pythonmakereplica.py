"""
test_pythonmakereplica.py

Module for testing the python implementation of make replica

"""

from copy import deepcopy

import numpy as np
from pandas.testing import assert_frame_equal, assert_series_equal
import pytest

from validphys.api import API
from validphys.tests.conftest import DATA
from validphys.tests.test_covmats import CORR_DATA

SEED = 123456

# Datasets to be tested
SINGLE_SYS_DATASETS = [
    {"dataset": "DYE866_Z0_800GEV_DW_RATIO_PDXSECRATIO", "variant": "legacy"},
    {"dataset": "D0_Z0_1P96TEV_ZRAP", "variant": "legacy"},
    {"dataset": "NMC_NC_NOTFIXED_P_EM-SIGMARED", "variant": "legacy"},
    {"dataset": "NMC_NC_NOTFIXED_EM-F2", "variant": "legacy"},
    {"dataset": "ATLAS_Z0J_8TEV_PT-M", "variant": "legacy"},
    {"dataset": "ATLAS_DY_7TEV_36PB_ETA", "variant": "legacy"},
    {"dataset": "ATLAS_Z0_7TEV_49FB_HIMASS", "variant": "legacy"},
    {"dataset": "CMS_WPWM_7TEV_ELECTRON_ASY"},
    {"dataset": "CMS_WPWM_7TEV_MUON_ASY", "variant": "legacy"},
]


@pytest.mark.parametrize("use_cuts", ["nocuts", "internal"])
@pytest.mark.parametrize("dataset_inputs", [DATA, CORR_DATA, SINGLE_SYS_DATASETS])
def test_commondata_unchanged(data_config, dataset_inputs, use_cuts):
    """Load the commondata, then generate some pseudodata using make replica
    Check that the following attributes of the commondata have not been
    modified: central_values, commondata_table, systematics_table.

    """
    config = dict(data_config)
    config["dataset_inputs"] = dataset_inputs
    config["use_cuts"] = use_cuts
    config["replica_mcseed"] = SEED
    ld_cds = API.dataset_inputs_loaded_cd_with_cuts(**config)

    # keep a copy of all dataframes/series pre make replica
    pre_mkrep_cvs = [deepcopy(cd.central_values) for cd in ld_cds]
    pre_mkrep_sys_tabs = [deepcopy(cd.systematics_table) for cd in ld_cds]
    pre_mkrep_cd_tabs = [deepcopy(cd.commondata_table) for cd in ld_cds]
    API.make_replica(**config)

    for post_mkrep_cd, pre_mkrep_cv in zip(ld_cds, pre_mkrep_cvs):
        assert_series_equal(post_mkrep_cd.central_values, pre_mkrep_cv)

    for post_mkrep_cd, pre_mkrep_sys_tab in zip(ld_cds, pre_mkrep_sys_tabs):
        assert_frame_equal(post_mkrep_cd.systematics_table, pre_mkrep_sys_tab)

    for post_mkrep_cd, pre_mkrep_cd_tab in zip(ld_cds, pre_mkrep_cd_tabs):
        assert_frame_equal(post_mkrep_cd.commondata_table, pre_mkrep_cd_tab)


@pytest.mark.parametrize("use_cuts", ["nocuts", "internal"])
@pytest.mark.parametrize("dataset_inputs", [DATA, CORR_DATA, SINGLE_SYS_DATASETS])
def test_pseudodata_seeding(data_config, dataset_inputs, use_cuts):
    """Check that using a seed reproduces the pseudodata. Note that this also
    will check that the commondata hasn't been modified since reproducing
    the same commondata requires that the commondata is unchanged and that
    the random numbers are generated and used identically.

    """
    config = dict(data_config)
    config["dataset_inputs"] = dataset_inputs
    config["use_cuts"] = use_cuts
    config["replica_mcseed"] = SEED
    rep_1 = API.make_replica(**config)
    rep_2 = API.make_replica(**config)
    np.testing.assert_allclose(rep_1, rep_2)


@pytest.mark.parametrize("use_cuts", ["nocuts", "internal"])
@pytest.mark.parametrize("dataset_inputs", [DATA, CORR_DATA, SINGLE_SYS_DATASETS])
def test_pseudodata_has_correct_ndata(data_config, dataset_inputs, use_cuts):
    """Check that we get the correct ndata when generating pseudodata"""
    config = dict(data_config)
    config["dataset_inputs"] = dataset_inputs
    config["use_cuts"] = use_cuts
    config["replica_mcseed"] = SEED
    ld_cds = API.dataset_inputs_loaded_cd_with_cuts(**config)
    rep = API.make_replica(**config)
    ndata = np.sum([cd.ndata for cd in ld_cds])
    assert len(rep) == ndata


@pytest.mark.parametrize("use_cuts", ["nocuts", "internal"])
@pytest.mark.parametrize("dataset_inputs", [DATA, CORR_DATA, SINGLE_SYS_DATASETS])
def test_genrep_off(data_config, dataset_inputs, use_cuts):
    """Check that when genrep is set to False replicas are not generated."""
    config = dict(data_config)
    config["dataset_inputs"] = dataset_inputs
    config["use_cuts"] = use_cuts
    config["replica_mcseed"] = SEED
    config["genrep"] = False

    # `make_replica` reorders the data according to the groups in the exp covmat
    not_replica = API.make_replica(**config)

    # And so it needs to be tested against the central values loaded by group
    ld_cds = API.groups_dataset_inputs_loaded_cd_with_cuts(**config)

    central_data = np.concatenate([d.central_values for d in ld_cds])
    np.testing.assert_allclose(not_replica, central_data)
