"""
test_pythonmakereplica.py

Module for testing the python implementation of make replica

"""
import numpy as np
import pytest

from validphys.api import API
from validphys.pseudodata import make_replica
from validphys.tests.conftest import DATA
from validphys.tests.test_covmats import CORR_DATA


SINGLE_SYS_DATASETS = [
    {"dataset": "DYE886R"},
    {"dataset": "D0ZRAP", "cfac": ["QCD"]},
    {"dataset": "ATLAS_SINGLETOP_TCH_R_13TEV", "cfac": ["QCD"]},
    {"dataset": "CMS_SINGLETOP_TCH_R_8TEV", "cfac": ["QCD"]},
    {"dataset": "CMS_SINGLETOP_TCH_R_13TEV", "cfac": ["QCD"]},
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
    ld_cds = API.dataset_inputs_loaded_cd_with_cuts(**config)
    # keep a copy of all central values and all systematics
    original_cvs = [np.array(cd.central_values.to_numpy(copy=True)) for cd in ld_cds]
    original_sys = [np.array(cd.systematics_table.to_numpy(copy=True)) for cd in ld_cds]
    make_replica(ld_cds)
    for post_mkrep_cd, pre_mkrep_cv, pre_mkrep_sys in zip(ld_cds, original_cvs, original_sys):
        np.testing.assert_allclose(
            post_mkrep_cd.central_values.to_numpy(), pre_mkrep_cv
        )
        np.testing.assert_allclose(
            post_mkrep_cd.systematics_table.to_numpy(), pre_mkrep_sys
        )


@pytest.mark.parametrize("use_cuts", ["nocuts", "internal"])
@pytest.mark.parametrize("dataset_inputs", [DATA, CORR_DATA, SINGLE_SYS_DATASETS])
def test_pseudodata_seeding(data_config, dataset_inputs, use_cuts):
    """Check that using a seed reproduces the pseudodata. Note that this also
    will check that the commondata hasn't been modified since reproducing
    the same commondata requires that the commondata is unchanged and that
    the random numbers are generated and used identically.

    """
    seed = 123456
    config = dict(data_config)
    config["dataset_inputs"] = dataset_inputs
    config["use_cuts"] = use_cuts
    ld_cds = API.dataset_inputs_loaded_cd_with_cuts(**config)
    rep_1 = make_replica(ld_cds, seed=seed)
    rep_2 = make_replica(ld_cds, seed=seed)
    np.testing.assert_allclose(rep_1, rep_2)


@pytest.mark.parametrize("use_cuts", ["nocuts", "internal"])
@pytest.mark.parametrize("dataset_inputs", [DATA, CORR_DATA, SINGLE_SYS_DATASETS])
def test_pseudodata_has_correct_ndata(data_config, dataset_inputs, use_cuts):
    """Check that we get the correct ndata when generating pseudodata"""
    config = dict(data_config)
    config["dataset_inputs"] = dataset_inputs
    config["use_cuts"] = use_cuts
    ld_cds = API.dataset_inputs_loaded_cd_with_cuts(**config)
    rep = make_replica(ld_cds)
    ndata = np.sum([cd.ndata for cd in ld_cds])
    assert len(rep) == ndata
