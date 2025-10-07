"""
test_loader.py

Test loading utilities.
"""

import os
from pathlib import Path
import shutil
import subprocess as sp
import sys

from hypothesis import given, settings
from hypothesis.strategies import composite, sampled_from, sets
import numpy as np
import pytest

from nnpdf_data import legacy_to_new_map
from nnpdf_data.utils import DEFAULT_PATH_VPDATA, NNPDF_DIR, yaml_fast
from reportengine.configparser import ConfigError
from validphys.api import API
from validphys.loader import FallbackLoader, FitNotFound
from validphys.plotoptions.core import get_info, kitable
from validphys.tests.conftest import FIT, FIT_3REPLICAS, THEORYID

l = FallbackLoader()
# The sorted is to appease hypothesis
dss = sorted(l.available_datasets - {"PDFEVOLTEST"})


class MockCuts:
    def __init__(self, arr):
        self.arr = arr

    def load(self):
        return self.arr


@composite
def commondata_and_cuts(draw):
    old_name = draw(sampled_from(dss))
    new_name, variant = legacy_to_new_map(old_name)
    cd = l.check_commondata(new_name, variant=variant)
    ndata = cd.metadata.ndata
    # Get a cut mask with at least one selected datapoint
    masks = sets(sampled_from(range(ndata)), min_size=1)
    mask = sorted(draw(masks))
    return cd, mask


@given(inp=commondata_and_cuts())
@settings(deadline=None)
def test_kitable_with_cuts(inp):
    cd, cuts = inp
    info = get_info(cd, cuts=cuts)
    tb = kitable(cd, info, cuts=MockCuts(cuts))
    assert len(tb) == len(cuts)


def test_load_fit():
    assert l.check_fit(FIT)
    with pytest.raises(FitNotFound):
        l.check_fit(f"{FIT}/")


def test_load_old_error(data_config):
    """Checks that when loading a dataset with an old name, it errors out
    if ``allow_legacy_names`` is set to False"""
    idict = {**data_config, "dataset_input": {"dataset": "LHCBZ940PB"}}
    # This works
    _ = API.dataset(**idict, allow_legacy_names=True)
    # This doesn't
    with pytest.raises(ConfigError):
        _ = API.dataset(**idict, allow_legacy_names=False)


### nnprofile testing
def _check_download_resource(results_dir, profile=None, res_type="fit", resource=FIT_3REPLICAS):
    """Downloads a resource (by default the fit FIT_3REPLICAS and
    checks whether it can indeed be found in the given ``results_dir``
    Accepts a ``profile`` pointing to a custom nnprofile which will be set with NNPDF_PROFILE_PATH
    """
    custom_environ = dict(os.environ)
    if profile is not None:
        custom_environ["NNPDF_PROFILE_PATH"] = Path(profile).as_posix()

    sp.run(["vp-get", res_type, resource], env=custom_environ)

    if res_type == "theoryID":
        res_location = results_dir / f"theory_{resource}"
    else:
        res_location = results_dir / resource

    assert res_location.exists()
    assert res_location.is_dir()


def test_custom_profile_nnpdf_share(tmp_path):
    """Creates a custom profile with a key ``nnpdf_share` and checks whether it is used"""
    profile_path = tmp_path / "nnprofile.yaml"
    profile_path.write_text(f"nnpdf_share: {tmp_path}")

    # Check a fit
    _check_download_resource(tmp_path / "results", profile=profile_path)
    # Check a theory
    _check_download_resource(
        tmp_path / "theories", profile=profile_path, res_type="theoryID", resource=str(THEORYID)
    )


def test_custom_profile_explicit(tmp_path):
    """Creates a custom profile with explicit paths that should take precedence"""
    profile_path = tmp_path / "nnprofile.yaml"
    theory_path = tmp_path / "test_theory"
    results_path = tmp_path / "test_fit"

    profile_path.write_text(
        f"""
theories_path: {theory_path.as_posix()}
results_path: {results_path.as_posix()}"""
    )

    # Check a fit
    _check_download_resource(results_path, profile=profile_path)
    # Check a theory
    _check_download_resource(
        theory_path, profile=profile_path, res_type="theoryID", resource=str(THEORYID)
    )


def test_home_profile(tmp_path):
    """Check that {XDG_CONFIG_HOME} profile takes precedence"""
    original_xdg = os.environ.get("XDG_CONFIG_HOME")

    nnpdf_path = tmp_path / NNPDF_DIR
    nnpdf_path.mkdir(exist_ok=True, parents=True)

    profile_path = nnpdf_path / "nnprofile.yaml"
    profile_path.write_text(f"nnpdf_share: {tmp_path}")

    os.environ["XDG_CONFIG_HOME"] = tmp_path.as_posix()

    try:
        _check_download_resource(tmp_path / "results")
    finally:
        if original_xdg is None or not original_xdg:
            os.environ.pop("XDG_CONFIG_HOME")
        else:
            os.environ["XDG_CONFIG_HOME"] = original_xdg


def test_profile_relative_to_python(tmp_path):
    """Check that when the key RELATIVE_TO_PYTHON is used, the sys.prefix is used"""
    profile_path = tmp_path / "nnprofile.yaml"
    profile_path.write_text("nnpdf_share: RELATIVE_TO_PYTHON")

    results_sys_prefix = Path(sys.prefix) / "share" / NNPDF_DIR / "results"
    _check_download_resource(results_sys_prefix, profile=profile_path)


def test_extra_data_sources(tmp_path):
    """Creates a custom profile with a ``data_path`` key and checks whether it is used
    and whether it took precedence.
    """
    # Target dataset (and name of the data file), both need to be changed toghether
    dataset = "LHCB_Z0_8TEV_MUON_Y"
    faketaset = "LHCB_Z0_8TEV_FAKE_Y"
    data_file = "data.yaml"

    # Create the loader with a profile pointing to the temporary folder for data
    profile_path = tmp_path / "nnprofile.yaml"
    profile_path.write_text(f"data_path: [{tmp_path}]")
    l = FallbackLoader(profile=profile_path)

    original_set = dataset.rsplit("_", 1)[0]
    original_folder = DEFAULT_PATH_VPDATA / original_set
    fake_folder = tmp_path / faketaset.rsplit("_", 1)[0]
    # Copy the dataset, first with the fake name
    shutil.copytree(original_folder, fake_folder)

    # And load it
    fake_ds = l.check_dataset(faketaset, theoryid=40_000_000)
    assert fake_ds.commondata.metadata._parent.folder == fake_folder

    # But we can still load the original one, right?
    original_ds = l.check_dataset(dataset, theoryid=40_000_000)
    assert original_ds.commondata.metadata._parent.folder == DEFAULT_PATH_VPDATA / original_folder

    # And if we substitute it now with random data but with the same name?
    shutil.copytree(original_folder, tmp_path / original_set)

    target_data = tmp_path / original_set / data_file
    data_dict = yaml_fast.load(target_data.read_text(encoding="utf-8"))
    random_data = np.random.rand(len(data_dict["data_central"]))
    data_dict["data_central"] = random_data.tolist()
    with target_data.open("w", encoding="utf-8") as sf:
        yaml_fast.dump(data_dict, sf)

    new_original_ds = l.check_dataset(dataset, theoryid=40_000_000)
    read_data = new_original_ds.load_commondata().central_values
    read_cuts = new_original_ds.cuts.load()
    np.testing.assert_array_equal(read_data, random_data[read_cuts])
