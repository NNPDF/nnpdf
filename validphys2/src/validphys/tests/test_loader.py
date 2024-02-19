"""
test_loader.py

Test loading utilities.
"""
import os
from pathlib import Path
import subprocess as sp
import sys

from hypothesis import given, settings
from hypothesis.strategies import composite, sampled_from, sets
import pytest

from validphys.loader import FallbackLoader, FitNotFound, NNPDF_DIR
from validphys.plotoptions.core import kitable, get_info
from validphys.tests.conftest import FIT, FIT_3REPLICAS, THEORYID_NEW

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
    cd = l.check_commondata(old_name, force_old_format=True)
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
        tmp_path / "theories", profile=profile_path, res_type="theoryID", resource=str(THEORYID_NEW)
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
        theory_path, profile=profile_path, res_type="theoryID", resource=str(THEORYID_NEW)
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
