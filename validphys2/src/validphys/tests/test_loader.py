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
import numpy as np
import pytest

from validphys.core import CommonDataSpec, Cuts
from validphys.loader import FallbackLoader, FitNotFound, rebuild_commondata_without_cuts
from validphys.plotoptions import get_info, kitable
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
    cd = l.check_commondata(draw(sampled_from(dss)))
    ndata = cd.metadata.ndata
    # Get a cut mask with at least one selected datapoint
    masks = sets(sampled_from(range(ndata)), min_size=1)
    mask = sorted(draw(masks))
    return cd, mask


@given(arg=commondata_and_cuts())
@settings(deadline=None)
def test_rebuild_commondata_without_cuts(tmp_path_factory, arg):
    # We need to create a new directory for each call of the test
    # otherwise we get files mixed together
    tmp = tmp_path_factory.mktemp("test_loader")

    cd, cuts = arg
    lcd = cd.load()
    cutspec = None
    if cuts:
        cutpath = tmp / "cuts.txt"
        np.savetxt(cutpath, np.asarray(cuts, dtype=int), fmt="%u")
        cutspec = Cuts(cd, cutpath)
        lcd = lcd.with_cuts(cuts)
    lcd.export(tmp)
    # We have to reconstruct the name here...
    with_cuts = tmp / f"DATA_{cd.name}.dat"
    newpath = tmp / "commondata.dat"
    rebuild_commondata_without_cuts(with_cuts, cutspec, cd.datafile, newpath)
    newcd = CommonDataSpec(newpath, cd.sysfile, cd.plotfiles)
    # Note this one is without cuts
    t1 = kitable(cd, get_info(cd))
    t2 = kitable(newcd, get_info(newcd))
    assert (t1 == t2).all
    lncd = newcd.load()
    if cuts:
        assert np.allclose(lncd.get_cv()[cuts], lcd.get_cv())
        nocuts = np.ones(cd.ndata, dtype=bool)
        nocuts[cuts] = False
        assert (lncd.get_cv()[nocuts] == 0).all()


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

    nnpdf_path = tmp_path / "NNPDF"
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

    results_sys_prefix = Path(sys.prefix) / "share" / "NNPDF" / "results"
    _check_download_resource(results_sys_prefix, profile=profile_path)
