"""
Regression tests for full fits.
"""

import pathlib
import shutil

import pytest

from n3fit.tests.helpers import run_n3fit
from n3fit.tests.test_fit import check_fit_results
from validphys.utils import yaml_safe

REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regression_fits")
SETUPFIT_FOLDER = REGRESSION_FOLDER / "setupfits"

# Avoid always round-number replicas or 1/2
RUNCARD_AND_REPLICAS = {
    "thcovmat": 43,
    "no_positivity": 316,
    "normal_fit": 72,
    "central": 16,
    "no_diagonal": 45,
    "feature_scaling": 81,
    "flavour": 29,
    "no_msr": 92,
    "no_sumrules": 18,
    "no_vsr": 54,
    "trainable_prepro": 61,
    "no_lagrange": 27,
    "no_csr": 613,
    "polarized_evol": 34,
    "t0theoryid": 100,
    "no_t0_sampling": 430,
    "hyperopt_sampling": 4,
}

# Some runcards need to be a bit more lenient with the tolerances
# WARNING: 0.01 is too big of a difference, we should also check # of mismatched elements
extra_tolerances_exportgrid = {"hyperopt_sampling": 0.01, "t0theoryid": 1e-3, "no_diagonal": 0.01}
extra_tolerances_rel = {"hyperopt_sampling": 3e-2}


def prepare_runcard(tmp_path, runcard_filename):
    """Copies the the runcard file to the path in which the test will run.
    If there are weights to be loaded, they will also be copied to the same folder.
    """
    # Copy the runcard to the run folder
    runcard_file = REGRESSION_FOLDER / runcard_filename
    shutil.copy(runcard_file, tmp_path)

    # If weights have to be loaded, copy also the weights
    runcard_info = yaml_safe.load(runcard_file.read_text())
    if (wname := runcard_info.get("load")) is not None:
        shutil.copy(REGRESSION_FOLDER / wname, tmp_path)


@pytest.mark.parametrize("runcard,replica", RUNCARD_AND_REPLICAS.items())
def test_regression_fit(tmp_path, runcard, replica, regenerate):
    """Runs the runcard <runcard> for <replica> in the <tmp_path> folder.
    This test starts from an already ran setupfit and, often, from set weights.

    If regenerate is set to `True`, then setupfit will be also be run.
    """
    runcard_filename = f"{runcard}.yml"
    prepare_runcard(tmp_path, runcard_filename)

    # Copy setupfit, then run n3fit
    setupfit_files = SETUPFIT_FOLDER / runcard
    if not regenerate:
        if not setupfit_files.exists():
            raise FileNotFoundError(f"The setupfit folder {setupfit_files} could not be found")
        shutil.copytree(setupfit_files, tmp_path / runcard)

    run_n3fit(runcard_filename, f"{replica}", cwd=tmp_path, check=True, setupfit=regenerate)
    old_json_file = REGRESSION_FOLDER / f"{runcard}_{replica}.json"

    rel_error = extra_tolerances_rel.get(runcard, 1e-2)
    exportgrid_error = extra_tolerances_exportgrid.get(runcard, 5e-5)

    check_fit_results(
        tmp_path,
        runcard,
        replica,
        old_json_file,
        regenerate=regenerate,
        rel_error=rel_error,
        atol_exportgrid=exportgrid_error,
    )
    # When regenerating, move the setupfit folder to SETUPFIT_FOLDER to have them organized
    if regenerate:
        shutil.copytree(REGRESSION_FOLDER / runcard, setupfit_files, dirs_exist_ok=True)
