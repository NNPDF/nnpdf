"""
Regression tests
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
runcard_and_replicas = {
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


@pytest.mark.parametrize("runcard,replica", runcard_and_replicas.items())
def test_regression_fit(tmp_path, runcard, replica, regenerate):
    """Runs the runcard <runcard> for <replica> in the <tmp_path> folder.
    This test starts from an already ran setupfit and, often, from set weights.

    If regenerate is set to `True`, then setupfit will be also be run.
    """
    # Copy the runcard to the run folder
    runcard_name = f"{runcard}.yml"
    runcard_file = REGRESSION_FOLDER / runcard_name
    shutil.copy(runcard_file, tmp_path)

    # If weights have to be loaded, copy also the weights
    runcard_info = yaml_safe.load(runcard_file.read_text())
    if (wname := runcard_info.get("load")) is not None:
        shutil.copy(REGRESSION_FOLDER / wname, tmp_path)

    # Copy setupfit, then run n3fit
    setupfit_files = SETUPFIT_FOLDER / runcard
    if not regenerate:
        if not setupfit_files.exists():
            raise FileNotFoundError(f"The setupfit folder {setupfit_files} could not be found")
        shutil.copytree(setupfit_files, tmp_path / runcard)

    run_n3fit(runcard_name, f"{replica}", cwd=tmp_path, check=True, setupfit=regenerate)
    old_json_file = REGRESSION_FOLDER / f"{runcard}_{replica}.json"

    check_fit_results(
        tmp_path, runcard, replica, old_json_file, regenerate=regenerate, rel_error=1e-2
    )
    # When regenerating, move the setupfit folder to SETUPFIT_FOLDER to have them organized
    if regenerate:
        shutil.copytree(REGRESSION_FOLDER / runcard, setupfit_files, dirs_exist_ok=True)
