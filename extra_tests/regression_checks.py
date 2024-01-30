"""
    Regression tests
"""

from n3fit.tests.test_fit import check_fit_results, EXE
import pathlib
from reportengine.compat import yaml
import shutil
import subprocess as sp

import pytest



REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regression_fits")

# Avoid always round-number replicas or 1/2
runcard_and_replicas = {
        "normal_fit": 72,
        "central": 16,
        "diagonal": 45,
        "feature_scaling": 81,
        "flavour": 29,
        "no_msr": 92,
        "no_sumrules": 18,
        "no_vsr": 54,
        "trainable_prepro": 61,
        "no_lagrange": 27
        }

@pytest.mark.parametrize("runcard,replica", runcard_and_replicas.items())
def test_regression_fit(tmp_path, runcard, replica):
    runcard_name = f"{runcard}.yml"
    runcard_file = REGRESSION_FOLDER / runcard_name
    shutil.copy(runcard_file, tmp_path)

    runcard_info = yaml.load(runcard_file.read_text())
    if (wname := runcard_info.get("load")) is not None:
        shutil.copy(REGRESSION_FOLDER / wname, tmp_path)

    sp.run(f"{EXE} {runcard_name} {replica}".split(), cwd=tmp_path, check=True)
    old_json_file = REGRESSION_FOLDER / f"{runcard}_{replica}.json"

    check_fit_results(tmp_path, runcard, replica, old_json_file, regenerate=False, rel_error=1e-2)
