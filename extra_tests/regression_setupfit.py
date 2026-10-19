"""
Regression tests for vp-setupfit.
The thresholds set in this file are based on the empirical numerical variability
we see among runs.
"""

import numpy as np
import pandas as pd
import pytest

from n3fit.tests.helpers import run_setupfit

from .regression_checks import RUNCARD_AND_REPLICAS, SETUPFIT_FOLDER, prepare_runcard

# Some runcards need to be a bit more lenient with the tolerances
# WARNING: the tolerance on no_diagonal is way too high
extra_tolerances_exportgrid = {"hyperopt_sampling": 1e-4, "t0theoryid": 1e-3, "no_diagonal": 0.01}
extra_tolerances_rel = {"hyperopt_sampling": 3e-2}


@pytest.mark.parametrize("runcard", RUNCARD_AND_REPLICAS.keys())
def test_setupfit(tmp_path, runcard, rtol=1e-3):
    """Runs vp-setupfit for runcard <runcard>.
    Checks the exact same files are generated as the ones currently saved in SETUPFIT_FOLDER.

    This test do not accept --regenerate as all regeneration is done by the fit's regression
    for consistency (and to ensure they run in order).

    The relative tolerance allowance is empirically chosen.
    """
    runcard_filename = f"{runcard}.yml"
    prepare_runcard(tmp_path, runcard_filename)
    run_setupfit(runcard_filename, cwd=tmp_path, no_eko_download=True)

    # We are going to check _every_ file in the tracked folder
    runfolder = SETUPFIT_FOLDER / runcard
    for tracked in runfolder.rglob("*"):
        # First check that whatever you find there it exists also in the new one
        newfile = runfolder / tracked.relative_to(runfolder)
        if not newfile.exists():
            raise FileNotFoundError(
                f"Was expecting {newfile} for {runcard} but it has not been generated"
            )
        # Then check the files have the same content unless they are _not_ files or are .csv
        if not tracked.is_file():
            continue
        # Then check whether the _content_ is identical
        new_content = newfile.read_text()
        old_content = tracked.read_text()
        if new_content == old_content:
            continue
        # If it is not, we can accept differences but _only_ for CSV files
        if tracked.suffix == ".csv":
            # TODO read the two csv files and check whether the content is ok up to some tolerance
            if tracked.suffix == ".csv":
                old_df = pd.read_csv(tracked, header=None)
                new_df = pd.read_csv(newfile, header=None)
                if not np.allclose(old_df, new_df, atol=1e-5, rtol=rtol):
                    raise ValueError(f"The content of {newfile} has changed beyond tolerance")
        else:
            raise ValueError(f"The content of {newfile} has changed")
