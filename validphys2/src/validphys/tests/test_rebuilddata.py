"""
test_rebuilddata.py

filter some simple closure data and then check that rebuild data runs on it
and produces expected results

"""
import pathlib
import shutil
import subprocess as sp

import pytest

from n3fit.scripts.n3fit_exec import N3FitConfig
from reportengine import api
from validphys.app import providers
from validphys.config import Environment
from validphys.scripts.vp_rebuild_data import REBUILD_CONFIG
from validphys.tableloader import sane_load
from validphys.tests.test_regressions import make_table_comp

FIT_NAME = "dummy_closure_runcard"

REGRESSION_FOLDER = pathlib.Path(__file__).with_name("regressions")


def parse_test_output(filename):
    """Parse the output of groups_data_values."""
    df = sane_load(filename, header=0, index_col=[0, 1, 2])
    return df


@pytest.mark.linux
@make_table_comp(parse_test_output)
def test_filter_rebuild_closure_data(tmp):
    """
    Takes a closure test runcard from the regressions directory
    and then runs ``vp-setupfit`` in a temp directory and then
    ``vp-rebuild-data`` on the resulting fit folder.

    The test then loads the filtered and rebuilt data and checks that the
    experimental central values (which generated during ``vp-setupfit``)
    take on the expected values.

    """
    runcard_name = FIT_NAME + ".yaml"
    runcard = REGRESSION_FOLDER / runcard_name

    # cp runcard to tmp folder
    shutil.copy(runcard, tmp)
    # filter the runcard
    sp.run(f"vp-setupfit {runcard_name}".split(), cwd=tmp, check=True)

    sp.run(f"vp-rebuild-data {FIT_NAME}".split(), cwd=tmp, check=True)

    API = api.API(providers, N3FitConfig, Environment, output=str(tmp / FIT_NAME))
    # to use groups_data_values we need to do some gymnastics with data spec
    # because taking data_input from fitinputcontext overwrites any subsequent
    # grouping.
    input_params = dict(REBUILD_CONFIG)
    input_params.pop("data_input")
    input_params["dataset_inputs"] = {"from_": "fit"}
    df = API.groups_data_values(**input_params, pdf="NNPDF31_nnlo_as_0118")
    return df.to_frame()
