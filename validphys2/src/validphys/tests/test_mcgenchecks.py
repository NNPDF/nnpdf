"""
Tests for the mg_gen_checks module which deals with 
C++ replica generation

"""
import pandas as pd
from validphys.tests.test_covmats import CORR_DATA
from validphys.tests.test_regressions import make_table_comp
from validphys.tableloader import sane_load
from validphys.api import API

@make_table_comp(sane_load)
def test_art_rep_generation(data_config):
    config = dict(data_config)
    config["dataset_inputs"] = CORR_DATA
    config["fitting"] = {"seed": 123456}
    config["nreplica"] = 1
    _, art_replicas, _,_ = API.art_rep_generation(**config)
    return pd.DataFrame(art_replicas.T, columns=['rep0'])