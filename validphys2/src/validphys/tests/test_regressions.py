"""
test_regression.py

Write files with data to disk and assert they are the same upon
updates.
"""
import pathlib
import logging
import functools

import numpy as np
import scipy.linalg as la
import pandas as pd
from pandas.testing import assert_frame_equal

from reportengine.table import savetable

import NNPDF
from validphys import results
from validphys.api import API
from validphys.tableloader import (parse_exp_mat, load_perreplica_chi2_table,
                                   sane_load, load_fits_chi2_table)



log = logging.getLogger(__name__)

REGRESSION_FOLDER = pathlib.Path(__file__).with_name('regressions')

#TODO: Move these to a library
def compare_tables(produced_table, storage_path, loader_func):
    """Test that the ``produced_table`` is equal (as in allclose) to
    the one loaded from the `storage_path` using the `loader_func`"""
    if not storage_path.exists():
        savetable(produced_table, storage_path)
        #Fail test
        assert False, "Storage path does not exist"
    stored_table = loader_func(storage_path)
    assert_frame_equal(produced_table, stored_table)

def make_table_comp(loader_func):
    """Compare the dataframe that the decorated function outputs with
    a file with the same name as the function and extension csv, loaded
    using the provided `loader_func`"""
    def decorator(f):
        @functools.wraps(f)
        def f_(*args, **kwargs):
            filename = f'{f.__name__}.csv'
            produced_table = f(*args, **kwargs)
            compare_tables(produced_table, REGRESSION_FOLDER/filename, loader_func)
        return f_
    return decorator

@make_table_comp(parse_exp_mat)
def test_expcovmat(data_config):
    mat = API.experiments_covmat_no_table(**data_config)
    covmats = []
    for exp in API.experiments(**data_config):
        cd = exp.datasets[0].commondata.load()
        covmats.append(NNPDF.ComputeCovMat(cd, cd.get_cv()))
    othermat = la.block_diag(*covmats)
    assert np.allclose(mat.values, othermat)
    return mat

@make_table_comp(parse_exp_mat)
def test_t0covmat(data_witht0_config):
    return API.experiments_covmat_no_table(**data_witht0_config)

@make_table_comp(parse_exp_mat)
def test_expsqrtcovmat(data_config):
    return API.experiments_sqrtcovmat(**data_config)

@make_table_comp(parse_exp_mat)
def test_t0sqrtcovmat(data_witht0_config):
    return API.experiments_sqrtcovmat(**data_witht0_config)


@make_table_comp(sane_load)
def test_predictions(data_config):
    # TODO: ideally we would change the baseline to just be corresponding columns
    # of `experiment_result_table`, however sane_load expects just a single level
    # of column and index - if we use a different format like parquet this could
    # be changed.
    exp_res_tab = API.experiment_result_table_no_table(**data_config)
    th = exp_res_tab.iloc[:, 2:].values
    return pd.DataFrame(th, columns=map(str, range(th.shape[1])))

@make_table_comp(sane_load)
def test_dataset_t0_predictions(data_witht0_config):
    # TODO: As in `test_predictions`
    exp_res_tab = API.experiment_result_table_no_table(**data_witht0_config)
    th = exp_res_tab.iloc[:, 2:].values
    return pd.DataFrame(th, columns=map(str, range(th.shape[1])))

@make_table_comp(sane_load)
def test_cv(data_config):
    # TODO: As in `test_predictions`
    exp_res_tab = API.experiment_result_table_no_table(**data_config)
    data_values = exp_res_tab.iloc[:, 0].values[:, np.newaxis]
    return pd.DataFrame(data_values, columns=['CV'])

@make_table_comp(load_perreplica_chi2_table)
def test_replicachi2data(data_witht0_config):
    return API.perreplica_chi2_table(**data_witht0_config)

@make_table_comp(load_fits_chi2_table)
def test_datasetchi2(data_singleexp_witht0_config):
    # This is a bit hacky but avoids requiring a fit
    exps = API.experiments(**data_singleexp_witht0_config)
    chi2s = API.each_dataset_chi2(**data_singleexp_witht0_config)
    return results.fits_datasets_chi2_table(['test'], [exps], chi2s)
