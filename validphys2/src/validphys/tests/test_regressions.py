"""
test_regression.py

Write files with data to disk and assert they are the same upon
updates.
"""
import pathlib
import logging
import functools

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from reportengine.table import savetable

from validphys import results
from validphys.api import API
from validphys.tests.test_covmats import CORR_DATA
from validphys.tableloader import (parse_data_cv, parse_exp_mat, load_perreplica_chi2_table,
                                   sane_load, load_fits_chi2_table)
from validphys.tests.test_covmats import CORR_DATA



log = logging.getLogger(__name__)

REGRESSION_FOLDER = pathlib.Path(__file__).with_name('regressions')

#TODO: Move these to a library
def compare_tables(produced_table, storage_path, loader_func, tolerance=1e-8):
    """Test that the ``produced_table`` is equal (as in allclose) to
    the one loaded from the `storage_path` using the `loader_func`"""
    if not storage_path.exists():
        savetable(produced_table, storage_path)
        #Fail test
        assert False, "Storage path does not exist"
    stored_table = loader_func(storage_path)
    assert_frame_equal(produced_table, stored_table, atol=tolerance)

def make_table_comp(loader_func, tolerance=1e-8):
    """Compare the dataframe that the decorated function outputs with
    a file with the same name as the function and extension csv, loaded
    using the provided `loader_func`"""
    def decorator(f):
        @functools.wraps(f)
        def f_(*args, **kwargs):
            filename = f'{f.__name__}.csv'
            produced_table = f(*args, **kwargs)
            compare_tables(produced_table, REGRESSION_FOLDER/filename, loader_func, tolerance=tolerance)
        return f_
    return decorator


@make_table_comp(parse_data_cv)
def test_mcreplica(data_config):
    config = dict(data_config)
    config["dataset_inputs"] = CORR_DATA
    seed = 123456
    # Use no cuts because if filter rules change in the
    # future then this test will end up failing
    rep = API.indexed_make_replica(**config, replica_mcseed=seed)
    return rep


@make_table_comp(parse_exp_mat)
def test_expcovmat(data_config):
    return API.groups_covmat_no_table(**data_config)


@make_table_comp(parse_exp_mat)
def test_t0covmat(data_witht0_internal_cuts_config):
    return API.groups_covmat_no_table(**data_witht0_internal_cuts_config)


@make_table_comp(parse_exp_mat)
def test_expsqrtcovmat(data_config):
    return API.groups_sqrtcovmat(**data_config)


@make_table_comp(parse_exp_mat)
def test_t0sqrtcovmat(data_witht0_internal_cuts_config):
    return API.groups_sqrtcovmat(**data_witht0_internal_cuts_config)


@make_table_comp(parse_exp_mat)
def test_pdf_plus_exp_covmat(data_internal_cuts_config):
    return API.groups_covmat_no_table(use_pdferr=True, **data_internal_cuts_config)


@make_table_comp(parse_exp_mat)
def test_hessian_pdf_plus_exp_covmat(hessian_data_internal_cuts_config):
    return API.groups_covmat_no_table(use_pdferr=True, **hessian_data_internal_cuts_config)


@make_table_comp(sane_load)
def test_predictions(data_internal_cuts_config):
    # TODO: ideally we would change the baseline to just be corresponding columns
    # of `experiment_result_table`, however sane_load expects just a single level
    # of column and index - if we use a different format like parquet this could
    # be changed.
    res_tab = API.group_result_table_no_table(**data_internal_cuts_config)
    th = res_tab.iloc[:, 2:].values
    return pd.DataFrame(th, columns=map(str, range(th.shape[1])), dtype='float64')

@make_table_comp(sane_load)
def test_thprediction_results(single_data_internal_cuts_config):
    """Test the central prediction and the resulting std deviation for a MC PDF"""
    pdf = API.pdf(**single_data_internal_cuts_config)
    dataset = API.dataset(**single_data_internal_cuts_config)
    res = results.ThPredictionsResult.from_convolution(pdf, dataset)
    tp = np.stack([res.central_value, res.std_error]).T
    return pd.DataFrame(tp, columns=map(str, range(tp.shape[1])), dtype='float64')

@make_table_comp(sane_load)
def test_thprediction_results_hessian(hessian_single_data_internal_cuts_config):
    """Test the central prediction and the resulting std deviation for a hessian PDF"""
    pdf = API.pdf(**hessian_single_data_internal_cuts_config)
    dataset = API.dataset(**hessian_single_data_internal_cuts_config)
    res = results.ThPredictionsResult.from_convolution(pdf, dataset)
    tp = np.stack([res.central_value, res.std_error]).T
    return pd.DataFrame(tp, columns=map(str, range(tp.shape[1])), dtype='float64')

@make_table_comp(sane_load)
def test_dataset_t0_predictions(data_witht0_internal_cuts_config):
    # TODO: As in `test_predictions`
    res_tab = API.group_result_table_no_table(**data_witht0_internal_cuts_config)
    th = res_tab.iloc[:, 2:].values
    return pd.DataFrame(th, columns=map(str, range(th.shape[1])), dtype='float64')

@make_table_comp(sane_load)
def test_cv(data_internal_cuts_config):
    # TODO: As in `test_predictions`
    res_tab = API.group_result_table_no_table(**data_internal_cuts_config)
    data_values = res_tab.iloc[:, 0].values[:, np.newaxis]
    return pd.DataFrame(data_values, columns=['CV'])

@make_table_comp(load_perreplica_chi2_table)
def test_replicachi2data(data_witht0_internal_cuts_config):
    return API.perreplica_chi2_table(**data_witht0_internal_cuts_config)

@make_table_comp(load_fits_chi2_table)
def test_datasetchi2(data_singleexp_witht0_config):
    # This is a bit hacky but avoids requiring a fit
    exps = API.groups_data(**data_singleexp_witht0_config)
    chi2s = API.groups_datasets_chi2_data(**data_singleexp_witht0_config)
    return results.fits_datasets_chi2_table(['test'], [exps], [chi2s])

@make_table_comp(sane_load)
def test_art_rep_generation(data_config):
    config = dict(data_config)
    config["dataset_inputs"] = CORR_DATA
    config["mcseed"] = 123456
    config["genrep"] = True
    config["nreplica"] = 1
    _, art_replicas, _,_ = API.art_rep_generation(**config)
    return pd.DataFrame(art_replicas.T, columns=['rep0'])
