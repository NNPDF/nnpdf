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
from validphys.tableloader import (parse_exp_mat, load_perreplica_chi2_table,
                                   sane_load)



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
def test_expcovmat(data):
    pdf, exps = data
    eindex = results.experiments_index(exps)
    mat = results.experiments_covmat(exps, eindex, t0set=None)
    covmats = []
    for exp in exps:
        cd = exp.datasets[0].commondata.load()
        covmats.append(NNPDF.ComputeCovMat(cd, cd.get_cv()))
    othermat = la.block_diag(*covmats)
    assert np.allclose(mat.values, othermat)
    return mat

@make_table_comp(parse_exp_mat)
def test_t0covmat(data):
    pdf, exps = data
    eindex = results.experiments_index(exps)
    return results.experiments_covmat(exps, eindex, pdf)

@make_table_comp(parse_exp_mat)
def test_expsqrtcovmat(data):
    pdf, exps = data
    eindex = results.experiments_index(exps)
    return results.experiments_sqrtcovmat(exps, eindex, t0set=None)

@make_table_comp(parse_exp_mat)
def test_t0sqrtcovmat(data):
    pdf, exps = data
    eindex = results.experiments_index(exps)
    return results.experiments_sqrtcovmat(exps, eindex, pdf)


@make_table_comp(sane_load)
def test_predictions(convolution_results):
    ths = []
    for convolution_result in convolution_results:
        dt, th = convolution_result
        ths.append(th._rawdata.astype(float))
    th = np.concatenate(ths)
    return pd.DataFrame(th,
        columns=map(str,
        range(th.shape[1])))

@make_table_comp(sane_load)
def test_dataset_t0_predictions(dataset_t0_convolution_results):
    ths = []
    for convolution_result in dataset_t0_convolution_results:
        dt, th = convolution_result
        ths.append(th._rawdata.astype(float))
    th = np.concatenate(ths)
    return pd.DataFrame(th,
        columns=map(str,
        range(th.shape[1])))

@make_table_comp(sane_load)
def test_cv(convolution_results):
    cvs = []
    for convolution_result in convolution_results:
        dt, _ = convolution_result
        cvs.append(dt.central_value)
    data_values = np.concatenate(cvs)
    return pd.DataFrame(data_values, columns=['CV'])

@make_table_comp(load_perreplica_chi2_table)
def test_replicachi2data(data, chi2data):
    pdf, exps = data
    return results.perreplica_chi2_table(exps, chi2data)

@make_table_comp(sane_load)
def test_datasetchi2(data, chi2data, dataset_chi2data):
    pdf, exps = data
    return results.experiments_chi2_table(exps, pdf, chi2data, dataset_chi2data)
