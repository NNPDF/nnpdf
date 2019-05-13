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
from validphys import theorycovariance as tc
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
def test_expcovmat(data):
    pdf, exps = data
    eindex = results.experiments_index(exps)
    mat = results.experiments_covmat_no_table(exps, eindex, t0set=None)
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
    return results.experiments_covmat_no_table(exps, eindex, pdf)

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

@make_table_comp(parse_exp_mat)
def test_theorycovmat(theory_data):
    pdf, exps_by_theoryid, exps_central_theory, theoryids = theory_data
    eindex = results.experiments_index(exps_central_theory)

    data_theory_centrals_1 = [results.experiment_results(exp, pdf) for exp in exps_by_theoryid[0]]
    data_theory_centrals_2 = [results.experiment_results(exp, pdf) for exp in exps_by_theoryid[1]]
    each_dataset_results_bytheory = [data_theory_centrals_1, data_theory_centrals_2]

    commondata_experiments = [ds.commondata for exp in exps_by_theoryid for ds in exp[0].datasets]

    dataset_names = tc.dataset_names(commondata_experiments)
    process_lookup = tc.process_lookup(commondata_experiments)
    combine_by_type = tc.combine_by_type(process_lookup, each_dataset_results_bytheory, dataset_names)
    covmap = tc.covmap(combine_by_type, dataset_names)
    process_starting_points = tc.process_starting_points(combine_by_type)
    covs_pt_prescrip = tc.covs_pt_prescrip(combine_by_type, process_starting_points, theoryids)

    return tc.theory_covmat_custom(covs_pt_prescrip, covmap, eindex)

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

@make_table_comp(load_fits_chi2_table)
def test_datasetchi2(single_exp_data, dataset_chi2data):
    _, exp = single_exp_data
    return results.fits_datasets_chi2_table(['test'], [[exp]], dataset_chi2data)
