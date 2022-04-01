"""
test_covmatreg

functions to test the covariance matrix regularization
"""
import numpy as np
import pandas as pd

from validphys.tests.test_regressions import make_table_comp
from validphys.tableloader import parse_exp_mat
from validphys.calcutils import regularize_covmat, regularize_l2
from validphys.api import API

def test_withidentity():
    identity_mat = np.eye(3)
    np.testing.assert_allclose(identity_mat, regularize_l2(identity_mat, 1))
    np.testing.assert_allclose(identity_mat, regularize_covmat(identity_mat, 1))
    np.testing.assert_allclose(2*identity_mat, regularize_l2(identity_mat, 0.5))

@make_table_comp(parse_exp_mat)
def test_regularize_expcov(data_config):
    """Test if the higher level covmat is regularized by procedure"""
    inp = dict(**data_config, norm_threshol=3)
    df1 = API.groups_covmat(**inp)
    df2 = API.groups_covmat(**data_config)
    # check here that regularization occured
    assert ~np.allclose(df1.values, df2.values)
    # check that square of sqrt matches
    sqrt_df1 = API.groups_sqrtcovmat(**inp)
    np.testing.assert_allclose(df1.values, sqrt_df1.values@sqrt_df1.values.T)
    # check that same result obtained
    return df1

@make_table_comp(parse_exp_mat)
def test_regularized_covmat_generation(data_config):
    covmat = API.dataset_inputs_covariance_matrix(**data_config, norm_threshold=3)
    index = API.groups_index(**data_config)
    return pd.DataFrame(covmat, index=index, columns=index)

def test_regularization_matches_sane():
    """Check that regularizing the sqrt cov produces same result as regularizing
    on the covariance matrix with a fairly sane matrix
    """
    # choose some sensible matrix to be regularized
    a = np.ones((3, 3)) + 0.5*np.diag(np.ones(3))
    cov = a@a.T
    a_reg = regularize_l2(a, 3)
    np.testing.assert_allclose(regularize_covmat(cov, 3), a_reg@a_reg.T)

def test_regularization_matches_zero_eig():
    """Check that regularizing the sqrt cov produces the same result as regularizing
    on the covmat with a matrix with almost zero eigenvalue
    """
    # choose matrix with ~zero eigenvalue
    a = np.arange(9).reshape(3, 3)
    cov = a@a.T
    a_reg = regularize_l2(a, 3)
    np.testing.assert_allclose(regularize_covmat(cov, 3), a_reg@a_reg.T)


@make_table_comp(parse_exp_mat)
def test_no_regularization(data_config):
    """Test if the higher level covmat is regularized by procedure"""
    inp = dict(**data_config, norm_threshol=10)
    df1 = API.groups_covmat(**inp)
    df2 = API.groups_covmat(**data_config)
    # check here that regularization occured
    np.testing.assert_allclose(df1.values, df2.values)
    # check that same result obtained
    return df1
