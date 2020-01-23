"""
closuretest/multiclosure.py

containing all of the statistical estimators which are
averaged across multiple fits or a single replica proxy fit
"""

import numpy as np
import pandas as pd

from reportengine import collect
from reportengine.table import table

fits_bias_experiment = collect("bias_experiment", ("fits", "fitcontext"))

def expected_bias(fits_bias_experiment):
    bias_centrals = [
        biasdata.bias*biasdata.ndata for biasdata in fits_bias_experiment]
    #ndata = fits_bias_experiment[0].ndata # ndata should be same for all fits
    bias_mean = np.mean(bias_centrals)
    return bias_mean

exps_expected_bias = collect("expected_bias", ("experiments",))

fits_phi_experiment = collect("phi_data_experiment", ("fits", "fitcontext"))

def expected_variance(fits_phi_experiment):
    var_centrals = [
        phi[0]**2 * phi[1] for phi in fits_phi_experiment]
    var_mean = np.mean(var_centrals)
    return var_mean

exps_expected_var = collect("expected_variance", ("experiments",))

@table
def experiments_bias_variance_ratio(
    exps_expected_bias, exps_expected_var, experiments,):
    records = []
    bias_tot = 0
    var_tot = 0
    for exp, bias, var in zip(experiments, exps_expected_bias, exps_expected_var):
        records.append(dict(
            experiment=str(exp),
            ratio=bias/var
        ))
        bias_tot += bias
        var_tot += var
    records.append(dict(
            experiment="total",
            ratio=bias_tot/var_tot
    ))
    df = pd.DataFrame.from_records(
        records,
        index="experiment",
        columns=("experiment", "ratio")
    )
    df.columns = ["bias/variance"]
    return df
