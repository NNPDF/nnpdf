import numpy as np
import pandas as pd
from scipy import stats

from reportengine.table import table

from validphys.closuretest.inconsistent_closuretest.multiclosure_inconsistent import principal_components_dataset

@table
def table_bias_variance_datasets(principal_components_bias_variance_datasets, each_dataset):
    """
    TODO
    """
    records = []
    for pc_bias_var_dataset, ds in zip(principal_components_bias_variance_datasets, each_dataset):
        biases, variances, n_comp = pc_bias_var_dataset
        
        try:
            bias = np.mean(biases)
            variance = np.mean(variances)
            rbv = bias / variance
            sqrt_rbv = np.sqrt(bias / variance)
            records.append(
            dict(
                dataset=str(ds),
                dof=n_comp,
                bias=bias,
                variance=variance,
                ratio=rbv,
                ratio_sqrt=sqrt_rbv,
            )
        )
        
        except:
            records.append(
            dict(
                dataset=str(ds),
                dof=n_comp,
                bias=bias,
                variance=variance,
                ratio=np.nan,
                ratio_sqrt=np.nan,
            ))
    
        
            

    df = pd.DataFrame.from_records(
            records,
            index="dataset",
            columns=("dataset", "dof", "bias", "variance", "ratio", "ratio_sqrt"),
        )
    df.columns = ["dof", "bias", "variance", "ratio", "sqrt(ratio)"]
    return df