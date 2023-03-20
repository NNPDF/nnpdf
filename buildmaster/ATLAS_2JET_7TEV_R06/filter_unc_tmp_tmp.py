# Filter for ATLAS_2JET_7TEV_R06
"""
Created on Mar  2023

@author: Mark N. Costantini
"""

import yaml
import numpy as np
import pandas as pd
from scipy.linalg import block_diag

from validphys.covmats import covmat_from_systematics
from validphys.loader import Loader
import copy
import warnings
warnings.filterwarnings('ignore')



def process_err(error,cv):
    """
    Converts an error given in 
    percentage of the central data
    value in absolute value.
    Follows the experimental prescription

    Parameters
    ----------

    error : dictionary
            e.g. {'label': 'sys', 'symerror': 0.1%}
    
    cv : float
        central value

    Returns
    -------


    """
    if error['label'] == 'sys':

        if 'asymerror' in error:
            d_p = float(error['asymerror']['plus'].strip('%')) / 100. * cv
            d_m = float(error['asymerror']['minus'].strip('%')) / 100. * cv

            return d_p / np.sqrt(2.), d_m / np.sqrt(2.)
        
        else:
            sigma = float(error['symerror'].strip('%')) / 100. * cv
            return sigma / np.sqrt(2.), -sigma / np.sqrt(2.)

    

def HEP_table_to_df(heptable,scenario=0):
    """
    Given hep data table return a pandas
    DataFrame with index given by Ndata,
    columns by the uncertainties and 
    np.nan entries

    Parameters
    ----------
    heptable : str
            path to hepdata table 

    scenario : 0, 1, 2
            0 = nominal, 1 = stronger, 2 = weaker
    """

    with open(heptable, 'r') as file:
        card = yaml.safe_load(file)
    
    # list of len ndata. Each entry is dict with
    # keys errors and value
    card = card['dependent_variables'][scenario]['values']
    df = pd.DataFrame(index = range(1,len(card)+1))
    
    errors = card[0]['errors']
    for i, err in enumerate(errors):
        if i==67:
            df["lum_plus"]=np.nan
            df["lum_minus"]=np.nan
        elif err['label'] == 'sys':
            df[f"{err['label']}_plus_{i}"]=np.nan
            df[f"{err['label']}_minus_{i}"]=np.nan

    return df


def fill_df(heptable,scenario=0):
    """
    Fill a data frame with index 
    corresponding to dijet mass bins
    and columns to different uncertainties
    Each df is for a fixed rapidity bin.
    """

    with open(heptable, 'r') as file:
        card = yaml.safe_load(file)
    
    # list of len ndata. Each entry is dict with
    # keys errors and value
    card = card['dependent_variables'][scenario]['values']
    df_nan = HEP_table_to_df(heptable,scenario)
    
    for i,dat in enumerate(card):
        cv = dat['value']
    
        for j, err in enumerate(dat['errors']):

            if j == 67:
                d_p, d_m = process_err(err,cv)
                df_nan.loc[df_nan.index == i+1,"lum_plus"] = d_p
                df_nan.loc[df_nan.index == i+1,"lum_minus"] = d_m

            elif err['label'] == 'sys':
                d_p, d_m = process_err(err,cv)
                df_nan.loc[df_nan.index == i+1,f"{err['label']}_plus_{j}"] = d_p
                df_nan.loc[df_nan.index == i+1,f"{err['label']}_minus_{j}"] = d_m
                
    return df_nan



def filter_ATLAS_2JET_7TEV_R06_test():
    """
    
    """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']

    dfs = []
    for table in tables:
        hepdata_tables = f"rawdata/HEPData-ins1268975-v{version}-Table_{table}.yaml"
        df = fill_df(hepdata_tables,scenario=0)
        dfs.append(df)
    
    # Construct Covariance matrix for Systematics
    Asys = dfs[0].drop(['lum_plus','lum_minus'],axis =1).to_numpy()
    for i in range(len(dfs)-1):
        Asys = np.concatenate((Asys, dfs[i+1].drop(['lum_plus','lum_minus'], axis = 1).to_numpy()), axis = 0)
    Csys = np.einsum('ij,kj->ik',Asys,Asys)

    # Construct Special Sys (Lum) Cov matrix
    Alum = dfs[0][['lum_plus','lum_minus']].to_numpy()
    for i in range(len(dfs)-1):
        Alum = np.concatenate((Alum, dfs[i+1][['lum_plus','lum_minus']].to_numpy()), axis = 0)
    Clum = np.einsum('ij,kj->ik',Alum,Alum)

    # construct Block diagonal Covariance matrix
    # for statistical errors
    ndata = [21, 21, 19, 17, 8, 4]
    stat = pd.read_csv(f"rawdata/dijet_statcov/hepcov_R06_Eta0.txt", sep = " ", header = None)
    BD_stat = stat.loc[13:,1:ndata[0]].to_numpy().astype(float)
    
    for i in range(1,6):
        stat = pd.read_csv(f"rawdata/dijet_statcov/hepcov_R06_Eta{i}.txt", sep = " ", header = None)
        stat = stat.loc[13:,1:ndata[i]].to_numpy().astype(float)    
        BD_stat = block_diag(BD_stat,stat)

    covmat = Csys + BD_stat + Clum
    return covmat



if __name__ == "__main__":
    

    setname = "ATLAS_2JET_7TEV_R06"
    l = Loader()
    cd = l.check_commondata(setname = setname).load_commondata_instance()
    dataset_input = l.check_dataset(setname,theoryid=200)

    cov = covmat_from_systematics(
        loaded_commondata_with_cuts = cd,
        dataset_input = dataset_input,
        use_weights_in_covmat=False,
        norm_threshold=None,
        _central_values=None
    )

    covmat = filter_ATLAS_2JET_7TEV_R06_test()
    
    ones = cov / covmat
    print(ones[0,:])
    print(f"min covmat/cov = {np.min(ones)}")
    print(f"max covmat/cov = {np.max(ones)}")
    print(f"mean covmat / cov = {np.mean(ones)}")
    print(f"np.allclose: {np.allclose(covmat,cov, rtol = 0.001, atol =1e-5)}")
    
