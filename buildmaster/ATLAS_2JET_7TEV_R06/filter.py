# Filter for ATLAS_2JET_7TEV_R06
"""
Created on Mar  2023

@author: Mark N. Costantini
"""

import yaml
import numpy as np
import pandas as pd
from scipy.linalg import block_diag, cholesky


from validphys.covmats import dataset_inputs_covmat_from_systematics
from validphys.loader import Loader
import warnings
warnings.filterwarnings('ignore')



def filter_ATLAS_2JET_7TEV_R06(scenario=0,kin_data=True):
    """
    Parameters
    ----------

    scenario : 0, 1, 2
            0 = nominal, 1 = stronger, 2 = weaker
            default is 0
    
    kin_data : bool
            whether to save kin and data
    
    """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']

    data_central = []
    kin = []
    dfs = []
    for table in tables:
        hepdata_tables = f"rawdata/HEPData-ins1268975-v{version}-Table_{table}.yaml"
        
        # uncertainties dataframe
        df = fill_df(hepdata_tables,scenario=scenario)
        dfs.append(df)

        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
        
        ystar = range_str_to_floats(input['dependent_variables'][0]['qualifiers'][8]['value'])
        sqrts = float(input['dependent_variables'][0]['qualifiers'][7]['value'])

        # measurements of several dijet mass values
        values = input['dependent_variables'][0]['values']

        # loop over different diijet mass bins
        for i,value in enumerate(values):

            data_central_value = value['value']
            data_central.append(data_central_value)
              
            m12 = input['independent_variables'][0]['values'][i]
            m12['low'], m12['high'] = 1e3 * m12['low'], 1e3 * m12['high']
            m12['mid'] = float(f"{0.5 * (m12['low']+m12['high']):.3f}")

            kin_value = {'ystar' : {'min': ystar['min'], 'mid': ystar['mid'] , 'max': ystar['max']}, 
                        'm12' : {'min': m12['low'], 'mid': m12['mid'] , 'max': m12['high']} ,
                         'sqrt_s' : {'min': None, 'mid': sqrts , 'max': None}}

            kin.append(kin_value)
    
    
    # Construct Covariance matrix for Systematics
    Asys = pd.concat([df.drop(['lum_plus','lum_minus'],axis=1) for df in dfs], axis = 0).to_numpy()
    Csys = np.einsum('ij,kj->ik',Asys,Asys)

    # Construct Special Sys (Lum) Cov matrix
    Alum = pd.concat([df[['lum_plus','lum_minus']] for df in dfs], axis = 0).to_numpy()
    Clum = np.einsum('ij,kj->ik',Alum,Alum)

    # construct Block diagonal statistical Covariance matrix
    ndata = [21, 21, 19, 17, 8, 4]
    stat = pd.read_csv(f"rawdata/dijet_statcov/hepcov_R06_Eta0.txt", sep = " ", header = None)
    BD_stat = stat.loc[13:,1:ndata[0]].to_numpy().astype(float)
    
    for i in range(1,6):
        stat = pd.read_csv(f"rawdata/dijet_statcov/hepcov_R06_Eta{i}.txt", sep = " ", header = None)
        stat = stat.loc[13:,1:ndata[i]].to_numpy().astype(float)    
        BD_stat = block_diag(BD_stat,stat)

    # full covariance matrix
    covmat = Csys + BD_stat + Clum

    # generate artificial systematics 
    A_art_sys = decompose_covmat(covmat=covmat)
    
    # error definition
    error_definition = {f"art_sys_{i}" : {"description":  f"artificial systematic {i}",
                                        "treatment": "ADD", "type": "CORR"}
                        for i in range(1,A_art_sys.shape[0]+1)}

    # store error in dict
    error = []
    for n in range(A_art_sys.shape[0]):
        error_value={}
        for m in range(A_art_sys.shape[1]):
            error_value[f"art_sys_{m+1}"] = float(A_art_sys[n,m])
        error.append(error_value)

    
    uncertainties_yaml = {"definition": error_definition, "bins": error}
    data_central_yaml  = { 'data_central' : data_central }
    kinematics_yaml    = { 'bins' : kin }

    if kin_data:
        with open('data.yaml', 'w') as file:
            yaml.dump(data_central_yaml, file, sort_keys=False)

        with open('kinematics.yaml', 'w') as file:
            yaml.dump(kinematics_yaml, file, sort_keys=False)

    if scenario==0:
        with open("uncertainties.yaml",'w') as file:
            yaml.dump(uncertainties_yaml,file, sort_keys=False)
    elif scenario==1:
        with open("uncertainties_stronger.yaml",'w') as file:
            yaml.dump(uncertainties_yaml,file, sort_keys=False)
    elif scenario==2:
        with open("uncertainties_weaker.yaml",'w') as file:
            yaml.dump(uncertainties_yaml,file, sort_keys=False)

    return covmat


def range_str_to_floats(str_range):
    """
    converts a string range to a list,
    e.g. "0.5 - 1.0" --> [0.5,1.0]
    """
    # Split the range string into two parts
    str_nums = str_range.split('-')
    # Convert the parts to floats
    min = float(str_nums[0])
    max = float(str_nums[1])
    mid = float(f"{0.5 * (min + max):.3f}")
    # Return a dict
    return {"min": min, "mid": mid, "max": max}

def decompose_covmat(covmat):
     """Given a covmat it return an array sys with shape (ndat,ndat)
     giving ndat correlated systematics for each of the ndat point.
     The original covmat is obtained by doing sys@sys.T"""

     lamb, mat = np.linalg.eig(covmat)
     sys = np.multiply(np.sqrt(lamb), mat)
     return sys

def process_err(error,cv):
    """
    Converts an error given in percentage
    of the central data value in absolute value.

    Note: the d'Agostini prescription for the
    symmetrization of the error does not hold here.
    We follow here the experimental prescription

    Parameters
    ----------

    error : dictionary
            e.g. {'label': 'sys', 'symerror': 0.1%}
    
    cv : float
        central value

    Returns
    -------
    tuple
        tuple containing two floats

    """
    if error['label'] == 'sys':

        if 'asymerror' in error:
            d_p = float(error['asymerror']['plus'].strip('%')) / 100. * cv
            d_m = float(error['asymerror']['minus'].strip('%')) / 100. * cv

            
            tmp1 = d_p
            tmp2 = d_m
            # case 1: d_p and d_m are both negative
            if (tmp1<0.0 and tmp2<0.0):
                if(tmp2<tmp1):    
                    d_p = 0.0
                    d_m = tmp2
                else:    
                    d_p = 0.0
                    d_m = tmp1
            
            # case 2: d_p and d_m are both positive
            if (tmp1>0.0 and tmp2>0.0):
                if(tmp1>tmp2):    
                    d_p = tmp1
                    d_m = 0.0
                else:
                    d_p = tmp2
                    d_m = 0.0
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
        if (scenario == 0 and i == 67) or (scenario == 1 and i == 57) or (scenario == 2 and i == 69):
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

            if (scenario == 0 and j == 67) or (scenario == 1 and j == 57) or (scenario == 2 and j == 69):
                d_p, d_m = process_err(err,cv)
                df_nan.loc[df_nan.index == i+1,"lum_plus"] = d_p
                df_nan.loc[df_nan.index == i+1,"lum_minus"] = d_m

            elif err['label'] == 'sys':
                d_p, d_m = process_err(err,cv)
                df_nan.loc[df_nan.index == i+1,f"{err['label']}_plus_{j}"] = d_p
                df_nan.loc[df_nan.index == i+1,f"{err['label']}_minus_{j}"] = d_m
                
    return df_nan







if __name__ == "__main__":
    
    # generate all uncertainties scenarios
    filter_ATLAS_2JET_7TEV_R06(scenario=0,kin_data=True)
    filter_ATLAS_2JET_7TEV_R06(scenario=1, kin_data=False)
    filter_ATLAS_2JET_7TEV_R06(scenario=2, kin_data=False)
    setname = "ATLAS_2JET_7TEV_R06"
    l = Loader()
    cd = l.check_commondata(setname = setname).load_commondata_instance()
    dataset_input = l.check_dataset(setname,theoryid=200)
    from validphys.commondataparser import parse_commondata
    dat_file = '/Users/markcostantini/codes/nnpdfgit/nnpdf/buildmaster/results/DATA_ATLAS_2JET_7TEV_R06.dat'
    sys_file = '/Users/markcostantini/codes/nnpdfgit/nnpdf/buildmaster/results/systypes/SYSTYPE_ATLAS_2JET_7TEV_R06_DEFAULT.dat'
    
    cd = parse_commondata(dat_file,sys_file,setname)

    cmat = dataset_inputs_covmat_from_systematics(
        dataset_inputs_loaded_cd_with_cuts=[cd],
        data_input = [dataset_input],
        use_weights_in_covmat=False,
        norm_threshold=None,
        _list_of_central_values=None,
        _only_additive=False,
             )


    covmat = filter_ATLAS_2JET_7TEV_R06(scenario=0,kin_data=False)
    
    ones = cmat / covmat
    print(ones[0,:])
    print(f"min covmat/cov = {np.min(ones)}")
    print(f"max covmat/cov = {np.max(ones)}")
    print(f"mean covmat / cov = {np.mean(ones)}")
    print(f"np.allclose: {np.allclose(covmat,cmat, rtol = 1e05, atol =1e-5)}")
    
