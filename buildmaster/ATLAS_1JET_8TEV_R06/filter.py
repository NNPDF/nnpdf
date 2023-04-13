import yaml
import numpy as np
import pandas as pd
from scipy.linalg import block_diag
from filter_utils import (
                            get_data_values, get_kinematics,
                            get_stat_errors, get_lumi_errors,
                            fill_df, decompose_covmat
                        )

# ignore pandas warning 
import warnings
warnings.filterwarnings('ignore')


def filter_ATLAS_1JET_8TEV_data_kinetic():
    """
    write kinematics and central data values
    in kinematics.yaml and data.yaml files
    respectively.
    """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']

    # get kinematics from hepdata tables
    kin = get_kinematics(tables,version)

    # get central values from hepdata tables
    data_central = get_data_values(tables,version)
    
    data_central_yaml  = { 'data_central' : data_central }
    kinematics_yaml    = { 'bins' : kin }

    # write central values and kinematics to yaml file
    with open('data.yaml', 'w') as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

def filter_ATLAS_1JET_8TEV_uncertainties():
    """
    write uncertainties to uncertainties.yaml
    file.
    There are three types of uncertainties:

    1. Statistical Uncertainties: CORR
       -> correlated over the full dataset

    2. Artificial Uncertainties: CORR, these
       are obtained by following the steps:

       - Construct an Error matrix in which 
         each part of an asymmetric error is considered
         as as separate error (hence dividing by sqrt(2))
         see also filter_utils/process_error and 
         filter_utils/HEP_table_to_df

       - Construct covariance matrix from the Error matrix

       - Decompose covariance matrix so as to get ndat unc

    3. Luminosity Uncertainty: ATLASLUMI12
       this uncertainty is correlated with all 
       the other ATLASLUMI12 datasets
    """
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']
    
    # get df of uncertainties
    dfs = []
    for table in tables:
        # uncertainties dataframe
        df = fill_df(table,version)
        dfs.append(df)

    df_unc = pd.concat([df for df in dfs], axis = 0)
    
    # statistical errors fully uncorrelated
    stat_errors = df_unc['stat'].to_numpy()

    # BD_stat = np.einsum('i,j->ij',dfs[0]['stat'],dfs[0]['stat'])
    # for i in range(1,len(dfs)):
    #     cov = np.einsum('i,j->ij',dfs[i]['stat'].to_numpy(),dfs[i]['stat'].to_numpy())
    #     BD_stat = block_diag(BD_stat,cov)

    # luminosity errors
    lum_errors = df_unc["syst_lumi"].to_numpy()
    

    # col_to_drop =  ("syst_JES_EtaIntercalibration_Stat98", "syst_JES_MJB_Alpha",
    #                 "syst_JES_Pileup_MuOffset")
    # df_unc = df_unc.loc[:, ~df.columns.str.startswith(col_to_drop)]
    
    A_corr = df_unc.drop(["stat","syst_lumi"],axis =1).to_numpy() / np.sqrt(2.)
    cov_corr = np.einsum("ij,kj->ik",A_corr,A_corr)
    A_art_corr = decompose_covmat(covmat=cov_corr)
    cov_corr1 = np.einsum("ij,kj->ik",A_art_corr,A_art_corr)
    print(f"same cov_corr mat? :{np.allclose(cov_corr,cov_corr1)}")

     # error definition
    error_definition = {f"art_sys_corr_{i}" : {"description":  f"artificial systematic {i}",
                                        "treatment": "ADD", "type": "CORR"}
                        for i in range(1,A_art_corr.shape[0]+1)}

    error_definition["luminosity_uncertainty"] = {"description": "luminosity uncertainty",
                                                    "treatment": "ADD", "type": "ATLASLUMI12"}
    
    error_definition["statistical_uncertainty"] = {"description": "statistical uncertainty",
                                                    "treatment": "ADD", "type": "CORR"}
    
    # store error in dict
    error = []
    for n in range(A_art_corr.shape[0]):
        error_value={}
        for m in range(A_art_corr.shape[1]):
            error_value[f"art_sys_corr_{m+1}"] = float(A_art_corr[n,m])
        
        error_value["luminosity_uncertainty"] = float(lum_errors[n])
        error_value["statistical_uncertainty"] = float(stat_errors[n])
        error.append(error_value)

    uncertainties_yaml = {"definition": error_definition, "bins": error}
    
    with open(f"uncertainties.yaml",'w') as file:
        yaml.dump(uncertainties_yaml,file, sort_keys=False)
    

    cov_lum = np.einsum("i,j->ij",lum_errors,lum_errors)
    cov_stat = np.diag(stat_errors**2)
    
    

    covmat = cov_corr + cov_stat + cov_lum 
    covmat1 = cov_corr1 + cov_stat + cov_lum 
    print(f"same covmat ?: {np.allclose(covmat,covmat1)}")
    return np.real(covmat1)



if __name__ == "__main__" :

    # write kinematics and central data values
    # filter_ATLAS_1JET_8TEV_data_kinetic()

    covmat = filter_ATLAS_1JET_8TEV_uncertainties()
    

    from validphys.loader import Loader
    from validphys.covmats import dataset_inputs_covmat_from_systematics
    from validphys.api import API
    from validphys.commondataparser import parse_commondata

    setname = "ATLAS_1JET_8TEV_R06"
    dsinps = [
         {'dataset': setname},
            ]
    inp = dict(dataset_inputs=dsinps, theoryid=200, use_cuts="internal")
    cov = API.dataset_inputs_covmat_from_systematics(**inp)

    l = Loader()
    cd = l.check_commondata(setname = setname).load_commondata_instance()
    dataset_input = l.check_dataset(setname,theoryid=200)
    
    dat_file = '/Users/markcostantini/codes/nnpdfgit/nnpdf/buildmaster/results/DATA_ATLAS_1JET_8TEV_R06.dat'
    # dat_file = '/Users/markcostantini/codes/nnpdfgit/nnpdf/nnpdfcpp/data/commondata/DATA_ATLAS_1JET_8TEV_R06.dat'

    sys_file = '/Users/markcostantini/codes/nnpdfgit/nnpdf/buildmaster/results/systypes/SYSTYPE_ATLAS_1JET_8TEV_R06_DEFAULT.dat'
    sys_file = '/Users/markcostantini/codes/nnpdfgit/nnpdf/nnpdfcpp/data/commondata/systypes/SYSTYPE_ATLAS_1JET_8TEV_R06_DEFAULT.dat'

    cd = parse_commondata(dat_file,sys_file,setname)

    cmat = dataset_inputs_covmat_from_systematics(
        dataset_inputs_loaded_cd_with_cuts=[cd],
        data_input = [dataset_input],
        use_weights_in_covmat=False,
        norm_threshold=None,
        _list_of_central_values=None,
        _only_additive=False,
             )

    # print(cov / cmat)
    # print(np.allclose(cov / cmat, np.ones(cmat.shape)))

    ones = cov / covmat
    print(ones)
    print(np.max(ones), np.min(ones))
    print(np.allclose(ones, np.ones(cov.shape), rtol = 1e-5))