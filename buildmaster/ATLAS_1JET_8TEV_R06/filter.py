import yaml
import numpy as np
import pandas as pd
from scipy.linalg import block_diag
from filter_utils import (
                            get_data_values, get_kinematics,
                            get_stat_errors, get_lumi_errors,
                            fill_df
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
    
    # statistical errors fully uncorrelated ?
    stat_errors = df_unc['stat'].to_numpy()
    cov_stat = np.diag(stat_errors**2)
    

    BD_stat = np.einsum('i,j->ij',dfs[0]['stat'],dfs[0]['stat'])
    for i in range(1,len(dfs)):
        cov = np.einsum('i,j->ij',dfs[i]['stat'].to_numpy(),dfs[i]['stat'].to_numpy())
        BD_stat = block_diag(BD_stat,cov)

    # luminosity errors
    lum_errors = df_unc["syst_lumi"].to_numpy()
    cov_lum = np.einsum("i,j->ij",lum_errors,lum_errors)

    col_to_drop =  ("syst_JES_EtaIntercalibration_Stat98", "syst_JES_MJB_Alpha",
                    "syst_JES_Pileup_MuOffset")

    # df_unc = df_unc.loc[:, ~df.columns.str.startswith(col_to_drop)]
    
    A_corr = df_unc.drop(["stat","syst_lumi"],axis =1).to_numpy() / np.sqrt(2.)

    cov_corr = np.einsum("ij,kj->ik",A_corr,A_corr)

    covmat = cov_corr + BD_stat + cov_lum
    return covmat



if __name__ == "__main__" :

    # write kinematics and central data values
    # filter_ATLAS_1JET_8TEV_data_kinetic()

    covmat = filter_ATLAS_1JET_8TEV_uncertainties()
    # print(covmat)

    from validphys.loader import Loader
    from validphys.covmats import dataset_inputs_covmat_from_systematics
    setname = "ATLAS_1JET_8TEV_R06"
    l = Loader()
    cd = l.check_commondata(setname = setname).load_commondata_instance()
    dataset_input = l.check_dataset(setname,theoryid=200)
    from validphys.commondataparser import parse_commondata
    dat_file = '/Users/markcostantini/codes/nnpdfgit/nnpdf/buildmaster/results/DATA_ATLAS_1JET_8TEV_R06.dat'
    sys_file = '/Users/markcostantini/codes/nnpdfgit/nnpdf/buildmaster/results/systypes/SYSTYPE_ATLAS_1JET_8TEV_R06_DEFAULT.dat'
    
    cd = parse_commondata(dat_file,sys_file,setname)

    cmat = dataset_inputs_covmat_from_systematics(
        dataset_inputs_loaded_cd_with_cuts=[cd],
        data_input = [dataset_input],
        use_weights_in_covmat=False,
        norm_threshold=None,
        _list_of_central_values=None,
        _only_additive=False,
             )

    print(covmat / cmat)
    print(np.allclose(covmat / cmat, np.ones(cmat.shape)))