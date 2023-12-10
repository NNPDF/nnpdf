import yaml
import numpy as np
import pandas as pd
from filter_utils import get_data_values, get_kinematics, fill_df

# ignore pandas warning
import warnings

warnings.filterwarnings("ignore")


def filter_ATLAS_1JET_8TEV_data_kinetic():
    """
    write kinematics and central data values
    in kinematics.yaml and data.yaml files
    respectively.
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    tables = metadata["hepdata"]["tables"]

    # get kinematics from hepdata tables
    kin = get_kinematics(tables, version)

    # get central values from hepdata tables
    data_central = get_data_values(tables, version)

    data_central_yaml = {"data_central": data_central}
    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_ATLAS_1JET_8TEV_uncertainties(variant='nominal'):
    """
    Writes the uncertainties to a .yaml file.
    Two possible variants are implemented: nominal and decorrelated
    
    There are three types of uncertainties:

    1. Statistical Uncertainties: ADD, UNCORR
       

    2. Systematic Uncertainties: ADD, CORR
       Constructed following the exp. prescription:

       - Construct an Error matrix in which
         each part of an asymmetric error is considered
         as as separate error (hence dividing by sqrt(2))
         see also filter_utils/process_error and
         filter_utils/HEP_table_to_df


    3. Luminosity Uncertainty: ATLASLUMI12
       this uncertainty is correlated with all
       the other ATLASLUMI12 datasets
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    tables = metadata["hepdata"]["tables"]

    # get df of uncertainties
    dfs = []
    for table in tables:
        # uncertainties dataframe
        df = fill_df(table, version, variant)
        dfs.append(df)

    df_unc = pd.concat([df for df in dfs], axis=0)

    # statistical errors fully uncorrelated
    stat_errors = df_unc["stat"].to_numpy()

    # luminosity errors
    lum_errors = df_unc["syst_lumi"].to_numpy()

    A_corr = df_unc.drop(["stat", "syst_lumi"], axis=1).to_numpy() / np.sqrt(2.0)
    cov_corr = np.einsum("ij,kj->ik", A_corr, A_corr)

    # error definition
    error_definition = {
        f"{col}": {
            "description": f"correlated systematic {col}",
            "treatment": "MULT",
            "type": "CORR",
        }
        for col in df_unc.drop(["stat", "syst_lumi"], axis=1).columns
    }

    error_definition["luminosity_uncertainty"] = {
        "description": "luminosity uncertainty",
        "treatment": "MULT",
        "type": "ATLASLUMI12",
    }

    error_definition["statistical_uncertainty"] = {
        "description": "statistical uncertainty",
        "treatment": "MULT",
        "type": "UNCORR",
    }

    # store error in dict
    error = []
    for n in range(A_corr.shape[0]):
        error_value = {}
        for col, m in zip(df_unc.drop(["stat", "syst_lumi"], axis=1).columns, range(A_corr.shape[1])):
            error_value[f"{col}"] = float(A_corr[n, m])

        error_value["luminosity_uncertainty"] = float(lum_errors[n])
        error_value["statistical_uncertainty"] = float(stat_errors[n])
        error.append(error_value)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    # write uncertainties to file
    if variant=='nominal':
        with open(f"uncertainties.yaml", "w") as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)
    else:
        with open(f"uncertainties_{variant}.yaml", "w") as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)

    # @@@@@@@@@@@ code below for testing only, should be removed at some point @@@@@@@@@@@@#
    cov_lum = np.einsum("i,j->ij", lum_errors, lum_errors)
    cov_stat = np.diag(stat_errors**2)

    covmat = cov_corr + cov_stat + cov_lum

    return np.real(covmat)


if __name__ == "__main__":
    # write kinematics and central data values
    filter_ATLAS_1JET_8TEV_data_kinetic()

    # write uncertainties file
    filter_ATLAS_1JET_8TEV_uncertainties(variant='nominal')

    # write decorrelated uncertainties file
    filter_ATLAS_1JET_8TEV_uncertainties(variant='decorrelated')

    ## 

    # # code below for testing only. Should be removed at some point
    # covmat = filter_ATLAS_1JET_8TEV_uncertainties()

    # from validphys.api import API

    # setname = "ATLAS_1JET_8TEV_R06"
    # dsinps = [
    #     {"dataset": setname},
    # ]
    # inp = dict(dataset_inputs=dsinps, theoryid=200, use_cuts="internal")
    # cov = API.dataset_inputs_covmat_from_systematics(**inp)

    # ones = cov / covmat
    # print(ones)
    # print(np.max(ones), np.min(ones))
    # print(np.allclose(ones, np.ones(cov.shape), rtol=1e-5))
