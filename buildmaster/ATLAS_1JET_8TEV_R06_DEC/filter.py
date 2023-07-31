import yaml
import numpy as np
import pandas as pd

import sys

sys.path.insert(0, ".")
from filter_utils import get_data_values, get_kinematics, fill_df, decompose_covmat

# ignore pandas warning
import warnings

warnings.filterwarnings("ignore")


def filter_ATLAS_1JET_8TEV_R06_DEC_data_kinetic():
    """
    write kinematics and central data values
    in kinematics.yaml and data.yaml files
    respectively.
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    tables = metadata["hepdata"]["tables"]

    # get kinematics from hepdata tables
    kin = get_kinematics(tables)

    # get central values from hepdata tables
    data_central = get_data_values(tables)

    data_central_yaml = {"data_central": data_central}
    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_ATLAS_1JET_8TEV_R06_DEC_uncertainties():
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

       - Decompose covariance matrix so as to get ndat art unc

    3. Luminosity Uncertainty: ATLASLUMI12
       this uncertainty is correlated with all
       the other ATLASLUMI12 datasets
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    tables = metadata["hepdata"]["tables"]

    # get df of uncertainties
    dfs = []
    for table in tables:
        # uncertainties dataframe
        df = fill_df(table)
        dfs.append(df)

    df_unc = pd.concat([df for df in dfs], axis=0)

    # statistical errors fully uncorrelated
    stat_errors = df_unc["stat"].to_numpy()

    # luminosity errors
    lum_errors = df_unc["syst_lumi"].to_numpy()

    A_corr = df_unc.drop(["stat", "syst_lumi"], axis=1).to_numpy() / np.sqrt(2.0)
    cov_corr = np.einsum("ij,kj->ik", A_corr, A_corr)
    A_art_corr = decompose_covmat(covmat=cov_corr)

    # error definition
    error_definition = {
        f"art_sys_corr_{i}": {
            "description": f"artificial systematic {i}",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(1, A_art_corr.shape[0] + 1)
    }

    error_definition["luminosity_uncertainty"] = {
        "description": "luminosity uncertainty",
        "treatment": "ADD",
        "type": "ATLASLUMI12",
    }

    error_definition["statistical_uncertainty"] = {
        "description": "statistical uncertainty",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    # store error in dict
    error = []
    for n in range(A_art_corr.shape[0]):
        error_value = {}
        for m in range(A_art_corr.shape[1]):
            error_value[f"art_sys_corr_{m+1}"] = float(A_art_corr[n, m])

        error_value["luminosity_uncertainty"] = float(lum_errors[n])
        error_value["statistical_uncertainty"] = float(stat_errors[n])
        error.append(error_value)

    uncertainties_yaml = {"definition": error_definition, "bins": error}

    # write uncertainties to file
    with open(f"uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)

    # @@@@@@@@@@@@@@@ code below is only needed for testing and can be removed at some point @@@@@@@@@@@@ #
    cov_lum = np.einsum("i,j->ij", lum_errors, lum_errors)
    cov_stat = np.diag(stat_errors**2)

    covmat = cov_corr + cov_stat + cov_lum
    # @@@@@@@@@@@@@@@ code above is only needed for testing and can be removed at some point @@@@@@@@@@@@ #
    return np.real(covmat)


if __name__ == "__main__":
    # write kinematics and central data values
    # filter_ATLAS_1JET_8TEV_R06_DEC_data_kinetic()

    # # write uncertainties 
    # filter_ATLAS_1JET_8TEV_R06_DEC_uncertainties()

    # # @@@@@@@@@@@@@@@ CODE BELOW ONLY NEEDED TO TEST AND SHOULD BE REMOVED
    covmat = filter_ATLAS_1JET_8TEV_R06_DEC_uncertainties()

    from validphys.api import API
    
    setname = "ATLAS_1JET_8TEV_R06_DEC"
    dsinps = [
         {'dataset': setname},
            ]
    inp = dict(dataset_inputs=dsinps, theoryid=200, use_cuts="internal")
    cov = API.dataset_inputs_covmat_from_systematics(**inp)


    ones = cov / covmat
    print(ones)
    print(np.max(ones), np.min(ones))
    print(np.allclose(ones, np.ones(cov.shape), rtol=1e-5))
