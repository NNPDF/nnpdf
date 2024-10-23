import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

from nnpdf_data.filter_utils.legacy_jets_utils import (
    fill_df_ATLAS_1JET_8TEV_R06,
    get_data_values_ATLAS_1JET_8TEV_R06,
    get_kinematics_ATLAS_1JET_8TEV_R06,
)


def filter_ATLAS_1JET_8TEV_data_kinetic():
    """
    write kinematics and central data values
    in kinematics.yaml and data.yaml files
    respectively.
    """

    with open("metadata.yaml") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    tables = metadata["hepdata"]["tables"]

    # get kinematics from hepdata tables
    kin = get_kinematics_ATLAS_1JET_8TEV_R06(tables, version)

    # get central values from hepdata tables
    data_central = get_data_values_ATLAS_1JET_8TEV_R06(tables, version)

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

    with open("metadata.yaml") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    tables = metadata["hepdata"]["tables"]

    # get df of uncertainties
    dfs = []
    for table in tables:
        # uncertainties dataframe
        df = fill_df_ATLAS_1JET_8TEV_R06(table, version, variant)
        dfs.append(df)

    df_unc = pd.concat([df for df in dfs], axis=0)

    # statistical errors fully uncorrelated
    stat_errors = df_unc["stat"].to_numpy()

    # luminosity errors
    lum_errors = df_unc["syst_lumi"].to_numpy()

    A_corr = df_unc.drop(["stat", "syst_lumi"], axis=1).to_numpy() / np.sqrt(2.0)

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
        for col, m in zip(
            df_unc.drop(["stat", "syst_lumi"], axis=1).columns, range(A_corr.shape[1])
        ):
            error_value[f"{col}"] = float(A_corr[n, m])

        error_value["luminosity_uncertainty"] = float(lum_errors[n])
        error_value["statistical_uncertainty"] = float(stat_errors[n])
        error.append(error_value)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    # write uncertainties to file
    if variant == 'nominal':
        with open(f"uncertainties.yaml", "w") as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)
    else:
        with open(f"uncertainties_{variant}.yaml", "w") as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    # write kinematics and central data values
    filter_ATLAS_1JET_8TEV_data_kinetic()

    # write uncertainties file
    filter_ATLAS_1JET_8TEV_uncertainties(variant='nominal')

    # write decorrelated uncertainties file
    filter_ATLAS_1JET_8TEV_uncertainties(variant='decorrelated')
