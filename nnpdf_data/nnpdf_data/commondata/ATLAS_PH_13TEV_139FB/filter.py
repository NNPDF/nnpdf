"""
This file contains the piece of code needed to implement the ATLAS direct
photon measurement at 13 TeV. Distributions differential in the transverse
energy of the photon in various bins of photon rapidity are consisdered.
Systematic uncertainties are implemented starting from the breakdown
available on HepData. The correlation treatment follows the approach mentioned
in the paper (see Sect. 7.7): "There are bin-to-bin correlations of the 
systematic uncertainties for each source. [...] The following uncertainties are 
considered as uncorrelated bin-to-bin: photon-identification efficiency, choice 
of background control regions, E_t^iso modelling and MC statistical 
uncertainties.
"""

import pathlib

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

from nnpdf_data.filter_utils.utils import symmetrize_errors as se


def load_yaml(table_id: int, version: int = 1) -> dict:
    """Load the HEP data table in yaml format.

    Parameters
    ----------
    table_id: int
        table ID number

    Returns
    -------
    dict:
        ditionary containing the table contents

    """
    filename = f"HEPData-ins2628741-v{version}-Table_{table_id}"
    table = pathlib.Path(f"./rawdata/{filename}.yaml")

    return yaml.safe_load(table.read_text())


def get_kinematics(hepdata: dict, bin_index: list = [], indx: int = 0, min_rap=None, max_rap=None) -> list:
    """Read the version and list of tables from metadata.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
    bin_index: list
        list of Non-empty bin index
    indx: int
        Column index from which to read, default=0

    Returns
    -------
    kinematics: list
        kinematic info

    """
    bins = hepdata["independent_variables"][indx]["values"]

    kinematics = []
    for i in bin_index:
        min_et, max_et = bins[i]["low"], bins[i]["high"]

        kin_value = {
            "eta": {"min": min_rap, "mid": ((min_rap + max_rap) / 2) , "max": max_rap},
            "ET": {"min": min_et, "mid": ((min_et + max_et) / 2), "max": max_et},
            "sqrts": {"min": None, "mid": 13000, "max": None},
        }
        kinematics.append(kin_value)

    return kinematics


def get_data_values(hepdata: dict, bin_index: list, indx: int = 0) -> list:
    """Extract the central values from the HepData yaml file.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
    bin_index: list
        Bin indices that must be parsed
    indx: int
        Column index from which to read the central value, default=0

    Returns
    -------
    list:
        list of dictionaries whose contents are the central values

    """
    central = hepdata["dependent_variables"][indx]["values"]
    return np.array([central[i]["value"] for i in bin_index])


def get_errors(hepdata: dict, bin_index: list) -> dict:
    """
    Extract the uncertainties from hepdata and computes the shift of the central value in case of
    asymmetric uncertainties

    Parameters
    ----------
    hepdata: dict
        Hepdata yaml file loaded as dictionary
    bin_index: list
        Bin indices that must be parsed

    Returns
    -------
    dict:
        Dictionary containing the errors (as pandas DataFrame) and shifts of central values
    """
    # parse the systematics
    central_values = []  # relevant for asymmetric uncertainties
    df_errors = pd.DataFrame()
    for i, bin in enumerate(hepdata["dependent_variables"][0]["values"]):

        error_sources = []
        shift_cv = 0
        error_names = []
        for source in bin["errors"]:
            error_names.append(source["label"])
            if source["label"] == "stat":
                error_sources.append(source["symerror"])
            elif "asymerror" in source:
                delta_min = source["asymerror"]["minus"]
                delta_plus = source["asymerror"]["plus"]
                se_delta, se_sigma = se(delta_plus, delta_min)
                error_sources.append(se_sigma)
                shift_cv += se_delta
            elif "symerror" in source:
                se_sigma = source["symerror"]
                error_sources.append(se_sigma)
        df_bin = pd.DataFrame([error_sources], columns=error_names, index=[f"bin {i}"])
        df_errors = pd.concat([df_errors, df_bin])
        cv_i = bin["value"] + shift_cv
        central_values.append(cv_i)

    return central_values, df_errors


def format_uncertainties(uncs: dict) -> list:
    """Format the uncertainties to be dumped into the yaml file.

    Parameters
    ----------
    uncs: dict
        Dictionary containing the various source of uncertainties

    Returns
    -------
    list:
        list of dictionaries whose elements are the various errors

    """

    combined_errors = []
    n_bins = uncs["systematics"].index.str.startswith("bin").sum()
    for i in range(n_bins):
        errors = {}
        if "statistics" in uncs:
            errors["stat"] = uncs["statistics"].loc[f"bin {i}"].values.item()
        for j, unc in enumerate(uncs["systematics"].loc[f"bin {i}"].values):
            errors[f"sys_corr_{j + 1}"] = float(unc)

        combined_errors.append(errors)

    return combined_errors


def dump_commondata(kinematics: list, data: list, errors: dict, obs: str) -> None:
    """Function that generates and writes the commondata files.

    Parameters
    ----------
    kinematics: list
        list containing the kinematic values
    data: list
        list containing the central values
    errors: dict
        Dictionary containing the different errors
    obs: str
        Name to append to the file names
    """

    if "statistics" in errors:
        error_definition = {
            "stat": {
                "description": "Uncorrelated statistical uncertainties",
                "treatment": errors["statistics"].loc["treatment"].iloc[0],
                "type": errors["statistics"].loc["type"].iloc[0],
            }
        }
    else:
        error_definition = {}

    n_sys = errors["systematics"].shape[1]
    for i in range(n_sys):

        error_definition[f"sys_corr_{i + 1}"] = {
            "description": errors["systematics"].columns[i],
            "treatment": errors["systematics"].loc["treatment"].iloc[i],
            "type": errors["systematics"].loc["type"].iloc[i],
        }

    errors_formatted = format_uncertainties(errors)
    with open(f"data_{obs}.yaml", "w") as file:
        yaml.dump({"data_central": data.tolist()}, file, sort_keys=False)

    with open(f"kinematics_{obs}.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open(f"uncertainties_{obs}.yaml", "w") as file:
        yaml.dump(
            {"definitions": error_definition, "bins": errors_formatted}, file, sort_keys=False
        )


def main_filter(obs=None) -> None:
    """
    Main function that reads the HepData yaml files and generates the commondata files
    """
    unc_err_lbl = ["sysPhotonID",
                   "sysBackgroundIsolation",
                   "sysBackgroundIsolationUpperLimit",
                   "sysBackgroundID",
                   "sysBackground",
                   "sysIsolationMC",
                   "sysMCstats",]
    uncertainties_all = pd.DataFrame()
    central_values_all = np.array([])
    kinematics_all = []
    n_datapoints = [12, 11, 11, 10, 9, 9]
    min_rapidities = [0.0, 0.6, 0.8, 1.56, 1.81, 2.01]
    max_rapidities = [0.6, 0.8, 1.37, 1.81, 2.01, 2.37]

    if obs=="ET-ETA-R04":
        yaml_content_data = [load_yaml(table_id=i, version=1) for i in range(1,7)]
    elif obs == "ET-ETA-R02":
        yaml_content_data = [load_yaml(table_id=i, version=1) for i in range(7,13)]
    else:
        print("Wrong observable.")
        print("Available observables are:")
        print("- ET-ETA-R04")
        print("- ET-ETA-R02")
        exit()
    
    for i, yaml_content in enumerate(yaml_content_data):
        kinematics = get_kinematics(
            yaml_content, bin_index=range(n_datapoints[i]), min_rap=min_rapidities[i], max_rap=max_rapidities[i]
        )
        central_values, uncertainties = get_errors(yaml_content, bin_index=range(n_datapoints[i]))
        uncertainties_all = pd.concat([uncertainties_all, uncertainties])
        central_values_all = np.concatenate([central_values_all, central_values])
        kinematics_all += kinematics
        
    uncertainties_all.index = [f"bin {i}" for i in range(uncertainties_all.shape[0])]
    
    n_sources = uncertainties_all.shape[1]
    sys_types = {
        "treatment": ["ADD"] + ["MULT"] * (n_sources - 1),
        "type": ["UNCORR"] + ["CORR"] * (n_sources - 2) + ["ATLASLUMIRUNII"],
    }
    sys_types_df = pd.DataFrame(sys_types, index=uncertainties_all.columns).T
    df_errors = pd.concat([sys_types_df, uncertainties_all])
    
    errors = {"statistics": df_errors.iloc[:, [0]], "systematics": df_errors.iloc[:, 1:]}

    for lbl in unc_err_lbl:
        #errors["systematics"][lbl] = errors["systematics"][lbl].replace('CORR','UNCORR')
        #errors["systematics"][lbl] = errors["systematics"][lbl].replace('MULT','ADD')
        errors["systematics"].loc["treatment",lbl] = "ADD"
        errors["systematics"].loc["type", lbl] = "UNCORR"
        
    dump_commondata(kinematics_all, central_values_all, errors, obs=obs)
    
    return


if __name__ == "__main__":
    main_filter("ET-ETA-R02")
    main_filter("ET-ETA-R04")
