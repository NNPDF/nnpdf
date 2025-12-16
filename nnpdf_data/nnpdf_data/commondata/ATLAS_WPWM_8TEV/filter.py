"""
This file contains the piece of code needed to implement the ATLAS Measurement
of the cross-section and charge asymmetry of W bosons at 8TeV.
Systematic uncertainties are implemented starting from the breakdown
available on HepData.
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
    version: int
        version number

    Returns
    -------
    dict:
        dictionary containing the table contents
    """
    filename = f"HEPData-ins1729240-v{version}-Table_{table_id}"
    table = pathlib.Path(f"./rawdata/{filename}.yaml")
    return yaml.safe_load(table.read_text())


def get_kinematics(hepdata: dict) -> list:
    """
    Returns the kinematics in the form of a list of dictionaries.
    
    Parameters
    ----------
    hepdata: dict
        HepData yaml content loaded as dictionary
    
    Returns
    -------
    list:
        List of kinematic dictionaries
    """
    kin = []
    
    for i, M in enumerate(hepdata["independent_variables"][0]['values']):
        kin_value = {
            'abs_eta': {
                'min': None,
                'mid': 0.5 * (M['low'] + M['high']),
                'max': None,
            },
            'm_W2': {'min': None, 'mid': 6463.838404, 'max': None},
            'sqrts': {'min': None, 'mid': 8000.0, 'max': None},
        }
        kin.append(kin_value)
    
    return kin


def get_data_values(hepdata: dict) -> list:
    """
    Returns the central data values in the form of a list.
    
    Parameters
    ----------
    hepdata: dict
        HepData yaml content loaded as dictionary
    
    Returns
    -------
    list:
        List of central data values
    """
    data_central = []
    values = hepdata['dependent_variables'][0]['values']
    
    for value in values:
        # Store data central and convert units (pb to fb)
        data_central.append(value['value'] * 1000)
    
    return data_central


def get_errors(hepdata: dict) -> tuple:
    """
    Extract the uncertainties from hepdata and compute the shift of the central value 
    in case of asymmetric uncertainties.

    Parameters
    ----------
    hepdata: dict
        HepData yaml file loaded as dictionary

    Returns
    -------
    tuple:
        (central_values, df_errors) where df_errors is a pandas DataFrame
    """
    central_values = []
    df_errors = pd.DataFrame()
    
    for i, bin in enumerate(hepdata["dependent_variables"][0]["values"]):
        error_sources = []
        shift_cv = 0
        error_names = []
        
        for source in bin["errors"]:
            error_names.append(source["label"])
            
            if "asymerror" in source:
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
        central_values.append(cv_i * 1000)  # Convert to fb
    
    return central_values, df_errors


def format_uncertainties(uncs: dict) -> list:
    """Format the uncertainties to be dumped into the yaml file.

    Parameters
    ----------
    uncs: dict
        Dictionary containing the various sources of uncertainties

    Returns
    -------
    list:
        List of dictionaries whose elements are the various errors
    """
    combined_errors = []
    n_bins = uncs["systematics"].index.str.startswith("bin").sum()
    
    for i in range(n_bins):
        errors = {}
        if "statistics" in uncs:
            errors["stat"] = float(uncs["statistics"].loc[f"bin {i}"].values.item())
        
        for j, unc in enumerate(uncs["systematics"].loc[f"bin {i}"].values):
            errors[f"sys_corr_{j + 1}"] = float(unc)
        
        combined_errors.append(errors)
    
    return combined_errors


def dump_commondata(kinematics: list, data: np.ndarray, errors: dict) -> None:
    """Function that generates and writes the commondata files.

    Parameters
    ----------
    kinematics: list
        List containing the kinematic values
    data: np.ndarray
        Array containing the central values
    errors: dict
        Dictionary containing the different errors
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
    
    with open("data.yaml", "w") as file:
        yaml.dump({"data_central": data.tolist()}, file, sort_keys=False)
    
    with open("kinematics.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)
    
    with open("uncertainties.yaml", "w") as file:
        yaml.dump(
            {"definitions": error_definition, "bins": errors_formatted}, 
            file, 
            sort_keys=False
        )


def main_filter() -> None:
    """
    Main function that reads the HepData yaml files and generates the commondata files.
    Combines both W+ and W- data into single files.
    """
    # Load both W+ and W- data files
    hepdata_file_wplus = "rawdata/HEPData-ins1729240-v1-Table_5a.yaml"
    hepdata_file_wminus = "rawdata/HEPData-ins1729240-v1-Table_5b.yaml"
    
    with open(hepdata_file_wplus, 'r') as file:
        yaml_content_wplus = yaml.safe_load(file)
    
    with open(hepdata_file_wminus, 'r') as file:
        yaml_content_wminus = yaml.safe_load(file)
    
    # Process W+ data
    kinematics_wplus = get_kinematics(yaml_content_wplus)
    central_values_wplus, uncertainties_wplus = get_errors(yaml_content_wplus)
    
    # Process W- data
    kinematics_wminus = get_kinematics(yaml_content_wminus)
    central_values_wminus, uncertainties_wminus = get_errors(yaml_content_wminus)
    
    # Combine W+ and W- data
    kinematics_all = kinematics_wplus + kinematics_wminus
    central_values_all = np.concatenate([central_values_wplus, central_values_wminus])
    
    # Combine uncertainties DataFrames
    n_wplus = len(central_values_wplus)
    n_wminus = len(central_values_wminus)
    
    # Re-index before concatenating
    uncertainties_wplus.index = [f"bin {i}" for i in range(n_wplus)]
    uncertainties_wminus.index = [f"bin {i}" for i in range(n_wplus, n_wplus + n_wminus)]
    
    uncertainties_all = pd.concat([uncertainties_wplus, uncertainties_wminus])
    
    # Determine types for each systematic
    # First column is stat, rest are systematics
    treatment_list = ["ADD"]  # stat is additive
    type_list = ["UNCORR"]  # stat is uncorrelated
    
    for col in uncertainties_all.columns[1:]:
        # Most systematics are multiplicative
        treatment_list.append("MULT")
        
        # Determine correlation type
        if col == "sys,Modelling":
            type_list.append("CORR")
        elif col.startswith("sys,MJ"):
            type_list.append("CORR")
        elif col.startswith("sys,Muon"):
            type_list.append("CORR")
        elif col.startswith("sys,MET"):
            type_list.append("CORR")
        elif col.startswith("sys,Soft"):
            type_list.append("CORR")
        elif col.startswith("sys,PDF"):
            type_list.append("CORR")
        elif "MC stat" in col:
            type_list.append("UNCORR")
        elif col == "sys,Pileup":
            type_list.append("CORR")
        elif col.startswith("sys,WW") or col.startswith("sys,WZ") or col.startswith("sys,ZZ"):
            type_list.append("CORR")
        elif col.startswith("sys,Wtaunu") or col.startswith("sys,Ztautau") or col.startswith("sys,Zmumu"):
            type_list.append("CORR")
        elif col.startswith("sys,ttbar") or col == "sys,photon" or col == "sys,single top":
            type_list.append("CORR")
        elif col == "sys,Z-vertex":
            type_list.append("CORR")
        elif col == "sys,F/B Compatibility":
            type_list.append("CORR")
        else:
            type_list.append("CORR")
    
    sys_types = {
        "treatment": treatment_list,
        "type": type_list,
    }
    
    sys_types_df = pd.DataFrame(sys_types, index=uncertainties_all.columns).T
    df_errors = pd.concat([sys_types_df, uncertainties_all])
    
    # Split into statistics and systematics
    errors = {
        "statistics": df_errors.iloc[:, [0]], 
        "systematics": df_errors.iloc[:, 1:]
    }
    
    # Dump the commondata files
    dump_commondata(kinematics_all, central_values_all, errors)
    
    print("Generated data.yaml, kinematics.yaml, and uncertainties.yaml")
    print(f"Total bins: {len(central_values_all)} (W+: {n_wplus}, W-: {n_wminus})")


if __name__ == "__main__":
    # Process combined W+ and W- data
    main_filter()