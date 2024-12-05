import pathlib

import numpy as np
import pandas
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

NB_POINTS = 10
MZ_VALUE = 91.1876  # GeV
MW_VALUE = 80.398  # GeV
SQRT_S = 1_960.0


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
    filename = f"HEPData-ins1253555-v{version}-Table_{table_id}"
    table = pathlib.Path(f"./rawdata/{filename}.yaml")

    return yaml.safe_load(table.read_text())


def get_kinematics(hepdata: dict) -> list:
    """Read the version and list of tables from metadata.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info

    Returns
    -------
    tuple(int, list):
        data version and list of hepdata tables

    """
    rapbins = hepdata["independent_variables"][0]["values"]

    kinematics = []
    for bins in range(NB_POINTS):
        ymin = float(rapbins[bins]["low"])
        ymax = float(rapbins[bins]["high"])
        kin_value = {
            "abs_eta": {"min": ymin, "mid": (ymin + ymax) / 2, "max": ymax},
            "m_W2": {"min": None, "mid": MW_VALUE**2, "max": None},
            "sqrts": {"min": None, "mid": SQRT_S, "max": None},
        }
        kinematics.append(kin_value)

    return kinematics


def get_data_values(hepdata: dict, indx: int = 0) -> list:
    """Extract the central values from the HepData yaml file.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
    idx: int
        index from which to read the central value, default=0

    Returns
    -------
    list:
        list of dictionaries whose contents are the central values

    """
    central = hepdata["dependent_variables"][indx]["values"]
    return [central[i]["value"] * 1e-2 for i in range(NB_POINTS)]


def get_errors(hepdata: list) -> dict:
    """Extract the error values from the HepData yaml file.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info

    Returns
    -------
    list:
        list of dictionaries whose contents are the various
        source of uncertainties

    """
    # extract the systematics
    systematics = []
    for row in hepdata[1]["dependent_variables"]:
        systematics.append([value["value"] * 1e-2 for value in row["values"]])
    systematics = np.array(systematics).T

    # append the statistical uncertainties per datapoint
    stat = []
    for row in hepdata[0]["dependent_variables"][0]["values"]:
        stat.append(row["errors"][0]["symerror"] * 1e-2)

    return {"stat": stat, "sys_corr": systematics}


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
    for i in range(NB_POINTS):
        error_value = {}
        error_value["stat"] = uncs["stat"][i]
        for j, sys in enumerate(uncs["sys_corr"][i]):
            error_value[f"sys_corr_{j+1}"] = float(sys)
        combined_errors.append(error_value)

    return combined_errors


def dump_commondata(kinematics: list, data: list, errors: dict) -> None:
    """Function that generates and writes the commondata files.

    Parameters
    ----------
    kinematics: list
        list containing the kinematic values
    data: list
        list containing the central values
    errors: dict
        Dictionary containing the different errors

    """

    error_definition = {
        "stat": {
            "description": "Uncorrelated statistical uncertainties",
            "treatment": "ADD",
            "type": "UNCORR",
        }
    }

    n_sys = errors["sys_corr"].shape[1]

    for i in range(n_sys):
        error_definition[f"sys_corr_{i + 1}"] = {
            "description": f"Systematic uncertainty {i + 1}",
            "treatment": "ADD",
            "type": "CORR",
        }

    # update lumi entry
    error_definition[f'sys_corr_{n_sys}']['type'] = "UNCORR"

    error_definition["stat"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    errors_formatted = format_uncertainties(errors)

    with open("data_ASY.yaml", "w") as file:
        yaml.dump({"data_central": data}, file, sort_keys=False)

    with open("kinematics_ASY.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open("uncertainties_ASY.yaml", "w") as file:
        yaml.dump(
            {"definitions": error_definition, "bins": errors_formatted}, file, sort_keys=False
        )


def main_filter() -> None:
    """Main driver of the filter that produces commmondata.

    There are four main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR

    2. Correlated Systematic uncertainties: ADD, CORR

    3. Uncorrelated Systematic uncertainties: ADD, UNCORR

    """

    yaml_content_data = load_yaml(table_id=1, version=1)
    yaml_content_uncertainties = load_yaml(table_id=3, version=1)
    kinematics = get_kinematics(yaml_content_data)
    data_central = get_data_values(yaml_content_data)
    uncertainties = get_errors([yaml_content_data, yaml_content_uncertainties])

    # Generate all the necessary files
    dump_commondata(kinematics, data_central, uncertainties)

    return


if __name__ == "__main__":
    main_filter()
