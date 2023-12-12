import numpy as np
import pandas as pd
import pathlib
import yaml

MW_VALUE = 80.398  # GeV
SQRT_S = 8_000.0  # GeV
NORM_FACTOR = 1_000.0  # from pb -> fb
OBSERVABLE = ['R', 'A']
MAP_STATE = {'R': 4, 'A': 5}
MAP_METADATA = {'R': 0, 'A':1}

def load_yaml(observable: str) -> dict:
    """Load the HEP data table in yaml format.

    Parameters
    ----------
    table_id: int
        table ID number
    observable: str
        type of observable (R/A)

    Returns
    -------
    dict:
        ditionary containing the table contents

    """
    filename = f"./rawdata/Table{MAP_STATE[observable]}.yaml"
    table = pathlib.Path(filename)

    return yaml.safe_load(table.read_text())


def read_metadata(observable: str) -> tuple[int, int, list]:
    """Read the version and list of tables from metadata.

    Returns
    -------
    tuple(int, list):
        data version and list of hepdata tables
    observable: str
        type of observable (R/A)

    """
    metadata = pathlib.Path("./metadata.yaml")
    content = yaml.safe_load(metadata.read_text())
    idx = MAP_METADATA[observable]

    version = content["hepdata"]["version"]
    nb_datapoints = sum(content["implemented_observables"][idx]["npoints"])
    tables = content["implemented_observables"][idx]["tables"]

    return version, nb_datapoints, tables

def get_kinematics(hepdata: dict, bin_index: list) -> list:
    """Read the version and list of tables from metadata.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
    bin_index: list
        list of Non-empty bin index

    Returns
    -------
    tuple(int, list):
        data version and list of hepdata tables

    """
    etabins = hepdata["independent_variables"][1]["values"]

    kinematics = []
    for bins in bin_index:
        etamin = float(etabins[bins]["low"])
        etamax = float(etabins[bins]["high"])
        kin_value = {
            "eta": {"min": etamin, "mid": 0.5 * (etamin + etamax), "max": etamax},
            "M2": {"min": None, "mid": MW_VALUE ** 2, "max": None},
            "sqrt_s": {"min": None, "mid": SQRT_S, "max": None},
        }
        kinematics.append(kin_value)

    return kinematics

def get_data_values(hepdata: dict, bin_index: list) -> list:
    """Extract the central values from the HepData yaml file.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
    bin_index: list
        list of Non-empty bin index

    Returns
    -------
    list:
        list of dictionaries whose contents are the central values

    """
    central = hepdata["dependent_variables"][0]["values"]

    return [NORM_FACTOR * central[i]["value"] for i in bin_index]


def get_errors(hepdata: dict, bin_index: list) -> dict:
    """Extract the error values from the HepData yaml file.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
    bin_index: list
        list of Non-empty bin index

    Returns
    -------
    list:
        list of dictionaries whose contents are the various
        source of uncertainties

    """
    errors = hepdata["dependent_variables"][0]["values"]

    stat, sys_uncorr, sys_beam = [], [], []
    for idx in bin_index:
        stat.append(NORM_FACTOR * errors[idx]["errors"][0]["symerror"])
        sys_uncorr.append(NORM_FACTOR * errors[idx]["errors"][1]["symerror"])
        sys_beam.append(NORM_FACTOR * errors[idx]["errors"][2]["symerror"])

    return {
        "stat": stat,
        "sys_uncorr": sys_uncorr,
        "sys_beam": sys_beam,
    }
    
def concatenate_dicts(multidict: list[dict]) -> dict:
    """Given a list of dictionaries combined the values.

    Parameters
    ----------
    multidict: list[dict]
        list of dictionaries with the same keys

    Returns
    -------
    dict:
        dictionary whose keys are combined
    """
    new_dict = {}
    for key in multidict[0].keys():
        new_dict[key] = []
        for element in multidict:
            new_dict[key] += element[key]

    return new_dict


def format_uncertainties(uncs: dict, bin_index: np.ndarray) -> list:
    """Format the uncertainties to be dumped into the yaml file.

    Parameters
    ----------
    uncs: dict
        dictionary containing the various source of uncertainties

    Returns
    -------
    list:
        list of ditionaries whose elements are the various errors

    """
    combined_errors = []

    for idat in range(len(bin_index)):
        error_value = {}
        # for jdat in range(len(artunc[idat])):
        #     error_value[f"sys_corr_{jdat + 1}"] = artunc[idat][jdat]

        error_value["stat"] = uncs["stat"][idat]
        error_value["sys_uncorr"] = uncs["sys_uncorr"][idat]
        error_value["sys_beam"] = uncs["sys_beam"][idat]
        combined_errors.append(error_value)

    return combined_errors


def dump_commondata(kinematics: list, data: list, errors: list, observable: str) -> None:
    """Function that generates and writes the commondata files.

    Parameters
    ----------
    kinematics: list
        list containing the kinematic values
    data: list
        list containing the central values
    errors: list
        list containing the different errors
    observable: str
        type of final state (R/A)

    """
    error_definition["sys_uncorr"] = {
            "description": "Uncorrelated systematic uncertainties",
            "treatment": "MULT",
            "type": "UNCORR",
    }

    error_definition["stat"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    error_definition["sys_beam"] = {
        "description": "Systematic Beam uncertainties",
        "treatment": "MULT",
        "type": "LHCBBEAM8TEV",
    }

    with open(f"data_{observable}.yaml", "w") as file:
        yaml.dump({"data_central": data}, file, sort_keys=False)

    with open(f"kinematics_{observable}.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open(f"uncertainties_{observable}.yaml", "w") as file:
        yaml.dump(
            {"definitions": error_definition, "bins": errors},
            file,
            sort_keys=False,
        )


def main_filter():
    """Main driver of the filter that produces commmondata.

    There are three main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR

    2. Uncorrelated Systematic uncertainties: MULT, UNCORR

    3. Beam Systematic uncertainties: MULT, LHCBBEAM8TEV

    """
    for obs in OBSERVABLE:
        _, nbpoints, tables = read_metadata(observable=obs)
        bin_index = [i for i in range(nbpoints)]  # Non-empty Bins

        comb_kins, comb_data = [], []
        combined_errors = []
        for tabid in tables:
            yaml_content = load_yaml(observable=obs)
            
            kinematics = get_kinematics(yaml_content, bin_index)

            # Extract the data, and uncertainties
            data_central = get_data_values(yaml_content, bin_index)
            uncertainties = get_errors(yaml_content, bin_index)

            # Collect all the results from different tables
            comb_kins += kinematics
            comb_data += data_central
            combined_errors.append(uncertainties)

        errors_combined = concatenate_dicts(combined_errors)
        
        errors = format_uncertainties(errors_combined, bin_index)

        # Generate all the necessary files
        dump_commondata(comb_kins, comb_data, errors, obs)

    return


if __name__ == "__main__":
    main_filter()
