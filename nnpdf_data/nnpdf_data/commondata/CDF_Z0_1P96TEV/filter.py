import pathlib

import numpy as np
import pandas as pd
import yaml

NB_POINTS = 28
MZ_VALUE = 91.1876  # GeV
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
    filename = f"HEPData-ins856131-v{version}-Table_{table_id}"
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
            "k1": {"min": ymin, "mid": 0.5 * (ymin + ymax), "max": ymax},
            "k2": {"min": None, "mid": MZ_VALUE ** 2, "max": None},
            "k3": {"min": None, "mid": SQRT_S, "max": None},
        }
        kinematics.append(kin_value)

    return kinematics

def get_data_values(hepdata: dict, indx: int = 0) -> list:
    """Extract the central values from the HepData yaml file.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
    bin_index: list
        list of Non-empty bin index
    idx: int
        index from which to read the central value, default=0

    Returns
    -------
    list:
        list of dictionaries whose contents are the central values

    """
    central = hepdata["dependent_variables"][indx]["values"]
    return [central[i]["value"] for i in range(NB_POINTS)]

def get_errors(hepdata: dict, indx: int = 0) -> dict:
    """Extract the error values from the HepData yaml file.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
    indx: int
        index from which the errors will be read

    Returns
    -------
    list:
        list of dictionaries whose contents are the various
        source of uncertainties

    """
    errors = hepdata["dependent_variables"][indx]["values"]

    stat, sys_corr = [], []
    for idx in range(NB_POINTS):
        stat.append(errors[idx]["errors"][0]["symerror"])
        sys_corr.append(errors[idx]["errors"][1]["symerror"])

    return {"stat": stat, "sys_corr": sys_corr}

#
# kinematics = get_kinematics(yaml_content, bin_index, MAP_TABLE[tabid])
# data_central = get_data_values(yaml_content, bin_index, indx=idx)
# uncertainties = get_errors(yaml_content, bin_index, indx=idx)

def read_metadata() -> tuple[int, int, list]:
    """Read the version and list of tables from metadata.

    Returns
    -------
    tuple(int, list):
        data version and list of hepdata tables

    """
    metadata = pathlib.Path("./metadata.yaml")
    content = yaml.safe_load(metadata.read_text())

    version = content["hepdata"]["version"]
    nb_datapoints = sum(content["implemented_observables"][0]["npoints"])
    tables = content["implemented_observables"][0]["tables"]

    return version, nb_datapoints, tables

def main_filter() -> None:
    """Main driver of the filter that produces commmondata.

    There are four main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR

    2. Correlated Systematic uncertainties: MULT, CORR:
        constructed by symmetrizing the correlation matrix and extracting
        the artificial systematic uncertainties from the latter

    3. Beam Systematic uncertainties: MULT, LHCBBEAM7TEV

    4. Luminosity Systematic uncertainties: MULT, LHCBLUMI7TEV

    """

    yaml_content = load_yaml(table_id=2, version=1)

    kinematics = get_kinematics(yaml_content)
    data_central = get_data_values(yaml_content)
    uncertainties = get_errors(yaml_content)

    # correlations from https://inspirehep.net/literature/806697
    # compile using
    # g++ -c error_propagator_g++_032610.C
    # g++ error_propagator_g++_032610.o -o systematics

    import pdb; pdb.set_trace()


    nbp_idx = 0
    comb_kins, comb_data = [], []
    combined_errors = []
    for tabid in TABLES.keys():  # Loop over tables
        for idx in TABLES[tabid]:  # Loop over Bosons [Z, W+, W-]
            bin_index = [i for i in range(NB_POINTS[nbp_idx])]
            yaml_content = load_yaml(table_id=tabid, version=version)

            # Extract the kinematic, data, and uncertainties
            kinematics = get_kinematics(yaml_content, bin_index, MAP_TABLE[tabid])
            data_central = get_data_values(yaml_content, bin_index, indx=idx)
            uncertainties = get_errors(yaml_content, bin_index, indx=idx)

            # Collect all the results from different tables
            comb_kins += kinematics
            comb_data += data_central
            combined_errors.append(uncertainties)

            nbp_idx += 1

    errors_combined = concatenate_dicts(combined_errors)
    # Compute the Artifical Systematics from CovMat
    corrmat = read_corrmatrix(nb_datapoints=nbpoints)
    covmat = multiply_syst(corrmat, errors_combined["sys_corr"])
    artunc = generate_artificial_unc(ndata=nbpoints, covmat_list=covmat.tolist(), no_of_norm_mat=0)
    errors = format_uncertainties(errors_combined, artunc)

    # Generate all the necessary files
    dump_commondata(comb_kins, comb_data, errors)

    return


# with open("data.yaml", "w") as file:
#     yaml.dump({"data_central": data}, file, sort_keys=False)
#
# with open("kinematics.yaml", "w") as file:
#     yaml.dump({"bins": kinematics}, file, sort_keys=False)
#
# with open("uncertainties.yaml", "w") as file:
#     yaml.dump({"definitions": error_definition, "bins": errors}, file, sort_keys=False)

if __name__ == "__main__":
    main_filter()
