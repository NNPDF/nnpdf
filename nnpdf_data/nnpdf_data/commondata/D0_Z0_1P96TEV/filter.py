import pathlib

import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

NB_POINTS = 28
MZ_VALUE = 91.1876  # GeV
MW_VALUE = 80.398  # GeV
SQRT_S = 1_960.0

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
    filename = f"HEPData-ins744624-v{version}-Table_{table_id}"
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
            "y": {"min": ymin, "mid": (ymin + ymax) / 2, "max": ymax},
            "m_Z2": {"min": None, "mid": MZ_VALUE**2, "max": None},
            "sqrts": {"min": None, "mid": SQRT_S, "max": None},
        }
        kinematics.append(kin_value)

    return kinematics


def get_errors(hepdata: dict) -> dict:
    """Extract the error values from the HepData yaml file.

    Parameters
    ----------
    hepdata: Dict
        Dictionary containing all data info

    Returns
    -------
    list:
        list of dictionaries whose contents are the various
        source of uncertainties

    """

    stat = []
    systematics = []
    central_values = []
    for data_i in hepdata["dependent_variables"][0]["values"]:

        stat_i = data_i["errors"][0]["symerror"]
        stat.append(stat_i)

        if "asymerror" in data_i["errors"][1]:
            delta_min = data_i["errors"][1]["asymerror"]["minus"]
            delta_plus = data_i["errors"][1]["asymerror"]["plus"]
            se_delta, se_sigma = se(delta_plus, delta_min)
        else:
            se_delta = 0
            se_sigma = data_i["errors"][1]["symerror"]

        cv_i = data_i["value"] + se_delta
        systematics.append(se_sigma)
        central_values.append(cv_i)

    return central_values, {"stat": stat, "sys_corr": systematics}


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
        error_value[f"sys_corr_1"] = uncs["sys_corr"][i]
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
        },
        "sys_corr_1": {
            "description": f"Systematic uncertainty 1",
            "treatment": "MULT",
            "type": "UNCORR",
        },
    }

    errors_formatted = format_uncertainties(errors)

    with open("data_ZRAP.yaml", "w") as file:
        yaml.dump({"data_central": data}, file, sort_keys=False)

    with open("kinematics_ZRAP.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open("uncertainties_ZRAP.yaml", "w") as file:
        yaml.dump(
            {"definitions": error_definition, "bins": errors_formatted}, file, sort_keys=False
        )


def main_filter() -> None:
    """Main driver of the filter that produces commmondata.

    There are four main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR

    2. Uncorrelated Systematic uncertainties: MULT, UNCORR


    """

    yaml_content = load_yaml(table_id=1, version=1)
    kinematics = get_kinematics(yaml_content)
    data_central, uncertainties = get_errors(yaml_content)

    # Generate all the necessary files
    dump_commondata(kinematics, data_central, uncertainties)

    return


if __name__ == "__main__":
    main_filter()
