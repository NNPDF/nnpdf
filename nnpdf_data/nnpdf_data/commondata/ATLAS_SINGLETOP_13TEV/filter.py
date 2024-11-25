import pathlib

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

NB_POINTS = 4
MT_VALUE = 172.5
SQRT_S = 13_000.0

from nnpdf_data.filter_utils.utils import cormat_to_covmat, covmat_to_artunc
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
    filename = f"HEPData-ins1486394-v{version}-Table_{table_id}"
    table = pathlib.Path(f"./rawdata/{filename}.yaml")

    return yaml.safe_load(table.read_text())


def get_kinematics(hepdata: dict, bin_index: list = [], indx: int = 0) -> list:
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
    if len(hepdata["independent_variables"]) > 0:
        bins = hepdata["independent_variables"][indx]["values"]
        if len(bin_index) > 0:
            bins = [bins[i] for i in bin_index]
    else:
        bins = []

    kinematics = []
    if len(bins) > 1:  # differential case
        for i in bin_index:
            ymin, ymax = [float(value) for value in bins[i]["value"].split('-')]
            kin_value = {
                "y_t": {"min": ymin, "mid": (ymin + ymax) / 2, "max": ymax},
                "m_t2": {"min": None, "mid": MT_VALUE**2, "max": None},
                "sqrts": {"min": None, "mid": SQRT_S, "max": None},
            }
            kinematics.append(kin_value)
    else:  # inclusive case
        kin_value = {"m_t2": {"min": None, "mid": MT_VALUE**2, "max": None}}
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
    errors = []
    shifts_cv = []  # relevant for asymmetric uncertainties
    for i in bin_index:
        bin = hepdata["dependent_variables"][i]

        # loop over sources of uncertainty
        error_sources = []
        shift_cv = 0
        for source in bin["values"]:
            if "errors" not in source:
                error_sources.append(None)
            elif "symerror" in source["errors"][0]:
                error_sources.append(source["errors"][0]["symerror"])
            elif "asymerror" in source["errors"][0]:

                delta_min = source["errors"][0]["asymerror"]["minus"]
                delta_plus = source["errors"][0]["asymerror"]["plus"]
                se_delta, se_sigma = se(delta_plus, delta_min)

                # shift central value
                shift_cv += se_delta

                error_sources.append(se_sigma)
        shifts_cv.append(shift_cv)
        errors.append(error_sources)
    errors = np.array(errors)

    # get the description of the uncertainty from hepdata
    error_names = [
        error_name["value"] for error_name in hepdata["independent_variables"][0]["values"]
    ]
    errors = pd.DataFrame(
        errors, columns=error_names, index=[f"bin {i}" for i in range(errors.shape[0])]
    )

    return {"errors": errors, "shifts_cv": np.array(shifts_cv)}


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
            errors["stat"] = float(uncs["statistics"].iloc[i, 0])
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
                "treatment": "ADD",
                "type": "UNCORR",
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


def main_filter() -> None:
    """
    This filter produces the commondata for the following three observables:
        1) T-Y-NORM
        2) TBAR-Y-NORM
        3) TCHANNEL-XSEC (ratio tq/tqbar)
    1) and 2) follow the same procedure.
    """

    # TCHANNEL-XSEC RATIO
    yaml_content_data = load_yaml(table_id=3, version=1)
    yaml_uncertainties = load_yaml(table_id=4, version=1)

    data_central = get_data_values(yaml_content_data, bin_index=[0], indx=2)
    kinematics = get_kinematics(yaml_content_data)
    uncertainties = get_errors(yaml_uncertainties, bin_index=[2])

    uncertainties_absolute = uncertainties["errors"] * data_central.reshape(1, -1) * 1e-2
    stat_abs = uncertainties_absolute.iloc[[0], [0]]
    sys_abs = uncertainties_absolute.iloc[[0], 1:-2]
    errors = {"statistics": stat_abs, "systematics": sys_abs}

    n_sys = sys_abs.shape[1]
    sys_types = {"treatment": ["MULT"] * n_sys, "type": ["UNCORR"] * (n_sys - 1) + ["ATLASLUMI15"]}
    sys_types_df = pd.DataFrame(sys_types, index=sys_abs.columns).T
    errors["systematics"] = pd.concat([sys_types_df, sys_abs])

    dump_commondata(kinematics, data_central, errors, obs="TCHANNEL-XSEC")

    return


if __name__ == "__main__":
    main_filter()
