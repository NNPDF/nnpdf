import pathlib
import yaml
import pandas as pd
import numpy as np

from validphys.commondata_utils import percentage_to_absolute

SQRT_S = 8_000.0  # GeV
NORM_FACTOR = 1_000.0  # include conversion from pb -> fb

#############################################################
### Unormalized Distributions. Low and High Mass regions, ###
### inclusive in RAP, mass below & above Z peak           ###
### Fixed Rapidity: 0.0 < y_{ll} < 2.4                   ###
#############################################################
MZ1_VALUE = 16.0  # GeV
MZ2_VALUE = 25.0  # GeV
MZ3_VALUE = 38.0  # GeV
MZ4_VALUE = 56.0  # GeV
MZ5_VALUE = 138.0  # GeV

# Correct tables to read values of the different Z mass ranges
#   - [MZ1] 12 GeV  <  M_{ll} < 20 GeV
#   - [MZ2] 20 GeV  <  M_{ll} < 30 GeV
#   - [MZ3] 30 GeV  <  M_{ll} < 46 GeV
#   - [MZ4] 46 GeV  <  M_{ll} < 66 GeV
#   - [MZ5] 116 GeV <  M_{ll} < 150 GeV
TABLES = {
    35: [5],
    36: [5],
    37: [5],
    38: [5],
    40: [5],
}  # {table_id: [Combined Born]}
NB_POINTS = [8, 8, 8, 20, 20]

#############################################################
### Unnormalized Distributions for fixed Invariant Mass   ###
### Fixed Invariant Mass: 66 GeV < M_{ll} < 116 GeV       ###
#############################################################
Y1_VALUE = 0.2
Y2_VALUE = 0.6
Y3_VALUE = 1.0
Y4_VALUE = 1.4
Y5_VALUE = 1.8
Y6_VALUE = 2.2

# Correct tables to read values of the different Rapidity ranges
#   - [YZ1] 0.0 < y_{ll} < 0.4
#   - [YZ2] 0.4 < y_{ll} < 0.8
#   - [YZ3] 0.8 < y_{ll} < 1.2
#   - [YZ4] 1.2 < y_{ll} < 1.6
#   - [YZ5] 1.6 < y_{ll} < 2.0
#   - [YZ6] 2.0 < y_{ll} < 0.4
TABLES_RAP = {
    29: [5],
    30: [5],
    31: [5],
    32: [5],
    33: [5],
    34: [5],
}  # {table_id: [Combined Born]}
NB_POINTS_RAP = [20, 20, 20, 20, 20, 20]

### Various Dictionary Mappings ###
MAP_TABLE = {
    29: "Y1",
    30: "Y2",
    31: "Y3",
    32: "Y4",
    33: "Y5",
    34: "Y6",
    35: "Z1",
    36: "Z2",
    37: "Z3",
    38: "Z4",
    40: "Z5",
}

MAP_BOSON = {
    "Y1": Y1_VALUE,
    "Y2": Y2_VALUE,
    "Y3": Y4_VALUE,
    "Y4": Y4_VALUE,
    "Y5": Y5_VALUE,
    "Y6": Y6_VALUE,
    "Z1": MZ1_VALUE,
    "Z2": MZ2_VALUE,
    "Z3": MZ4_VALUE,
    "Z4": MZ4_VALUE,
    "Z5": MZ5_VALUE,
}

MAP_TABID_UNCS = {
    29: "ZcombPt_born_m66116_y0004",
    30: "ZcombPt_born_m66116_y0408",
    31: "ZcombPt_born_m66116_y0812",
    32: "ZcombPt_born_m66116_y1216",
    33: "ZcombPt_born_m66116_y1620",
    34: "ZcombPt_born_m66116_y2024",
    35: "ZcombPt_born_m1220_y0024",
    36: "ZcombPt_born_m2030_y0024",
    37: "ZcombPt_born_m3046_y0024",
    38: "ZcombPt_born_m4666_y0024",
    40: "ZcombPt_born_m116150_y0024",
}


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
    filename = f"HEPData-ins1408516-v{version}-Table_{table_id}"
    table = pathlib.Path(f"./rawdata/{filename}.yaml")

    return yaml.safe_load(table.read_text())


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


def get_kinematics(hepdata: dict, bin_index: list, boson: str = "Z") -> list:
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
    rapbins = hepdata["independent_variables"][0]["values"]

    kinematics = []
    for bins in bin_index:
        ymin = float(rapbins[bins]["low"])
        ymax = float(rapbins[bins]["high"])

        if boson.startswith("Z"):
            kin_value = {
                "p_T": {"min": ymin, "mid": 0.5 * (ymin + ymax), "max": ymax},
                "M2": {"min": None, "mid": MAP_BOSON[boson] ** 2, "max": None},
                "sqrts": {"min": None, "mid": SQRT_S, "max": None},
            }
        elif boson.startswith("Y"):
            kin_value = {
                "y": {"min": None, "mid": MAP_BOSON[boson], "max": None},
                "p_T": {"min": ymin, "mid": 0.5 * (ymin + ymax), "max": ymax},
                "sqrts": {"min": None, "mid": SQRT_S, "max": None},
            }
        else:
            raise ValueError(f"Type {boson} is not recognised!")

        kinematics.append(kin_value)

    return kinematics


def get_data_values(hepdata: dict, bin_index: list, indx: int = 0) -> list:
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

    return [NORM_FACTOR * central[i]["value"] for i in bin_index]


def get_errors(
    hepdata: dict, central: list, bin_index: list, indx: int = 0
) -> dict:
    """Extract the error values from the HepData yaml file.

    NOTE: Everything is expressed as percentage.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
    central: list
        list containing the central values
    bin_index: list
        list of Non-empty bin index
    indx: int
        index from which the errors will be read

    Returns
    -------
    list:
        list of dictionaries whose contents are the various
        source of uncertainties

    """
    errors = hepdata["dependent_variables"][indx]["values"]

    stat, sys_corr, sys_uncorr, sys_lumi = [], [], [], []
    for idx in bin_index:
        stat_perc = errors[idx]["errors"][0]["symerror"]
        uncr_perc = errors[idx]["errors"][1]["symerror"]
        corr_perc = errors[idx]["errors"][2]["symerror"]
        lumi_perc = "2.8%"  # Normalisation Luminosity `ATLASLUMI12`

        stat.append(percentage_to_absolute(stat_perc, central[idx]))
        sys_uncorr.append(percentage_to_absolute(uncr_perc, central[idx]))
        sys_corr.append(percentage_to_absolute(corr_perc, central[idx]))
        sys_lumi.append(percentage_to_absolute(lumi_perc, central[idx]))

    return {
        "stat": stat,
        "sys_uncorr": sys_uncorr,
        "sys_corr": sys_corr,
        "sys_lumi": sys_lumi,
    }


def get_breakdown_syst(
    tablename: str, obs_type: str = "unnormalized"
) -> pd.DataFrame:
    """Extract the breakdown of the correlated systematic uncertainties
    from the resource files.

    Parameters
    ----------
    tablename: str
        name of the corresponding table
    obs_type: str
        type of the observable, either normalized or unnormalized

    Returns
    -------
    pd.DataFrame:
        pandas table containing the source of systematic

    """
    df_file = pd.read_csv(
        f"./rawdata/resources/{obs_type}/{tablename}/tab.dat",
        skiprows=4,
        delim_whitespace=True,
        header=None,
    )

    return df_file.dropna()


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


def format_uncertainties(
    uncs: dict, corr_systs: np.ndarray, central: list
) -> list:
    """Format the uncertainties to be dumped into the yaml file.

    Parameters
    ----------
    uncs: dict
        dictionary containing the various source of uncertainties
    corr_systs: np.ndarray
        list of correlated systematics whose first dimension refers
        to the datapoints while the second refers to the source of
        systematics

    Returns
    -------
    list:
        list of ditionaries whose elements are the various errors

    """
    combined_errors = []

    assert corr_systs.shape[0] == len(central)
    for idat in range(len(central)):
        error_value = {}
        for jdat in range(corr_systs.shape[-1]):
            error_value[f"sys_corr_{jdat + 1}"] = float(corr_systs[idat][jdat])
            # Convert Values into Absolute
            error_value[f"sys_corr_{jdat + 1}"] *= 1e-2 * central[idat]

        error_value["stat"] = uncs["stat"][idat]
        error_value["sys_uncorr"] = uncs["sys_uncorr"][idat]
        error_value["sys_luminosity"] = uncs["sys_lumi"][idat]
        error_value["sys_skipcorr"] = 1e-2 * central[idat]  # [1%]

        combined_errors.append(error_value)

    return combined_errors


def dump_commondata(
    kinematics: list,
    data: list,
    errors: list,
    number_systematics: int,
    rap_meas: str,
) -> None:
    """Function that generates and writes the commondata files.

    Parameters
    ----------
    kinematics: list
        list containing the kinematic values
    data: list
        list containing the central values
    errors: list
        list containing the different errors
    number_systematics: int
        number of systematics

    """
    error_definition = {
        f"sys_corr_{i + 1}": {
            "description": "Correlated systematic uncertainties",
            "treatment": "MULT",
            "type": "CORR",
        }
        for i in range(number_systematics)
    }

    error_definition["stat"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    error_definition["sys_uncorr"] = {
        "description": "Uncorrelated systematic uncertainties",
        "treatment": "MULT",
        "type": "UNCORR",
    }

    error_definition["sys_luminosity"] = {
        "description": "Systematic Luminosity uncertainties",
        "treatment": "MULT",
        "type": "ATLASLUMI12",
    }

    error_definition["sys_skipcorr"] = {
        "description": "Additional systematic uncertainties",
        "treatment": "MULT",
        "type": "SKIP",
    }

    suffix = "inv_dist" if rap_meas == "mdist_unnorm" else "rap_dist"

    with open(f"data_{suffix}.yaml", "w") as file:
        yaml.dump({"data_central": data}, file, sort_keys=False)

    with open(f"kinematics_{suffix}.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open(f"uncertainties_{suffix}.yaml", "w") as file:
        yaml.dump(
            {"definitions": error_definition, "bins": errors},
            file,
            sort_keys=False,
        )


def main_filter(rap_meas: str = "mdist_unnorm") -> None:
    """Main driver of the filter that produces commmondata.

    There are five main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR

    2. Correlated Systematic uncertainties: MULT, CORR

    3. Uncorrelated Systematic uncertainties: MULT, UNCORR

    4. Luminosity Systematic uncertainties: MULT, ATLASLUMI12

    5. Additional Systematic uncertainties: MULT, SKIP

    """
    version, _, _ = read_metadata()

    corr_nbpots = NB_POINTS if rap_meas == "mdist_unnorm" else NB_POINTS_RAP
    corr_tables = TABLES if rap_meas == "mdist_unnorm" else TABLES_RAP

    nbp_idx = 0
    comb_kins, comb_data = [], []
    combined_errors, corr_syserrors = [], []
    for tabid in corr_tables.keys():  # Loop over tables
        syscorr_df = get_breakdown_syst(MAP_TABID_UNCS[tabid])
        for idx in corr_tables[tabid]:
            bin_index = [i for i in range(corr_nbpots[nbp_idx])]
            yaml_content = load_yaml(table_id=tabid, version=version)

            # Extract the kinematic, data, and uncertainties
            kinematics = get_kinematics(
                yaml_content, bin_index, MAP_TABLE[tabid]
            )
            data_central = get_data_values(yaml_content, bin_index, indx=idx)
            uncertainties = get_errors(
                yaml_content, data_central, bin_index, indx=idx
            )

            # Collect all the results from different tables
            comb_kins += kinematics
            comb_data += data_central
            combined_errors.append(uncertainties)
            corr_syserrors.append(syscorr_df.values)

            nbp_idx += 1

    combine_syscorrs = np.concatenate(corr_syserrors, axis=0)
    errors_combined = concatenate_dicts(combined_errors)
    errors = format_uncertainties(errors_combined, combine_syscorrs, comb_data)

    # Generate all the necessary files
    dump_commondata(
        comb_kins, comb_data, errors, combine_syscorrs.shape[-1], rap_meas
    )
    print(f"Generation of {rap_meas} completed!!!")

    return


if __name__ == "__main__":
    main_filter(rap_meas="mdist_unnorm")  # Unnormalized Distr. in Inv. Mass
    main_filter(rap_meas="ydist_unnorm")  # Unnormalized Distr. in Rapidity
