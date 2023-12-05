import pathlib
import yaml
import pandas as pd

from validphys.commondata_utils import percentage_to_absolute

MZ_LOW = 56.0  # GeV
MZ_PEAK = 91.0  # GeV
MZ_HIGH = 133.0  # GeV
MW_VALUE = 80.398  # GeV
SQRT_S = 7_000.0  # GeV
NORM_FACTOR = 1_000.0  # include conversion from pb -> fb

# Correct tables to read values [[W+], [W-], [Z_low], [Z_peak], [Z_high]]
TABLES = {9: [0], 10: [0], 11: [0], 12: [0], 13: [0]}  # {table_id: [indexes]}
NB_POINTS = [11, 11, 6, 12, 6]  # [_W+, N_W-, N_Zlow, N_Zpeak, N_Zhigh]

# Correct tables to read values for FORWARD [[Z_peak], [Z_high]]
TABLES_FWD = {14: [0], 15: [0]}  # {table_id: [indexes]}
NB_POINTS_FWD = [9, 6]  # [N_Zpeak, N_Zhigh]

MAP_TAB_UNC = {
    9: "wplus",
    10: "wminus",
    11: "zylow_cc",
    12: "zypeak_cc",
    13: "zyhigh_cc",
    14: "zypeak_cf",
    15: "zyhigh_cf",
}  # Map Data table to Unc. tables

# MAP Boson string to M values
MAP_BOSON = {
    "Z_low": MZ_LOW,
    "Z_peak": MZ_PEAK,
    "Z_high": MZ_HIGH,
    "W": MW_VALUE,
}
MAP_TABLE = {
    9: "W",
    10: "W",
    11: "Z_low",
    12: "Z_peak",
    13: "Z_high",
    14: "Z_peak",
    15: "Z_high",
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
    filename = f"HEPData-ins1502620-v{version}-Table_{table_id}"
    table = pathlib.Path(f"./rawdata/{filename}.yaml")

    return yaml.safe_load(table.read_text())


def load_fulluncs(table_name: str) -> pd.DataFrame:
    """Load the complete breakdown uncertainties into Pandas table.

    Parameters
    ----------
    table_name: str
        name of the .dat file

    Returns
    -------
    pd.DataFrame:
        pandas table containing all the uncertainty breakdown

    """
    return pd.read_csv(
        f"./rawdata/{table_name}.dat",
        skiprows=28,
        delim_whitespace=True,
        header=None,
    )


def get_adduncorr(full_table: pd.DataFrame) -> dict:
    """Extract the additional Uncorrelated [SYST] uncertainties.

    Parameters
    ----------
    full_table: pd.DataFrame
        complete table of uncertainties

    Returns
    -------
    dict:
        dictionary containing all the values of the uncertainty

    """
    # TODO: double-check the origin of this uncertainty
    return {"sys_adduncorr": full_table[full_table.columns[-1]].values.tolist()}


def get_corrsyst(full_table: pd.DataFrame) -> list:
    """Extract the Correlated Systematic uncertainties form the breakdown.

    Parameters
    ----------
    full_table: pd.DataFrame
        complete table of uncertainties

    Returns
    -------
    list:
        array containing all the values of the correlated systematics

    """
    syscorr = full_table.iloc[:, list(range(6, full_table.shape[1] - 1))]
    return syscorr.values.tolist()


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
        kin_value = {
            "y": {"min": ymin, "mid": 0.5 * (ymin + ymax), "max": ymax},
            "M2": {"min": None, "mid": MAP_BOSON[boson] ** 2, "max": None},
            "sqrt_s": {"min": None, "mid": SQRT_S, "max": None},
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
        lumi_perc = errors[idx]["errors"][3]["symerror"]

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


def format_uncertainties(uncs: dict, corr_systs: list, central: list) -> list:
    """Format the uncertainties to be dumped into the yaml file.

    Parameters
    ----------
    uncs: dict
        dictionary containing the various source of uncertainties
    corr_systs: list[list]
        list of correlated systematics whose first dimension refers
        to the datapoints while the second refers to the source of
        systematics

    Returns
    -------
    list:
        list of ditionaries whose elements are the various errors

    """
    combined_errors = []

    for idat in range(len(corr_systs)):
        error_value = {}
        for jdat in range(len(corr_systs[idat])):
            error_value[f"sys_corr_{jdat + 1}"] = corr_systs[idat][jdat]
            # Convert Values into Absolute
            error_value[f"sys_corr_{jdat + 1}"] *= 1e-2 * central[idat]

        error_value["stat"] = uncs["stat"][idat]
        error_value["sys_uncorr"] = uncs["sys_uncorr"][idat]
        error_value["sys_luminosity"] = uncs["sys_lumi"][idat]

        # TODO: Add here the missing correlated systematics
        # error_value["sys_adduncorr"] = uncs["sys_adduncorr"][idat]
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
            "type": f"ATLASWZRAP11_{i + 1001}",
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
        "type": "ATLASLUMI11",
    }

    # error_definition["sys_adduncorr"] = {
    #     "description": "Additional Uncorrelated uncertainty",
    #     "treatment": "MULT",
    #     "type": "UNCORR",
    # }

    suffix = "crap" if rap_meas == "central" else "frap"

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


def main_filter(rap_meas: str = "central") -> None:
    """Main driver of the filter that produces commmondata.

    There are five main different sources of uncertainties.

    1. Statistical uncertainties: MULT, UNCORR

    2. Correlated Systematic uncertainties: MULT, CORR

    3. Uncorrelated Systematic uncertainties: MULT, UNCORR

    4. Luminosity Systematic uncertainties: MULT, ATLASLUMI11

    5. Uncorrelated [Systematic] uncertainties: MULT, UNCORR

    """
    version, _, _ = read_metadata()

    corr_nbpots: list = NB_POINTS if rap_meas == "central" else NB_POINTS_FWD
    corr_tables: dict = TABLES if rap_meas == "central" else TABLES_FWD

    nbp_idx = 0
    comb_kins, comb_data = [], []
    combined_errors, corr_syserrors = [], []
    for tabid in corr_tables.keys():  # Loop over tables
        fulluncs = load_fulluncs(table_name=MAP_TAB_UNC[tabid])
        for idx in corr_tables[
            tabid
        ]:  # Loop over Bosons [Z_low, Z_peak, Z_high, W+, W-]
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
            # TODO: Activate when addition correlated is sorted out
            # uncertainties.update(get_adduncorr(full_table=fulluncs))

            # Collect all the results from different tables
            comb_kins += kinematics
            comb_data += data_central
            combined_errors.append(uncertainties)
            # Combine all the correlated systematic uncertainties
            corr_syserrors += get_corrsyst(full_table=fulluncs)

            nbp_idx += 1

    errors_combined = concatenate_dicts(combined_errors)
    errors = format_uncertainties(errors_combined, corr_syserrors, comb_data)

    # Generate all the necessary files
    dump_commondata(
        comb_kins, comb_data, errors, len(corr_syserrors[0]), rap_meas
    )
    print(f"Generation of {rap_meas} completed!!!")

    return


if __name__ == "__main__":
    main_filter(rap_meas="central")
    main_filter(rap_meas="forward")
