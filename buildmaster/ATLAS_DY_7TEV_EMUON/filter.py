import pathlib
import yaml

from validphys.commondata_utils import percentage_to_absolute


LUMI_CORR = 1.0187  # Correction Factor due to Luminosity Upgrade

MZ_VALUE = 91.1876  # GeV
MW_VALUE = 80.398  # GeV
SQRT_S = 7_000.0  # GeV
NORM_FACTOR = 1_000.0 * LUMI_CORR  # include conversion from pb -> fb

# Correct tables to read values [[W+], [W-], [Z]]
TABLES = {5: [0], 3: [0], 1: [0]}  # {table_id: [indexes]}
NB_POINTS = [11, 11, 8]  # [N_Z, N_W+, N_W-]
MAP_TAB_UNC = {1: 2, 3: 4, 5: 6}  # Map Data table to Unc. tables

# MAP Boson string to M values
MAP_BOSON = {"Z": MZ_VALUE, "W": MW_VALUE}
MAP_TABLE = {1: "Z", 3: "W", 5: "W"}


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
    filename = f"HEPData-ins928289-v{version}-Table_{table_id}"
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
        stat.append(NORM_FACTOR * errors[idx]["errors"][0]["symerror"])
        sys_uncorr.append(NORM_FACTOR * errors[idx]["errors"][1]["symerror"])
        sys_corr.append(NORM_FACTOR * errors[idx]["errors"][2]["symerror"])

        syslumi_abs = percentage_to_absolute(
            errors[idx]["errors"][3]["symerror"], central[idx]
        )
        sys_lumi.append(syslumi_abs)  # TODO: Why OLD=3.5???

    return {
        "stat": stat,
        "sys_uncorr": sys_uncorr,
        "sys_corr": sys_corr,
        "sys_lumi": sys_lumi,
    }


def get_corrsyst(tabid_corrsyst: int, central: list) -> list:
    """Load the Hepdata table that contains the correlated Systematics.

    NOTE: The correlated systematics are given as percentage

    Parameters
    ----------
    tabid_corrsyst: int
        table ID corresponding to the given Boson
    central: list
        list containing the central values

    Returns
    -------
    list[list]:
        list of lists whose first dimension represents the different
        datapoints while the second represent the number of systs.

    """
    uncfile = load_yaml(table_id=tabid_corrsyst)
    datapts = []
    for idx, dpts in enumerate(uncfile["dependent_variables"]):
        source_syst = []
        for sys in dpts["values"]:
            # Convert values into Absolute
            source_syst.append(1e-2 * sys["value"] * central[idx])
        datapts.append(source_syst)
    return datapts


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


def format_uncertainties(uncs: dict, corr_systs: list) -> list:
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

        error_value["stat"] = uncs["stat"][idat]
        error_value["sys_uncorr"] = uncs["sys_uncorr"][idat]
        error_value["sys_luminosity"] = uncs["sys_lumi"][idat]
        combined_errors.append(error_value)

    return combined_errors


def dump_commondata(kinematics: list, data: list, errors: list, number_systematics: int) -> None:
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
        "type": "ATLASLUMI10",
    }

    with open("data.yaml", "w") as file:
        yaml.dump({"data_central": data}, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(
            {"definitions": error_definition, "bins": errors},
            file,
            sort_keys=False,
        )


def main_filter() -> None:
    """Main driver of the filter that produces commmondata.

    There are four main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR

    2. Correlated Systematic uncertainties: MULT, CORR:

    3. Uncorrelated Systematic uncertainties: MULT, UNCORR

    4. Luminosity Systematic uncertainties: MULT, ATLASLUMI10

    """
    version, _, _ = read_metadata()

    nbp_idx = 0
    comb_kins, comb_data = [], []
    combined_errors, corr_syserrors = [], []
    for tabid in TABLES.keys():  # Loop over tables
        for idx in TABLES[tabid]:  # Loop over Bosons [Z, W+, W-]
            bin_index = [i for i in range(NB_POINTS[nbp_idx])]
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
            # Combine all the correlated systematic uncertainties
            corr_syserrors += get_corrsyst(MAP_TAB_UNC[tabid], data_central)

            nbp_idx += 1

    assert sum(NB_POINTS) == len(corr_syserrors)
    errors_combined = concatenate_dicts(combined_errors)
    errors = format_uncertainties(errors_combined, corr_syserrors)

    # Generate all the necessary files
    dump_commondata(comb_kins, comb_data, errors, len(corr_syserrors))

    return


if __name__ == "__main__":
    main_filter()
