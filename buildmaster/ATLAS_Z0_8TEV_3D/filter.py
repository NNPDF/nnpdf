import pathlib
import yaml


MZ_VALUE = 91.1876  # GeV
MW_VALUE = 80.398  # GeV
SQRT_S = 8_000.0  # GeV
NORM_FACTOR = 1_000.0  # include conversion from pb -> fb

# Correct tables to read values [[Z]]
TABLES = {4: [0]}  # {table_id: [indexes]}
NB_POINTS = [123]  # [N_Z]

# MAP Boson string to M values
MAP_BOSON = {"Z": MZ_VALUE, "W": MW_VALUE}
MAP_TABLE = {4: "W"}


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
    filename = f"HEPData-ins1630886-v{version}-Table_{table_id}"
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
    del boson
    cotheta = hepdata["independent_variables"][0]["values"]  # Cos(Theta)
    rapbins = hepdata["independent_variables"][1]["values"]  # Rapidity
    massbin = hepdata["independent_variables"][2]["values"]  # Invariant Masses

    kinematics = []
    for bins in bin_index:
        tmin, tmax = float(cotheta[bins]["low"]), float(cotheta[bins]["high"])
        ymin, ymax = float(rapbins[bins]["low"]), float(rapbins[bins]["high"])
        mmin, mmax = float(massbin[bins]["low"]), float(massbin[bins]["high"])

        kin_value = {
            "y": {"min": ymin, "mid": 0.5 * (ymin + ymax), "max": ymax},
            "M": {"min": mmin, "mid": 0.5 * (mmin + mmax), "max": mmax},
            "cos_theta": {"min": tmin, "mid": 0.5 * (tmin + tmax), "max": tmax},
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
    massbin = hepdata["independent_variables"][-1]["values"]

    # Normalization to match APPLgrid Predictions
    return [
        0.4
        * (massbin[i]["high"] - massbin[i]["low"])
        * NORM_FACTOR
        * central[i]["value"]
        for i in bin_index
    ]


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
    masbin = hepdata["independent_variables"][-1]["values"]  # Lepton Mass

    stat, sys_lumi, sys_uncorr = [], [], []

    sysdicts: dict = {
        f"sys_corr_{n}": [] for n in range(1, len(errors[0]["errors"]) - 1)
    }

    for idc, idx in enumerate(bin_index):
        # Prefactor that includes APPLgrid Normalization and Unit conversion
        prefac = 0.4 * NORM_FACTOR * (masbin[idx]["high"] - masbin[idx]["low"])

        stat.append(prefac * errors[idx]["errors"][0]["symerror"])

        for sysid in range(1, len(errors[idx]["errors"]) - 1):
            sysdicts[f"sys_corr_{sysid}"].append(
                prefac * errors[idx]["errors"][sysid]["symerror"]
            )

        sys_uncorr.append(prefac * errors[idx]["errors"][-1]["symerror"])

        # Overall Luminosity uncertainty of 1.9%
        sys_lumi.append(1.9 * central[idc] * 1e-2)

    sysdicts.update(
        {"stat": stat, "sys_lumi": sys_lumi, "sys_uncorr": sys_uncorr}
    )

    return sysdicts


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


def format_uncertainties(uncs: dict) -> list:
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

    for idat in range(len(uncs['stat'])):
        error_value = {k: uncs[k][idat] for k in uncs.keys()}
        combined_errors.append(error_value)

    return combined_errors


def dump_commondata(
    kinematics: list, data: list, errors: list, rap_meas: str
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
    rap_meas: str
        measurement on the rapidity: central or forward

    """
    nb_syscorr = [1 for n in errors[0].keys() if "sys_corr_" in n]
    suffix = "crap" if rap_meas == "central" else "frap"

    error_definition = {
        f"sys_corr_{i + 1}": {
            "description": "Correlated systematic uncertainties",
            "treatment": "MULT",
            "type": "CORR",
        }
        for i in range(sum(nb_syscorr))
    }

    error_definition["stat"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    error_definition["sys_uncorr"] = {
        "description": "Additive Uncorrelated Systematic uncertainties",
        "treatment": "MULT",
        "type": "UNCORR",
    }

    error_definition["sys_lumi"] = {
        "description": "Luminosity systematic uncertainties",
        "treatment": "MULT",
        "type": "ATLASLUMI12",
    }

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

    There are four main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR

    2. Correlated Systematic uncertainties: MULT, CORR:

    3. Uncorrelated Systematic uncertainties: MULT, UNCORR

    4. Luminosity Systematic uncertainties: MULT, ATLASLUMI12

    """
    version, _, _ = read_metadata()
    start_dpts = 0 if rap_meas == "central" else 458

    nbp_idx = 0
    comb_kins, comb_data = [], []
    combined_errors = []
    for tabid in TABLES.keys():  # Loop over tables
        for idx in TABLES[tabid]:  # Loop over Bosons [Z, W+, W-]
            bin_index = [i + start_dpts for i in range(NB_POINTS[nbp_idx])]
            yaml_content = load_yaml(table_id=tabid, version=version)

            # Extract the kinematic, data, and uncertainties
            kinematics = get_kinematics(
                yaml_content, bin_index, MAP_TABLE[tabid]
            )
            data_central = get_data_values(yaml_content, bin_index, indx=idx)
            uncertainties = get_errors(
                yaml_content, data_central, bin_index, idx
            )

            # Collect all the results from different tables
            comb_kins += kinematics
            comb_data += data_central
            combined_errors.append(uncertainties)

            nbp_idx += 1

    errors_combined = concatenate_dicts(combined_errors)
    errors = format_uncertainties(errors_combined)

    # Generate all the necessary files
    dump_commondata(comb_kins, comb_data, errors, rap_meas)
    print(f"Finished generating Z 3D for {rap_meas}!")

    return


if __name__ == "__main__":
    main_filter(rap_meas="central")
    main_filter(rap_meas="forward")
