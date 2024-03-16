import pathlib
import yaml
import numpy as np

from validphys.commondata_utils import covmat_to_artunc


MZ_VALUE = 91.1876  # GeV
MW_VALUE = 80.385  # GeV
SQRT_S = 13_000.0  # GeV
NORM_FACTOR = 1e6  # include conversion from pb -> fb

# Correct tables to read values [[W-], [W+], [Z]]
TABLES = {9: [0], 8: [0], 11: [0]}  # {table_id: [indexes]}
NB_POINTS = [1, 1, 1]  # [N_W-, N_W+, N_Z]

# MAP Boson string to M values
MAP_BOSON = {"Z": MZ_VALUE, "W": MW_VALUE}
MAP_TABLE = {8: "W", 9: "W", 11: "W"}
CORRMAT_NDIM = 3


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
    filename = f"HEPData-ins1436497-v{version}-Table_{table_id}"
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
    del hepdata
    kinematics = []
    for _ in bin_index:
        kin_value = {
            "mu2": {
                "min": None,
                "mid": NORM_FACTOR * MAP_BOSON[boson],
                "max": None,
            },
            "sqrts": {"min": None, "mid": 1e3 * SQRT_S, "max": None},
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
    del central
    errors = hepdata["dependent_variables"][indx]["values"]

    stat, sys_lumi, sys_corr = [], [], []
    for idx in bin_index:
        stat.append(NORM_FACTOR * errors[idx]["errors"][0]["symerror"])
        sys_corr.append(NORM_FACTOR * errors[idx]["errors"][1]["symerror"])
        sys_lumi.append(NORM_FACTOR * errors[idx]["errors"][2]["symerror"])

    return {"stat": stat, "sys_corr": sys_corr, "sys_lumi": sys_lumi}


def compute_correlation_ceoffs(table_id: int) -> np.ndarray:
    """Extract the correlation coefficients from HepData and compute
    the correlation matrix.

    Parameters
    ----------
    table_id: int
        table ID that contains the correlation matrix

    Returns
    -------
    np.ndarray: (n_corr, n_corr)
        a square symmetrized matrix containing the correlations

    """
    hepdata = load_yaml(table_id=table_id)
    correlation = hepdata["dependent_variables"][0]["values"]

    corrmatrix = np.identity(n=3)
    corrmatrix[0][1] = correlation[-3]["value"]
    corrmatrix[0][2] = correlation[-2]["value"]
    corrmatrix[1][2] = correlation[-1]["value"]

    # Symmetrize the correlation matrix
    return corrmatrix.T + np.triu(corrmatrix, 1)


def multiply_syst(covmat: np.ndarray, sys: list) -> np.ndarray:
    """Multiply the covmat with the total systematics.

    Parameters
    ----------
    covmat: np.ndarray
        a squared matrix representing the covmat
    sys: list
        a vector containing the correlatdd systematics

    Returns
    -------
    np.ndarray:
        covariance matrix multiplied by the total systematics

    """
    covmat_sys = np.zeros(covmat.shape)

    for i in range(covmat.shape[0]):
        for j in range(covmat.shape[-1]):
            covmat_sys[i][j] = covmat[i][j] * sys[i] * sys[j]

    return covmat_sys.flatten()


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


def generate_artificial_unc(**kwargs) -> np.ndarray:
    """A wrapper around `covmat_to_artunc` to return an array.

    Returns
    -------
    np.ndarray:
        a squared matrix with artificial uncertainties
    """
    artunc = covmat_to_artunc(**kwargs)
    return np.array(artunc)


def format_uncertainties(uncs: dict, artunc: np.ndarray) -> list:
    """Format the uncertainties to be dumped into the yaml file.

    Parameters
    ----------
    uncs: dict
        dictionary containing the various source of uncertainties
    artunc: no.ndarray
        array of the artificial systematic uncertainties

    Returns
    -------
    list:
        list of ditionaries whose elements are the various errors

    """
    combined_errors = []

    artunc = artunc.tolist()
    for idat in range(len(artunc)):
        error_value = {}
        for jdat in range(len(artunc[idat])):
            error_value[f"sys_corr_{jdat + 1}"] = artunc[idat][jdat]

        error_value["stat"] = uncs["stat"][idat]
        error_value["sys_lumi"] = uncs["sys_lumi"][idat]
        combined_errors.append(error_value)

    return combined_errors


def dump_commondata(kinematics: list, data: list, errors: list) -> None:
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
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(CORRMAT_NDIM)
    }

    error_definition["stat"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    error_definition["sys_lumi"] = {
        "description": "Luminosity systematic uncertainties",
        "treatment": "MULT",
        "type": "ATLASLUMI13",
    }

    with open(f"data.yaml", "w") as file:
        yaml.dump({"data_central": data}, file, sort_keys=False)

    with open(f"kinematics.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open(f"uncertainties.yaml", "w") as file:
        yaml.dump(
            {"definitions": error_definition, "bins": errors},
            file,
            sort_keys=False,
        )


def main_filter() -> None:
    """Main driver of the filter that produces commmondata.

    There are three main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR

    2. Correlated Systematic uncertainties: ADD, CORR:

    3. Luminosity Systematic uncertainties: MULT, ATLASLUMI13

    """
    version, _, _ = read_metadata()
    nbpoints = sum(NB_POINTS)

    nbp_idx = 0
    comb_kins, comb_data = [], []
    combined_errors = []
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
                yaml_content, data_central, bin_index, idx
            )

            # Collect all the results from different tables
            comb_kins += kinematics
            comb_data += data_central
            combined_errors.append(uncertainties)

            nbp_idx += 1

    errors_combined = concatenate_dicts(combined_errors)
    # Compute the Artifical Systematics from Correlation Matrix
    corrmat = compute_correlation_ceoffs(table_id=24)
    covmat = multiply_syst(corrmat, errors_combined["sys_corr"])
    artunc = generate_artificial_unc(
        ndata=nbpoints, covmat_list=covmat.tolist(), no_of_norm_mat=0
    )
    errors = format_uncertainties(errors_combined, artunc)

    # Generate all the necessary files
    dump_commondata(comb_kins, comb_data, errors)

    return


if __name__ == "__main__":
    main_filter()
