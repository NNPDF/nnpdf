import pathlib

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.new_commondata.ATLAS_TTBAR_13TEV_HADR_DIF.utils import covmat_to_artunc

# MZ_VALUE = 91.1876  # GeV
# MW_VALUE = 80.398  # GeV
# NORM_FACTOR = 1e-2  # Correct for Units
#
TABLES = {1: [0]}  # {table_id: [indexes]}
NB_POINTS = [132]
SQRT_S = 7_000.0  # GeV
#
# # MAP Boson string to M values
# MAP_BOSON = {"Z": MZ_VALUE, "W": MW_VALUE}
# MAP_TABLE = {1: "W"}


def load_rawdata() -> pd.DataFrame:
    """Load the raw data into a Panda table.

    Returns
    -------
    pd.DataFrame:
        table containing the information on the dataset


    """
    return pd.read_csv(
        "./rawdata/CMS-DY2D11-ABS.data", delim_whitespace=True, names=['y', 'M', 'sigma']
    )


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


def get_kinematics(hepdata: pd.DataFrame, bin_index: list, boson: str = "Z") -> list:
    """Read the version and list of tables from metadata.

    Parameters
    ----------
    hepdata: pd.DataFrame
        table containing all data info
    bin_index: list
        list of Non-empty bin index

    Returns
    -------
    tuple(int, list):
        data version and list of hepdata tables

    """
    del boson
    y = hepdata["y"].values
    M = hepdata["M"].values

    kinematics = []
    for bins in bin_index:
        kin_value = {
            "y": {"min": None, "mid": float(y[bins]), "max": None},
            "M2": {"min": None, "mid": float(M[bins]) ** 2, "max": None},
            "sqrts": {"min": None, "mid": SQRT_S, "max": None},
        }
        kinematics.append(kin_value)

    return kinematics


def get_data_values(hepdata: pd.DataFrame, bin_index: list, indx: int = 0) -> list:
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
    del indx
    sigma = hepdata["sigma"].values

    return [float(sigma[i]) for i in bin_index]


def read_corrmatrix(nb_datapoints: int) -> np.ndarray:
    """Read the covariance matrix from a table.

    Parameters
    ----------
    nb_datapoints: int
        total number of datapoints

    Returns
    -------
    np.ndarray: (nb_datapoints, nb_datapoints)
        entries of the corr/cov-mat as an array

    """
    df_corrmat = pd.read_csv("./rawdata/covmat.corr", delim_whitespace=True, header=None)
    corrmat = df_corrmat.iloc[:, 2].values
    return corrmat.reshape(nb_datapoints, nb_datapoints)


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

    return covmat_sys.T.flatten()


def generate_artificial_unc(**kwargs) -> np.ndarray:
    """A wrapper around `covmat_to_artunc` to return an array.

    Returns
    -------
    np.ndarray:
        a squared matrix with artificial uncertainties
    """
    artunc = covmat_to_artunc(**kwargs)
    return np.array(artunc)


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


def format_uncertainties(artunc: np.ndarray, central: list) -> list:
    """Format the uncertainties to be dumped into the yaml file.

    Parameters
    ----------
    artunc: no.ndarray
        array of the artificial systematic uncertainties
    central: list
        list containing the complete central values

    Returns
    -------
    list:
        list of ditionaries whose elements are the various errors

    """
    combined_errors = []

    artunc = artunc.tolist()
    for idat in range(len(artunc)):
        error_value = {}
        for jdat, value in enumerate(artunc[idat]):
            error_value[f"sys_corr_{jdat + 1}"] = value

        error_value["stat"] = 0.0  # Stat. Errs contained in CovMat
        error_value["luminosity"] = 2.2 * 1e-2 * central[idat]
        combined_errors.append(error_value)

    return combined_errors


def dump_commondata(kinematics: list, data: list, errors: list, nb_syscorr: int) -> None:
    """Function that generates and writes the commondata files.

    Parameters
    ----------
    kinematics: list
        list containing the kinematic values
    data: list
        list containing the central values
    errors: list
        list containing the different errors
    nb_syscorr: int
        number of source of Systematics

    """
    error_definition = {
        f"sys_corr_{i + 1}": {
            "description": "Correlated systematic uncertainties",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(nb_syscorr)
    }

    error_definition["stat"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    error_definition["luminosity"] = {
        "description": "Correlated luminosity Normalization",
        "treatment": "MULT",
        "type": "CMSLUMI11",
    }

    with open("data.yaml", "w") as file:
        yaml.dump({"data_central": data}, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open("uncertainties.yaml", "w") as file:
        yaml.dump({"definitions": error_definition, "bins": errors}, file, sort_keys=False)


def main_filter() -> None:
    """Main driver of the filter that produces commmondata.

    There are two main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR == ZERO

    2. Correlated Systematic uncertainties: ADD, CORR:
        constructed from the covariance matrix provided with the
        measurements

    3. Normalization Luminosity uncertainties: MULT, CMSLUMI11

    """
    nbpoints = sum(NB_POINTS)

    nbp_idx = 0
    comb_kins, comb_data = [], []
    for tabid in TABLES.keys():  # Loop over tables
        for _ in TABLES[tabid]:  # Loop over Bosons
            bin_index = [i for i in range(NB_POINTS[nbp_idx])]
            df_content = load_rawdata()

            # Extract the kinematic, data, and uncertainties
            kinematics = get_kinematics(df_content, bin_index)
            data_central = get_data_values(df_content, bin_index)

            # Collect all the results from different tables
            comb_kins += kinematics
            comb_data += data_central
            nbp_idx += 1

    # Compute the Artifical Systematics from CovMat
    corrmat = read_corrmatrix(nb_datapoints=nbpoints)
    artunc = generate_artificial_unc(
        ndata=nbpoints, covmat_list=corrmat.T.flatten().tolist(), no_of_norm_mat=0
    )
    errors = format_uncertainties(artunc, comb_data)

    # Generate all the necessary files
    dump_commondata(comb_kins, comb_data, errors, artunc.shape[-1])

    return


if __name__ == "__main__":
    main_filter()
