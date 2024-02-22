import pathlib

import numpy as np
import pandas as pd
import yaml

from validphys.commondata_utils import covmat_to_artunc

MZ_VALUE = 91.1876  # GeV
MW_VALUE = 80.398  # GeV
SQRT_S = 7_000.0  # GeV
NORM_FACTOR = 1_000.0  # from pb -> fb

# Correct tables to read values [[Z], [W+, W-]]
TABLES = {1: [0], 4: [0, 2]}  # {table_id: [indexes]}
NB_POINTS = [17, 8, 8]  # [N_Z, N_W+, N_W-]

# MAP Boson string to M values
MAP_BOSON = {"Z": MZ_VALUE, "W": MW_VALUE}
MAP_TABLE = {1: "Z", 4: "W"}


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
    filename = f"HEPData-ins1373300-v{version}-Table_{table_id}"
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
            "sqrts": {"min": None, "mid": SQRT_S, "max": None},
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


def get_errors(hepdata: dict, bin_index: list, indx: int = 0) -> dict:
    """Extract the error values from the HepData yaml file.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
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

    stat, sys_corr, sys_beam, sys_lumi = [], [], [], []
    for idx in bin_index:
        stat.append(NORM_FACTOR * errors[idx]["errors"][0]["symerror"])
        sys_corr.append(NORM_FACTOR * errors[idx]["errors"][1]["symerror"])
        sys_beam.append(NORM_FACTOR * errors[idx]["errors"][2]["symerror"])
        sys_lumi.append(NORM_FACTOR * errors[idx]["errors"][3]["symerror"])

    return {"stat": stat, "sys_corr": sys_corr, "sys_beam": sys_beam, "sys_lumi": sys_lumi}


def read_corrmatrix(nb_datapoints: int) -> np.ndarray:
    """Read the matrix and returns a symmetrized verions.

    Parameters
    ----------
    nb_datapoints: int
        total number of datapoints

    Returns
    -------
    np.ndarray: (nb_datapoints, nb_datapoints)
        a symmetrized square matrix

    """
    corrmat = pd.read_csv(
        "./rawdata/corrmat.corr",
        names=[f'{i}' for i in range(nb_datapoints)],
        delim_whitespace=True,
    )
    corrmat = corrmat.iloc[:, :].values

    # Symmetrize the Correlation Matrix
    symmetrised = np.triu(corrmat.T) + np.triu(corrmat.T, 1).T
    assert np.allclose(symmetrised, symmetrised.T)

    return symmetrised


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


def format_uncertainties(uncs: dict, artunc: np.ndarray, bslice: slice) -> list:
    """Format the uncertainties to be dumped into the yaml file.

    Parameters
    ----------
    uncs: dict
        dictionary containing the various source of uncertainties
    artunc: no.ndarray
        array of the artificial systematic uncertainties

    bslice: slice
        slice from which the subtables will be selected

    Returns
    -------
    list:
        list of ditionaries whose elements are the various errors

    """
    del bslice
    combined_errors = []

    artunc = artunc.tolist()
    for idat in range(len(artunc)):
        error_value = {}
        for jdat, value in enumerate(artunc[idat]):
            error_value[f"sys_corr_{jdat + 1}"] = value

        error_value["stat"] = uncs["stat"][idat]
        error_value["sys_beam"] = uncs["sys_beam"][idat]
        error_value["sys_luminosity"] = uncs["sys_lumi"][idat]
        combined_errors.append(error_value)

    return combined_errors


def dump_commondata(kinematics: list, data: list, errors: list, nbpoints: int) -> None:
    """Function that generates and writes the commondata files.

    Parameters
    ----------
    kinematics: list
        list containing the kinematic values
    data: list
        list containing the central values
    errors: list
        list containing the different errors
    nbpoints: int
        total number of points including Z, W+/-

    """
    error_definition = {
        f"sys_corr_{i + 1}": {
            "description": "Correlated systematic uncertainties",
            "treatment": "MULT",
            "type": f"LHCBWZMU7TEV_{i}",
        }
        for i in range(nbpoints)
    }

    error_definition["stat"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    error_definition["sys_beam"] = {
        "description": "Systematic Beam uncertainties",
        "treatment": "MULT",
        "type": "LHCBBEAM7TEV",
    }

    error_definition["sys_luminosity"] = {
        "description": "Systematic Luminosity uncertainties",
        "treatment": "MULT",
        "type": "LHCBLUMI7TEV",
    }

    with open("data.yaml", "w") as file:
        yaml.dump({"data_central": data}, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open("uncertainties.yaml", "w") as file:
        yaml.dump({"definitions": error_definition, "bins": errors}, file, sort_keys=False)


def main_filter(boson: str = "Z") -> None:
    """Main driver of the filter that produces commmondata.

    There are four main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR

    2. Correlated Systematic uncertainties: MULT, CORR:
        constructed by symmetrizing the correlation matrix and extracting
        the artificial systematic uncertainties from the latter

    3. Beam Systematic uncertainties: MULT, LHCBBEAM7TEV

    4. Luminosity Systematic uncertainties: MULT, LHCBLUMI7TEV

    """
    version, _, _ = read_metadata()
    nbpoints = sum(NB_POINTS)

    # Select only the bins corresponding to the correct boson
    bslice = slice(17) if boson == "Z" else slice(17, nbpoints + 1)

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
    errors = format_uncertainties(errors_combined, artunc, bslice)

    # Generate all the necessary files
    dump_commondata(comb_kins[bslice], comb_data[bslice], errors[bslice], nbpoints)

    return


if __name__ == "__main__":
    main_filter(boson="Z")
