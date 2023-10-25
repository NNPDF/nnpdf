import numpy as np
import pandas as pd
import pathlib
import yaml

from validphys.commondata_utils import covmat_to_artunc, percentage_to_absolute


MZ_VALUE = 91.1876  # GeV
SQRT_S = 7_000.0  # GeV
NORM_FACTOR = 1_000.0  # from pb -> fb


def load_yaml(table_id: int) -> dict:
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
    table = pathlib.Path(f"./rawdata/Table{table_id}.yaml")

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


def get_kinematics(hepdata: dict) -> list:
    """Read the version and list of tables from metadata.

    Returns
    -------
    tuple(int, list):
        data version and list of hepdata tables

    """
    rapbins = hepdata["independent_variables"][0]["values"]

    kinematics = []
    for bins in rapbins:
        kin_value = {
            "y": {
                "min": bins["low"],
                "mid": 0.5 * (bins["low"] + bins["high"]),
                "max": bins["high"],
            },
            "M2": {"min": None, "mid": MZ_VALUE**2, "max": None},
            "sqrt_s": {"min": None, "mid": SQRT_S, "max": None},
        }
        kinematics.append(kin_value)

    return kinematics


def get_data_values(hepdata: dict) -> list:
    """Extract the central values from the HepData yaml file.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info

    Returns
    -------
    list:
        list of dictionaries whose contents are the central values

    """
    central = hepdata["dependent_variables"][0]["values"]

    return [NORM_FACTOR * v["value"] for v in central]


def get_errors(hepdata: dict, central: list) -> dict:
    """Extract the error values from the HepData yaml file.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
    central: list
        list containing the central values

    Returns
    -------
    list:
        list of dictionaries whose contents are the various
        source of uncertainties

    """
    errors = hepdata["dependent_variables"][0]["values"]

    stat, sys_uncorr, sys_corr, sys_lumi = [], [], [], []
    for idx, err in enumerate(errors):
        stat.append(NORM_FACTOR * err["errors"][0]["symerror"])
        sys_uncorr.append(NORM_FACTOR * err["errors"][1]["symerror"])

        # Corr. syst and FSR need to be added in Quadrature
        sys_cor = err["errors"][2]["symerror"]
        sys_fsr = err["errors"][3]["symerror"]
        sys_corr.append(NORM_FACTOR * np.sqrt(sys_cor**2 + sys_fsr**2))

        # Convert systematic Luminosity into absolute
        syslumi = percentage_to_absolute(
            err["errors"][4]["symerror"],  # [%]
            central[idx],
        )
        sys_lumi.append(NORM_FACTOR * syslumi)

    return {
        "stat": stat,
        "sys_uncorr": sys_uncorr,
        "sys_corr": sys_corr,
        "sys_lumi": sys_lumi,
    }


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
        error_value["sys_uncorr"] = uncs["sys_uncorr"][idat]
        error_value["sys_luminosity"] = uncs["sys_lumi"][idat]
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

    """
    error_definition = {
        f"sys_corr_{i + 1}": {
            "description": "Correlated systematic uncertainties",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(len(data))
    }

    error_definition["stat"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    error_definition["sys_uncorr"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    error_definition["sys_luminosity"] = {
        "description": "Systematic Luminosity uncertainties",
        "treatment": "ADD",
        "type": "LHCBLUMI10",
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


def main_filter():
    """Main driver of the LHCB_ZEE_8TEV filter that produces commmondata.

    There are three main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR

    2. Uncorrelated Systematic uncertainties: ADD, UNCORR

    3. Correlated Systematic uncertainties: ADD, CORR:
        constructed by symmetrizing the correlation matrix and extracting
        the artificial systematic uncertainties from the latter

    4. Luminosity Systematic uncertainties: ADD, LHCBLUMI8TEV

    """
    _, nbpoints, tables = read_metadata()

    comb_kins, comb_data = [], []
    combined_errors = []
    for tabid in tables:
        yaml_content = load_yaml(table_id=tabid)

        # Extract the kinematic, data, and uncertainties
        kinematics = get_kinematics(hepdata=yaml_content)
        data_central = get_data_values(hepdata=yaml_content)
        uncertainties = get_errors(yaml_content, data_central)

        # Collect all the results from different tables
        comb_kins += kinematics
        comb_data += data_central
        combined_errors.append(uncertainties)

    errors_combined = concatenate_dicts(combined_errors)
    # Compute the Artifical Systematics from CovMat
    corrmat = read_corrmatrix(nb_datapoints=nbpoints)
    covmat = multiply_syst(corrmat, errors_combined["sys_corr"])
    artunc = generate_artificial_unc(
        ndata=nbpoints,
        covmat_list=covmat.tolist(),
        no_of_norm_mat=0,
    )
    errors = format_uncertainties(errors_combined, artunc)

    # Generate all the necessary files
    dump_commondata(comb_kins, comb_data, errors)

    return


if __name__ == "__main__":
    main_filter()
