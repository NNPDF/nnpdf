import pathlib

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import covmat_to_artunc, prettify_float

yaml.add_representer(float, prettify_float)

MW_VALUE = 80.398  # GeV
SQRT_S = 5_020.0  # GeV


def load_yaml(W_sign: str, version: int = 1) -> dict:
    """Load the HEP data table in yaml format.

    Parameters
    ----------
    W_sign: sign of cc dy
    version: hepdata version

    Returns
    -------
    dict:
        ditionary containing the table contents

    """
    filename = f"HEPData-ins2972386-v{version}-W^{W_sign}_to_mu^{W_sign}_nu_{{mu}}_dsigma_dp_T"
    table = pathlib.Path(f"rawdata/{filename}.yaml")

    return yaml.safe_load(table.read_text())

def read_metadata() -> tuple[int, int, list]:
    """Read the version, number of datapoints, table identificators
    from metadata.

    Returns
    -------
    tuple(int, int, list):
        data version, number of datapoints (per table), table names

    """
    metadata = pathlib.Path("metadata.yaml")
    content = yaml.safe_load(metadata.read_text())

    version = content["hepdata"]["version"]
    nb_datapoints = content["implemented_observables"][0]["ndata"]
    tables = []
    for i in range(2):
        tables.append(content['implemented_observables'][i]["tables"][0])
    # tables = content["implemented_observables"][0]["tables"]

    return version, nb_datapoints, tables

def get_kinematics(hepdata: dict, bin_index: list) -> list:
    """Read the kinematics (bins) for a chosen charged current sign from hepdata
    W_sign determined implicitly through load_yaml

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info. We get it from load_yaml
    bin_index: list
        list of Non-empty bin index. We get its length from read_metadata

    Returns
    -------
    list of kinematics bins, i.e. list of dics of form [{p_T: {min: , max: , mid: }, m_W2: {...}]

    """
    ptbins = hepdata["independent_variables"][0]["values"]

    kinematics = []
    for bins in bin_index:
        ptmin = float(ptbins[bins]["low"])
        ptmax = float(ptbins[bins]["high"])
        kin_value = {
            "p_T": {"min": ptmin, "mid": 0.5 * (ptmin + ptmax), "max": ptmax},
            "m_W2": {"min": None, "mid": MW_VALUE ** 2, "max": None},
            "sqrts": {"min": None, "mid": SQRT_S, "max": None},
        }
        kinematics.append(kin_value)

    return kinematics

def get_data_values(hepdata: dict, bin_index: list, indx: int = 0) -> list:
    """Extract the central values from the HepData yaml file. W_sign determined
    implicitly in load_yaml

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

    return [central[i]["value"] for i in bin_index]

def get_errors(hepdata: dict, bin_index: list, indx: int = 0) -> dict:
    """Extract the error values from the HepData yaml file.
    W_sign determined implicitly via load_yaml

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info: get it from load_yaml
    bin_index: list
        list of Non-empty bin index: for indexing datapoints [0, ..., ndata-1]
    indx: int
        index from which the errors will be read

    Returns
    -------
    list:
        list of dictionaries whose contents are the various
        source of uncertainties

    """
    errors = hepdata["dependent_variables"][indx]["values"]

    stat, sys = [], []
    for idx in bin_index:
        stat.append(errors[idx]["errors"][0]["symerror"]) # 0 for stat
        sys.append(errors[idx]["errors"][1]["symerror"]) # 1 for syst

    return {"stat": stat, "sys": sys}

def read_stat_corrmatrix(ndata: int, bin_index: list, W_sign: str = "+", version: int = 1) -> np.ndarray:
    """Reads the statistical correlation matrix for a chosen cc dy sign
    no need to symmetrise as it's already symmetric

    Parameters
    ----------
    nb_datapoints: int
        total number of datapoints

    Returns
    -------
    np.ndarray: (nb_datapoints, nb_datapoints)
        a symmetrized square matrix

    """

    filename = f"HEPData-ins2972386-v{version}-Statistical_correlation_matrix_for_W^{W_sign}"
    table = pathlib.Path(f"rawdata/{filename}.yaml")
    corr_vals = yaml.safe_load(table.read_text())['dependent_variables'][0]['values']
    corrmat = np.array([[float(0)]*ndata]*ndata)
    for i in bin_index:
        for j in bin_index:
            corrmat[i,j] = corr_vals[i*ndata+j]['value']

    return corrmat

def read_sys_corrmatrix(ndata: int, bin_index: list, version: int = 1) -> np.ndarray:
    """Reads the total systematic correlation matrix
    no need to symmetrise as it's already symmetric

    Parameters
    ----------
    nb_datapoints: int
        total number of datapoints (24, for both + and -)
        bin_index: [0, 1, ..., 23]

    Returns
    -------
    np.ndarray: (nb_datapoints, nb_datapoints)
        a symmetrized square matrix

    """
    filename = f"HEPData-ins2972386-v{version}-Systematic_correlation_matrix_for_W^{{pm}}_(combined)"
    table = pathlib.Path(f"rawdata/{filename}.yaml")
    corr_vals = yaml.safe_load(table.read_text())['dependent_variables'][0]['values']
    corrmat = np.array([[float(0)]*ndata]*ndata)
    for i in bin_index:
        for j in bin_index:
            corrmat[i,j] = corr_vals[i*ndata+j]['value']

    return corrmat

def corr_to_cov(corrmat: np.ndarray, unc: list) -> np.ndarray:
    """Multiply the corrmat to get covmat.

    Parameters
    ----------
    corrmat: np.ndarray
        a squared matrix representing the corrmat
    unc: list
        a vector containing the correlated stat uncertainties

    Returns
    -------
    np.ndarray:
        covariance matrix multiplied by the total systematics

    """
    covmat = np.zeros(corrmat.shape)

    for i in range(corrmat.shape[0]):
        for j in range(corrmat.shape[-1]):
            covmat[i][j] = corrmat[i][j] * unc[i] * unc[j]

    return covmat

def generate_artificial_unc(**kwargs) -> np.ndarray:
    """A wrapper around `covmat_to_artunc` to return an array.

    Returns
    -------
    np.ndarray:
        a squared matrix with artificial uncertainties
    """
    artunc = covmat_to_artunc(**kwargs)
    return np.array(artunc)

def generate_lum_uncertainties(centrals: list, rel_unc: float = 0.02) -> list:
    '''
    generates luminosity uncertainty on central value for each bin. As per paper, 
    the luminosity uncertainty is 2%

    Parameters
    ----------
    centrals: list
        a list of central values in each bin
    rel_unc: float
        relative luminosity uncertainty, as above, specified as 0.02

    Returns
    -------
    list:
        list of luminosity uncertainties

    '''
    lum_errors = []
    for value in centrals:
        lum_errors.append(value*rel_unc)
    return lum_errors

def aggregate_uncertainties(stat_art_unc: dict, sys_art_unc: np.ndarray, lum_unc: dict) -> dict:
    '''
    take statistical, systematic artificial uncertainties and luminosity uncertainties and aggregate them into a yaml-able lists

    Parameters
    ----------
    stat_art_unc: dict
        dictionary with two keys, + and -
        each value of the dictionary is a 12x12 array of artificial uncertainties
    sys_art_unc: np.ndarray
        a 24x24 array of systematic artificial uncertainties 
    lum_unc: dict
        a dictionary with two keys, + and -
        each value of the dictionary is a list of luminosity uncertainties for each bin
    Returns
    -------
    dict:
        dictionary with two keys, + and -
        each value of a dictionary is a list with each entry corresponding to one bin.
        each entry is a dictionary listing all uncertainties for a given bin
    '''
    errors = {}
    for W_sign in ['+', '-']:
        errors[W_sign]=[]
        for row in stat_art_unc[W_sign]:
            errors[W_sign].append({f"stat_corr_{i + 1}": float(row[i]) for i in range(len(row))})
    
    
    for i, row in enumerate(sys_art_unc):
        if i < 12:
            errors["+"][i] = errors["+"][i] | {f"sys_corr_{j+1}": float(row[j]) for j in range(len(row))}
        else:
            errors["-"][i-12] = errors["-"][i-12] | {f"sys_corr_{j+1}": float(row[j]) for j in range(len(row))}

    for W_sign in lum_unc.keys():
        for i, val in enumerate(lum_unc[W_sign]):
            errors[W_sign][i]['luminosity'] = val

    return errors

def dump_commondata(kinematics: dict, data: dict, errors: dict, nbpoints: int) -> None:
    """Function that generates and writes the commondata files.

    Parameters
    ----------
    kinematics: dict
        dictionary with two keys, + and -
        each value is a list of bin kinematics, obtained via get_kinematics
    data: dict
        dictionary with two keys, + and -
        each value is a list of bin central values, obtained via get_data_values
    errors: dict
        dictionary with two keys, + and -
        contains uncertainties for each bin in a form of dictionary, obtained via aggregate_uncertainties


    """
    error_definition = {
        f"stat_corr_{i + 1}": {
            "description": "Correlated statistical uncertainties",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(nbpoints)
    }

    error_definition.update({
        f"sys_corr_{i+1}": {
            "description": "Correlated systematic uncertainties",
            "treatment": "ADD",
            "type": "CORR"
        }
        for i in range(nbpoints*2)
    })

    error_definition["luminosity"] = {
        "description": "Luminosity uncertainty",
        "treatment": "MULT",
        "type": "LHCBLUMI5P02TEV"
    }

    errors_p = {"definitions": error_definition, "bins": errors["+"]}

    errors_m = {"definitions": error_definition, "bins": errors["-"]}

    with open("data_wm.yaml", "w") as file:
        yaml.dump({"data_central": data['-']}, file, sort_keys=False)

    with open("data_wp.yaml", "w") as file:
        yaml.dump({"data_central": data['+']}, file, sort_keys=False)

    with open("kinematics_wp.yaml", "w") as file:
        yaml.dump({"bins": kinematics['+']}, file, sort_keys=False)

    with open("kinematics_wm.yaml", "w") as file:
        yaml.dump({"bins": kinematics['-']}, file, sort_keys=False)

    with open("uncertainties_wp.yaml", "w") as file:
        yaml.dump(errors_p, file, sort_keys=False)

    with open("uncertainties_wm.yaml", "w") as file:
        yaml.dump(errors_m, file, sort_keys=False)

def main_filter() -> None:
    """Main driver of the filter that produces commmondata.

    There are three main different sources of uncertainties.

    1. Correlated Statistical uncertainties: ADD, CORR
        constructed by extracting the artificial statistical uncertainties
        from the covariance matrix

    2. Correlated Systematic uncertainties: ADD, CORR:
        constructed by extracting the artificial statistical uncertainties
        from the covariance matrix

    3. Luminosity Systematic uncertainties: MULT, LHCBLUMI5P02TEV
        as stated in the paper, 2% of each central value
    """
    version, nb_data, tables = read_metadata()

    kinematics = {} # bins of p_T, m_W2, sqrts, read from data
    data_central = {} # cetral values of diff XS, read from data
    uncertainties = {} # sys and stat uncertainties for each bin
    lum_unc = {} # luminosity uncertainty, calculated as 2% of each central value
    stat_art_unc = {} #statistical artificial uncertainties, obtained from correlation matrix
    for W_sign in tables:  # Loop over tables
            bin_index = list(range(nb_data))
            hepdata = load_yaml(W_sign)

            # Extract the kinematics, data, and uncertainties
            kinematics[W_sign] = get_kinematics(hepdata, bin_index)
            data_central[W_sign] = get_data_values(hepdata, bin_index)
            uncertainties[W_sign] = get_errors(hepdata, bin_index)

            # get luminosity uncertainties
            lum_unc[W_sign] = generate_lum_uncertainties(data_central[W_sign])
            
            # get statistical artificial uncertainties
            stat_corrmat = read_stat_corrmatrix(nb_data, bin_index, W_sign)
            stat_covmat_list = corr_to_cov(stat_corrmat, uncertainties[W_sign]['stat']).flatten()
            stat_art_unc[W_sign] = generate_artificial_unc(ndata = nb_data, covmat_list = stat_covmat_list)

    # get systematic artificial uncertainties
    sys_corrmat = read_sys_corrmatrix(ndata = 2*nb_data, bin_index = list(range(2*nb_data)))
    sys_uncertainties_all = uncertainties['+']['sys'] + uncertainties['-']['sys']
    sys_covmat_list = corr_to_cov(corrmat = sys_corrmat, unc = sys_uncertainties_all).flatten()
    sys_art_unc = generate_artificial_unc(ndata = 2*nb_data, covmat_list = sys_covmat_list)

    # aggregate all uncertainties
    errors = aggregate_uncertainties(stat_art_unc, sys_art_unc, lum_unc)

    # create commondata files
    dump_commondata(kinematics = kinematics, data = data_central, errors = errors, nbpoints = nb_data)

if __name__ == "__main__":
    main_filter()
