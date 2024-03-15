import yaml
import uproot
import numpy as np


def get_kinematics(version, figure):
    """
    returns the relevant kinematics values.

    Parameters
    ----------
    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    figure : str
        eg. 17a or 17b

    Returns
    -------
    list
        list containing the kinematic values for all
        hepdata tables
    """
    kin = []

    hepdata_table = f"rawdata/HEPData-ins1810913-v{version}-Figure_{figure}.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    for eta in input["independent_variables"][0]['values']:
        kin_value = {
            'eta': {'min': eta['low'], 'mid': 0.5 * (eta['low'] + eta['high']), 'max': eta['high']},
            'm_W2': {'min': None, 'mid': 6460.5, 'max': None},
            'sqrts': {'min': None, 'mid': 13000.0, 'max': None},
        }

        kin.append(kin_value)

    return kin


def get_data_values(version, figure):
    """
    returns the central data.

    Parameters
    ----------
    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    figure : str
        eg. 17a or 17b

    Returns
    -------
    list
        list containing the central values for all
        hepdata tables

    """

    data_central = []

    hepdata_table = f"rawdata/HEPData-ins1810913-v{version}-Figure_{figure}.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values = input['dependent_variables'][0]['values']

    for value in values:
        data_central.append(value['value'])

    return data_central


def decompose_covmat(covmat):
    """Given a covmat it return an array sys with shape (ndat,ndat)
    giving ndat correlated systematics for each of the ndat point.
    The original covmat is obtained by doing sys@sys.T"""

    lamb, mat = np.linalg.eig(covmat)
    sys = np.multiply(np.sqrt(lamb), mat)
    return sys


def get_systematics(observable, version, figure):
    """
    Following the CMS advice we take the covariance matrix from
    https://cms-results.web.cern.ch/cms-results/public-results/publications/SMP-18-012/index.html

    The root file sumpois 2dxsec contains the needed covariance matrix.

    Parameters
    ----------
    observable : str
        The observable for which the covariance matrix is needed
        can be either "W+" or "W-"

    """

    # Open the ROOT file
    file = uproot.open("rawdata/covariance_coefficients_sumpois_2dxsec.root")

    # access the TH2D histogram
    histogram = file["sub_cov_matrix"]

    # Get the labels along x-axis
    x_labels = histogram.axis(0).labels()

    # Select rows whose label starts with "POI" and contains the observable
    if observable == "W+":
        poi_indices = [
            i for i, label in enumerate(x_labels) if "POI, $W^{+}$, $|\\eta^{l}|" in label
        ]

    elif observable == "W-":
        poi_indices = [
            i for i, label in enumerate(x_labels) if "POI, $W^{-}$, $|\\eta^{l}|" in label
        ]

    # Extract the submatrix
    submatrix = histogram.values()[poi_indices][:, poi_indices]

    # Convert submatrix to numpy array
    submatrix_array = np.array(submatrix)

    # Get Luminosity covariance matrix
    if observable=="W+":
        with open("rawdata/HEPData-ins1810913-v1-Impacts_Figure_A23a.yaml", "r") as file:
            impacts = yaml.safe_load(file)
    elif observable=="W-":
        with open("rawdata/HEPData-ins1810913-v1-Impacts_Figure_A23b.yaml", "r") as file:
            impacts = yaml.safe_load(file)

        
    
    lumi_unc = np.array([val['value'] for val in impacts['dependent_variables'][2]['values']])

    hepdata_table = f"rawdata/HEPData-ins1810913-v{version}-Figure_{figure}.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values = input['dependent_variables'][0]['values']
    values = np.array([val['value'] for val in values])

    lumi_unc *= values / 100
    lumi_covmat = lumi_unc[:, np.newaxis] @ lumi_unc[:, np.newaxis].T

    artificial_uncertainties = np.real(decompose_covmat(lumi_covmat+submatrix_array))
    
    uncertainties = []

    for i, unc in enumerate(artificial_uncertainties.T):
        
        name = f"artificial_uncertainty_{i}"
        values = [unc[i] for i in range(len(unc))]
        uncertainties.append([{"name": name, "values": values}])

    return uncertainties


if __name__ == "__main__":
    get_systematics(observable="W+", version=1, figure='17a')
