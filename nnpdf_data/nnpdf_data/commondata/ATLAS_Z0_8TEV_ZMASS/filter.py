from filter_utils import get_data_values, get_kinematics, get_systematics
import numpy as np
import yaml


def filter_ATLAS_Z0_8TEV_data_kinetic():
    """
    writes data central values and kinematics
    to respective .yaml file
    """
    kin = get_kinematics()
    central_values = get_data_values()

    data_central_yaml = {"data_central": central_values}
    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_ATLAS_Z0_8TEV_uncertainties():
    """
    writes uncertainties to respective .yaml file
    """
    systematics = get_systematics()

    # load correlation matrix from .txt file
    corr_matrix = np.loadtxt("rawdata/zy.txt")

    # generate covariance matrix from correlation matrix
    tot_systematics = np.array([syst[0]['value'] for syst in systematics['tot']])

    # TODO: this should be done with utils.correlation_to_covariance once that is merged in master
    cov_matrix_no_lumi = np.outer(tot_systematics, tot_systematics) * corr_matrix

    # add lumi uncertainty
    lumi_unc = np.array([syst[0]['value'] for syst in systematics['lumi']])
    lumi_cov = lumi_unc[:, None] @ lumi_unc[:, None].T

    # add covariances
    cov_matrix = cov_matrix_no_lumi + lumi_cov

    # compute decomposition of covariance matrix so as to get artificial systematics
    # TODO: use utils once merged in master
    lamb, mat = np.linalg.eig(cov_matrix)
    art_sys = np.multiply(np.sqrt(lamb), mat)

    uncertainties = []

    for i, unc in enumerate(art_sys.T):

        name = f"artificial_uncertainty_{i+1}"
        values = [unc[i] for i in range(len(unc))]
        uncertainties.append([{"name": name, "values": values}])

    # error definition
    error_definitions = {}
    errors = []

    for sys in uncertainties:

        error_definitions[sys[0]['name']] = {
            "description": f"{sys[0]['name']}",
            "treatment": "ADD",
            "type": "CORR",
        }

    for i in range(cov_matrix.shape[0]):
        error_value = {}

        for sys in uncertainties:
            error_value[sys[0]['name']] = float(sys[0]['values'][i])

        errors.append(error_value)

    uncertainties_yaml = {"definitions": error_definitions, "bins": errors}

    # write uncertainties
    with open(f"uncertainties.yaml", 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_Z0_8TEV_data_kinetic()
    filter_ATLAS_Z0_8TEV_uncertainties()
