# Filter for E906R

import numpy as np
import yaml


def decompose_covmat(covmat):
    """Given a covmat it return an array sys with shape (ndat,ndat)
    giving ndat correlated systematics for each of the ndat point.
    The original covmat is obtained by doing sys@sys.T"""

    lamb, mat = np.linalg.eig(covmat)
    sys = np.multiply(np.sqrt(lamb), mat)
    return sys


def filter_E906R():

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    ndata = metadata["ndata"]

    Eb = 120  # beam energy
    mp = 0.938  # proton mass
    s = 2 * mp * mp + 2 * Eb * mp  # com energy square

    xt, xb, M, data, sys = np.loadtxt(
        "rawdata/data_paper.dat",
        skiprows=1,
        usecols=(2, 3, 4, 6, 8),
        unpack=True,
        dtype=float,
    )
    Y = 0.5 * np.log(xb / xt)
    s = np.sqrt(s)

    covmat = np.loadtxt("rawdata/cov.dat", usecols=range(6))
    art_sys = decompose_covmat(covmat)

    data_central = []
    kin = []
    error = []

    for n in range(ndata):

        data_central.append(float(data[n]))

        kin_value = {
            "Y": {"min": None, "mid": float(Y[n]), "max": None},
            "M2": {"min": None, "mid": float(M[n]), "max": None},
            "s": {"min": None, "mid": float(s), "max": None},
        }
        kin.append(kin_value)

        error_value = {
            "syst_1": float(sys[n]),
        }

        for m in range(ndata):
            error_value["syst_" + str(m + 2)] = float(art_sys[n, m])

        error.append(error_value)

    error_definition = {
        "sys_1": {
            "description": "systemtaic uncertainty",
            "treatment": "ADD",
            "type": "CORR",
        },
        "sys_2": {
            "description": "artificial systematic 1 ",
            "treatment": "ADD",
            "type": "CORR",
        },
        "sys_3": {
            "description": "artificial systematic 2 ",
            "treatment": "ADD",
            "type": "CORR",
        },
        "sys_4": {
            "description": "artificial systematic 3 ",
            "treatment": "ADD",
            "type": "CORR",
        },
        "sys_5": {
            "description": "artificial systematic 4 ",
            "treatment": "ADD",
            "type": "CORR",
        },
        "sys_6": {
            "description": "artificial systematic 5 ",
            "treatment": "ADD",
            "type": "CORR",
        },
        "sys_7": {
            "description": "artificial systematic 6 ",
            "treatment": "ADD",
            "type": "CORR",
        },
    }

    data_central_yaml = {"data_central": data_central}
    kinematics_yaml = {"bins": kin}
    uncertainties_yaml = {"definition": error_definition, "bins": error}

    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


filter_E906R()
