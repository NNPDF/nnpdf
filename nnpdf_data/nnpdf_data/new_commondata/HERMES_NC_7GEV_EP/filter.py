import pandas as pd
import yaml
import glob
import numpy as np
import pathlib
import pandas as pd
from numpy.linalg import eig


def read_data(fnames):
    df = pd.DataFrame()
    for fname in fnames:
        with open(fname, "r") as file:
            data = yaml.safe_load(file)

        xsub = data["independent_variables"][0]["values"]
        Qsub = data["independent_variables"][1]["values"]
        Gsub = data["dependent_variables"][0]["values"]

        for i in range(len(xsub)):
            df = pd.concat(
                [
                    df,
                    pd.DataFrame(
                        {
                            "x": [xsub[i]["value"]],
                            "Q2": [Qsub[i]["value"]],
                            "G": [Gsub[i]["value"]],
                            "stat": [Gsub[i]["errors"][0]["symerror"]],
                            "exp": [Gsub[i]["errors"][1]["symerror"]],
                            "param": [Gsub[i]["errors"][2]["symerror"]],
                            "evol": [Gsub[i]["errors"][3]["symerror"]],
                        }
                    ),
                ],
                ignore_index=True,
            )

    return df


def read_corrmatrix(nb_datapoints: int = 15) -> np.ndarray:
    """Load the correlation Matrix in Table 22."""
    file = pathlib.Path("./rawdata/HEPData-ins726689-v1-Table_22.yaml")
    loaded_file = yaml.safe_load(file.read_text())

    corrs = loaded_file['dependent_variables'][0]['values']
    df_corrs = pd.DataFrame(corrs)

    return df_corrs.value.values.reshape((nb_datapoints, nb_datapoints))


def covmat_to_artunc(ndata, covmat_list, no_of_norm_mat=0):
    r"""Convert the covariance matrix to a matrix of
    artificial uncertainties.

    NOTE: This function has been taken from validphys.newcommondata_utils.
    If those utils get merged in the future, we can replace this.

    Parameters
    ----------
    ndata : integer
        Number of data points
    covmat_list : list
        A one dimensional list which contains the elements of
        the covariance matrix row by row. Since experimental
        datasets provide these matrices in a list form, this
        simplifies the implementation for the user.
    no_of_norm_mat : int
        Normalized covariance matrices may have an eigenvalue
        of 0 due to the last data point not being linearly
        independent. To allow for this, the user should input
        the number of normalized matrices that are being treated
        in an instance. For example, if a single covariance matrix
        of a normalized distribution is being processed, the input
        would be 1. If a covariance matrix contains pertains to
        3 normalized datasets (i.e. cross covmat for 3
        distributions), the input would be 3. The default value is
        0 for when the covariance matrix pertains to an absolute
        distribution.

    Returns
    -------
    artunc : list
        A two dimensional matrix (given as a list of lists)
        which contains artificial uncertainties to be added
        to the commondata. i^th row (or list) contains the
        artificial uncertainties of the i^th data point.

    """
    epsilon = -0.0000000001
    neg_eval_count = 0
    psd_check = True
    covmat = np.zeros((ndata, ndata))
    artunc = np.zeros((ndata, ndata))
    for i in range(len(covmat_list)):
        a = i // ndata
        b = i % ndata
        covmat[a][b] = covmat_list[i]
    eigval, eigvec = eig(covmat)
    for j in range(len(eigval)):
        if eigval[j] < epsilon:
            psd_check = False
        elif eigval[j] > epsilon and eigval[j] <= 0:
            neg_eval_count = neg_eval_count + 1
            if neg_eval_count == (no_of_norm_mat + 1):
                psd_check = False
        elif eigval[j] > 0:
            continue
    if psd_check == False:
        raise ValueError('The covariance matrix is not positive-semidefinite')
    else:
        for i in range(ndata):
            for j in range(ndata):
                if eigval[j] < 0:
                    continue
                else:
                    artunc[i][j] = eigvec[i][j] * np.sqrt(eigval[j])
    return artunc.tolist()


def compute_covmat(corrmat: np.ndarray, df: pd.DataFrame, ndata: int) -> list:
    """Compute the covariance matrix with the artificial stat uncertanties."""
    # multiply by stat err
    stat = df["stat"]
    cov_mat = np.einsum("i,ij,j->ij", stat, corrmat, stat)
    return covmat_to_artunc(ndata, cov_mat.flatten().tolist())


def write_data(df):
    data_central = []
    for i in range(len(df["G"])):
        data_central.append(float(df.loc[i, "G"]))

    data_central_yaml = {"data_central": data_central}
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    # Write kin file
    kin = []
    for i in range(len(df["G"])):
        kin_value = {
            "x": {"min": None, "mid": float(df.loc[i, "x"]), "max": None},
            "Q2": {"min": None, "mid": float(df.loc[i, "Q2"]), "max": None},
        }
        kin.append(kin_value)

    kinematics_yaml = {"bins": kin}

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # Extract the correlation matrix and compute artificial systematics
    ndata_points = len(data_central)
    corrmatrix = read_corrmatrix(nb_datapoints=ndata_points)
    # Compute the covariance matrix
    compute_covmat(corrmatrix, df, ndata_points)

    # Compute the covariance matrix
    art_sys = compute_covmat(corrmatrix, df, ndata_points)

    error = []
    for i in range(ndata_points):
        e = {}
        # add the art sys
        for j in range(ndata_points):
            e[f"sys_{j}"] = art_sys[i][j]

        e[
            "stat"
        ] = 0  # This is set to 0 as the stat unc is correlated and reported in sys_0
        e["exp"] = float(
            df.loc[i, "exp"]
        )  # experimental including normalization
        e["param"] = float(df.loc[i, "param"])
        e["evol"] = float(df.loc[i, "evol"])
        error.append(e)

    error_definition = {
        f"sys_{i}": {
            "description": f"{i} artificial correlated statistical uncertainty",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(ndata_points)
    }

    error_definition.update(
        {
            "stat": {
                "description": "statistical uncertainty",
                "treatment": "ADD",
                "type": "UNCORR",
            },
            "exp": {
                "description": "experimental systematic uncertainty",
                "treatment": "ADD",
                "type": "UNCORR",
            },
            "param": {
                "description": "parametrization systematic uncertainty",
                "treatment": "ADD",
                "type": "CORR",
            },
            "evol": {
                "description": "evolution systematic uncertainty",
                "treatment": "ADD",
                "type": "CORR",
            },
        }
    )

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    fnames = glob.glob("./rawdata/HEPData-ins726689-v1-Table_13.yaml")
    nnames = sorted([i for i in fnames])
    df = read_data(nnames)
    write_data(df)
