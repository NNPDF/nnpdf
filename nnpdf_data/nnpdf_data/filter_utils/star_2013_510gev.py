"""This script provides the common filer to the jet and dijet STAR 2013 datasets.
Files need to be parsed all together as there are correlations provided. 
"""
import pandas as pd
import numpy as np
import pathlib
import yaml

from correlations import compute_covmat

# values from the paper https://arxiv.org/pdf/2110.11020.pdf
SQRTS = 510
ETA_ABS = 0.9
TOPOPLOGY_LIST = ["1JET", "A", "B", "C", "D"]

HERE = pathlib.Path(__file__).parents[1]
RAWDATA_PATH = HERE / "new_commondata/STAR_2013_1JET_510GEV/rawdata/"


def read_1jet_data():
    data_table = pathlib.Path(RAWDATA_PATH / "Figure3.csv")

    with open(data_table, "r", encoding="utf-8") as file:
        parton_jet_data = pd.read_csv(
            file, delimiter=",", skiprows=lambda x: (x <= 21 or x >= 38)
        )
    with open(data_table, "r", encoding="utf-8") as file:
        all_data = pd.read_csv(file, delimiter=",", skiprows=37)

    df = pd.DataFrame()
    df["pT"] = parton_jet_data[r"Parton Jet $p_{T}$ (GeV/$c$)"]
    df["pT_min"] = (
        parton_jet_data[r"Parton Jet $p_{T}$ (GeV/$c$)"] + parton_jet_data["syst -"]
    )
    df["pT_max"] = (
        parton_jet_data[r"Parton Jet $p_{T}$ (GeV/$c$)"] + parton_jet_data["syst +"]
    )
    df["eta"] = 0.0
    df["eta_min"] = -ETA_ABS
    df["eta_max"] = +ETA_ABS
    df["sqrts"] = SQRTS
    df["ALL"] = all_data[r"Inclusive Jet $A_{LL}$"]
    df["stat"] = all_data[r"stat +"]
    df["syst"] = all_data[r"syst +"]

    print("1JET data loaded. Npoint: ", len(df))
    return df


def read_2jet_data(topology):
    data_table = RAWDATA_PATH / f"Figure5topology{topology}.csv"
    with open(data_table, "r", encoding="utf-8") as file:
        mjj_data = pd.read_csv(
            file, delimiter=",", skiprows=lambda x: (x <= 5 or x >= 20)
        )
    with open(data_table, "r", encoding="utf-8") as file:
        all_data = pd.read_csv(file, delimiter=",", skiprows=20)

    df = pd.DataFrame()
    df["mjj"] = mjj_data[r"Parton Dijet $M_{inv}$ (GeV/$c^{2}$)"]
    df["mjj_min"] = (
        mjj_data[r"Parton Dijet $M_{inv}$ (GeV/$c^{2}$)"] + mjj_data["syst -"]
    )
    df["mjj_max"] = (
        mjj_data[r"Parton Dijet $M_{inv}$ (GeV/$c^{2}$)"] + mjj_data["syst +"]
    )
    df["sqrts"] = SQRTS
    df["ALL"] = all_data[r"Dijet $A_{LL}$, topology " + topology]
    df["stat"] = all_data[r"stat +"]
    df["syst"] = all_data[r"syst +"]

    print(f"2JET {topology} data loaded. Npoint: ", len(df))
    return df


def get_correlation_label(a):
    if a == "1JET":
        return "Inclusivejet"
    return f"Dijettopology{a}"


def upper_triangular_to_symmetric(ut, dim):
    """Build a symmetric matrix from the upper diagonal"""
    corr = np.zeros((dim, dim))
    last = dim
    first = 0
    for i in range(dim):
        corr[i, i:] = ut[first:last]
        last += dim - i - 1
        first += dim - i
    return corr + corr.T - np.eye(dim)


def read_correlations(ndata_dict):
    """Read the correlation files and build a big matix"""
    corr_rows = []
    # loop on block rows
    for a, ndata_a in ndata_dict.items():
        label_a = get_correlation_label(a)
        la = [a for _ in range(ndata_a)]
        corr_row = pd.DataFrame()
        # loop on block columns
        for b, ndata_b in ndata_dict.items():
            label_b = get_correlation_label(b)
            lb = [b for _ in range(ndata_b)]

            # build the block
            try:
                with open(
                    RAWDATA_PATH / f"{label_a}-{label_b}correlation.csv",
                    encoding="utf-8",
                ) as file:
                    corr_df = pd.read_csv(file, delimiter=",", skiprows=6)
                if a == b:
                    corr = upper_triangular_to_symmetric(corr_df.values[:, 2], ndata_a)
                else:
                    corr = corr_df.values[:, 2].reshape((ndata_a, ndata_b))
            except FileNotFoundError:
                corr = pd.DataFrame(np.zeros((ndata_a, ndata_b)), index=la, columns=lb)

            corr = pd.DataFrame(corr, index=la, columns=lb)
            corr_row = pd.concat([corr_row, corr], axis=1)
        corr_rows.append(corr_row)

    tot_corr = pd.concat(corr_rows)
    return tot_corr + tot_corr.T - np.eye(np.sum((*ndata_dict.values(),)))


def write_1jet_data(df, art_sys, ndata):
    STORE_PATH = HERE / "new_commondata/STAR_2013_1JET_510GEV/"

    # Write central data
    data_central_yaml = {"data_central": list(df["ALL"])}
    with open(STORE_PATH / "data.yaml", "w", encoding="utf-8") as file:
        yaml.dump(data_central_yaml, file)

    # Write kin file
    kin = []
    for i in range(len(df)):
        kin_value = {
            "pT": {
                "min": float(df.loc[i, "pT_min"]),
                "mid": float(df.loc[i, "pT"]),
                "max": float(df.loc[i, "pT_max"]),
            },
            "sqrts": {"min": None, "mid": float(df.loc[i, "sqrts"]), "max": None},
            "eta": {
                "min": float(df.loc[i, "eta_min"]),
                "mid": float(df.loc[i, "eta"]),
                "max": float(df.loc[i, "eta_max"]),
            },
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}
    with open(STORE_PATH / "kinematics.yaml", "w", encoding="utf-8") as file:
        yaml.dump(kinematics_yaml, file)

    # Write unc file
    error = []
    for i in range(ndata):
        e = {}
        # add the art sys
        for j in range(ndata):
            e[f"sys_{j}"] = art_sys[i][j]
        e[
            "stat"
        ] = 0  # This is set to 0 as the stat unc is correlated and reported in sys_0
        error.append(e)

    error_definition = {
        f"sys_{i}": {
            "description": f"{i} artificial correlated statistical + systematics uncertainty",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(ndata)
    }
    uncertainties_yaml = {"definitions": error_definition, "bins": error}
    with open(STORE_PATH / "uncertainties.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


def write_2jet_data(df, topology, art_sys, ndata):
    STORE_PATH = HERE / f"new_commondata/STAR_2013_2JET_{topology}_510GEV/"
    # Write central data
    data_central_yaml = {"data_central": list(df["ALL"])}
    with open(STORE_PATH / "data.yaml", "w", encoding="utf-8") as file:
        yaml.dump(data_central_yaml, file)

    # Write kin file
    kin = []
    for i in range(len(df)):
        kin_value = {
            "mjj": {
                "min": float(df.loc[i, "mjj_min"]),
                "mid": float(df.loc[i, "mjj"]),
                "max": float(df.loc[i, "mjj_max"]),
            },
            "sqrts": {"min": None, "mid": float(df.loc[i, "sqrts"]), "max": None},
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}
    with open(STORE_PATH / "kinematics.yaml", "w", encoding="utf-8") as file:
        yaml.dump(kinematics_yaml, file)

    # Write unc file
    error = []
    for i in range(ndata):
        e = {}
        # add the art sys
        for j in range(ndata):
            e[f"sys_{j}"] = art_sys[i][j]
        e["stat"] = float(df.loc[i, "stat"])
        error.append(e)

    error_definition = {
        f"sys_{i}": {
            "description": f"{i} artificial correlated systematics uncertainty",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(ndata)
    }
    error_definition.update(
        {
            "stat": {
                "description": "statistical uncertainty",
                "treatment": "ADD",
                "type": "UNCORR",
            }
        }
    )
    uncertainties_yaml = {"definitions": error_definition, "bins": error}
    with open(STORE_PATH / "uncertainties.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    # load all the data
    dfs = {"1JET": read_1jet_data()}
    for topo in TOPOPLOGY_LIST[1:]:
        dfs[topo] = read_2jet_data(topo)

    # load correlations
    ndata_dict = {a: len(b) for a, b in dfs.items()}
    correlation_df = read_correlations(ndata_dict)
    # from the paper we understand that stat dijet is not correlated
    #    I-I (stat + sys) | I-D (stat + sys)
    #    D-I (stat + sys) | D-D (sys)
    correlated_unc = np.sqrt(
        dfs["1JET"]["syst"] ** 2 + dfs["1JET"]["stat"] ** 2
    ).values.tolist()
    for a in TOPOPLOGY_LIST[1:]:
        correlated_unc.extend(dfs[a]["syst"].values)
    ndata_points = np.sum((*ndata_dict.values(),))
    # decompose uncertainties
    # TODO: how come this is not yet block diagonal??
    art_sys = np.array(compute_covmat(correlation_df, correlated_unc, ndata_points))

    # write data
    cnt = 0
    for topo, df in dfs.items():
        ndata = ndata_dict[topo]
        syst = art_sys[cnt:cnt+ndata, cnt:cnt+ndata].tolist()
        if topo == "1JET":
            write_1jet_data(df, syst, ndata)
        else:
            write_2jet_data(df, topo, syst, ndata)
        cnt += ndata
