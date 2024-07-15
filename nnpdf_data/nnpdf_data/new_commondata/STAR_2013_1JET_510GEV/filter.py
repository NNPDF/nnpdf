"""This script provides the common filer to the jet and dijet STAR 2013 datasets.
Files need to be parsed all together as there are correlations provided. 
"""
import pathlib

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.correlations import (
    compute_covmat,
    upper_triangular_to_symmetric,
)
from nnpdf_data.new_commondata.STAR_2012_1JET_510GEV.filter import TOPO_DEF

# values from the paper https://arxiv.org/pdf/2110.11020.pdf
SQRTS = 510
ETA_ABS = 0.9
POL_UNC = 0.064
LUMI_UNC = 0.00047
YEAR = 2013
TOPOPLOGY_LIST = ["I", "A", "B", "C", "D"]

HERE = pathlib.Path(__file__).parent
RAWDATA_PATH = HERE / "rawdata/"


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
    df["eta_min"] = -TOPO_DEF["I"]["abs_eta_max"]
    df["eta_max"] = +TOPO_DEF["I"]["abs_eta_max"]
    df["sqrts"] = SQRTS
    df["ALL"] = all_data[r"Inclusive Jet $A_{LL}$"]
    df["stat"] = all_data[r"stat +"]
    df["syst"] = all_data[r"syst +"]
    df["pol"] = POL_UNC * abs(df["ALL"])
    df["lumi"] = LUMI_UNC

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

    for p in ["1", "2"]:
        df[f"abs_eta{p}_min"] = TOPO_DEF[topology][f"abs_eta{p}_min"]
        df[f"abs_eta{p}_max"] = TOPO_DEF[topology][f"abs_eta{p}_max"]
        df[f"abs_eta{p}"] = (df[f"abs_eta{p}_min"] + df[f"abs_eta{p}_max"]) / 2

    df["sqrts"] = SQRTS
    df["ALL"] = all_data[r"Dijet $A_{LL}$, topology " + topology]
    df["stat"] = all_data[r"stat +"]
    df["syst"] = all_data[r"syst +"]
    df["pol"] = POL_UNC * abs(df["ALL"])
    df["lumi"] = LUMI_UNC

    print(f"2JET {topology} data loaded. Npoint: ", len(df))
    return df


def get_correlation_label(a):
    if a == "I":
        return "Inclusivejet"
    return f"Dijettopology{a}"


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
    if not np.allclose(tot_corr, np.triu(tot_corr)):
        raise ValueError("Correlation matrix not read correctly")
    return tot_corr + tot_corr.T - np.eye(np.sum((*ndata_dict.values(),)))


def write_1jet_data(df, art_sys):
    STORE_PATH = HERE

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
    error_definition = {
        "pol": {
            "description": "beam polarization uncertainty",
            "treatment": "MULT",
            "type": f"STAR{YEAR}POL",
        },
        "lumi": {
            "description": "luminosity uncertainty",
            "treatment": "ADD",
            "type": f"STAR{YEAR}LUMI",
        },
    }
    # loop on data points
    for i, sys_i in enumerate(art_sys):
        e = {"pol": float(df.loc[i, "pol"]), "lumi": float(df.loc[i, "lumi"])}
        # loop on art sys
        for j, val in enumerate(sys_i):
            e[f"sys_{j}"] = val

        if i == 0:
            error_definition.update(
                {
                    f"sys_{j}": {
                        "description": f"{j} artificial correlated statistical + systematics uncertainty",
                        "treatment": "ADD",
                        "type": f"STAR{YEAR}JETunc{j}",
                    }
                    for j in range(len(sys_i))
                }
            )
        error.append(e)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}
    with open(STORE_PATH / "uncertainties.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


def write_2jet_data(df, topology, art_sys):
    STORE_PATH = HERE / f"../STAR_{YEAR}_2JET_510GEV/"
    # Write central data
    data_central_yaml = {"data_central": list(df["ALL"])}
    with open(STORE_PATH / f"data_{topology}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(data_central_yaml, file)

    # Write kin file
    kin = []
    for i in range(len(df)):
        kin_value = {
            "m_jj": {
                "min": float(df.loc[i, "mjj_min"]),
                "mid": float(df.loc[i, "mjj"]),
                "max": float(df.loc[i, "mjj_max"]),
            },
            # "sqrts": {"min": None, "mid": float(df.loc[i, "sqrts"]), "max": None},
            "abs_eta_1": {
                "min": float(df.loc[i, "abs_eta1_min"]),
                "mid": float(df.loc[i, "abs_eta1"]),
                "max": float(df.loc[i, "abs_eta1_max"]),
            },
            "abs_eta_2": {
                "min": float(df.loc[i, "abs_eta2_min"]),
                "mid": float(df.loc[i, "abs_eta2"]),
                "max": float(df.loc[i, "abs_eta2_max"]),
            },
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}
    with open(
        STORE_PATH / f"kinematics_{topology}.yaml", "w", encoding="utf-8"
    ) as file:
        yaml.dump(kinematics_yaml, file)

    # Write unc file
    error = []
    error_definition = {
        "pol": {
            "description": "beam polarization uncertainty",
            "treatment": "MULT",
            "type": f"STAR{YEAR}POL",
        },
        "lumi": {
            "description": "luminosity uncertainty",
            "treatment": "ADD",
            "type": f"STAR{YEAR}LUMI",
        },
    }
    # loop on data points
    for i, sys_i in enumerate(art_sys):
        e = {
            "pol": float(df.loc[i, "pol"]),
            "lumi": float(df.loc[i, "lumi"]),
        }
        # loop on art sys
        for j, val in enumerate(sys_i):
            e[f"sys_{j}"] = val
        error.append(e)

        if i == 0:
            error_definition.update(
                {
                    f"sys_{j}": {
                        "description": f"{j} artificial correlated statistical + systematics uncertainty",
                        "treatment": "ADD",
                        "type": f"STAR{YEAR}JETunc{j}",
                    }
                    for j in range(len(sys_i))
                }
            )

    uncertainties_yaml = {"definitions": error_definition, "bins": error}
    with open(
        STORE_PATH / f"uncertainties_{topology}.yaml", "w", encoding="utf-8"
    ) as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    # load all the data
    dfs = {"I": read_1jet_data()}
    for topo in TOPOPLOGY_LIST[1:]:
        dfs[topo] = read_2jet_data(topo)

    # load correlations
    ndata_dict = {a: len(b) for a, b in dfs.items()}
    correlation_df = read_correlations(ndata_dict)
    # sum stat and syst for both jet and dijets as recommended
    # by E.Aschenauer, see https://github.com/NNPDF/nnpdf/pull/2035#issuecomment-2201979662
    correlated_unc = []
    for a in TOPOPLOGY_LIST:
        correlated_unc.extend(
            np.sqrt(dfs[a]["syst"] ** 2 + dfs[a]["stat"] ** 2).values.tolist()
        )
    ndata_points = np.sum((*ndata_dict.values(),))
    # decompose uncertainties
    art_sys = np.array(compute_covmat(correlation_df, correlated_unc, ndata_points))

    # write data
    cnt = 0
    for topo, df in dfs.items():
        ndata = ndata_dict[topo]
        syst = art_sys[cnt : cnt + ndata, :].tolist()
        if topo == "I":
            write_1jet_data(df, syst)
        else:
            write_2jet_data(df, topo, syst)
        cnt += ndata
