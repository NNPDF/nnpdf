"""This script provides the common filer to the jet and dijet STAR 2012 datasets.
Files need to be parsed all together as there are correlations provided. 
"""
import pathlib
import math

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.correlations import (
    compute_covmat,
    upper_triangular_to_symmetric,
)

# values from the paper https://arxiv.org/abs/1906.02740
SQRTS = 510
YEAR = 2012
ETA_ABS = 0.9
POL_UNC = 0.066
TOPOPLOGY_LIST = ["I", "A", "B", "C", "D"]

HERE = pathlib.Path(__file__).parent
RAWDATA_PATH = HERE / "rawdata/"


def read_1jet_data():
    data_table = pathlib.Path(RAWDATA_PATH / "Figure12(Run12).csv")

    with open(data_table, "r", encoding="utf-8") as file:
        all_data = pd.read_csv(file, delimiter=",", skiprows=6)

    df = pd.DataFrame()
    df["pT"] = all_data["Parton Jet $p_{T}$ [GeV/c]"]
    df["pT_min"] = all_data["Parton Jet $p_{T}$ [GeV/c] LOW"]
    df["pT_max"] = all_data["Parton Jet $p_{T}$ [GeV/c] HIGH"]
    df["eta"] = 0.0
    df["eta_min"] = -ETA_ABS
    df["eta_max"] = +ETA_ABS
    df["sqrts"] = SQRTS
    df["ALL"] = all_data["Inclusive Jet $A_{LL}$"]
    df["stat"] = all_data["stat +"]
    df["syst"] = all_data["sys +"]
    df["pol"] = POL_UNC * abs(df["ALL"])
    df["lumi_ue"] = all_data["UE/RL sys +"]

    print("1JET data loaded. Npoint: ", len(df))
    return df


def read_2jet_data(topology):
    data_table = RAWDATA_PATH / f"Figure14(b,{topology},data).csv"

    with open(data_table, "r", encoding="utf-8") as file:
        all_data = pd.read_csv(file, delimiter=",", skiprows=6)

    df = pd.DataFrame()
    df["mjj"] = all_data["Parton DiJet $M_{inv}$ [$GeV/c^{2}$]"]
    df["mjj_min"] = all_data["Parton DiJet $M_{inv}$ [$GeV/c^{2}$] LOW"]
    df["mjj_max"] = all_data["Parton DiJet $M_{inv}$ [$GeV/c^{2}$] HIGH"]
    df["sqrts"] = SQRTS
    df["ALL"] = all_data[r"DiJet $A_{LL}$"]
    df["stat"] = all_data[r"stat +"]
    df["syst"] = all_data[r"sys +"]
    df["pol"] = POL_UNC * abs(df["ALL"])
    df["lumi_ue"] = all_data["UE/RL sys +"]
    print(f"2JET {topology} data loaded. Npoint: ", len(df))
    return df


def read_correlations(ndata_dict):
    """Read the correlation files and build a big matix"""
    corr_rows = []
    # loop on block rows
    for a, ndata_a in ndata_dict.items():
        la = [a for _ in range(ndata_a)]
        corr_row = pd.DataFrame()
        # loop on block columns
        for b, ndata_b in ndata_dict.items():
            lb = [b for _ in range(ndata_b)]
            # build the block
            try:
                with open(RAWDATA_PATH / f"corr{a}Vs{b}.tex", encoding="utf-8") as file:
                    corr_df = pd.read_csv(
                        file, sep="&", skiprows=3, skipfooter=3, engine="python"
                    )

                    # add some parsing
                    corr_vals = []
                    for val in corr_df.values.flatten():
                        if val is None:
                            continue
                        if isinstance(val, str):
                            val = val.replace("\\", "")
                            try:
                                val = float(val)
                            except ValueError:
                                continue
                        if not math.isnan(val):
                            corr_vals.append(val)

                if a == b:
                    corr = upper_triangular_to_symmetric(corr_vals, ndata_a)
                else:
                    corr = np.array(corr_vals).reshape((ndata_a, ndata_b))
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
        "lumi_ue": {
            "description": "underlying event and relative luminosity uncertainty",
            "treatment": "ADD",
            "type": f"STAR{YEAR}LUMIUE",
        },
        "pol": {
            "description": "beam polarization uncertainty",
            "treatment": "MULT",
            "type": f"STAR{YEAR}POL",
        },
    }
    # loop on data points
    for i, sys_i in enumerate(art_sys):
        e = {
            "lumi_ue": float(df.loc[i, "lumi_ue"]),
            "pol": float(df.loc[i, "pol"]),
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
    with open(STORE_PATH / "uncertainties.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


def write_2jet_data(df, topology, art_sys):
    STORE_PATH = HERE / f"../STAR_{YEAR}_2JET_{SQRTS}GEV/"
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
            "sqrts": {"min": None, "mid": float(df.loc[i, "sqrts"]), "max": None},
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
        "stat": {
            "description": "statistical uncertainty",
            "treatment": "ADD",
            "type": "UNCORR",
        },
        "lumi_ue": {
            "description": "underlying event and relative luminosity uncertainty",
            "treatment": "ADD",
            "type": f"STAR{YEAR}LUMIUE",
        },
        "pol": {
            "description": "beam polarization uncertainty",
            "treatment": "MULT",
            "type": f"STAR{YEAR}POL",
        },
    }
    # loop on data points
    for i, sys_i in enumerate(art_sys):
        e = {
            "stat": float(df.loc[i, "stat"]),
            "lumi_ue": float(df.loc[i, "lumi_ue"]),
            "pol": float(df.loc[i, "pol"]),
        }
        # loop on art sys
        for j, val in enumerate(sys_i):
            e[f"sys_{j}"] = val
        error.append(e)

        if i == 0:
            error_definition.update(
                {
                    f"sys_{j}": {
                        "description": f"{j} artificial correlated systematics uncertainty",
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
    # from the paper we understand that stat dijet is not correlated
    #    I-I (stat + sys) | I-D (stat + sys)
    #    D-I (stat + sys) | D-D (sys)
    correlated_unc = np.sqrt(
        dfs["I"]["syst"] ** 2 + dfs["I"]["stat"] ** 2
    ).values.tolist()
    for a in TOPOPLOGY_LIST[1:]:
        correlated_unc.extend(dfs[a]["syst"].values)
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
