"""This script provides the common filer to the jet and dijet STAR 2015 datasets.
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

# values from the paper
SQRTS = 200
YEAR = 2015
# mapping between topologies, tables and abs_eta values
TOPOPLOGY_LIST = {
    "CC": ("bottom", 0, 0.5),
    "CF": ("top", 0.5, 1.0),
    "OS": "bottom",
    "SS": "top",
}

# mapping between correlations blocks and tables
MAP_CORR_TABLE = {
    ("CC", "CC"): 2,
    ("CC", "CF"): 4,
    ("CC", "OS"): 8,
    ("CC", "SS"): 7,
    ("CF", "CF"): 1,
    ("CF", "OS"): 6,
    ("CF", "SS"): 5,
    ("OS", "OS"): 12,
    ("OS", "SS"): 13,
    ("SS", "SS"): 11,
}

HERE = pathlib.Path(__file__).parent
RAWDATA_PATH = HERE / "rawdata/"


def read_1jet_data(topology):
    table_label = TOPOPLOGY_LIST[topology][0]
    max_eta = TOPOPLOGY_LIST[topology][2]
    min_eta = TOPOPLOGY_LIST[topology][1]
    data_table = pathlib.Path(RAWDATA_PATH / f"Table1{table_label}.csv")

    with open(data_table, "r", encoding="utf-8") as file:
        parton_jet_data = pd.read_csv(
            file, delimiter=",", skiprows=lambda x: (x <= 6 or x >= 20)
        )
    with open(data_table, "r", encoding="utf-8") as file:
        all_data = pd.read_csv(file, delimiter=",", skiprows=21)

    df = pd.DataFrame()
    df["pT"] = parton_jet_data[
        r"Inclusive jet transverse momentum $p_T$ at the parton level [GeV/$c$]"
    ]
    df["pT_min"] = df["pT"] + parton_jet_data["Syst -"]
    df["pT_max"] = df["pT"] + parton_jet_data["Syst +"]
    df["abs_eta"] = (max_eta + min_eta) / 2
    df["abs_eta_min"] = min_eta
    df["abs_eta_max"] = max_eta
    df["sqrts"] = SQRTS
    df["ALL"] = all_data[r"Double spin asymmetry $A_{LL}$"]
    df["stat"] = all_data["Stat +"]
    df["syst"] = all_data["Syst +"]
    df["lumi"] = all_data["Lumi +"]
    df["pol"] = [float(a[:-1]) for a in all_data["Pol +"]] * abs(df["ALL"]) / 100

    print(f"1JET {topology} data loaded. Npoint: ", len(df))
    return df


def read_2jet_data(topology):
    table_label = TOPOPLOGY_LIST[topology]
    data_table = RAWDATA_PATH / f"Table2{table_label}.csv"
    with open(data_table, "r", encoding="utf-8") as file:
        mjj_data = pd.read_csv(
            file, delimiter=",", skiprows=lambda x: (x <= 6 or x >= 16)
        )
    with open(data_table, "r", encoding="utf-8") as file:
        all_data = pd.read_csv(file, delimiter=",", skiprows=27)

    df = pd.DataFrame()
    df["mjj"] = mjj_data["Dijet invariant mass $M_{inv}$ at the parton level [GeV/$c$]"]
    df["mjj_min"] = df["mjj"] + mjj_data["Syst -"]
    df["mjj_max"] = df["mjj"] + mjj_data["Syst +"]
    df["sqrts"] = SQRTS
    df["ALL"] = all_data[r"Double spin asymmetry $A_{LL}$"]
    df["stat"] = all_data["Stat +"]
    df["syst"] = all_data["Syst +"]
    df["lumi"] = all_data["Lumi +"]
    df["pol"] = [float(a[:-1]) for a in all_data["Pol +"]] * abs(df["ALL"]) / 100

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
                idx = MAP_CORR_TABLE[(a, b)]
                with open(
                    RAWDATA_PATH / f"Table{idx}SupplementalMaterial.csv",
                    encoding="utf-8",
                ) as file:
                    corr_df = pd.read_csv(file, delimiter=",", skiprows=7)

                if a == b:
                    corr = upper_triangular_to_symmetric(corr_df.values[:, 2], ndata_a)
                else:
                    corr = corr_df.values[:, 2].reshape((ndata_a, ndata_b))
            except (FileNotFoundError, KeyError) as _:
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
            "abs_eta": {
                "min": float(df.loc[i, "abs_eta_min"]),
                "mid": float(df.loc[i, "abs_eta"]),
                "max": float(df.loc[i, "abs_eta_max"]),
            },
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}
    with open(STORE_PATH / "kinematics.yaml", "w", encoding="utf-8") as file:
        yaml.dump(kinematics_yaml, file)

    # Write unc file
    error = []
    error_definition = {
        "lumi": {
            "description": "luminosity uncertainty",
            "treatment": "ADD",
            "type": f"STAR{YEAR}LUMI",
        },
        "pol": {
            "description": "beam polarization uncertainty",
            "treatment": "MULT",
            "type": f"STAR{YEAR}POL",
        },
    }
    # loop on data points
    for i, sys_i in enumerate(art_sys):
        e = {"lumi": float(df.loc[i, "lumi"]), "pol": float(df.loc[i, "pol"])}
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
    STORE_PATH = HERE / f"../STAR_{YEAR}_2JET_MIDRAP_{SQRTS}GEV/"
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
        "lumi": {
            "description": "luminosity uncertainty",
            "treatment": "ADD",
            "type": f"STAR{YEAR}LUMI",
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
            "lumi": float(df.loc[i, "lumi"]),
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
    dfs = {}
    for topo in TOPOPLOGY_LIST:
        fcn = read_1jet_data if "C" in topo else read_2jet_data
        dfs[topo] = fcn(topo)

    # load correlations
    ndata_dict = {a: len(b) for a, b in dfs.items()}
    correlation_df = read_correlations(ndata_dict)

    # merge the 1JET topologies
    new_dfs = {
        "1JET": pd.concat([dfs["CC"], dfs["CF"]], ignore_index=True),
        "OS": dfs["OS"],
        "SS": dfs["SS"],
    }

    # from the paper we understand that stat dijet is not correlated
    #    I-I (stat + sys) | I-D (stat + sys)
    #    D-I (stat + sys) | D-D (sys)
    correlated_unc = np.sqrt(
        new_dfs["1JET"]["syst"] ** 2 + new_dfs["1JET"]["stat"] ** 2
    ).values.tolist()
    for a in (*TOPOPLOGY_LIST,)[2:]:
        correlated_unc.extend(new_dfs[a]["syst"].values)

    # decompose uncertainties
    ndata_dict = {a: len(b) for a, b in new_dfs.items()}
    ndata_points = np.sum((*ndata_dict.values(),))
    art_sys = np.array(compute_covmat(correlation_df, correlated_unc, ndata_points))

    # write data
    cnt = 0
    for topo, df in new_dfs.items():
        ndata = ndata_dict[topo]
        syst = art_sys[cnt : cnt + ndata, :].tolist()
        if topo == "1JET":
            write_1jet_data(df, syst)
        else:
            write_2jet_data(df, topo, syst)
        cnt += ndata
