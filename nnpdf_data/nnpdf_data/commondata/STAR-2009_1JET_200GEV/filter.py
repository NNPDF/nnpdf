"""This script provides the common filer to the jet and dijet STAR 2009 datasets.
Files need to be parsed all together as there are correlations provided.
"""

import pathlib

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.correlations import compute_covmat

TOPOPLOGY_LIST = ["CC", "CF", "SS", "OS", "A", "B", "C"]
POL_UNC = 0.065
LUMI_UNC = 0.0005
YEAR = 2009
HERE = pathlib.Path(__file__).parent

# NOTE: this is not the full relevant as the observable is symmetric
# for jet1 and jet2, so 1 and 2 are not ordered in pT and the
# information about the sign in missing.
TOPO_DEF = {
    "B": {"eta1_min": -0.8, "eta1_max": 0, "eta2_min": 0.8, "eta2_max": 1.8},
    "C": {"eta1_min": 0, "eta1_max": 0.8, "eta2_min": 0.8, "eta2_max": 1.8},
    "A": {"eta1_min": 0.8, "eta1_max": 1.8, "eta2_min": 0.8, "eta2_max": 1.8},
    "SS": {"abs_eta_min": 0, "abs_eta_max": 0.8},
    "OS": {"abs_eta_min": 0, "abs_eta_max": 0.8},
    "CC": {"abs_eta_min": 0, "abs_eta_max": 0.5},
    "CF": {"abs_eta_min": 0.5, "abs_eta_max": 1.0},
}


def read_1jet_data(topology):
    tab_number = 3 if "CC" in topology else 4
    fnames = [
        f"rawdata/Table_{tab_number}_ALL.csv",
        f"rawdata/Table_{tab_number}_pT.csv",
        f"rawdata/Table_5.csv",
    ]

    df = pd.DataFrame()
    for fname in fnames:
        dfi = pd.read_csv(fname)
        if "pT" in fname:
            df["pT"] = dfi["Parton Jet $p_T$ [GeV/c]"]
            df["pT_min"] = df["pT"] - dfi["sys +"]
            df["pT_max"] = df["pT"] + dfi["sys +"]

        elif "ALL" in fname:
            df["ALL"] = dfi["$A_{LL}$"]
            df["stat_min"] = dfi["stat -"]
            df["stat_max"] = dfi["stat +"]
            df["sys_min"] = dfi["sys -"]
            df["sys_max"] = dfi["sys +"]
        # elif "5" in fname:
        #     dfc_col = pd.read_csv(fname)

        #     dfc = pd.DataFrame()
        #     biny = 1
        #     for i in range(len(dfc_col)):
        #         if dfc_col.loc[i, "binx"] == "-":
        #             biny = float(dfc_col.loc[i, "val"])
        #         else:
        #             binx = float(dfc_col.loc[i, "binx"])
        #             dfc.loc[binx, biny] = dfc_col.loc[i, "val"]
        #     dfc = dfc.astype("float")

    df["abs_eta_min"] = TOPO_DEF[topology]["abs_eta_min"]
    df["abs_eta_max"] = TOPO_DEF[topology]["abs_eta_max"]
    df["abs_eta"] = (df["abs_eta_min"] + df["abs_eta_max"]) / 2
    df["stat"] = df["stat_max"]
    df["sys"] = df["sys_max"]
    df["pol"] = POL_UNC * abs(df["ALL"])
    df["lumi"] = LUMI_UNC
    df["sqrts"] = 200
    return df


def read_2jet_data(topology):
    if "S" in topology:
        fname = f"../STAR-{YEAR}_2JET_200GEV_MIDRAP/rawdata/Table_{topology}.yaml"
    else:
        fname = f"../STAR-{YEAR}_2JET_200GEV/rawdata/Table_{topology}.yaml"
    df = pd.DataFrame()
    with open(fname, "r", encoding="utf-8") as file:
        data = yaml.safe_load(file)

    if topology in ["A", "B", "C"]:
        Msub = data["independent_variables"][0]["values"]
        Gsub = data["dependent_variables"][0]["values"]

        for i in range(len(Msub)):
            df = pd.concat(
                [
                    df,
                    pd.DataFrame(
                        {
                            "m_low": [Msub[i]["low"]],
                            "m_high": [Msub[i]["high"]],
                            "ALL": [Gsub[i]["value"]],
                            "stat": [Gsub[i]["errors"][0]["symerror"]],
                            "sys": [Gsub[i]["errors"][1]["symerror"]],
                        }
                    ),
                ],
                ignore_index=True,
            )
        df["m"] = (df["m_low"] + df["m_high"]) / 2
        df["eta1_min"] = TOPO_DEF[topology]["eta1_min"]
        df["eta1_max"] = TOPO_DEF[topology]["eta1_max"]
        df["eta2_min"] = TOPO_DEF[topology]["eta2_min"]
        df["eta2_max"] = TOPO_DEF[topology]["eta2_max"]
        df["eta1"] = (df["eta1_min"] + df["eta1_max"]) / 2
        df["eta2"] = (df["eta2_min"] + df["eta2_max"]) / 2
    else:
        Msub = data["dependent_variables"][0]["values"]
        Gsub = data["dependent_variables"][1]["values"]

        for i in range(len(Msub)):
            df = pd.concat(
                [
                    df,
                    pd.DataFrame(
                        {
                            "m": [Msub[i]["value"]],
                            "ALL": [Gsub[i]["value"]],
                            "stat": [Gsub[i]["errors"][0]["symerror"]],
                            "sys": [Gsub[i]["errors"][1]["symerror"]],
                            "abs_eta_min": TOPO_DEF[topology]["abs_eta_min"],
                            "abs_eta_max": TOPO_DEF[topology]["abs_eta_max"],
                        }
                    ),
                ],
                ignore_index=True,
            )

    df["pol"] = POL_UNC * abs(df["ALL"])
    df["lumi"] = LUMI_UNC
    df["sqrts"] = 200
    return df


def write_1jet_data(df, topology, art_sys):
    STORE_PATH = HERE

    # Write central data
    data_central_yaml = {"data_central": list(df["ALL"])}
    with open(STORE_PATH / f"data_{topology}.yaml", "w", encoding="utf-8") as file:
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
    with open(STORE_PATH / f"kinematics_{topology}.yaml", "w", encoding="utf-8") as file:
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
    with open(STORE_PATH / f"uncertainties_{topology}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


def write_2jet_data(df, topology, art_sys):
    STORE_PATH = f"../STAR-{YEAR}_2JET_200GEV"
    if "S" in topology:
        STORE_PATH += "_MIDRAP"
    STORE_PATH = pathlib.Path(STORE_PATH)
    # Write central data
    data_central_yaml = {"data_central": list(df["ALL"])}
    with open(STORE_PATH / f"data_{topology}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(data_central_yaml, file)

    # Write kin file
    kin = []
    for i in range(len(df)):
        try:
            kin_value = {
                "m_jj": {
                    "min": float(df.loc[i, "m_low"]),
                    "mid": float(df.loc[i, "m"]),
                    "max": float(df.loc[i, "m_high"]),
                },
                # "sqrts": {"min": None, "mid": float(df.loc[i, "sqrts"]), "max": None},
                "eta_1": {
                    "min": float(df.loc[i, "eta1_min"]),
                    "mid": float(df.loc[i, "eta1"]),
                    "max": float(df.loc[i, "eta1_max"]),
                },
                "eta_2": {
                    "min": float(df.loc[i, "eta2_min"]),
                    "mid": float(df.loc[i, "eta2"]),
                    "max": float(df.loc[i, "eta2_max"]),
                },
            }
        except KeyError:
            df["abs_eta"] = (df["abs_eta_min"] + df["abs_eta_max"]) / 2
            kin_value = {
                "m_jj": {"min": None, "mid": float(df.loc[i, "m"]), "max": None},
                # "sqrts": {"min": None, "mid": float(df.loc[i, "sqrts"]), "max": None},
                "abs_eta_1": {
                    "min": float(df.loc[i, "abs_eta_min"]),
                    "mid": float(df.loc[i, "abs_eta"]),
                    "max": float(df.loc[i, "abs_eta_max"]),
                },
                "abs_eta_2": {
                    "min": float(df.loc[i, "abs_eta_min"]),
                    "mid": float(df.loc[i, "abs_eta"]),
                    "max": float(df.loc[i, "abs_eta_max"]),
                },
            }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}
    with open(STORE_PATH / f"kinematics_{topology}.yaml", "w", encoding="utf-8") as file:
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
    with open(STORE_PATH / f"uncertainties_{topology}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    # load all the data
    dfs = {}
    for topo in TOPOPLOGY_LIST[:2]:
        dfs[topo] = read_1jet_data(topo)
    for topo in TOPOPLOGY_LIST[2:]:
        dfs[topo] = read_2jet_data(topo)

    # load correlations
    ndata_dict = {a: len(b) for a, b in dfs.items()}
    correlation_df = pd.read_csv("rawdata/correlation.csv", index_col=0)
    # from the supplement material:
    # https://journals.aps.org/prd/supplemental/10.1103/PhysRevD.98.032011/Supplementalmaterial.pdf
    # we understand that stat jet and dijet are correlated,
    # see also https://github.com/NNPDF/nnpdf/pull/2035#issuecomment-2201979662
    correlated_unc = []
    for a in TOPOPLOGY_LIST:
        correlated_unc.extend(np.sqrt(dfs[a]["sys"] ** 2 + dfs[a]["stat"] ** 2).values.tolist())
    ndata_points = np.sum((*ndata_dict.values(),))
    # decompose uncertainties
    art_sys = np.array(compute_covmat(correlation_df, correlated_unc, ndata_points))

    # write data
    cnt = 0
    for topo, df in dfs.items():
        ndata = ndata_dict[topo]
        syst = art_sys[cnt : cnt + ndata, :].tolist()
        if topo in ["CC", "CF"]:
            write_1jet_data(df, topo, syst)
        else:
            write_2jet_data(df, topo, syst)
        cnt += ndata
