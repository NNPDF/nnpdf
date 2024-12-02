r"""Implement data form Hepdata reference, data from Fig 4,5.
The last points for (nu and nb) for Table121 have been copied manually from the 
legacy version of the dataset, as they are not present in the original Hepdata tables.
"""

import pathlib

import pandas as pd
import yaml

HERE = pathlib.Path(__file__).parent

# Nucleon Mass determined using scikithep/particle for Pb208 in GeV
M_NEUTRON = 939.565346 * 0.001
M_PROTON = 938.272013 * 0.001
# A = 208  # A(Pb): Atomic Mass
# Z = 82  # Z(Pb): Atomic Number
# MN = 193.729 / (Z * M_PROTON + (A - Z) * M_NEUTRON)
# The old implementation used the isoscalar average.
MN = (M_PROTON + M_NEUTRON) / 2

# Normalization factor, to match old implementation
NORM_FACT = 10


def read_legacy_uncertainties(legacy_unc_file):
    """Parse the CHORUSISOTARGCOR and CHORUSQEDRADCOR from the legacy variant."""
    with open(legacy_unc_file, encoding="utf-8") as f:
        legacy_unc_nb = yaml.safe_load(f)

    isotargt_err = []
    qed_err = []
    for p in legacy_unc_nb['bins']:
        isotargt_err.append(p["sys_corr_1"])
        qed_err.append(p["sys_corr_2"])
    return isotargt_err, qed_err


def read_tables():
    """Parse Tables."""
    dfs_nu = pd.DataFrame()
    dfs_nub = pd.DataFrame()
    columns = ["y", "sigma", "stat+", "stat-", "sys+", "sys-"]

    # here we want to keep the same ordering as
    # the hepdata tables
    for tab_idx in range(23, 122):
        file = HERE / "rawdata" / f"Table{tab_idx}.csv"
        with open(file, "r", encoding="utf-8") as f:
            # skip intro
            lines = f.readlines()[9:]
            # determine df size
            size = int((len(lines) - 2) / 2)

            # read values
            df_nu = pd.DataFrame([l.split(",") for l in lines[5:size]], columns=columns)
            df_nub = pd.DataFrame([l.split(",") for l in lines[size + 6 : -1]], columns=columns)
            df_nu["x"] = float(lines[3].split(",")[1])
            df_nub["x"] = float(lines[size + 4].split(",")[1])
            df_nu["E"] = float(lines[0].split(",")[1])
            df_nub["E"] = float(lines[size + 1].split(",")[1])

            dfs_nu = pd.concat([dfs_nu, df_nu], ignore_index=True) if not dfs_nu.empty else df_nu
            dfs_nub = (
                pd.concat([dfs_nub, df_nub], ignore_index=True) if not dfs_nub.empty else df_nub
            )

    # here we need to load also 2 uncertainties from the legacy file,
    # as these are not present in hepdata but wes still want to keep the,
    isotargt_err_nu, qed_err_nu = read_legacy_uncertainties("uncertainties_legacy_NU-SIGMARED.yaml")
    isotargt_err_nub, qed_err_nub = read_legacy_uncertainties(
        "uncertainties_legacy_NB-SIGMARED.yaml"
    )

    dfs_nu["isotargt_err"] = isotargt_err_nu
    dfs_nub["isotargt_err"] = isotargt_err_nub
    dfs_nu["qed_err"] = qed_err_nu
    dfs_nub["qed_err"] = qed_err_nub

    dfs_nu = dfs_nu.astype(float)
    dfs_nub = dfs_nub.astype(float)

    # compute Q2
    dfs_nu["Q2"] = 2 * MN * dfs_nu.x * dfs_nu.y * dfs_nu.E
    dfs_nub["Q2"] = 2 * MN * dfs_nub.x * dfs_nub.y * dfs_nub.E

    return dfs_nu, dfs_nub


def write_files(df, variant):
    """Write kinematics, central value and uncertainties files."""

    # Write central data
    data_central_yaml = {"data_central": [round(float(x), 4) for x in df["sigma"] / NORM_FACT]}
    with open(HERE / f"data_{variant}.yaml", "w", encoding="utf-8") as file:
        yaml.dump(data_central_yaml, file)

    # Write kin file
    kin = []
    for _, row in df.iterrows():
        kin_value = {
            "x": {"min": None, "mid": float(row.x), "max": None},
            "Q2": {"min": None, "mid": float(row.Q2), "max": None},
            "y": {"min": None, "mid": float(row.y), "max": None},
        }
        kin.append(kin_value)
    kinematics_yaml = {"bins": kin}
    # note kinematics is the same for nu and nb
    with open(HERE / "kinematics.yaml", "w", encoding="utf-8") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    # loop on data points
    error_definition = {
        "stat": {"description": "statistical uncertainty", "treatment": "ADD", "type": "UNCORR"},
        "sys": {"description": "systematical uncertainty", "treatment": "MULT", "type": "CORR"},
        # NOTE: these uncertainties are directly coming from the legacy variant
        "CHORUSISOTARGCOR": {
            "description": "Isoscalar target correction interpreted as uncertainty",
            "treatment": "MULT",
            "type": "CORR",
        },
        "CHORUSQEDRADCOR": {
            "description": "QED radiation correction interpreted as uncertainty",
            "treatment": "MULT",
            "type": "CORR",
        },
    }
    error = []
    for _, row in df.iterrows():
        e = {
            "stat": float(row["stat+"]) / NORM_FACT,
            "sys": float(row["sys+"]) / NORM_FACT,
            "CHORUSISOTARGCOR": float(row["isotargt_err"]),
            "CHORUSQEDRADCOR": float(row["qed_err"]),
        }
        error.append(e)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}
    with open(HERE / f"uncertainties_{variant}_hepdata.yaml", "w", encoding="utf-8") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    dfs_nu, dfs_nub = read_tables()
    write_files(dfs_nu, "nu")
    write_files(dfs_nub, "nb")
