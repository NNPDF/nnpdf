"""
This script is used to filter the data in:
rawdata/HEPData-ins1810913-v1-Figure_A*.yaml
with uncertainty breakdowns in:
rawdata/HEPData-ins1810913-v1-Impacts_Figure_A14*.yaml

This is for the "boson-rap" variant, which is a variant with respect
to the rapidity of the W boson.

The HepData table for WP and WM only contains relative uncertainties (in %) as a function
of |y_W|, so absolute uncertainties are computed as  (rel_unc / 100) * central_value.

For the Asymmetry error the absolute uncertainty is quoted in the HEP data table directly.
"""

import numpy as np
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

# -------- Filter WP and WM - See below for ASY

IMPACTS_TABLE = "rawdata/HEPData-ins1810913-v1-Impacts_Figure_A14b.yaml"
CENTRAL_TABLE = "rawdata/HEPData-ins1810913-v1-Figure_A3.yaml"  
UNIT_FACTOR = 1  # convert to fb
M_W2 = 6460.5


UNC_MAP = {
    "Relative uncertainty due to statistical uncertainty in data": (
        "Statistical uncertainty in data", "ADD", "UNCORR"),
    "Relative uncertainty due to statistical uncertainty on lepton efficiency": (
        "Statistical uncertainty on lepton efficiency", "MULT", "UNCORR"),
    "Relative uncertainty due to statistical uncertainty in simulation": (
        "Statistical uncertainty in simulation", "MULT", "UNCORR"),
    "Relative uncertainty due to QCD background": (
        "QCD background", "MULT", "CORR"),
    "Relative uncertainty due to QCD factorization and normalization scales": (
        "QCD factorization and normalization scales", "MULT", "CORR"),
    r"Relative uncertainty due to PDFs and $\alpha_{S}$": (
        "PDFs and alpha_S", "MULT", "CORR"),
    "Relative uncertainty due to luminosity": (
        "Luminosity", "MULT", "CMSLUMI13"),
}
SKIP = {"Relative uncertainty due to total uncertainty"} # This is the total uncertainty calculated by experimentalists taking errors into account. We are using the breakdown instead.


def load_yaml(path: str) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def get_kinematics(hepdata: dict) -> list:
    """Return the |y_W| bins as a list of kinematic dictionaries."""
    kin = []
    for b in hepdata["independent_variables"][0]["values"]:
        kin.append(
            {
                "abs_y": {
                    "min": b["low"],
                    "mid": 0.5 * (b["low"] + b["high"]),
                    "max": b["high"],
                },
                "m_W2": {"min": None, "mid": M_W2, "max": None},
            }
        )
    return kin


def get_central_values(hepdata: dict) -> np.ndarray:
    """Central cross-section values from the cross-section table."""
    values = hepdata["dependent_variables"][0]["values"]
    return np.array([v["value"] * UNIT_FACTOR for v in values])


def get_relative_uncertainties(hepdata: dict) -> dict:
    """Return {header name: np.ndarray of relative uncertainties in %}."""
    rel = {}
    for dep in hepdata["dependent_variables"]:
        name = dep["header"]["name"]
        rel[name] = np.array([v["value"] for v in dep["values"]])
    return rel


def main_filter() -> None:
    impacts = load_yaml(IMPACTS_TABLE)
    central = get_central_values(load_yaml(CENTRAL_TABLE))
    kinematics = get_kinematics(impacts)
    rel = get_relative_uncertainties(impacts)

    n_bins = len(kinematics)
    if len(central) != n_bins:
        raise ValueError(
            f"Central values ({len(central)}) and impact bins ({n_bins}) differ in length"
        )

    # Sanity check: unknown headers
    for name in rel:
        if name not in UNC_MAP and name not in SKIP:
            print(f"WARNING: uncertainty '{name}' is not in UNC_MAP and will be ignored")

    # ---- Build definitions and absolute uncertainties ------------------
    definitions = {}
    abs_unc = {}  # key -> array of absolute uncertainties

    stat_name = "Relative uncertainty due to statistical uncertainty in data"
    definitions["stat"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }
    abs_unc["stat"] = rel[stat_name] / 100.0 * central

    i = 1
    for name, (descr, treatment, utype) in UNC_MAP.items():
        if name == stat_name or name not in rel:
            continue
        key = f"sys_corr_{i}"
        definitions[key] = {
            "description": f"sys,{descr}",
            "treatment": treatment,
            "type": utype,
        }
        abs_unc[key] = rel[name] / 100.0 * central
        i += 1

    # ---- Format bins ----------------------------------------------------
    bins = []
    for b in range(n_bins):
        bins.append({key: float(abs(vals[b])) for key, vals in abs_unc.items()})

    # ---- Dump ------------------------------------------------------------
    with open("data_ASY_y.yaml", "w") as f:
        yaml.dump({"data_central": central.tolist()}, f, sort_keys=False)

    with open("kinematics_ASY_y.yaml", "w") as f:
        yaml.dump({"bins": kinematics}, f, sort_keys=False)

    with open("uncertainties_ASY_y.yaml", "w") as f:
        yaml.dump({"definitions": definitions, "bins": bins}, f, sort_keys=False)

    print(f"Total bins: {n_bins}")


# ---------- Filter asymmetry data


IMPACTS_TABLE_ASY = "rawdata/HEPData-ins1810913-v1-Impacts_Figure_A14c.yaml"
CENTRAL_TABLE_ASY = "rawdata/HEPData-ins1810913-v1-Figure_A4.yaml"  


UNC_MAP_ASY = {
    "Absolute uncertainty due to statistical uncertainty in data": (
        "Statistical uncertainty in data", "ADD", "UNCORR"),
    "Absolute uncertainty due to statistical uncertainty on lepton efficiency": (
        "Statistical uncertainty on lepton efficiency", "MULT", "UNCORR"),
    "Absolute uncertainty due to statistical uncertainty in simulation": (
        "Statistical uncertainty in simulation", "MULT", "UNCORR"),
    "Absolute uncertainty due to QCD background": (
        "QCD background", "MULT", "CORR"),
    "Absolute uncertainty due to QCD factorization and normalization scales": (
        "QCD factorization and normalization scales", "MULT", "CORR"),
    r"Absolute uncertainty due to PDFs and $\alpha_{S}$": (
        "PDFs and alpha_S", "MULT", "CORR"),
    "Absolute uncertainty due to luminosity": (
        "Luminosity", "MULT", "CMSLUMI13"),
}
SKIP = {"Absolute uncertainty due to total uncertainty"} # This is the total uncertainty calculated by experimentalists taking errors into account. We are using the breakdown instead.


def load_yaml_asy(path: str) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def get_kinematics_asy(hepdata: dict) -> list:
    """Return the |y_W| bins as a list of kinematic dictionaries."""
    kin = []
    for b in hepdata["independent_variables"][0]["values"]:
        kin.append(
            {
                "abs_y": {
                    "min": b["low"],
                    "mid": 0.5 * (b["low"] + b["high"]),
                    "max": b["high"],
                },
                "m_W2": {"min": None, "mid": M_W2, "max": None},
            }
        )
    return kin


def get_central_values_asy(hepdata: dict) -> np.ndarray:
    """Central cross-section values from the cross-section table."""
    values = hepdata["dependent_variables"][0]["values"]
    return np.array([v["value"] * UNIT_FACTOR for v in values])


def get_uncertainties_asy(hepdata: dict) -> dict:
    """Return {header name: np.ndarray of relative uncertainties in %}."""
    error = {}
    for dep in hepdata["dependent_variables"]:
        name = dep["header"]["name"]
        error[name] = np.array([v["value"] for v in dep["values"]])
    return error


def main_filter_asy() -> None:
    impacts = load_yaml_asy(IMPACTS_TABLE)
    central = get_central_values_asy(load_yaml(CENTRAL_TABLE))
    kinematics = get_kinematics_asy(impacts)
    error = get_uncertainties_asy(impacts)

    n_bins = len(kinematics)
    if len(central) != n_bins:
        raise ValueError(
            f"Central values ({len(central)}) and impact bins ({n_bins}) differ in length"
        )

    # Sanity check: unknown headers
    for name in error:
        if name not in UNC_MAP_ASY and name not in SKIP:
            print(f"WARNING: uncertainty '{name}' is not in UNC_MAP_ASY and will be ignored")

    # ---- Build definitions and absolute uncertainties ------------------
    definitions = {}
    abs_unc = {}  # key -> array of absolute uncertainties

    stat_name = "Absolute uncertainty due to statistical uncertainty in data"
    definitions["stat"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }
    abs_unc["stat"] = error[stat_name]

    i = 1
    for name, (descr, treatment, utype) in UNC_MAP_ASY.items():
        if name == stat_name or name not in error:
            continue
        key = f"sys_corr_{i}"
        definitions[key] = {
            "description": f"sys,{descr}",
            "treatment": treatment,
            "type": utype,
        }
        abs_unc[key] = error[name]
        i += 1

    # ---- Format bins ----------------------------------------------------
    bins = []
    for b in range(n_bins):
        bins.append({key: float(abs(vals[b])) for key, vals in abs_unc.items()})

    # ---- Dump ------------------------------------------------------------
    with open("data_ASY_y.yaml", "w") as f:
        yaml.dump({"data_central": central.tolist()}, f, sort_keys=False)

    with open("kinematics_ASY_y.yaml", "w") as f:
        yaml.dump({"bins": kinematics}, f, sort_keys=False)

    with open("uncertainties_ASY_y.yaml", "w") as f:
        yaml.dump({"definitions": definitions, "bins": bins}, f, sort_keys=False)

    print(f"Total bins: {n_bins}")


if __name__ == "__main__":
    main_filter()