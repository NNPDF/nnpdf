import yaml
import numpy as np
import pandas as pd

TABLE_TO_RAPIDITY = {
    1: [0.0, 0.5],
    2: [0.5, 1.0],
    3: [1.0, 1.5],
    4: [1.5, 2],
    5: [2.0, 2.5],
    6: [2.5, 3.0],
}


def get_data_values(tables):
    """
    returns the central data

    Parameters
    ----------
    tables : list
            list that enumerates the table number

    Returns
    -------
    list
        list containing the central values for all
        hepdata tables

    """

    data_central = []
    for table in tables:
        hepdata_table = f"rawdata/atlas_inclusive_jet2012_r06_altcorr1_eta{table}.yaml"

        with open(hepdata_table, "r") as file:
            input = yaml.safe_load(file)

        values = input["dependent_variables"][0]["values"]

        for value in values:
            data_central.append(value["value"])

    return data_central


def get_kinematics(tables):
    """
    returns the relevant kinematics values

    Parameters
    ----------
    tables : list
            list that enumerates the table number

    Returns
    -------
    list
        list containing the kinematic values for all
        hepdata tables
    """
    kin = []

    for table in tables:
        hepdata_table = f"rawdata/atlas_inclusive_jet2012_r06_altcorr1_eta{table}.yaml"

        with open(hepdata_table, "r") as file:
            input = yaml.safe_load(file)

        # rapidity
        rapidity_interval = TABLE_TO_RAPIDITY[table]
        rap = {}
        rap["min"], rap["max"] = rapidity_interval[0], rapidity_interval[1]
        rap["mid"] = 0.5 * (rap["min"] + rap["max"])

        # center of mass energy
        sqrts = float(input["dependent_variables"][0]["qualifiers"][1]["value"])

        # transverse momentum
        jet_kt_bins = input["independent_variables"][0]["values"]
        KT = {}
        for kt in jet_kt_bins:
            KT["min"], KT["max"] = kt["low"], kt["high"]
            KT["mid"] = float(f"{0.5 * (kt['low'] + kt['high']):.3f}")

            kin_value = {
                "y": {"min": rap["min"], "mid": rap["mid"], "max": rap["max"]},
                "kt": {"min": KT["min"], "mid": KT["mid"], "max": KT["max"]},
                "sqrt_s": {"min": None, "mid": sqrts, "max": None},
            }

            kin.append(kin_value)

    return kin



def decompose_covmat(covmat):
    """Given a covmat it return an array sys with shape (ndat,ndat)
    giving ndat correlated systematics for each of the ndat point.
    The original covmat is obtained by doing sys@sys.T"""

    lamb, mat = np.linalg.eig(covmat)
    sys = np.multiply(np.sqrt(lamb), mat)
    return sys


def HEP_table_to_df(table):
    """
    Given hep data table return a pandas
    DataFrame with index given by Ndata,
    columns by the uncertainties and
    np.nan entries

    Parameters
    ----------
    table : int
            number of table

    """
    hepdata_table = f"rawdata/atlas_inclusive_jet2012_r06_altcorr1_eta{table}.yaml"
    with open(hepdata_table, "r") as file:
        card = yaml.safe_load(file)

    # list of len ndata. Each entry is dict with
    # keys errors and value
    card = card["dependent_variables"][0]["values"]
    df_nan = pd.DataFrame(index=range(1, len(card) + 1))

    errors = card[0]["errors"]
    for err in errors:
        # luminosity and stat uncertainty, always symmetric
        if err["label"] == "syst_lumi":
            df_nan[err["label"]] = np.nan
        elif err["label"] == "stat":
            df_nan[err["label"]] = np.nan
        else:
            df_nan[f"{err['label']}_plus"] = np.nan
            df_nan[f"{err['label']}_minus"] = np.nan

    return df_nan


def fill_df(table):
    """
    Fill a data frame with index
    corresponding to measured datapoints
    and columns to different uncertainties
    Each df is for a fixed rapidity bin.

    Parameters
    ----------
    table : int
            number of table

    """
    hepdata_table = f"rawdata/atlas_inclusive_jet2012_r06_altcorr1_eta{table}.yaml"
    with open(hepdata_table, "r") as file:
        card = yaml.safe_load(file)

    # list of len ndata. Each entry is dict with
    # keys errors and value
    card = card["dependent_variables"][0]["values"]
    df = HEP_table_to_df(table)

    for i, dat in enumerate(card):
        # cv = dat['value']

        for j, err in enumerate(dat["errors"]):
            # lum and stat uncertainties are always symmetric
            if err["label"] == "stat":
                df.loc[df.index == i + 1, err["label"]] = err["asymerror"]["plus"]
            elif err["label"] == "syst_lumi":
                df.loc[df.index == i + 1, err["label"]] = err["asymerror"]["plus"]
            else:
                dp, dm = process_err(err)
                df.loc[
                    df.index == i + 1, f"{err['label']}_plus"
                ] = dp  # err['asymerror']['plus']
                df.loc[
                    df.index == i + 1, f"{err['label']}_minus"
                ] = dm  # err['asymerror']['minus']
    
    return df


def process_err(error):
    """
    Note: the d'Agostini prescription for the
    symmetrization of the error does not hold here.
    We follow here the experimental prescription

    Parameters
    ----------

    error : dictionary
            e.g. {'label': 'sys', 'symerror': 0.1%}

    Returns
    -------
    tuple
        tuple containing two floats

    """
    d_p = float(error["asymerror"]["plus"])
    d_m = float(error["asymerror"]["minus"])

    tmp1 = d_p
    tmp2 = d_m
    # case 1: d_p and d_m are both negative
    if tmp1 < 0.0 and tmp2 < 0.0:
        if tmp2 < tmp1:
            d_p = 0.0
            d_m = tmp2
        else:
            d_p = 0.0
            d_m = tmp1

    # case 2: d_p and d_m are both positive
    if tmp1 > 0.0 and tmp2 > 0.0:
        if tmp1 > tmp2:
            d_p = tmp1
            d_m = 0.0
        else:
            d_p = tmp2
            d_m = 0.0
    return d_p, d_m