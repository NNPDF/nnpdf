# ignore pandas warning
import warnings

import numpy as np
import pandas as pd
import yaml

warnings.filterwarnings('ignore')

SCENARIO = {'nominal': 0, 'stronger': 1, 'weaker': 2}


def range_str_to_floats(str_range):
    """
    converts a string range to a list,
    e.g. "0.5 - 1.0" --> [0.5,1.0]
    """
    # Split the range string into two parts
    str_nums = str_range.split('-')
    # Convert the parts to floats
    min = float(str_nums[0])
    max = float(str_nums[1])
    mid = float(f"{0.5 * (min + max):.3f}")
    # Return a dict
    return {"min": min, "mid": mid, "max": max}


def decompose_covmat(covmat):
    """Given a covmat it return an array sys with shape (ndat,ndat)
    giving ndat correlated systematics for each of the ndat point.
    The original covmat is obtained by doing sys@sys.T"""

    lamb, mat = np.linalg.eig(covmat)
    sys = np.multiply(np.sqrt(lamb), mat)
    return sys


def process_err(error, cv):
    """
    Converts an error given in percentage
    of the central data value in absolute value.

    Note: the d'Agostini prescription for the
    symmetrization of the error does not hold here.
    We follow here the experimental prescription

    Parameters
    ----------

    error : dictionary
            e.g. {'label': 'sys', 'symerror': 0.1%}

    cv : float
        central value

    Returns
    -------
    tuple
        tuple containing two floats

    """
    if "lum" in error:
        # luminosity uncertainty is always symmetric
        sigma = float(error['symerror'].strip('%')) / 100.0 * cv
        return sigma

    elif error['label'] == 'sys':
        if 'asymerror' in error:
            d_p = float(error['asymerror']['plus'].strip('%')) / 100.0 * cv
            d_m = float(error['asymerror']['minus'].strip('%')) / 100.0 * cv

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
            return d_p / np.sqrt(2.0), d_m / np.sqrt(2.0)

        else:
            sigma = float(error['symerror'].strip('%')) / 100.0 * cv
            return sigma / np.sqrt(2.0), -sigma / np.sqrt(2.0)


def HEP_table_to_df(heptable, scenario='nominal'):
    """
    Given hep data table return a pandas
    DataFrame with index given by Ndata,
    columns by the uncertainties and
    np.nan entries

    Parameters
    ----------
    heptable : str
            path to hepdata table

    scenario : 0, 1, 2
            0 = nominal, 1 = stronger, 2 = weaker
    """

    with open(heptable, 'r') as file:
        card = yaml.safe_load(file)

    # list of len ndata. Each entry is dict with
    # keys errors and value
    card = card['dependent_variables'][SCENARIO[scenario]]['values']
    df = pd.DataFrame(index=range(1, len(card) + 1))

    errors = card[0]['errors']
    for i, err in enumerate(errors):
        # luminosity uncertainty, always symmetric
        if (
            (scenario == 'nominal' and i == 67)
            or (scenario == 'stronger' and i == 57)
            or (scenario == 'weaker' and i == 69)
        ):
            df["lum"] = np.nan

        elif err['label'] == 'sys':
            df[f"{err['label']}_plus_{i}"] = np.nan
            df[f"{err['label']}_minus_{i}"] = np.nan

    return df


def fill_df(heptable, scenario='nominal'):
    """
    Fill a data frame with index
    corresponding to dijet mass bins
    and columns to different uncertainties
    Each df is for a fixed rapidity bin.
    """

    with open(heptable, 'r') as file:
        card = yaml.safe_load(file)

    # list of len ndata. Each entry is dict with
    # keys errors and value
    card = card['dependent_variables'][SCENARIO[scenario]]['values']
    df_nan = HEP_table_to_df(heptable, scenario)

    for i, dat in enumerate(card):
        cv = dat['value']

        for j, err in enumerate(dat['errors']):
            if (
                (scenario == 'nominal' and j == 67)
                or (scenario == 'stronger' and j == 57)
                or (scenario == 'weaker' and j == 69)
            ):
                # d_p, d_m = process_err(err,cv)
                # df_nan.loc[df_nan.index == i+1,"lum_plus"] = d_p
                # df_nan.loc[df_nan.index == i+1,"lum_minus"] = d_m
                err["lum"] = "lum"
                sigma = process_err(err, cv)

                df_nan.loc[df_nan.index == i + 1, "lum"] = sigma

            elif err['label'] == 'sys':
                d_p, d_m = process_err(err, cv)
                df_nan.loc[df_nan.index == i + 1, f"{err['label']}_plus_{j}"] = d_p
                df_nan.loc[df_nan.index == i + 1, f"{err['label']}_minus_{j}"] = d_m

    return df_nan
