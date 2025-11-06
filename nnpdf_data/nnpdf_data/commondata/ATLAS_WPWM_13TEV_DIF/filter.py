"""
This file contains the piece of code needed to implement the ATLAS double-
differential CCDY at high transverse masses measurement at 13 TeV. 
"""

import yaml
import numpy as np

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def symmetrize_errors(delta_plus, delta_minus):
    r"""Compute the symmetrized uncertainty and the shift in data point.

    Parameters
    ----------
    delta_plus : float
        The top/plus uncertainty with sign
    delta_minus : float
        The bottom/minus uncertainty with sign

    Returns
    -------
    se_delta : float
        The value to be added to the data point
    se_sigma : float
        The symmetrized uncertainty to be used in commondata

    """
    semi_diff = (delta_plus - delta_minus) / 2
    average = (delta_plus + delta_minus) / 2
    se_delta = semi_diff
    se_sigma = np.sqrt(average * average + 2 * semi_diff * semi_diff)
    return se_delta, se_sigma


def get_tables(observable=None):
    """
    Get the Hepdata tables, given the tables and version specified in metadata
    """
    prefix = "rawdata/HEPData-ins2895869"
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]

    if observable == "WPWM_DIF_LEP":
        tables = ["lep_physical_plus_mtw", "lep_physical_minus_mtw"]
    elif observable == "WPWM_DIF_MUON":
        tables = ["muo_plus_mtw", "muo_minus_mtw"]
    elif observable == "WPWM_DDIF_LEP":
        tables = []
        for i in range(5):
            tables.append(f"lep_physical_plus_absetamtw_mtw{i}")
        for i in range(5):
            tables.append(f"lep_physical_minus_absetamtw_mtw{i}")
    elif observable == "WPWM_DDIF_MUON":
        tables = []
        for i in range(5):
            tables.append(f"muo_plus_absetamtw_mtw{i}")
        for i in range(5):
            tables.append(f"muo_minus_absetamtw_mtw{i}")
    else:
        print("Observable not implemented.")
        print("Choose one of the following observables:")
        print("- WPWM_DIF_LEP")
        print("- WPWM_DIF_MUON")
        print("- WPWM_DDIF_LEP")
        print("- WPWM_DDIF_MUON")

    hepdata_tables = []

    for table in tables:
        hepdata_tables.append(f"{prefix}-v{version}-{table}.yaml")
    return hepdata_tables


def get_data(observable=None):

    data_central = []
    kinematics = []
    uncertainties = []

    hepdata_tables = get_tables(observable)

    table_num = 0
    mt_bins = [
        [200, 300],
        [300, 425],
        [425, 600],
        [600, 900],
        [900, 2000],
        [200, 300],
        [300, 425],
        [425, 600],
        [600, 900],
        [900, 2000],
    ]

    for table in hepdata_tables:

        with open(table, "r") as f:
            input = yaml.safe_load(f)

        data_values = input["dependent_variables"][0]["values"]

        for data_value in data_values:
            data_central.append(data_value["value"])

        kin_values = input["independent_variables"][0]["values"]
        kin_name = input["independent_variables"][0]["header"]["name"]

        if kin_name == "$m_T^W$":  # Single differential
            kin_label = "m_T^W"
            for kin_value in kin_values:
                kin = {
                    'm_T^W': {
                        'min': kin_value['low'],
                        'mid': 0.5 * (kin_value['low'] + kin_value['high']),
                        'max': kin_value['high'],
                    },
                    'm_W2': {'min': None, 'mid': 6.46174823e03, 'max': None},
                }
                kinematics.append(kin)

        if kin_name == "$|\eta|$":  # Double differential
            kin_label = "abs_eta"
            mt_low, mt_high = mt_bins[table_num]
            for kin_value in kin_values:
                kin = {
                    'abs_eta': {
                        'min': kin_value['low'],
                        'mid': 0.5 * (kin_value['low'] + kin_value['high']),
                        'max': kin_value['high'],
                    },
                    'm_T^W': {'min': mt_low, 'mid': 0.5 * (mt_low + mt_high), 'max': mt_high},
                    'm_W2': {'min': None, 'mid': 6.46174823e03, 'max': None},
                }
                kinematics.append(kin)
            table_num += 1

    ndata = len(data_central)
    syst_dict = {}

    lumi_unc = 0.83
    lumi_uncs = []
    for data in data_central:
        lumi_uncs.append(data * lumi_unc / 100.0)

    syst_dict["LUMI"] = np.array(lumi_uncs)
    value_id = 0
    deltas = np.zeros(ndata)
    for table in hepdata_tables:
        with open(table, "r") as f:
            input = yaml.safe_load(f)

        values = input["dependent_variables"][0]["values"]

        for point_id, point in enumerate(values):
            for err in point["errors"]:
                label = err["label"]
                if 'asymerror' in err:
                    minus = err['asymerror']['minus']
                    plus = err['asymerror']['plus']

                elif 'symerror' in err:
                    minus = plus = err['symerror']
                else:
                    raise ValueError(f"Unknown error type in {hepdata_table} for point {point_idx}")

                if label not in syst_dict:
                    syst_dict[label] = np.zeros(ndata)

                delta, symm_err = symmetrize_errors(plus, minus)
                syst_dict[label][value_id] = symm_err
                deltas[value_id] += delta

            value_id += 1

    sys_list = []
    for label, values in syst_dict.items():
        sys_list.append({"name": label, "values": values.tolist()})

    data_central = (data_central + deltas).tolist()

    return data_central, kinematics, sys_list


def filter_ATLAS_WPWM_13TEV_DIF(observable=None):
    central_values, kinematics, uncertainties = get_data(observable)

    data_central_yaml = {"data_central": central_values}

    kinematics_yaml = {"bins": kinematics}

    if "MUON" in observable:
        # Muon channel
        treatment = {
            "LUMI": "MULT",
            "Data stat. unc.": "ADD",
            "Sig. stat. unc.": "ADD",
            "Bkg. stat. unc.": "ADD",
            "Alternative MC unf. unc.": "ADD",
            "Stat. unc.": "ADD",
            "Uncor. syst. unc.": "ADD",
            "Others": "MULT",
        }
        correlation = {
            "LUMI": "ATLASLUMIRUNII",
            "" "Data stat. unc.": "UNCORR",
            "Sig. stat. unc.": "UNCORR",
            "Bkg. stat. unc.": "UNCORR",
            "Alternative MC unf. unc.": "UNCORR",
            "Stat. unc.": "UNCORR",
            "Uncor. syst. unc.": "UNCORR",
            "Others": "CORR",
        }
    else:
        # Lepton channel
        treatment = {
            "LUMI": "MULT",
            "Data stat. unc.": "ADD",
            "Sig. stat. unc.": "ADD",
            "Bkg. stat. unc.": "ADD",
            "Alternative MC unf. unc.": "ADD",
            "Stat. unc.": "ADD",
            "Uncor. syst. unc.": "ADD",
            "Others": "ADD",
        }
        correlation = {
            "LUMI": "ATLASLUMIRUNII",
            "" "Data stat. unc.": "UNCORR",
            "Sig. stat. unc.": "UNCORR",
            "Bkg. stat. unc.": "UNCORR",
            "Alternative MC unf. unc.": "UNCORR",
            "Stat. unc.": "UNCORR",
            "Uncor. syst. unc.": "UNCORR",
            "Others": "CORR",
        }

    definitions = {}

    errors_yaml = []
    for unc in uncertainties:

        if unc["name"] in treatment:
            definitions[unc["name"]] = {
                "description": unc["name"],
                "treatment": treatment[unc["name"]],
                "correlation": correlation[unc["name"]],
            }
        else:
            definitions[unc["name"]] = {
                "description": unc["name"],
                "treatment": treatment["Others"],
                "correlation": correlation["Others"],
            }

    ndata = len(central_values)
    for i in range(ndata):
        errors = {}
        for unc in uncertainties:
            errors[unc["name"]] = float(unc["values"][i])

        errors_yaml.append(errors)

    uncertainties_yaml = {"definitions": definitions, "bins": errors_yaml}

    with open("data_" + observable + ".yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics_" + observable + ".yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    with open("uncertainties_" + observable + ".yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_WPWM_13TEV_DIF("WPWM_DIF_LEP")
    filter_ATLAS_WPWM_13TEV_DIF("WPWM_DIF_MUON")
    filter_ATLAS_WPWM_13TEV_DIF("WPWM_DDIF_LEP")
    filter_ATLAS_WPWM_13TEV_DIF("WPWM_DDIF_MUON")
