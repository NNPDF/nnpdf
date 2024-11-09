"""
This module contains helper functions that are used to extract the uncertainties, kinematics and data values 
from the rawdata files.
"""

import yaml


def get_kinematics():
    """
    returns the kinematics in the form of a list of dictionaries.
    """
    kin = []

    mw2 = 80.385**2

    for i in range(2):

        kin_value = {
            'k1': {'min': None, 'mid': 0.0, 'max': None},
            'M2': {'min': None, 'mid': mw2, 'max': None},
            'sqrts': {'min': None, 'mid': 13000.0, 'max': None},
        }

        kin.append(kin_value)

    return kin


def get_data_values():
    """
    returns the central data values in the form of a list.
    """
    hepdata_table_wp, hepdata_table_wm = (
        "rawdata/HEPData-ins1436497-v1-Table_8.yaml",
        "rawdata/HEPData-ins1436497-v1-Table_9.yaml",
    )

    with open(hepdata_table_wp, 'r') as file:
        input_wp = yaml.safe_load(file)

    with open(hepdata_table_wm, 'r') as file:
        input_wm = yaml.safe_load(file)

    values_wm = input_wm['dependent_variables'][0]['values']
    values_wp = input_wp['dependent_variables'][0]['values']

    data_central = [values_wm[0]['value'] * 1000000, values_wp[0]['value'] * 1000000]

    return data_central


def get_systematics(version=3):
    """ """

    uncertainties = []

    hepdata_table = f"rawdata/HEPData-ins1630886-v{version}-Table_5.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    # loop over systematics
    for unc_labels in input['dependent_variables'][0]['values'][0]['errors']:

        name = f"{unc_labels['label']}"
        values = []

        # loop over data points
        for unc in input['dependent_variables'][0]['values']:
            err = unc['errors']
            # convert unc from TeV to GeV
            for e in err:
                if e['label'] == name:
                    if name == 'Lumi:M':
                        values.append(e['symerror'] * unc['value'] * 1000)
                    else:
                        values.append(e['symerror'] * 1000)

        uncertainties.append([{"name": name, "values": values}])

    # # Luminosity uncertainty is 1.8 % of the central value (see https://inspirehep.net/literature/1630886)
    if version == 3:  # in version 1 Lumi is included in the hepdata file already
        name = "ATLAS_LUMI"
        values = [ATLAS_LUMI_UNC * val for val in get_data_values()]
        uncertainties.append([{"name": name, "values": values}])
    return uncertainties
