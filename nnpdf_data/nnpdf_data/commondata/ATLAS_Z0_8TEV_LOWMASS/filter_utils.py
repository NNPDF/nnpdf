"""
This module contains helper functions that are used to extract the uncertainties, kinematics and data values 
from the rawdata files.
"""

import yaml

ATLAS_LUMI_UNC = 0.018

def get_kinematics():
    """
    returns the kinematics in the form of a list of dictionaries.
    """
    kin = []

    hepdata_table = f"rawdata/HEPData-ins1630886-v3-Table_5.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    for indep_var1, indep_var2 in zip(
        input["independent_variables"][1]['values'], input["independent_variables"][2]['values']
    ):

        kin_value = {
            'abs_y': {
                'min': indep_var1['low'],
                'mid': 0.5 * (indep_var1['low'] + indep_var1['high']),
                'max': indep_var1['high'],
            },
            'm_Z2': {
                'min': indep_var2['low']**2,
                'mid': (0.5 * (indep_var2['low'] + indep_var2['high']))**2,
                'max': indep_var2['high']**2,
            },
            'sqrts': {'min': None, 'mid': 8000.0, 'max': None},
        }

        kin.append(kin_value)

    return kin


def get_data_values():
    """
    returns the central data values in the form of a list.
    """

    data_central = []

    hepdata_table = f"rawdata/HEPData-ins1630886-v3-Table_5.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values = input['dependent_variables'][0]['values']

    for value in values:
        # store data central and convert the units
        data_central.append(value['value'] * 1000)

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
