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

    hepdata_table = f"rawdata/HEPData-ins1467454-v1-Table_2.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    for indep_var1, indep_var2 in zip(
        input["independent_variables"][1]['values'], input["independent_variables"][0]['values']
    ):

        kin_value = {
            'y': {
                'min': indep_var1['low'],
                'mid': 0.5 * (indep_var1['low'] + indep_var1['high']),
                'max': indep_var1['high'],
            },
            'M2': {
                'min': indep_var2['low'] ** 2,
                'mid': (0.5 * (indep_var2['low'] + indep_var2['high'])) ** 2,
                'max': indep_var2['high'] ** 2,
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

    hepdata_table = f"rawdata/HEPData-ins1467454-v1-Table_2.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values = input['dependent_variables'][0]['values']

    for value, mass_bins in zip(values, input["independent_variables"][0]['values']):
        # store data central and normalize to match applgrid predictions
        data_central.append(value['value'] * 2 * (mass_bins['high'] - mass_bins['low']))

    return data_central


def get_systematics():
    """ """

    uncertainties = []

    hepdata_table = f"rawdata/HEPData-ins1467454-v1-Table_2.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    # loop over systematics
    for unc_labels in input['dependent_variables'][0]['values'][0]['errors']:

        name = f"{unc_labels['label']}"
        values = []

        # loop over data points
        for unc, mass_bins in zip(
            input['dependent_variables'][0]['values'], input["independent_variables"][0]['values']
        ):
            err = unc['errors']
            # normalize the central values
            cv = unc['value'] * 2 * (mass_bins['high'] - mass_bins['low'])

            # convert unc from TeV to GeV
            for e in err:
                if e['label'] == name:

                    if 'asymerror' in e:
                        # the errors are actually symmetric.
                        values.append(float(e['asymerror']['plus'][:-1]) * cv / 100.0)

                    else:
                        values.append(float(e['symerror'][:-1]) * cv / 100.0)

        uncertainties.append([{"name": name, "values": values}])

    return uncertainties
