import yaml


def get_kinematics():
    """
    returns the relevant kinematics values.
    Parameters
    ----------
    table : list
    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables
    Returns
    -------
    list
        list containing the kinematic values for all
        hepdata tables
    """
    kin = []

    hepdata_table = f"rawdata/HEPData-ins2698794-v1-Table_9.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    for yll in input["independent_variables"][0]['values']:
        kin_value = {
            'y': {'min': yll['low'], 'mid': 0.5 * (yll['low'] + yll['high']), 'max': yll['high']},
            'm_Z2': {'min': None, 'mid': 8317.44, 'max': None},
            'sqrts': {'min': None, 'mid': 8000.0, 'max': None},
        }

        kin.append(kin_value)

    return kin


def get_data_values():
    """
    returns the central data.
    Parameters
    ----------
    tables : list
            list that enumerates the table number
    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables
    Returns
    -------
    list
        list containing the central values for all
        hepdata tables
    """

    data_central = []

    hepdata_table = f"rawdata/HEPData-ins2698794-v1-Table_9.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values = input['dependent_variables'][0]['values']

    for value in values:
        data_central.append(value['value'])

    return data_central


def get_systematics():
    """ """
    tot_uncertainties = []
    lumi_uncertainties = []

    hepdata_table = f"rawdata/HEPData-ins2698794-v1-Table_9.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    dependent_vars = input['dependent_variables'][0]

    # skip 1st entry as these are central data values
    for err_values in dependent_vars['values']:

        tot_uncertainties.append(
            [
                {
                    "name": err_values['errors'][0]['label'],
                    "value": err_values['errors'][0]['symerror'],
                }
            ]
        )
        lumi_uncertainties.append(
            [
                {
                    "name": err_values['errors'][1]['label'],
                    "value": err_values['errors'][1]['symerror'],
                }
            ]
        )

    return {'tot': tot_uncertainties, 'lumi': lumi_uncertainties}
