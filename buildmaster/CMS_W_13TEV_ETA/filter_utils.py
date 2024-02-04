import yaml


def get_kinematics(version, figure):
    """
    returns the relevant kinematics values.

    Parameters
    ----------
    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    figure : str
        eg. 17a or 17b

    Returns
    -------
    list
        list containing the kinematic values for all
        hepdata tables
    """
    kin = []

    hepdata_table = f"rawdata/HEPData-ins1810913-v{version}-Figure_{figure}.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    for eta in input["independent_variables"][0]['values']:
        kin_value = {
            'eta': {'min': eta['low'], 'mid': 0.5 * (eta['low'] + eta['high']), 'max': eta['high']},
            'mZ_2': {'min': None, 'mid': 6460.5, 'max': None},
            'sqrt_s': {'min': None, 'mid': 13000.0, 'max': None},
        }

        kin.append(kin_value)

    return kin


def get_data_values(version, figure):
    """
    returns the central data.

    Parameters
    ----------
    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    figure : str
        eg. 17a or 17b

    Returns
    -------
    list
        list containing the central values for all
        hepdata tables

    """

    data_central = []

    hepdata_table = f"rawdata/HEPData-ins1810913-v{version}-Figure_{figure}.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values = input['dependent_variables'][0]['values']

    for value in values:
        data_central.append(value['value'])

    return data_central


def get_systematics(version, figure):
    """
    Get the systematics from the hepdata tables.

    Parameters
    ----------
    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    figure : str
        string indicating the figure number for the impacts
        can be A23a or A23b

    Returns
    -------
    list
        list containing the systematics
    """
    uncertainties = []

    hepdata_table = f"rawdata/HEPData-ins1810913-v{version}-Impacts_Figure_{figure}.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    dependent_vars = input['dependent_variables']

    for dep_var in dependent_vars:

        name = dep_var['header']['name']
        values = [d['value'] for d in dep_var['values']]
        uncertainties.append([{"name": name, "values": values}])

    return uncertainties


if __name__ == "__main__":
    get_systematics(version=1, figure="A23a")
