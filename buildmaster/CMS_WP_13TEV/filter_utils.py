import yaml


def get_kinematics(version):
    """
    returns the relevant kinematics values.

    Parameters
    ----------
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

    hepdata_table = f"rawdata/HEPData-ins1810913-v{version}-Figure_A2.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    for y in input["independent_variables"][0]['values']:
        kin_value = {
            'yW': {'min': y['low'], 'mid': 0.5 * (y['low'] + y['high']), 'max': y['high']},
            'mZ_2': {'min': None, 'mid': 6460.5, 'max': None},
            'sqrt_s': {'min': None, 'mid': 13000.0, 'max': None},
        }

        kin.append(kin_value)

    return kin



def get_data_values(version):
    """
    returns the central data.

    Parameters
    ----------
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

    hepdata_table = f"rawdata/HEPData-ins1810913-v{version}-Figure_A2.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values = input['dependent_variables'][0]['values']
      
    for value in values:
        data_central.append(value['value'])

    return data_central


def get_systematics(tables, version):
    """
    """
    uncertainties = []

    hepdata_table = f"rawdata/HEPData-ins1768911-v{version}-Table_{tables[0]}.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    dependent_vars = input['dependent_variables']

    # skip 1st entry as these are central data values
    for dep_var in dependent_vars[1:]:
    
        name = dep_var['header']['name'] 
        values = [d['value'] for d in dep_var['values']]
        uncertainties.append([{"name":name, "values":values}])
    
    return uncertainties

if __name__ == "__main__":
    get_kinematics(version=1)
    get_data_values(version=1)