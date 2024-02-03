import yaml


def get_kinematics(table, version):
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

    hepdata_table = f"rawdata/HEPData-ins1768911-v3-Tabulated_Figure_6.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    for pt in input["independent_variables"][0]['values']:
        kin_value = {
            'pT': {'min': pt['low'], 'mid': 0.5 * (pt['low'] + pt['high']), 'max': pt['high']},
            'mZ_2': {'min': None, 'mid': 8317.44, 'max': None},
            'sqrt_s': {'min': None, 'mid': 13000.0, 'max': None},
        }

        kin.append(kin_value)

    return kin


def get_data_values(tables, version):
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

    hepdata_table = f"rawdata/HEPData-ins1768911-v3-Tabulated_Figure_6.yaml"

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
        # skip first 4 and last entry as these central values are not contained
        # in Tabulated Figure 6
        values = [d['value'] for d in dep_var['values'][4:-1]]
        uncertainties.append([{"name":name, "values":values}])
    
    return uncertainties

if __name__ == "__main__":
    # get_kinematics(tables=[7], version=3)
    # get_data_values(tables=[7], version=3)
    get_systematics(tables=[5], version=3)
