import yaml

# paths to hepdata table
table = f"rawdata/HEPData-ins2771257-v2-Table_5.yaml"

# functions to filter the relevant information from the hepdata tables
def get_kinematics(table):
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

    hepdata_table = table

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    for pt in input["independent_variables"][0]['values']:
        kin_value = {
            'pT': {'min': pt['low'], 'mid': 0.5 * (pt['low'] + pt['high']), 'max': pt['high']},
            'm_Z2': {'min': None, 'mid': 8317.44, 'max': None},
            'sqrts': {'min': None, 'mid': 13000.0, 'max': None},
        }

        kin.append(kin_value)

    return kin

def get_data_values(table):
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

    hepdata_table = table

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values = input['dependent_variables'][0]['values']

    for value in values:
        data_central.append(value['value'])

    return data_central

'''def get_systematics():
    """ 
    extracts the systematic uncertainties from the HEPData tables
    """
    uncertainties = []

    hepdata_table = f"rawdata/HEPData-ins2771257-v2-Table_4.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    dependent_vars = input['dependent_variables']

    # skip 1st entry as these are central data values
    for dep_var in dependent_vars[1:]:

        name = dep_var['header']['name']
        # skip first 4 and last entry as these central values are not contained
        # in Tabulated Figure 6
        values = [d['value'] for d in dep_var['values'][4:-1]]
        uncertainties.append([{"name": name, "values": values}])

    return uncertainties'''

def get_uncertainties(table):

    hepdata_table = table
    with open(hepdata_table, 'r') as f:
        data = yaml.safe_load(f)
    
    dep_vars = data['dependent_variables'][0]  # Values and uncertainties
    values = dep_vars['values']

    # Collect values of all unique error labels
    labels = set()
    for v in values:
        for err in v['errors']:
            labels.add(err['label'])
    labels = sorted(labels)

    # Build definitions of uncertainties
    definitions = {}
    for label in labels:
        if label == "stat":
            definitions[label] = {
                "description": label,
                "treatment": "ADD",
                "type": "UNCORR"
            }
        elif "LUMI" in label or "lumi" in label.lower():
            definitions[label] = {
                "description": label,
                "treatment": "MULT",
                "type": "ATLASLUMI15"
            }
        else:
            definitions[label] = {
                "description": label,
                "treatment": "ADD",
                "type": "CORR"
            }

    # build bins
    bins = []
    for v in values:
        bin_entry = {}
        for err in v['errors']:
            bin_entry[err['label']] = err['symerror']
        bins.append(bin_entry)

    return {
        "definitions": definitions,
        "bins": bins
    }


def filter_ATLAS_Z1B_13TEV_PTZ_data_kinetic():
    """
    writes data central values and kinematics
    to respective .yaml file
    """
    '''with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]

    # tables for Z->l+l- observable
    tables = metadata["implemented_observables"][0]["tables"]'''

    kin = get_kinematics(table)
    central_values = get_data_values(table)

    data_central_yaml = {"data_central": central_values}
    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    with open("data_ATLAS_Z1B_13TEV_PTB.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics_ATLAS_Z1B_13TEV_PTB.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_ATLAS_Z_13TEV_PT_uncertainties():
    """
    writes uncertainties to respective .yaml file
    """
    uncertainties = get_uncertainties(table)

    # write uncertainties to yaml file
    with open("uncertainties_ATLAS_Z1B_13TEV_PTB.yaml", 'w') as file:
        yaml.dump(uncertainties, file, sort_keys=False)

if __name__ == "__main__":
    # Kinematics and central values:
    filter_ATLAS_Z1B_13TEV_PTZ_data_kinetic()
    # Uncertainties:
    filter_ATLAS_Z_13TEV_PT_uncertainties()



