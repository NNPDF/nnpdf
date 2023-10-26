import yaml
from validphys.commondata_utils import percentage_to_absolute as pta

def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    tables = metadata['implemented_observables'][0]['tables']
    tables_norm = metadata['implemented_observables'][1]['tables']

    ndata = metadata['implemented_observables'][0]['ndata']
    ndata_norm = metadata['implemented_observables'][1]['ndata']

    data_central = []
    kin = []
    error = []
    data_central_norm = []
    kin_norm = []
    error_norm = []

# jet data

    hepdata_tables="rawdata/Table"+str(tables[0])+".yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    sqrt_s = float(input['dependent_variables'][0]['qualifiers'][3]['value'])
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = values[i]['value']
        data_central.append(data_central_value)
        Q2_max = input['independent_variables'][0]['values'][i]['high']
        Q2_min = input['independent_variables'][0]['values'][i]['low']
        pT_max = input['independent_variables'][1]['values'][i]['high']
        pT_min = input['independent_variables'][1]['values'][i]['low']
        kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'Q2': {'min': Q2_min, 'mid': None, 'max': Q2_max},'pT': {'min': pT_min, 'mid': None, 'max': pT_max}}
        kin.append(kin_value)
        error_value = {}
        error_value['stat'] = pta(values[i]['errors'][0]['symerror'], data_central_value)
        error_value['sys'] = pta(values[i]['errors'][1]['symerror'], data_central_value)
        error.append(error_value)

    error_definition = {'stat':{'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}, 'sys':{'description': 'total systematic uncertainty', 'treatment':'MULT' , 'type': 'CORR'}}

    data_central_yaml = {'data_central': data_central}
    kinematics_yaml = {'bins': kin}
    uncertainties_yaml = {'definitions': error_definition, 'bins': error}

    with open('data.yaml', 'w') as file:
         yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
         yaml.dump(kinematics_yaml, file, sort_keys=False)

    with open('uncertainties.yaml', 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)

 # jet_norm data

    hepdata_tables="rawdata/Table"+str(tables_norm[0])+".yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    sqrt_s = float(input['dependent_variables'][0]['qualifiers'][3]['value'])
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = values[i]['value']
        data_central_norm.append(data_central_value)
        Q2_max = input['independent_variables'][0]['values'][i]['high']
        Q2_min = input['independent_variables'][0]['values'][i]['low']
        pT_max = input['independent_variables'][1]['values'][i]['high']
        pT_min = input['independent_variables'][1]['values'][i]['low']
        kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'Q2': {'min': Q2_min, 'mid': None, 'max': Q2_max},'pT': {'min': pT_min, 'mid': None, 'max': pT_max}}
        kin_norm.append(kin_value)
        error_value = {}
        error_value['stat'] = pta(values[i]['errors'][0]['symerror'], data_central_value)
        error_value['sys'] = pta(values[i]['errors'][1]['symerror'], data_central_value)
        error_norm.append(error_value)

    error_definition_norm = {'stat':{'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}, 'sys':{'description': 'total systematic uncertainty', 'treatment':'MULT' , 'type': 'CORR'}}

    data_central_norm_yaml = {'data_central': data_central_norm}
    kinematics_norm_yaml = {'bins': kin_norm}
    uncertainties_norm_yaml = {'definitions': error_definition_norm, 'bins': error_norm}

    with open('data_norm.yaml', 'w') as file:
         yaml.dump(data_central_norm_yaml, file, sort_keys=False)

    with open('kinematics_norm.yaml', 'w') as file:
         yaml.dump(kinematics_norm_yaml, file, sort_keys=False)

    with open('uncertainties_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_norm_yaml, file, sort_keys=False)

processData()