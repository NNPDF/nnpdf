import yaml

def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    tables_r04 = metadata['implemented_observables'][0]['tables']
    tables_r07 = metadata['implemented_observables'][1]['tables']

    data_central_r04 = []
    kin_r04 = []
    error_r04 = []
    data_central_r07 = []
    kin_r07 = []
    error_r07 = []

# r04 data

    for i in tables_r04:
        if i == 1:
            y_min = 0
            y_max = 0.5
        elif i == 2:
            y_min = 0.5
            y_max = 1
        elif i == 3:
            y_min = 1
            y_max = 1.5
        elif i == 4:
            y_min = 1.5
            y_max = 2
        hepdata_tables="rawdata/ak4_xsec_ybin"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
        
        values = input['dependent_variables'][0]['values']
        sqrts = 13000

        for j in range(len(values)):
            data_central_value = values[j]['value']
            data_central_r04.append(data_central_value)
            pT_min = input['independent_variables'][0]['values'][j]['low']
            pT_max = input['independent_variables'][0]['values'][j]['high']
            kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'pT': {'min': pT_min, 'mid': None, 'max': pT_max}, 'y': {'min': y_min, 'mid': None, 'max': y_max}}
            kin_r04.append(kin_value)
            error_value = {}
            error_value['all uncorr. unc.'] = values[j]['errors'][0]['symerror']
            for k in range(1, len(values[j]['errors'])):
                error_label = values[j]['errors'][k]['label']
                error_value[error_label] = values[j]['errors'][k]['symerror']
            error_r04.append(error_value)

    hepdata_tables="rawdata/ak4_xsec_ybin1.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    error_definition_r04 = {}
    error_definition_r04['all uncorr. unc.'] = {'description': 'all uncorrelated uncertainties', 'treatment': 'ADD', 'type': 'UNCORR'}
    for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
        error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
        error_definition_r04[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    
    data_central_r04_yaml = {'data_central': data_central_r04}
    kinematics_r04_yaml = {'bins': kin_r04}
    uncertainties_r04_yaml = {'definitions': error_definition_r04, 'bins': error_r04}

    with open('data_r04.yaml', 'w') as file:
         yaml.dump(data_central_r04_yaml, file, sort_keys=False)

    with open('kinematics_r04.yaml', 'w') as file:
         yaml.dump(kinematics_r04_yaml, file, sort_keys=False)

    with open('uncertainties_r04.yaml', 'w') as file:
        yaml.dump(uncertainties_r04_yaml, file, sort_keys=False)

# r07 data

    for i in tables_r07:
        if i == 1:
            y_min = 0
            y_max = 0.5
        elif i == 2:
            y_min = 0.5
            y_max = 1
        elif i == 3:
            y_min = 1
            y_max = 1.5
        elif i == 4:
            y_min = 1.5
            y_max = 2
        hepdata_tables="rawdata/ak7_xsec_ybin"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
        
        values = input['dependent_variables'][0]['values']
        sqrts = 13000

        for j in range(len(values)):
            data_central_value = values[j]['value']
            data_central_r07.append(data_central_value)
            pT_min = input['independent_variables'][0]['values'][j]['low']
            pT_max = input['independent_variables'][0]['values'][j]['high']
            kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'pT': {'min': pT_min, 'mid': None, 'max': pT_max}, 'y': {'min': y_min, 'mid': None, 'max': y_max}}
            kin_r07.append(kin_value)
            error_value = {}
            error_value['all uncorr. unc.'] = values[j]['errors'][0]['symerror']
            for k in range(1, len(values[j]['errors'])):
                error_label = values[j]['errors'][k]['label']
                error_value[error_label] = values[j]['errors'][k]['symerror']
            error_r07.append(error_value)

    hepdata_tables="rawdata/ak7_xsec_ybin1.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    error_definition_r07 = {}
    error_definition_r07['all uncorr. unc.'] = {'description': 'all uncorrelated uncertainties', 'treatment': 'ADD', 'type': 'UNCORR'}
    for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
        error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
        error_definition_r07[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    
    data_central_r07_yaml = {'data_central': data_central_r07}
    kinematics_r07_yaml = {'bins': kin_r07}
    uncertainties_r07_yaml = {'definitions': error_definition_r07, 'bins': error_r07}

    with open('data_r07.yaml', 'w') as file:
         yaml.dump(data_central_r07_yaml, file, sort_keys=False)

    with open('kinematics_r07.yaml', 'w') as file:
         yaml.dump(kinematics_r07_yaml, file, sort_keys=False)

    with open('uncertainties_r07.yaml', 'w') as file:
        yaml.dump(uncertainties_r07_yaml, file, sort_keys=False)

processData()
