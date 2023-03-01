# implemented by Tanishq Sharma

import yaml
from utils import covMat_to_artUnc as cta
from utils import percentage_to_absolute_num as ptan

def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    ndata_dSig_dpTt = metadata['implemented_observables'][0]['ndata']
    ndata_dSig_dyt = metadata['implemented_observables'][1]['ndata']
    ndata_dSig_dyttBar = metadata['implemented_observables'][2]['ndata']
    ndata_dSig_dmttBar = metadata['implemented_observables'][3]['ndata']

    data_central_dSig_dpTt = []
    kin_dSig_dpTt = []
    error_dSig_dpTt = []
    data_central_dSig_dyt = []
    kin_dSig_dyt = []
    error_dSig_dyt = []
    data_central_dSig_dyttBar = []
    kin_dSig_dyttBar = []
    error_dSig_dyttBar = []
    data_central_dSig_dmttBar = []
    kin_dSig_dmttBar = []
    error_dSig_dmttBar = []

    covMatArray_dSig_dpTt = []
    covMatArray_dSig_dyt = []
    covMatArray_dSig_dyttBar = []
    covMatArray_dSig_dmttBar = []

# dSig_dpTt data

    hepdata_tables="rawdata/Table15.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table16.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)

    systematics_breakdown="rawdata/Table17.yaml"
    with open(systematics_breakdown, 'r') as file3:
        input3 = yaml.safe_load(file3)

    for i in range(ndata_dSig_dpTt*ndata_dSig_dpTt):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dpTt.append(covMatEl)
    artUncMat_dSig_dpTt = cta(ndata_dSig_dpTt, covMatArray_dSig_dpTt)

    sqrt_s = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        pT_t_min = input['independent_variables'][0]['values'][i]['low']
        pT_t_mid = input['independent_variables'][1]['values'][i]['value']
        pT_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        error_value['stat'] = 0
        # error_value['sys'] = values[i]['errors'][1]['symerror']
        for j in range(ndata_dSig_dpTt):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dpTt[i][j])
        data_central_value = values[i]['value']
        for j in range(11):
            error_value[input3['independent_variables'][0]['values'][j]['value']] = ptan(input3['dependent_variables'][i]['values'][j]['value'], data_central_value)
        kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'pT_t': {'min': pT_t_min, 'mid': pT_t_mid, 'max': pT_t_max}}
        data_central_dSig_dpTt.append(data_central_value)
        kin_dSig_dpTt.append(kin_value)
        error_dSig_dpTt.append(error_value)

    error_definition_dSig_dpTt = {}
    error_definition_dSig_dpTt['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # error_definition_dSig_dpTt['sys'] = {'description': 'total systematic uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dpTt):
        error_definition_dSig_dpTt['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}
    for i in range(11):
        error_definition_dSig_dpTt[input3['independent_variables'][0]['values'][i]['value']] = {'definition': 'systematic uncertainty- '+str(input3['independent_variables'][0]['values'][i]['value']), 'treatment': 'MULT', 'type': 'CORR'}

    data_central_dSig_dpTt_yaml = {'data_central': data_central_dSig_dpTt}
    kinematics_dSig_dpTt_yaml = {'bins': kin_dSig_dpTt}
    uncertainties_dSig_dpTt_yaml = {'definitions': error_definition_dSig_dpTt, 'bins': error_dSig_dpTt}

    with open('data_dSig_dpTt.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dpTt_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dpTt.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dpTt_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dpTt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_yaml, file, sort_keys=False)

# dSig_dyt data

    hepdata_tables="rawdata/Table21.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table22.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)

    systematics_breakdown="rawdata/Table23.yaml"
    with open(systematics_breakdown, 'r') as file3:
        input3 = yaml.safe_load(file3)

    for i in range(ndata_dSig_dyt*ndata_dSig_dyt):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyt.append(covMatEl)
    artUncMat_dSig_dyt = cta(ndata_dSig_dyt, covMatArray_dSig_dyt)

    sqrt_s = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_mid = input['independent_variables'][1]['values'][i]['value']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        error_value['stat'] = 0
        # error_value['sys'] = values[i]['errors'][1]['symerror']
        for j in range(ndata_dSig_dyt):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dyt[i][j])
        data_central_value = values[i]['value']
        for j in range(11):
            error_value[input3['independent_variables'][0]['values'][j]['value']] = ptan(input3['dependent_variables'][i]['values'][j]['value'], data_central_value)
        kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'y_t': {'min': y_t_min, 'mid': y_t_mid, 'max': y_t_max}}
        data_central_dSig_dyt.append(data_central_value)
        kin_dSig_dyt.append(kin_value)
        error_dSig_dyt.append(error_value)

    error_definition_dSig_dyt = {}
    error_definition_dSig_dyt['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # error_definition_dSig_dyt['sys'] = {'description': 'total systematic uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dyt):
        error_definition_dSig_dyt['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}
    for i in range(11):
        error_definition_dSig_dyt[input3['independent_variables'][0]['values'][i]['value']] = {'definition': 'systematic uncertainty- '+str(input3['independent_variables'][0]['values'][i]['value']), 'treatment': 'MULT', 'type': 'CORR'}

    data_central_dSig_dyt_yaml = {'data_central': data_central_dSig_dyt}
    kinematics_dSig_dyt_yaml = {'bins': kin_dSig_dyt}
    uncertainties_dSig_dyt_yaml = {'definitions': error_definition_dSig_dyt, 'bins': error_dSig_dyt}

    with open('data_dSig_dyt.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dyt_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyt.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dyt_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_yaml, file, sort_keys=False)

# dSig_dyttBar data

    hepdata_tables="rawdata/Table36.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table37.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)

    systematics_breakdown="rawdata/Table38.yaml"
    with open(systematics_breakdown, 'r') as file3:
        input3 = yaml.safe_load(file3)

    for i in range(ndata_dSig_dyttBar*ndata_dSig_dyttBar):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyttBar.append(covMatEl)
    artUncMat_dSig_dyttBar = cta(ndata_dSig_dyttBar, covMatArray_dSig_dyttBar)

    sqrt_s = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_mid = input['independent_variables'][1]['values'][i]['value']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        error_value['stat'] = 0
        # error_value['sys'] = values[i]['errors'][1]['symerror']
        for j in range(ndata_dSig_dyttBar):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dyttBar[i][j])
        data_central_value = values[i]['value']
        for j in range(11):
            error_value[input3['independent_variables'][0]['values'][j]['value']] = ptan(input3['dependent_variables'][i]['values'][j]['value'], data_central_value)
        kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'y_ttBar': {'min': y_ttBar_min, 'mid': y_ttBar_mid, 'max': y_ttBar_max}}
        data_central_dSig_dyttBar.append(data_central_value)
        kin_dSig_dyttBar.append(kin_value)
        error_dSig_dyttBar.append(error_value)

    error_definition_dSig_dyttBar = {}
    error_definition_dSig_dyttBar['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # error_definition_dSig_dyttBar['sys'] = {'description': 'total systematic uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dyttBar):
        error_definition_dSig_dyttBar['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}
    for i in range(11):
        error_definition_dSig_dyttBar[input3['independent_variables'][0]['values'][i]['value']] = {'definition': 'systematic uncertainty- '+str(input3['independent_variables'][0]['values'][i]['value']), 'treatment': 'MULT', 'type': 'CORR'}

    data_central_dSig_dyttBar_yaml = {'data_central': data_central_dSig_dyttBar}
    kinematics_dSig_dyttBar_yaml = {'bins': kin_dSig_dyttBar}
    uncertainties_dSig_dyttBar_yaml = {'definitions': error_definition_dSig_dyttBar, 'bins': error_dSig_dyttBar}

    with open('data_dSig_dyttBar.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyttBar.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_yaml, file, sort_keys=False)

# dSig_dmttBar data

    hepdata_tables="rawdata/Table39.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table40.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)

    systematics_breakdown="rawdata/Table41.yaml"
    with open(systematics_breakdown, 'r') as file3:
        input3 = yaml.safe_load(file3)

    for i in range(ndata_dSig_dmttBar*ndata_dSig_dmttBar):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dmttBar.append(covMatEl)
    artUncMat_dSig_dmttBar = cta(ndata_dSig_dmttBar, covMatArray_dSig_dmttBar)

    sqrt_s = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_mid = input['independent_variables'][1]['values'][i]['value']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        error_value['stat'] = 0
        # error_value['sys'] = values[i]['errors'][1]['symerror']
        for j in range(ndata_dSig_dmttBar):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dmttBar[i][j])
        data_central_value = values[i]['value']
        for j in range(11):
            error_value[input3['independent_variables'][0]['values'][j]['value']] = ptan(input3['dependent_variables'][i]['values'][j]['value'], data_central_value)
        kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'm_ttBar': {'min': m_ttBar_min, 'mid': m_ttBar_mid, 'max': m_ttBar_max}}
        data_central_dSig_dmttBar.append(data_central_value)
        kin_dSig_dmttBar.append(kin_value)
        error_dSig_dmttBar.append(error_value)

    error_definition_dSig_dmttBar = {}
    error_definition_dSig_dmttBar['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # error_definition_dSig_dmttBar['sys'] = {'description': 'total systematic uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dmttBar):
        error_definition_dSig_dmttBar['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}
    for i in range(11):
        error_definition_dSig_dmttBar[input3['independent_variables'][0]['values'][i]['value']] = {'definition': 'systematic uncertainty- '+str(input3['independent_variables'][0]['values'][i]['value']), 'treatment': 'MULT', 'type': 'CORR'}

    data_central_dSig_dmttBar_yaml = {'data_central': data_central_dSig_dmttBar}
    kinematics_dSig_dmttBar_yaml = {'bins': kin_dSig_dmttBar}
    uncertainties_dSig_dmttBar_yaml = {'definitions': error_definition_dSig_dmttBar, 'bins': error_dSig_dmttBar}

    with open('data_dSig_dmttBar.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_yaml, file, sort_keys=False)

processData()
