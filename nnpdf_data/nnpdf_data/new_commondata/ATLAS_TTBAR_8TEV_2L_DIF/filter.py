import yaml
from math import sqrt
import numpy as np
from numpy.linalg import eig

def cta(ndata, covmat_list, no_of_norm_mat=0):
    
    epsilon = -0.0000000001
    neg_eval_count = 0
    psd_check = True
    covmat = np.zeros((ndata, ndata))
    artunc = np.zeros((ndata, ndata))
    for i in range(len(covmat_list)):
        a = i // ndata
        b = i % ndata
        covmat[a][b] = covmat_list[i]
    eigval, eigvec = eig(covmat)
    for j in range(len(eigval)):
        if eigval[j] < epsilon:
            psd_check = False
        elif eigval[j] > epsilon and eigval[j] <= 0:
            neg_eval_count = neg_eval_count + 1
            if neg_eval_count == (no_of_norm_mat + 1):
                psd_check = False
        elif eigval[j] > 0:
            continue
    if psd_check == False:
        raise ValueError('The covariance matrix is not positive-semidefinite')
    else:
        for i in range(ndata):
            for j in range(ndata):
                if eigval[j] < 0:
                    continue
                else:
                    artunc[i][j] = eigvec[i][j] * sqrt(eigval[j]) 
    return artunc.tolist()

def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    data_central_dSig_dmttBar = []
    kin_dSig_dmttBar = []
    error_dSig_dmttBar = []
    data_central_dSig_dmttBar_norm = []
    kin_dSig_dmttBar_norm = []
    error_dSig_dmttBar_norm = []
    data_central_dSig_dyttBar = []
    kin_dSig_dyttBar = []
    error_dSig_dyttBar = []
    data_central_dSig_dyttBar_norm = []
    kin_dSig_dyttBar_norm = []
    error_dSig_dyttBar_norm = []

    covMatArray_dSig_dmttBar = []
    covMatArray_dSig_dmttBar_norm = []
    covMatArray_dSig_dyttBar = []
    covMatArray_dSig_dyttBar_norm = []

# dSig_dmttBar data

    hepdata_tables="rawdata/Table_10.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table_22.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    
    sqrts = float(input['dependent_variables'][0]['qualifiers'][1]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']
    
    for i in range(len(values)*len(values)):
        covMatArray_dSig_dmttBar.append(input2['dependent_variables'][0]['values'][i]['value'])
    artUnc_dSig_dmttBar = cta(len(values), covMatArray_dSig_dmttBar)

    for i in range(len(values)):
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(len(values)):
            error_value['ArtUnc_'+str(j+1)] = artUnc_dSig_dmttBar[i][j]
        data_central_value = values[i]['value']
        data_central_dSig_dmttBar.append(data_central_value)
        error_dSig_dmttBar.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max}}
        kin_dSig_dmttBar.append(kin_value)

    error_definition_dSig_dmttBar = {}
    for i in range(len(values)):
        error_definition_dSig_dmttBar['ArtUnc_'+str(i+1)] = {'description': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_dSig_dmttBar_yaml = {'data_central': data_central_dSig_dmttBar}
    kinematics_dSig_dmttBar_yaml = {'bins': kin_dSig_dmttBar}
    uncertainties_dSig_dmttBar_yaml = {'definitions': error_definition_dSig_dmttBar, 'bins': error_dSig_dmttBar}

    with open('data_dSig_dmttBar.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_yaml, file, sort_keys=False)

# dSig_dmttBar_norm data

    hepdata_tables="rawdata/Table_4.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table_16.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    
    sqrts = float(input['dependent_variables'][0]['qualifiers'][1]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']
    
    for i in range(len(values)*len(values)):
        covMatArray_dSig_dmttBar_norm.append(input2['dependent_variables'][0]['values'][i]['value']*1e-6)
    artUnc_dSig_dmttBar_norm = cta(len(values), covMatArray_dSig_dmttBar_norm, 1)

    for i in range(len(values)):
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(len(values)):
            error_value['ArtUnc_'+str(j+1)] = artUnc_dSig_dmttBar_norm[i][j]
        data_central_value = values[i]['value']*1e-3
        data_central_dSig_dmttBar_norm.append(data_central_value)
        error_dSig_dmttBar_norm.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max}}
        kin_dSig_dmttBar_norm.append(kin_value)

    error_definition_dSig_dmttBar_norm = {}
    for i in range(len(values)):
        error_definition_dSig_dmttBar_norm['ArtUnc_'+str(i+1)] = {'description': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_dSig_dmttBar_norm_yaml = {'data_central': data_central_dSig_dmttBar_norm}
    kinematics_dSig_dmttBar_norm_yaml = {'bins': kin_dSig_dmttBar_norm}
    uncertainties_dSig_dmttBar_norm_yaml = {'definitions': error_definition_dSig_dmttBar_norm, 'bins': error_dSig_dmttBar_norm}

    with open('data_dSig_dmttBar_norm.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar_norm.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_norm_yaml, file, sort_keys=False)

# dSig_dyttBar data

    hepdata_tables="rawdata/Table_12.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table_24.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    
    sqrts = float(input['dependent_variables'][0]['qualifiers'][1]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']
    
    for i in range(len(values)*len(values)):
        covMatArray_dSig_dyttBar.append(input2['dependent_variables'][0]['values'][i]['value'])
    artUnc_dSig_dyttBar = cta(len(values), covMatArray_dSig_dyttBar)

    for i in range(len(values)):
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(len(values)):
            error_value['ArtUnc_'+str(j+1)] = artUnc_dSig_dyttBar[i][j]
        data_central_value = values[i]['value']
        data_central_dSig_dyttBar.append(data_central_value)
        error_dSig_dyttBar.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max}}
        kin_dSig_dyttBar.append(kin_value)

    error_definition_dSig_dyttBar = {}
    for i in range(len(values)):
        error_definition_dSig_dyttBar['ArtUnc_'+str(i+1)] = {'description': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_dSig_dyttBar_yaml = {'data_central': data_central_dSig_dyttBar}
    kinematics_dSig_dyttBar_yaml = {'bins': kin_dSig_dyttBar}
    uncertainties_dSig_dyttBar_yaml = {'definitions': error_definition_dSig_dyttBar, 'bins': error_dSig_dyttBar}

    with open('data_dSig_dyttBar.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyttBar.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_yaml, file, sort_keys=False)

# dSig_dyttBar_norm data

    hepdata_tables="rawdata/Table_6.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table_18.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    
    sqrts = float(input['dependent_variables'][0]['qualifiers'][1]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']
    
    for i in range(len(values)*len(values)):
        covMatArray_dSig_dyttBar_norm.append(input2['dependent_variables'][0]['values'][i]['value'])
    artUnc_dSig_dyttBar_norm = cta(len(values), covMatArray_dSig_dyttBar_norm, 1)

    for i in range(len(values)):
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(len(values)):
            error_value['ArtUnc_'+str(j+1)] = artUnc_dSig_dyttBar_norm[i][j]
        data_central_value = values[i]['value']
        data_central_dSig_dyttBar_norm.append(data_central_value)
        error_dSig_dyttBar_norm.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max}}
        kin_dSig_dyttBar_norm.append(kin_value)

    error_definition_dSig_dyttBar_norm = {}
    for i in range(len(values)):
        error_definition_dSig_dyttBar_norm['ArtUnc_'+str(i+1)] = {'description': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_dSig_dyttBar_norm_yaml = {'data_central': data_central_dSig_dyttBar_norm}
    kinematics_dSig_dyttBar_norm_yaml = {'bins': kin_dSig_dyttBar_norm}
    uncertainties_dSig_dyttBar_norm_yaml = {'definitions': error_definition_dSig_dyttBar_norm, 'bins': error_dSig_dyttBar_norm}

    with open('data_dSig_dyttBar_norm.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyttBar_norm.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_norm_yaml, file, sort_keys=False)

processData()
