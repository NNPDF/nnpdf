import yaml

import numpy as np

from math import sqrt
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

    data_central_dSig_dpTt = []
    kin_dSig_dpTt = []
    error_dSig_dpTt = []
    data_central_dSig_dpTt_norm = []
    kin_dSig_dpTt_norm = []
    error_dSig_dpTt_norm = []
    data_central_dSig_dmttBar = []
    kin_dSig_dmttBar = []
    error_dSig_dmttBar = []
    data_central_dSig_dmttBar_norm = []
    kin_dSig_dmttBar_norm = []
    error_dSig_dmttBar_norm = []
    data_central_dSig_dyt = []
    kin_dSig_dyt = []
    error_dSig_dyt = []
    data_central_dSig_dyt_norm = []
    kin_dSig_dyt_norm = []
    error_dSig_dyt_norm = []
    data_central_dSig_dyttBar = []
    kin_dSig_dyttBar = []
    error_dSig_dyttBar = []
    data_central_dSig_dyttBar_norm = []
    kin_dSig_dyttBar_norm = []
    error_dSig_dyttBar_norm = []

    covmat_dSig_dpTt = []
    covmat_dSig_dpTt_norm = []
    covmat_dSig_dmttBar = []
    covmat_dSig_dmttBar_norm = []
    covmat_dSig_dyt = []
    covmat_dSig_dyt_norm = []
    covmat_dSig_dyttBar = []
    covmat_dSig_dyttBar_norm = []

# dSig_dpTt

    hepdata_tables="rawdata/d01-x01-y01.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    covariance_matrix="rawdata/d01-x01-y01_cov.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    for i in range(36):
        covmat_dSig_dpTt.append(input2['dependent_variables'][0]['values'][i]['value'])
    artunc_dSig_dpTt = cta(6, covmat_dSig_dpTt, 0)


    sqrts = 13000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(6):
        data_central_value = values[i]['value']
        pT_t_min = input['independent_variables'][0]['values'][i]['low']
        pT_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(6):
            error_value['ArtUnc_'+str(j+1)] = artunc_dSig_dpTt[i][j]
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max}}
        data_central_dSig_dpTt.append(data_central_value)
        kin_dSig_dpTt.append(kin_value)
        error_dSig_dpTt.append(error_value)

    error_definition_dSig_dpTt = {}
    for i in range(6):
        error_definition_dSig_dpTt['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_dSig_dpTt_yaml = {'data_central': data_central_dSig_dpTt}
    kinematics_dSig_dpTt_yaml = {'bins': kin_dSig_dpTt}
    uncertainties_dSig_dpTt_yaml = {'definitions': error_definition_dSig_dpTt, 'bins': error_dSig_dpTt}

    with open('data_dSig_dpTt.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dpTt_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dpTt.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dpTt_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dpTt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_yaml, file, sort_keys=False)

# dSig_dpTt_norm

    hepdata_tables="rawdata/d02-x01-y01.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    covariance_matrix="rawdata/d02-x01-y01_cov.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    for i in range(25):
        covmat_dSig_dpTt_norm.append(input2['dependent_variables'][0]['values'][i]['value'])
    artunc_dSig_dpTt_norm = cta(5, covmat_dSig_dpTt_norm, 1)


    sqrts = 13000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(5):
        data_central_value = values[i]['value']
        pT_t_min = input['independent_variables'][0]['values'][i]['low']
        pT_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(5):
            error_value['ArtUnc_'+str(j+1)] = artunc_dSig_dpTt_norm[i][j]
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max}}
        data_central_dSig_dpTt_norm.append(data_central_value)
        kin_dSig_dpTt_norm.append(kin_value)
        error_dSig_dpTt_norm.append(error_value)

    error_definition_dSig_dpTt_norm = {}
    for i in range(5):
        error_definition_dSig_dpTt_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_dSig_dpTt_norm_yaml = {'data_central': data_central_dSig_dpTt_norm}
    kinematics_dSig_dpTt_norm_yaml = {'bins': kin_dSig_dpTt_norm}
    uncertainties_dSig_dpTt_norm_yaml = {'definitions': error_definition_dSig_dpTt_norm, 'bins': error_dSig_dpTt_norm}

    with open('data_dSig_dpTt_norm.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dpTt_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dpTt_norm.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dpTt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dpTt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_norm_yaml, file, sort_keys=False)

# dSig_dmttBar

    hepdata_tables="rawdata/d045-x01-y01.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    covariance_matrix="rawdata/d045-x01-y01_cov.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    for i in range(49):
        covmat_dSig_dmttBar.append(input2['dependent_variables'][0]['values'][i]['value'])
    artunc_dSig_dmttBar = cta(7, covmat_dSig_dmttBar, 0)


    sqrts = 13000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(7):
        data_central_value = values[i]['value']
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(7):
            error_value['ArtUnc_'+str(j+1)] = artunc_dSig_dmttBar[i][j]
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max}}
        data_central_dSig_dmttBar.append(data_central_value)
        kin_dSig_dmttBar.append(kin_value)
        error_dSig_dmttBar.append(error_value)

    error_definition_dSig_dmttBar = {}
    for i in range(7):
        error_definition_dSig_dmttBar['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_dSig_dmttBar_yaml = {'data_central': data_central_dSig_dmttBar}
    kinematics_dSig_dmttBar_yaml = {'bins': kin_dSig_dmttBar}
    uncertainties_dSig_dmttBar_yaml = {'definitions': error_definition_dSig_dmttBar, 'bins': error_dSig_dmttBar}

    with open('data_dSig_dmttBar.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_yaml, file, sort_keys=False)

# dSig_dmttBar_norm

    hepdata_tables="rawdata/d046-x01-y01.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    covariance_matrix="rawdata/d046-x01-y01_cov.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    for i in range(36):
        covmat_dSig_dmttBar_norm.append(input2['dependent_variables'][0]['values'][i]['value'])
    artunc_dSig_dmttBar_norm = cta(6, covmat_dSig_dmttBar_norm, 1)


    sqrts = 13000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(6):
        data_central_value = values[i]['value']
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(6):
            error_value['ArtUnc_'+str(j+1)] = artunc_dSig_dmttBar_norm[i][j]
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max}}
        data_central_dSig_dmttBar_norm.append(data_central_value)
        kin_dSig_dmttBar_norm.append(kin_value)
        error_dSig_dmttBar_norm.append(error_value)

    error_definition_dSig_dmttBar_norm = {}
    for i in range(6):
        error_definition_dSig_dmttBar_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_dSig_dmttBar_norm_yaml = {'data_central': data_central_dSig_dmttBar_norm}
    kinematics_dSig_dmttBar_norm_yaml = {'bins': kin_dSig_dmttBar_norm}
    uncertainties_dSig_dmttBar_norm_yaml = {'definitions': error_definition_dSig_dmttBar_norm, 'bins': error_dSig_dmttBar_norm}

    with open('data_dSig_dmttBar_norm.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar_norm.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_norm_yaml, file, sort_keys=False)

# dSig_dyt

    hepdata_tables="rawdata/d021-x01-y01.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    covariance_matrix="rawdata/d021-x01-y01_cov.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    for i in range(100):
        covmat_dSig_dyt.append(input2['dependent_variables'][0]['values'][i]['value'])
    artunc_dSig_dyt = cta(10, covmat_dSig_dyt, 0)


    sqrts = 13000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(10):
        data_central_value = values[i]['value']
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(10):
            error_value['ArtUnc_'+str(j+1)] = artunc_dSig_dyt[i][j]
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max}}
        data_central_dSig_dyt.append(data_central_value)
        kin_dSig_dyt.append(kin_value)
        error_dSig_dyt.append(error_value)

    error_definition_dSig_dyt = {}
    for i in range(10):
        error_definition_dSig_dyt['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_dSig_dyt_yaml = {'data_central': data_central_dSig_dyt}
    kinematics_dSig_dyt_yaml = {'bins': kin_dSig_dyt}
    uncertainties_dSig_dyt_yaml = {'definitions': error_definition_dSig_dyt, 'bins': error_dSig_dyt}

    with open('data_dSig_dyt.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dyt_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyt.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dyt_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_yaml, file, sort_keys=False)

# dSig_dyt_norm

    hepdata_tables="rawdata/d022-x01-y01.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    covariance_matrix="rawdata/d022-x01-y01_cov.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    for i in range(81):
        covmat_dSig_dyt_norm.append(input2['dependent_variables'][0]['values'][i]['value'])
    artunc_dSig_dyt_norm = cta(9, covmat_dSig_dyt_norm, 1)


    sqrts = 13000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(9):
        data_central_value = values[i]['value']
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(9):
            error_value['ArtUnc_'+str(j+1)] = artunc_dSig_dyt_norm[i][j]
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max}}
        data_central_dSig_dyt_norm.append(data_central_value)
        kin_dSig_dyt_norm.append(kin_value)
        error_dSig_dyt_norm.append(error_value)

    error_definition_dSig_dyt_norm = {}
    for i in range(9):
        error_definition_dSig_dyt_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_dSig_dyt_norm_yaml = {'data_central': data_central_dSig_dyt_norm}
    kinematics_dSig_dyt_norm_yaml = {'bins': kin_dSig_dyt_norm}
    uncertainties_dSig_dyt_norm_yaml = {'definitions': error_definition_dSig_dyt_norm, 'bins': error_dSig_dyt_norm}

    with open('data_dSig_dyt_norm.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dyt_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyt_norm.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dyt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_norm_yaml, file, sort_keys=False)

# dSig_dyttBar

    hepdata_tables="rawdata/d041-x01-y01.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    covariance_matrix="rawdata/d041-x01-y01_cov.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    for i in range(100):
        covmat_dSig_dyttBar.append(input2['dependent_variables'][0]['values'][i]['value'])
    artunc_dSig_dyttBar = cta(10, covmat_dSig_dyttBar, 0)


    sqrts = 13000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(10):
        data_central_value = values[i]['value']
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(10):
            error_value['ArtUnc_'+str(j+1)] = artunc_dSig_dyttBar[i][j]
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max}}
        data_central_dSig_dyttBar.append(data_central_value)
        kin_dSig_dyttBar.append(kin_value)
        error_dSig_dyttBar.append(error_value)

    error_definition_dSig_dyttBar = {}
    for i in range(10):
        error_definition_dSig_dyttBar['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_dSig_dyttBar_yaml = {'data_central': data_central_dSig_dyttBar}
    kinematics_dSig_dyttBar_yaml = {'bins': kin_dSig_dyttBar}
    uncertainties_dSig_dyttBar_yaml = {'definitions': error_definition_dSig_dyttBar, 'bins': error_dSig_dyttBar}

    with open('data_dSig_dyttBar.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyttBar.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_yaml, file, sort_keys=False)

# dSig_dyttBar_norm

    hepdata_tables="rawdata/d042-x01-y01.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    covariance_matrix="rawdata/d042-x01-y01_cov.yaml"
    with open(covariance_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
    for i in range(81):
        covmat_dSig_dyttBar_norm.append(input2['dependent_variables'][0]['values'][i]['value'])
    artunc_dSig_dyttBar_norm = cta(9, covmat_dSig_dyttBar_norm, 1)


    sqrts = 13000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(9):
        data_central_value = values[i]['value']
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(9):
            error_value['ArtUnc_'+str(j+1)] = artunc_dSig_dyttBar_norm[i][j]
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max}}
        data_central_dSig_dyttBar_norm.append(data_central_value)
        kin_dSig_dyttBar_norm.append(kin_value)
        error_dSig_dyttBar_norm.append(error_value)

    error_definition_dSig_dyttBar_norm = {}
    for i in range(9):
        error_definition_dSig_dyttBar_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

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
