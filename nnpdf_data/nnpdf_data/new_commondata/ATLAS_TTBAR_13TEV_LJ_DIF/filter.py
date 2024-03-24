import yaml
import numpy as np

from math import sqrt
from numpy.linalg import eig


def se(delta_plus, delta_minus):
    
    semi_diff = (delta_plus + delta_minus)/2
    average = (delta_plus - delta_minus)/2
    se_delta = semi_diff
    se_sigma = sqrt(average*average + 2*semi_diff*semi_diff)
    return se_delta, se_sigma

def pta(percentage, value):
    
    if type(percentage) is str:
        percentage = float(percentage.replace("%", ""))
        absolute = percentage * value * 0.01
        return absolute 
    else:
        absolute = percentage * value * 0.01
        return absolute

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

    ndata_dSig_dmttBar = metadata['implemented_observables'][0]['ndata']
    ndata_dSig_dmttBar_norm = metadata['implemented_observables'][1]['ndata']
    ndata_dSig_dpTt = metadata['implemented_observables'][2]['ndata']
    ndata_dSig_dpTt_norm = metadata['implemented_observables'][3]['ndata']
    ndata_dSig_dyt = metadata['implemented_observables'][4]['ndata']
    ndata_dSig_dyt_norm = metadata['implemented_observables'][5]['ndata']
    ndata_dSig_dyttBar = metadata['implemented_observables'][6]['ndata']
    ndata_dSig_dyttBar_norm = metadata['implemented_observables'][7]['ndata']

    data_central_dSig_dmttBar = []
    kin_dSig_dmttBar = []
    error_dSig_dmttBar = []
    data_central_dSig_dmttBar_norm = []
    kin_dSig_dmttBar_norm = []
    error_dSig_dmttBar_norm = []
    data_central_dSig_dpTt = []
    kin_dSig_dpTt = []
    error_dSig_dpTt = []
    data_central_dSig_dpTt_norm = []
    kin_dSig_dpTt_norm = []
    error_dSig_dpTt_norm = []
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

    covMatArray_dSig_dmttBar = []
    covMatArray_dSig_dmttBar_norm = []
    covMatArray_dSig_dpTt = []
    covMatArray_dSig_dpTt_norm = []
    covMatArray_dSig_dyt = []
    covMatArray_dSig_dyt_norm = []
    covMatArray_dSig_dyttBar = []
    covMatArray_dSig_dyttBar_norm = []

# dSig_dmttBar data
    hepdata_tables="rawdata/Table618.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    
    covariance_matrix="rawdata/Table619.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dmttBar*ndata_dSig_dmttBar):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dmttBar.append(covMatEl)
    artUncMat_dSig_dmttBar = cta(ndata_dSig_dmttBar, covMatArray_dSig_dmttBar)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value'] # + value_delta
        data_central_dSig_dmttBar.append(data_central_value)
        for j in range(ndata_dSig_dmttBar):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dmttBar[i][j])
        error_dSig_dmttBar.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max}}
        kin_dSig_dmttBar.append(kin_value)

    error_definition_dSig_dmttBar = {}
    # error_definition_dSig_dmttBar['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dmttBar[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dmttBar):
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

# dSig_dmttBar_norm data
    hepdata_tables="rawdata/Table616.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table617.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dmttBar_norm*ndata_dSig_dmttBar_norm):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dmttBar_norm.append(covMatEl)
    artUncMat_dSig_dmttBar_norm = cta(ndata_dSig_dmttBar_norm, covMatArray_dSig_dmttBar_norm)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value'] # + value_delta
        data_central_dSig_dmttBar_norm.append(data_central_value)
        for j in range(ndata_dSig_dmttBar_norm):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dmttBar_norm[i][j])
        error_dSig_dmttBar_norm.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max}}
        kin_dSig_dmttBar_norm.append(kin_value)

    error_definition_dSig_dmttBar_norm = {}
    # error_definition_dSig_dmttBar_norm['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dmttBar_norm[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dmttBar_norm):
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

# dSig_dpTt data
    hepdata_tables="rawdata/Table610.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table611.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dpTt*ndata_dSig_dpTt):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dpTt.append(covMatEl)
    artUncMat_dSig_dpTt = cta(ndata_dSig_dpTt, covMatArray_dSig_dpTt)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        pT_t_min = input['independent_variables'][0]['values'][i]['low']
        pT_t_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value'] # + value_delta
        data_central_dSig_dpTt.append(data_central_value)
        for j in range(ndata_dSig_dpTt):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dpTt[i][j])
        error_dSig_dpTt.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max}}
        kin_dSig_dpTt.append(kin_value)

    error_definition_dSig_dpTt = {}
    # error_definition_dSig_dpTt['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dpTt[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dpTt):
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

# dSig_dpTt_norm data
    hepdata_tables="rawdata/Table608.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table609.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dpTt_norm*ndata_dSig_dpTt_norm):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dpTt_norm.append(covMatEl)
    artUncMat_dSig_dpTt_norm = cta(ndata_dSig_dpTt_norm, covMatArray_dSig_dpTt_norm)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        pT_t_min = input['independent_variables'][0]['values'][i]['low']
        pT_t_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value'] # + value_delta
        data_central_dSig_dpTt_norm.append(data_central_value)
        for j in range(ndata_dSig_dpTt_norm):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dpTt_norm[i][j])
        error_dSig_dpTt_norm.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max}}
        kin_dSig_dpTt_norm.append(kin_value)

    error_definition_dSig_dpTt_norm = {}
    # error_definition_dSig_dpTt_norm['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dpTt_norm[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dpTt_norm):
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

# dSig_dyt data
    hepdata_tables="rawdata/Table614.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table615.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dyt*ndata_dSig_dyt):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyt.append(covMatEl)
    artUncMat_dSig_dyt = cta(ndata_dSig_dyt, covMatArray_dSig_dyt)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value'] # + value_delta
        data_central_dSig_dyt.append(data_central_value)
        for j in range(ndata_dSig_dyt):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dyt[i][j])
        error_dSig_dyt.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max}}
        kin_dSig_dyt.append(kin_value)

    error_definition_dSig_dyt = {}
    # error_definition_dSig_dyt['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dyt[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dyt):
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

# dSig_dyt_norm data
    hepdata_tables="rawdata/Table612.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table613.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dyt_norm*ndata_dSig_dyt_norm):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyt_norm.append(covMatEl)
    artUncMat_dSig_dyt_norm = cta(ndata_dSig_dyt_norm, covMatArray_dSig_dyt_norm)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value'] # + value_delta
        data_central_dSig_dyt_norm.append(data_central_value)
        for j in range(ndata_dSig_dyt_norm):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dyt_norm[i][j])
        error_dSig_dyt_norm.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max}}
        kin_dSig_dyt_norm.append(kin_value)

    error_definition_dSig_dyt_norm = {}
    # error_definition_dSig_dyt_norm['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dyt_norm[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dyt_norm):
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

# dSig_dyttBar data
    hepdata_tables="rawdata/Table626.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table627.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dyttBar*ndata_dSig_dyttBar):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyttBar.append(covMatEl)
    artUncMat_dSig_dyttBar = cta(ndata_dSig_dyttBar, covMatArray_dSig_dyttBar)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value'] # + value_delta
        data_central_dSig_dyttBar.append(data_central_value)
        for j in range(ndata_dSig_dyttBar):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dyttBar[i][j])
        error_dSig_dyttBar.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max}}
        kin_dSig_dyttBar.append(kin_value)

    error_definition_dSig_dyttBar = {}
    # error_definition_dSig_dyttBar['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dyttBar[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dyttBar):
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

# dSig_dyttBar_norm data
    hepdata_tables="rawdata/Table624.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix="rawdata/Table625.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dyttBar_norm*ndata_dSig_dyttBar_norm):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyttBar_norm.append(covMatEl)
    artUncMat_dSig_dyttBar_norm = cta(ndata_dSig_dyttBar_norm, covMatArray_dSig_dyttBar_norm)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value'] # + value_delta
        data_central_dSig_dyttBar_norm.append(data_central_value)
        for j in range(ndata_dSig_dyttBar_norm):
            error_value['ArtUnc_'+str(j+1)] = float(artUncMat_dSig_dyttBar_norm[i][j])
        error_dSig_dyttBar_norm.append(error_value)
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max}}
        kin_dSig_dyttBar_norm.append(kin_value)

    error_definition_dSig_dyttBar_norm = {}
    # error_definition_dSig_dyttBar_norm['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dyttBar_norm[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dyttBar_norm):
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
