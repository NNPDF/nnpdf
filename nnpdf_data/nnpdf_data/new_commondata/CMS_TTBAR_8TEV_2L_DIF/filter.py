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

def ctc(err_list, cormat_list):
    
    covmat_list = []
    for i in range(len(cormat_list)):
        a = i // len(err_list)
        b = i % len(err_list)
        covmat_list.append(cormat_list[i] * err_list[a] * err_list[b])
    return covmat_list

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

def ttf(mode, tri_mat_list):
    
    dim = int((np.sqrt(1 + 8*len(tri_mat_list)) - 1)/2)
    matrix = np.zeros((dim, dim))
    if mode == 0:
        for i in range(dim):
            for j in range(i + 1):
                list_el = len(tri_mat_list) - 1 - ((i*(i + 1))//2 + j)
                if i == j:
                    matrix[dim - 1 - i][dim - 1 - j] = tri_mat_list[list_el]
                else:
                    matrix[dim - 1 - i][dim - 1 - j] = tri_mat_list[list_el]
                    matrix[dim - 1 - j][dim - 1 - i] = tri_mat_list[list_el]
    elif mode == 1:
        for i in range(dim):
            for j in range(i + 1):
                list_el = (i*(i + 1))//2 + j
                if i == j:
                    matrix[i][j] = tri_mat_list[list_el]
                else:
                    matrix[i][j] = tri_mat_list[list_el]
                    matrix[j][i] = tri_mat_list[list_el]
    else:
        raise Exception('Mode should be 0 or 1')
    mat_list = []
    for i in range(dim):
        for j in range(dim):
            mat_list.append(matrix[i][j])
    return mat_list

def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    data_central_d2Sig_dyt_dpTt_norm = []
    kin_d2Sig_dyt_dpTt_norm = []
    error_d2Sig_dyt_dpTt_norm = []
    data_central_d2Sig_dyt_dmttBar_norm = []
    kin_d2Sig_dyt_dmttBar_norm = []
    error_d2Sig_dyt_dmttBar_norm = []
    data_central_d2Sig_dmttBar_dyttBar_norm = []
    kin_d2Sig_dmttBar_dyttBar_norm = []
    error_d2Sig_dmttBar_dyttBar_norm = []

# d2Sig_dyt_dpTt_norm

    hepdata_tables="rawdata/CMS_8TeV_ttbar_DoubleDiff_yt_ptt.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    correlation_matrix="rawdata/CMS_8TeV_ttbar_DoubleDiff_yt_ptt_statcorr.yaml"
    with open(correlation_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
#    systematics_breakdown="rawdata/CMS_8TeV_ttbar_DoubleDiff_yt_ptt_syst.yaml"
#    with open(systematics_breakdown, 'r') as file:
#        input3 = yaml.safe_load(file)

    sqrts = 8000
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']
    statlist1 = []
    for i in range(len(values)):
        statlist1.append(pta(str(values[i]['errors'][0]['symerror']), values[i]['value']))
    trimatlist1 = []
    for i in range(len(input2['dependent_variables'][0]['values'])):
        trimatlist1.append(input2['dependent_variables'][0]['values'][i]['value'])
    cormatlist1 = ttf(0, trimatlist1)
    covmatlist1 = ctc(statlist1, cormatlist1)
    artunc1 = cta(len(values), covmatlist1, 1)
    

    for i in range(len(values)):
        data_central_value = values[i]['value']
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        pT_t_min = input['independent_variables'][1]['values'][i]['low']
        pT_t_max = input['independent_variables'][1]['values'][i]['high']
        error_value = {}
        plus = pta(str(values[i]['errors'][1]['asymerror']['plus']), data_central_value)
        minus = pta(str(values[i]['errors'][1]['asymerror']['minus']), data_central_value)
        se_delta, se_sigma = se(plus, minus)
        data_central_value = data_central_value + se_delta
        error_value['sys'] = se_sigma
        for j in range(len(values)):
            error_value['ArtUnc_'+str(j+1)] = artunc1[i][j]
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max}, 'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max}}
        data_central_d2Sig_dyt_dpTt_norm.append(data_central_value)
        kin_d2Sig_dyt_dpTt_norm.append(kin_value)
        error_d2Sig_dyt_dpTt_norm.append(error_value)

    error_definition_d2Sig_dyt_dpTt_norm = {}
    error_definition_d2Sig_dyt_dpTt_norm['sys'] = {'definition': 'total systematic uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(16):
        error_definition_d2Sig_dyt_dpTt_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}
    
    data_central_d2Sig_dyt_dpTt_norm_yaml = {'data_central': data_central_d2Sig_dyt_dpTt_norm}
    kinematics_d2Sig_dyt_dpTt_norm_yaml = {'bins': kin_d2Sig_dyt_dpTt_norm}
    uncertainties_d2Sig_dyt_dpTt_norm_yaml = {'definitions': error_definition_d2Sig_dyt_dpTt_norm, 'bins': error_d2Sig_dyt_dpTt_norm}

    with open('data_d2Sig_dyt_dpTt_norm.yaml', 'w') as file:
         yaml.dump(data_central_d2Sig_dyt_dpTt_norm_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dyt_dpTt_norm.yaml', 'w') as file:
         yaml.dump(kinematics_d2Sig_dyt_dpTt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dyt_dpTt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dyt_dpTt_norm_yaml, file, sort_keys=False)

# d2Sig_dyt_dmttBar_norm

    hepdata_tables="rawdata/CMS_8TeV_ttbar_DoubleDiff_mtt_yt.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    correlation_matrix="rawdata/CMS_8TeV_ttbar_DoubleDiff_mtt_yt_statcorr.yaml"
    with open(correlation_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
#    systematics_breakdown="rawdata/CMS_8TeV_ttbar_DoubleDiff_mtt_yt_syst.yaml"
#    with open(systematics_breakdown, 'r') as file:
#        input3 = yaml.safe_load(file)

    sqrts = 8000
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']
    statlist2 = []
    for i in range(len(values)):
        statlist2.append(pta(str(values[i]['errors'][0]['symerror']), values[i]['value']))
    trimatlist2 = []
    for i in range(len(input2['dependent_variables'][0]['values'])):
        trimatlist2.append(input2['dependent_variables'][0]['values'][i]['value'])
    cormatlist2 = ttf(0, trimatlist2)
    covmatlist2 = ctc(statlist2, cormatlist2)
    artunc2 = cta(len(values), covmatlist2, 1)
    

    for i in range(len(values)):
        data_central_value = values[i]['value']
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        y_t_min = input['independent_variables'][1]['values'][i]['low']
        y_t_max = input['independent_variables'][1]['values'][i]['high']
        error_value = {}
        plus = pta(str(values[i]['errors'][1]['asymerror']['plus']), data_central_value)
        minus = pta(str(values[i]['errors'][1]['asymerror']['minus']), data_central_value)
        se_delta, se_sigma = se(plus, minus)
        data_central_value = data_central_value + se_delta
        error_value['sys'] = se_sigma
        for j in range(len(values)):
            error_value['ArtUnc_'+str(j+1)] = artunc2[i][j]
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max}, 'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max}}
        data_central_d2Sig_dyt_dmttBar_norm.append(data_central_value)
        kin_d2Sig_dyt_dmttBar_norm.append(kin_value)
        error_d2Sig_dyt_dmttBar_norm.append(error_value)

    error_definition_d2Sig_dyt_dmttBar_norm = {}
    error_definition_d2Sig_dyt_dmttBar_norm['sys'] = {'definition': 'total systematic uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(16):
        error_definition_d2Sig_dyt_dmttBar_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}
    
    data_central_d2Sig_dyt_dmttBar_norm_yaml = {'data_central': data_central_d2Sig_dyt_dmttBar_norm}
    kinematics_d2Sig_dyt_dmttBar_norm_yaml = {'bins': kin_d2Sig_dyt_dmttBar_norm}
    uncertainties_d2Sig_dyt_dmttBar_norm_yaml = {'definitions': error_definition_d2Sig_dyt_dmttBar_norm, 'bins': error_d2Sig_dyt_dmttBar_norm}

    with open('data_d2Sig_dyt_dmttBar_norm.yaml', 'w') as file:
         yaml.dump(data_central_d2Sig_dyt_dmttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dyt_dmttBar_norm.yaml', 'w') as file:
         yaml.dump(kinematics_d2Sig_dyt_dmttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dyt_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dyt_dmttBar_norm_yaml, file, sort_keys=False)

# d2Sig_dmttBar_dyttBar_norm

    hepdata_tables="rawdata/CMS_8TeV_ttbar_DoubleDiff_mtt_ytt.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    correlation_matrix="rawdata/CMS_8TeV_ttbar_DoubleDiff_mtt_ytt_statcorr.yaml"
    with open(correlation_matrix, 'r') as file:
        input2 = yaml.safe_load(file)
#    systematics_breakdown="rawdata/CMS_8TeV_ttbar_DoubleDiff_mtt_ytt_syst.yaml"
#    with open(systematics_breakdown, 'r') as file:
#        input3 = yaml.safe_load(file)

    sqrts = 8000
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']
    statlist3 = []
    for i in range(len(values)):
        statlist3.append(pta(str(values[i]['errors'][0]['symerror']), values[i]['value']))
    trimatlist3 = []
    for i in range(len(input2['dependent_variables'][0]['values'])):
        trimatlist3.append(input2['dependent_variables'][0]['values'][i]['value'])
    cormatlist3 = ttf(0, trimatlist3)
    covmatlist3 = ctc(statlist3, cormatlist3)
    artunc3 = cta(len(values), covmatlist3, 1)
    

    for i in range(len(values)):
        data_central_value = values[i]['value']
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        y_ttBar_min = input['independent_variables'][1]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][1]['values'][i]['high']
        error_value = {}
        plus = pta(str(values[i]['errors'][1]['asymerror']['plus']), data_central_value)
        minus = pta(str(values[i]['errors'][1]['asymerror']['minus']), data_central_value)
        se_delta, se_sigma = se(plus, minus)
        data_central_value = data_central_value + se_delta
        error_value['sys'] = se_sigma
        for j in range(len(values)):
            error_value['ArtUnc_'+str(j+1)] = artunc3[i][j]
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max}, 'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max}}
        data_central_d2Sig_dmttBar_dyttBar_norm.append(data_central_value)
        kin_d2Sig_dmttBar_dyttBar_norm.append(kin_value)
        error_d2Sig_dmttBar_dyttBar_norm.append(error_value)

    error_definition_d2Sig_dmttBar_dyttBar_norm = {}
    error_definition_d2Sig_dmttBar_dyttBar_norm['sys'] = {'definition': 'total systematic uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(16):
        error_definition_d2Sig_dmttBar_dyttBar_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}
    
    data_central_d2Sig_dmttBar_dyttBar_norm_yaml = {'data_central': data_central_d2Sig_dmttBar_dyttBar_norm}
    kinematics_d2Sig_dmttBar_dyttBar_norm_yaml = {'bins': kin_d2Sig_dmttBar_dyttBar_norm}
    uncertainties_d2Sig_dmttBar_dyttBar_norm_yaml = {'definitions': error_definition_d2Sig_dmttBar_dyttBar_norm, 'bins': error_d2Sig_dmttBar_dyttBar_norm}

    with open('data_d2Sig_dmttBar_dyttBar_norm.yaml', 'w') as file:
         yaml.dump(data_central_d2Sig_dmttBar_dyttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dmttBar_dyttBar_norm.yaml', 'w') as file:
         yaml.dump(kinematics_d2Sig_dmttBar_dyttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dmttBar_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dmttBar_dyttBar_norm_yaml, file, sort_keys=False)    

processData()
