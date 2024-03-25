from math import sqrt

import numpy as np
from numpy.linalg import eig
import yaml


def percentage_to_absolute(percentage, value):
    r"""Compute the absolute value of uncertainty from percentage.

    Parameters
    ----------
    percentage : string/float
        Experimental datasets can provide the percentage
        uncertainties with a % sign or without one.
        The function will autostrip % sign and convert to
        a float type in case the percentage uncertainty
        comes with a % sign. Else, it will directly perform
        the computation.
    value : float
        The data point

    Returns
    -------
    absolute : float
        The absolute value of the uncertainty

    """
    if type(percentage) is str:
        percentage = float(percentage.replace("%", ""))
        absolute = percentage * value * 0.01
        return absolute
    else:
        absolute = percentage * value * 0.01
        return absolute


def covmat_to_artunc(ndata, covmat_list, no_of_norm_mat=0):
    r"""Convert the covariance matrix to a matrix of
    artificial uncertainties.

    Parameters
    ----------
    ndata : integer
        Number of data points
    covmat_list : list
        A one dimensional list which contains the elements of
        the covariance matrix row by row. Since experimental
        datasets provide these matrices in a list form, this
        simplifies the implementation for the user.
    no_of_norm_mat : int
        Normalized covariance matrices may have an eigenvalue
        of 0 due to the last data point not being linearly
        independent. To allow for this, the user should input
        the number of normalized matrices that are being treated
        in an instance. For example, if a single covariance matrix
        of a normalized distribution is being processed, the input
        would be 1. If a covariance matrix contains pertains to
        3 normalized datasets (i.e. cross covmat for 3
        distributions), the input would be 3. The default value is
        0 for when the covariance matrix pertains to an absolute
        distribution.

    Returns
    -------
    artunc : list
        A two dimensional matrix (given as a list of lists)
        which contains artificial uncertainties to be added
        to the commondata. i^th row (or list) contains the
        artificial uncertainties of the i^th data point.

    """
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

    hepdata_tables = "rawdata/Table15.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix = "rawdata/Table16.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)

    systematics_breakdown = "rawdata/Table17.yaml"
    with open(systematics_breakdown, 'r') as file3:
        input3 = yaml.safe_load(file3)

    for i in range(ndata_dSig_dpTt * ndata_dSig_dpTt):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dpTt.append(covMatEl)
    artUncMat_dSig_dpTt = covmat_to_artunc(ndata_dSig_dpTt, covMatArray_dSig_dpTt, 1)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        pT_t_min = input['independent_variables'][0]['values'][i]['low']
        pT_t_mid = input['independent_variables'][1]['values'][i]['value']
        pT_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        error_value['stat'] = 0
        # error_value['sys'] = values[i]['errors'][1]['symerror']
        for j in range(ndata_dSig_dpTt):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dpTt[i][j])
        data_central_value = values[i]['value']
        for j in range(11):
            error_value[input3['independent_variables'][0]['values'][j]['value']] = (
                percentage_to_absolute(
                    str(input3['dependent_variables'][i]['values'][j]['value']), data_central_value
                )
            )
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': pT_t_mid, 'max': pT_t_max},
        }
        data_central_dSig_dpTt.append(data_central_value)
        kin_dSig_dpTt.append(kin_value)
        error_dSig_dpTt.append(error_value)

    error_definition_dSig_dpTt = {}
    error_definition_dSig_dpTt['stat'] = {
        'description': 'total statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    # error_definition_dSig_dpTt['sys'] = {'description': 'total systematic uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dpTt):
        error_definition_dSig_dpTt['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }
    for i in range(11):
        error_definition_dSig_dpTt[input3['independent_variables'][0]['values'][i]['value']] = {
            'definition': 'systematic uncertainty- '
            + str(input3['independent_variables'][0]['values'][i]['value']),
            'treatment': 'MULT',
            'type': 'CORR',
        }

    data_central_dSig_dpTt_norm_yaml = {'data_central': data_central_dSig_dpTt}
    kinematics_dSig_dpTt_norm_yaml = {'bins': kin_dSig_dpTt}
    uncertainties_dSig_dpTt_norm_yaml = {
        'definitions': error_definition_dSig_dpTt,
        'bins': error_dSig_dpTt,
    }

    with open('data_dSig_dpTt_norm.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dpTt_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dpTt_norm.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dpTt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dpTt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_norm_yaml, file, sort_keys=False)

    # dSig_dyt data

    hepdata_tables = "rawdata/Table21.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix = "rawdata/Table22.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)

    systematics_breakdown = "rawdata/Table23.yaml"
    with open(systematics_breakdown, 'r') as file3:
        input3 = yaml.safe_load(file3)

    for i in range(ndata_dSig_dyt * ndata_dSig_dyt):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyt.append(covMatEl)
    artUncMat_dSig_dyt = covmat_to_artunc(ndata_dSig_dyt, covMatArray_dSig_dyt, 1)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_mid = input['independent_variables'][1]['values'][i]['value']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        error_value['stat'] = 0
        # error_value['sys'] = values[i]['errors'][1]['symerror']
        for j in range(ndata_dSig_dyt):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dyt[i][j])
        data_central_value = values[i]['value']
        for j in range(11):
            error_value[input3['independent_variables'][0]['values'][j]['value']] = (
                percentage_to_absolute(
                    str(input3['dependent_variables'][i]['values'][j]['value']), data_central_value
                )
            )
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'y_t': {'min': y_t_min, 'mid': y_t_mid, 'max': y_t_max},
        }
        data_central_dSig_dyt.append(data_central_value)
        kin_dSig_dyt.append(kin_value)
        error_dSig_dyt.append(error_value)

    error_definition_dSig_dyt = {}
    error_definition_dSig_dyt['stat'] = {
        'description': 'total statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    # error_definition_dSig_dyt['sys'] = {'description': 'total systematic uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dyt):
        error_definition_dSig_dyt['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }
    for i in range(11):
        error_definition_dSig_dyt[input3['independent_variables'][0]['values'][i]['value']] = {
            'definition': 'systematic uncertainty- '
            + str(input3['independent_variables'][0]['values'][i]['value']),
            'treatment': 'MULT',
            'type': 'CORR',
        }

    data_central_dSig_dyt_norm_yaml = {'data_central': data_central_dSig_dyt}
    kinematics_dSig_dyt_norm_yaml = {'bins': kin_dSig_dyt}
    uncertainties_dSig_dyt_norm_yaml = {
        'definitions': error_definition_dSig_dyt,
        'bins': error_dSig_dyt,
    }

    with open('data_dSig_dyt_norm.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dyt_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyt_norm.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dyt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_norm_yaml, file, sort_keys=False)

    # dSig_dyttBar data

    hepdata_tables = "rawdata/Table36.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix = "rawdata/Table37.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)

    systematics_breakdown = "rawdata/Table38.yaml"
    with open(systematics_breakdown, 'r') as file3:
        input3 = yaml.safe_load(file3)

    for i in range(ndata_dSig_dyttBar * ndata_dSig_dyttBar):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyttBar.append(covMatEl)
    artUncMat_dSig_dyttBar = covmat_to_artunc(ndata_dSig_dyttBar, covMatArray_dSig_dyttBar, 1)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_mid = input['independent_variables'][1]['values'][i]['value']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        error_value['stat'] = 0
        # error_value['sys'] = values[i]['errors'][1]['symerror']
        for j in range(ndata_dSig_dyttBar):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dyttBar[i][j])
        data_central_value = values[i]['value']
        for j in range(11):
            error_value[input3['independent_variables'][0]['values'][j]['value']] = (
                percentage_to_absolute(
                    str(input3['dependent_variables'][i]['values'][j]['value']), data_central_value
                )
            )
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'y_ttBar': {'min': y_ttBar_min, 'mid': y_ttBar_mid, 'max': y_ttBar_max},
        }
        data_central_dSig_dyttBar.append(data_central_value)
        kin_dSig_dyttBar.append(kin_value)
        error_dSig_dyttBar.append(error_value)

    error_definition_dSig_dyttBar = {}
    error_definition_dSig_dyttBar['stat'] = {
        'description': 'total statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    # error_definition_dSig_dyttBar['sys'] = {'description': 'total systematic uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dyttBar):
        error_definition_dSig_dyttBar['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }
    for i in range(11):
        error_definition_dSig_dyttBar[input3['independent_variables'][0]['values'][i]['value']] = {
            'definition': 'systematic uncertainty- '
            + str(input3['independent_variables'][0]['values'][i]['value']),
            'treatment': 'MULT',
            'type': 'CORR',
        }

    data_central_dSig_dyttBar_norm_yaml = {'data_central': data_central_dSig_dyttBar}
    kinematics_dSig_dyttBar_norm_yaml = {'bins': kin_dSig_dyttBar}
    uncertainties_dSig_dyttBar_norm_yaml = {
        'definitions': error_definition_dSig_dyttBar,
        'bins': error_dSig_dyttBar,
    }

    with open('data_dSig_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    # dSig_dmttBar data

    hepdata_tables = "rawdata/Table39.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix = "rawdata/Table40.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)

    systematics_breakdown = "rawdata/Table41.yaml"
    with open(systematics_breakdown, 'r') as file3:
        input3 = yaml.safe_load(file3)

    for i in range(ndata_dSig_dmttBar * ndata_dSig_dmttBar):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dmttBar.append(covMatEl)
    artUncMat_dSig_dmttBar = covmat_to_artunc(ndata_dSig_dmttBar, covMatArray_dSig_dmttBar, 1)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_mid = input['independent_variables'][1]['values'][i]['value']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        error_value['stat'] = 0
        # error_value['sys'] = values[i]['errors'][1]['symerror']
        for j in range(ndata_dSig_dmttBar):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dmttBar[i][j])
        data_central_value = values[i]['value']
        for j in range(11):
            error_value[input3['independent_variables'][0]['values'][j]['value']] = (
                percentage_to_absolute(
                    str(input3['dependent_variables'][i]['values'][j]['value']), data_central_value
                )
            )
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'm_ttBar': {'min': m_ttBar_min, 'mid': m_ttBar_mid, 'max': m_ttBar_max},
        }
        data_central_dSig_dmttBar.append(data_central_value)
        kin_dSig_dmttBar.append(kin_value)
        error_dSig_dmttBar.append(error_value)

    error_definition_dSig_dmttBar = {}
    error_definition_dSig_dmttBar['stat'] = {
        'description': 'total statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    # error_definition_dSig_dmttBar['sys'] = {'description': 'total systematic uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dmttBar):
        error_definition_dSig_dmttBar['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }
    for i in range(11):
        error_definition_dSig_dmttBar[input3['independent_variables'][0]['values'][i]['value']] = {
            'definition': 'systematic uncertainty- '
            + str(input3['independent_variables'][0]['values'][i]['value']),
            'treatment': 'MULT',
            'type': 'CORR',
        }

    data_central_dSig_dmttBar_norm_yaml = {'data_central': data_central_dSig_dmttBar}
    kinematics_dSig_dmttBar_norm_yaml = {'bins': kin_dSig_dmttBar}
    uncertainties_dSig_dmttBar_norm_yaml = {
        'definitions': error_definition_dSig_dmttBar,
        'bins': error_dSig_dmttBar,
    }

    with open('data_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_norm_yaml, file, sort_keys=False)


processData()
