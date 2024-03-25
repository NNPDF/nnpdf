from math import sqrt

import numpy as np
from numpy.linalg import eig
import yaml


def symmetrize_errors(delta_plus, delta_minus):
    r"""Compute the symmterized uncertainty and the shift in data point.

    Parameters
    ----------
    delta_plus : float
        The top/plus uncertainty with sign
    delta_minus : float
        The bottom/minus uncertainty with sign

    Returns
    -------
    se_delta : float
        The value to be added to the data point
    se_sigma : float
        The symmetrized uncertainty to be used in commondata

    """
    semi_diff = (delta_plus + delta_minus) / 2
    average = (delta_plus - delta_minus) / 2
    se_delta = semi_diff
    se_sigma = sqrt(average * average + 2 * semi_diff * semi_diff)
    return se_delta, se_sigma


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


def cormat_to_covmat(err_list, cormat_list):
    r"""Convert correlation matrix elements to covariance
    matrix elements.

    Parameters
    ----------
    err_list : list
        A one dimensional list which contains the uncertainty
        associated to each data point in order.
    cormat_list : list
        A one dimensional list which contains the elements of
        the correlation matrix row by row. Since experimental
        datasets provide these matrices in a list form, this
        simplifies the implementation for the user.

    Returns
    -------
    covmat_list : list
        A one dimensional list which contains the elements of
        the covariance matrix row by row.

    """
    covmat_list = []
    for i in range(len(cormat_list)):
        a = i // len(err_list)
        b = i % len(err_list)
        covmat_list.append(cormat_list[i] * err_list[a] * err_list[b])
    return covmat_list


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


def trimat_to_fullmat(mode, tri_mat_list):
    r"""Convert a list of values of a triangular matrix
    to a symmetric matrix.

    Experimental datasets can provide the entries of
    correlation or covariance matrices as a triangular
    matrix, as these matrices are symmetric by their
    very nature. This function can convert these list to
    a complete symmetric matrix, that can be used for the
    dataset implementation.

    mode : bool
        Enter 0 or 1 based on the following scenarios:
        Use mode 0 if matrix entries are given row by
        row such as:
        0 1 2 3
          4 5 6
            7 8
              9
        Use mode 1 if the matrix entries are given column
        by column such as:
        0 1 3 6
          2 4 7
            5 8
              9
        Please note that the numbers above (0-9) are not
        entries of the matrix but rather the index of the
        entries of the list which contains the elements of
        the triangular matrix.
    tri_mat_list : list
        A list containing the elements of the triangular matrix,
        for example, for a 4*4 matrix, the list of
        triangular matrix entries could be:
        [a, b, c, d, e, f, g, h, i, j]

    Returns
    -------
    mat_list : list
        A one dimensional list which contains the elements of
        the fully populated, symmetric matrix row by row.

    """
    dim = int((np.sqrt(1 + 8 * len(tri_mat_list)) - 1) / 2)
    matrix = np.zeros((dim, dim))
    if mode == 0:
        for i in range(dim):
            for j in range(i + 1):
                list_el = len(tri_mat_list) - 1 - ((i * (i + 1)) // 2 + j)
                if i == j:
                    matrix[dim - 1 - i][dim - 1 - j] = tri_mat_list[list_el]
                else:
                    matrix[dim - 1 - i][dim - 1 - j] = tri_mat_list[list_el]
                    matrix[dim - 1 - j][dim - 1 - i] = tri_mat_list[list_el]
    elif mode == 1:
        for i in range(dim):
            for j in range(i + 1):
                list_el = (i * (i + 1)) // 2 + j
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

    hepdata_tables = "rawdata/CMS_8TeV_ttbar_DoubleDiff_yt_ptt.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    correlation_matrix = "rawdata/CMS_8TeV_ttbar_DoubleDiff_yt_ptt_statcorr.yaml"
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
        statlist1.append(
            percentage_to_absolute(str(values[i]['errors'][0]['symerror']), values[i]['value'])
        )
    trimatlist1 = []
    for i in range(len(input2['dependent_variables'][0]['values'])):
        trimatlist1.append(input2['dependent_variables'][0]['values'][i]['value'])
    cormatlist1 = trimat_to_fullmat(0, trimatlist1)
    covmatlist1 = cormat_to_covmat(statlist1, cormatlist1)
    artunc1 = covmat_to_artunc(len(values), covmatlist1, 1)

    for i in range(len(values)):
        data_central_value = values[i]['value']
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        pT_t_min = input['independent_variables'][1]['values'][i]['low']
        pT_t_max = input['independent_variables'][1]['values'][i]['high']
        error_value = {}
        plus = percentage_to_absolute(
            str(values[i]['errors'][1]['asymerror']['plus']), data_central_value
        )
        minus = percentage_to_absolute(
            str(values[i]['errors'][1]['asymerror']['minus']), data_central_value
        )
        se_delta, se_sigma = symmetrize_errors(plus, minus)
        data_central_value = data_central_value + se_delta
        error_value['sys'] = se_sigma
        for j in range(len(values)):
            error_value['ArtUnc_' + str(j + 1)] = artunc1[i][j]
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
        }
        data_central_d2Sig_dyt_dpTt_norm.append(data_central_value)
        kin_d2Sig_dyt_dpTt_norm.append(kin_value)
        error_d2Sig_dyt_dpTt_norm.append(error_value)

    error_definition_d2Sig_dyt_dpTt_norm = {}
    error_definition_d2Sig_dyt_dpTt_norm['sys'] = {
        'definition': 'total systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    }
    for i in range(16):
        error_definition_d2Sig_dyt_dpTt_norm['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }

    data_central_d2Sig_dyt_dpTt_norm_yaml = {'data_central': data_central_d2Sig_dyt_dpTt_norm}
    kinematics_d2Sig_dyt_dpTt_norm_yaml = {'bins': kin_d2Sig_dyt_dpTt_norm}
    uncertainties_d2Sig_dyt_dpTt_norm_yaml = {
        'definitions': error_definition_d2Sig_dyt_dpTt_norm,
        'bins': error_d2Sig_dyt_dpTt_norm,
    }

    with open('data_d2Sig_dyt_dpTt_norm.yaml', 'w') as file:
        yaml.dump(data_central_d2Sig_dyt_dpTt_norm_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dyt_dpTt_norm.yaml', 'w') as file:
        yaml.dump(kinematics_d2Sig_dyt_dpTt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dyt_dpTt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dyt_dpTt_norm_yaml, file, sort_keys=False)

    # d2Sig_dyt_dmttBar_norm

    hepdata_tables = "rawdata/CMS_8TeV_ttbar_DoubleDiff_mtt_yt.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    correlation_matrix = "rawdata/CMS_8TeV_ttbar_DoubleDiff_mtt_yt_statcorr.yaml"
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
        statlist2.append(
            percentage_to_absolute(str(values[i]['errors'][0]['symerror']), values[i]['value'])
        )
    trimatlist2 = []
    for i in range(len(input2['dependent_variables'][0]['values'])):
        trimatlist2.append(input2['dependent_variables'][0]['values'][i]['value'])
    cormatlist2 = trimat_to_fullmat(0, trimatlist2)
    covmatlist2 = cormat_to_covmat(statlist2, cormatlist2)
    artunc2 = covmat_to_artunc(len(values), covmatlist2, 1)

    for i in range(len(values)):
        data_central_value = values[i]['value']
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        y_t_min = input['independent_variables'][1]['values'][i]['low']
        y_t_max = input['independent_variables'][1]['values'][i]['high']
        error_value = {}
        plus = percentage_to_absolute(
            str(values[i]['errors'][1]['asymerror']['plus']), data_central_value
        )
        minus = percentage_to_absolute(
            str(values[i]['errors'][1]['asymerror']['minus']), data_central_value
        )
        se_delta, se_sigma = symmetrize_errors(plus, minus)
        data_central_value = data_central_value + se_delta
        error_value['sys'] = se_sigma
        for j in range(len(values)):
            error_value['ArtUnc_' + str(j + 1)] = artunc2[i][j]
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        data_central_d2Sig_dyt_dmttBar_norm.append(data_central_value)
        kin_d2Sig_dyt_dmttBar_norm.append(kin_value)
        error_d2Sig_dyt_dmttBar_norm.append(error_value)

    error_definition_d2Sig_dyt_dmttBar_norm = {}
    error_definition_d2Sig_dyt_dmttBar_norm['sys'] = {
        'definition': 'total systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    }
    for i in range(16):
        error_definition_d2Sig_dyt_dmttBar_norm['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }

    data_central_d2Sig_dyt_dmttBar_norm_yaml = {'data_central': data_central_d2Sig_dyt_dmttBar_norm}
    kinematics_d2Sig_dyt_dmttBar_norm_yaml = {'bins': kin_d2Sig_dyt_dmttBar_norm}
    uncertainties_d2Sig_dyt_dmttBar_norm_yaml = {
        'definitions': error_definition_d2Sig_dyt_dmttBar_norm,
        'bins': error_d2Sig_dyt_dmttBar_norm,
    }

    with open('data_d2Sig_dyt_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(data_central_d2Sig_dyt_dmttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dyt_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(kinematics_d2Sig_dyt_dmttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dyt_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dyt_dmttBar_norm_yaml, file, sort_keys=False)

    # d2Sig_dmttBar_dyttBar_norm

    hepdata_tables = "rawdata/CMS_8TeV_ttbar_DoubleDiff_mtt_ytt.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    correlation_matrix = "rawdata/CMS_8TeV_ttbar_DoubleDiff_mtt_ytt_statcorr.yaml"
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
        statlist3.append(
            percentage_to_absolute(str(values[i]['errors'][0]['symerror']), values[i]['value'])
        )
    trimatlist3 = []
    for i in range(len(input2['dependent_variables'][0]['values'])):
        trimatlist3.append(input2['dependent_variables'][0]['values'][i]['value'])
    cormatlist3 = trimat_to_fullmat(0, trimatlist3)
    covmatlist3 = cormat_to_covmat(statlist3, cormatlist3)
    artunc3 = covmat_to_artunc(len(values), covmatlist3, 1)

    for i in range(len(values)):
        data_central_value = values[i]['value']
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        y_ttBar_min = input['independent_variables'][1]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][1]['values'][i]['high']
        error_value = {}
        plus = percentage_to_absolute(
            str(values[i]['errors'][1]['asymerror']['plus']), data_central_value
        )
        minus = percentage_to_absolute(
            str(values[i]['errors'][1]['asymerror']['minus']), data_central_value
        )
        se_delta, se_sigma = symmetrize_errors(plus, minus)
        data_central_value = data_central_value + se_delta
        error_value['sys'] = se_sigma
        for j in range(len(values)):
            error_value['ArtUnc_' + str(j + 1)] = artunc3[i][j]
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        data_central_d2Sig_dmttBar_dyttBar_norm.append(data_central_value)
        kin_d2Sig_dmttBar_dyttBar_norm.append(kin_value)
        error_d2Sig_dmttBar_dyttBar_norm.append(error_value)

    error_definition_d2Sig_dmttBar_dyttBar_norm = {}
    error_definition_d2Sig_dmttBar_dyttBar_norm['sys'] = {
        'definition': 'total systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    }
    for i in range(16):
        error_definition_d2Sig_dmttBar_dyttBar_norm['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }

    data_central_d2Sig_dmttBar_dyttBar_norm_yaml = {
        'data_central': data_central_d2Sig_dmttBar_dyttBar_norm
    }
    kinematics_d2Sig_dmttBar_dyttBar_norm_yaml = {'bins': kin_d2Sig_dmttBar_dyttBar_norm}
    uncertainties_d2Sig_dmttBar_dyttBar_norm_yaml = {
        'definitions': error_definition_d2Sig_dmttBar_dyttBar_norm,
        'bins': error_d2Sig_dmttBar_dyttBar_norm,
    }

    with open('data_d2Sig_dmttBar_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(data_central_d2Sig_dmttBar_dyttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dmttBar_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(kinematics_d2Sig_dmttBar_dyttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dmttBar_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dmttBar_dyttBar_norm_yaml, file, sort_keys=False)


processData()
