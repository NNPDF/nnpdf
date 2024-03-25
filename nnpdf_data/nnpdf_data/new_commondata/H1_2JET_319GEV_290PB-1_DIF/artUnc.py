from math import sqrt

import numpy as np
from numpy.linalg import eig
import yaml


def percentage_to_absolute(percentage, value):

    if type(percentage) is str:
        percentage = float(percentage.replace("%", ""))
        absolute = percentage * value * 0.01
        return absolute
    else:
        absolute = percentage * value * 0.01
        return absolute


def covmat_to_artunc(ndata, covmat_list, no_of_norm_mat=0):

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


def artunc():

    with open('rawdata/data49.yaml', 'r') as file:
        corMatFile = yaml.safe_load(file)

    corMatHalfArr = []
    for i in range(4656):
        corMatHalfArr.append(float(corMatFile['dependent_variables'][0]['values'][i]['value']))

    errPercArr = []
    dataArr = []
    for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]:
        hepdata_tables = "rawdata/data" + str(i) + ".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
        values = input['dependent_variables'][0]['values']
        for j in range(len(values)):
            errPerc = values[j]['errors'][0]['symerror']
            errPercArr.append(errPerc)
            dataArr.append(float(values[j]['value']))

    errArr = []
    for i in range(96):
        errArr.append(percentage_to_absolute(errPercArr[i], dataArr[i]))

    covMat = np.zeros((96, 96))
    artUnc = np.zeros((96, 96))

    for i in range(96):
        for j in range(i + 1):
            cmhap = (i * (i + 1)) // 2 + j
            if i == j:
                covMat[i][j] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]
            else:
                covMat[i][j] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]
                covMat[j][i] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]

    covMatList = []
    for i in range(96):
        for j in range(96):
            covMatList.append(covMat[i][j])
    artUnc = covmat_to_artunc(96, covMatList, 0)

    return artUnc


def artunc_norm():

    with open('rawdata/data50.yaml', 'r') as file:
        corMatFile = yaml.safe_load(file)

    corMatHalfArr = []
    for i in range(4656):
        corMatHalfArr.append(float(corMatFile['dependent_variables'][0]['values'][i]['value']))

    errPercArr = []
    dataArr = []
    for i in [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]:
        hepdata_tables = "rawdata/data" + str(i) + ".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
        values = input['dependent_variables'][0]['values']
        for j in range(len(values)):
            errPerc = values[j]['errors'][0]['symerror']
            errPercArr.append(errPerc)
            dataArr.append(float(values[j]['value']))

    errArr = []
    for i in range(96):
        errArr.append(percentage_to_absolute(errPercArr[i], dataArr[i]))

    covMat = np.zeros((96, 96))
    artUnc = np.zeros((96, 96))

    for i in range(96):
        for j in range(i + 1):
            cmhap = (i * (i + 1)) // 2 + j
            if i == j:
                covMat[i][j] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]
            else:
                covMat[i][j] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]
                covMat[j][i] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]

    covMatList = []
    for i in range(96):
        for j in range(96):
            covMatList.append(covMat[i][j])
    artUnc = covmat_to_artunc(96, covMatList, 1)

    return artUnc
