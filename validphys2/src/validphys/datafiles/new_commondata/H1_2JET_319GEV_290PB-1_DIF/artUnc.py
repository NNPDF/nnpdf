import yaml
import numpy
from validphys.commondata_utils import covmat_to_artunc as cta
from validphys.commondata_utils import percentage_to_absolute as pta

def artunc():

    with open('rawdata/data49.yaml', 'r') as file:
        corMatFile = yaml.safe_load(file)

    corMatHalfArr = []
    for i in range(4656):
        corMatHalfArr.append(float(corMatFile['dependent_variables'][0]['values'][i]['value']))

    errPercArr = []
    dataArr = []
    for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]:
        hepdata_tables="rawdata/data"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
        values = input['dependent_variables'][0]['values']
        for j in range(len(values)):
            errPerc = values[j]['errors'][0]['symerror']
            errPercArr.append(errPerc)
            dataArr.append(float(values[j]['value']))


    errArr = []
    for i in range(96):
        errArr.append(pta(errPercArr[i], dataArr[i]))

    covMat = numpy.zeros((96, 96))
    artUnc = numpy.zeros((96, 96))

    for i in range(96):
        for j in range(i+1):
            cmhap = (i * (i+1)) // 2 + j
            if i == j:
                covMat[i][j] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]
            else:
                covMat[i][j] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]
                covMat[j][i] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]

    covMatList = []
    for i in range(96):
        for j in range(96):
            covMatList.append(covMat[i][j])
    artUnc = cta(96, covMatList, 0)

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
        hepdata_tables="rawdata/data"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
        values = input['dependent_variables'][0]['values']
        for j in range(len(values)):
            errPerc = values[j]['errors'][0]['symerror']
            errPercArr.append(errPerc)
            dataArr.append(float(values[j]['value']))

    errArr = []
    for i in range(96):
        errArr.append(pta(errPercArr[i], dataArr[i]))   

    covMat = numpy.zeros((96, 96))
    artUnc = numpy.zeros((96, 96))

    for i in range(96):
        for j in range(i+1):
            cmhap = (i * (i+1)) // 2 + j
            if i == j:
                covMat[i][j] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]
            else:
                covMat[i][j] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]
                covMat[j][i] = corMatHalfArr[cmhap] * errArr[i] * errArr[j]

    covMatList = []
    for i in range(96):
        for j in range(96):
            covMatList.append(covMat[i][j])
    artUnc = cta(96, covMatList, 1)

    return artUnc
