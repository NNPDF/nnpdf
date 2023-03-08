#    After adding a new util function, run the following command
#    in the terminal in buildmaster directory:
#    ....buildmaster$ for d in */; do cp utils.py "$d"; done
#    to ensure the utils file is consistent across all folders.
#    Existing functions should not be modified or removed as that might 
#    break already implemented experiments' filter scripts.

from math import sqrt
import numpy as np
from numpy.linalg import eig

def symmetrize_errors(delta_plus, delta_minus):
    semi_diff = (delta_plus + delta_minus)/2
    average = (delta_plus - delta_minus)/2
    se_delta = semi_diff
    se_sigma = sqrt(average*average + 2*semi_diff*semi_diff)
    return se_delta, se_sigma

def percentage_to_absolute(percentage, value):
    percentage = float(percentage.replace("%", ""))
    absolute = percentage * value * 0.01
    return absolute 

def corMat_to_covMat(ndata, errArray, corMatArray):
    covMatArray = []
    for i in range(len(corMatArray)):
        a = i // ndata
        b = i % ndata
        covMatArray.append(corMatArray[i] * errArray[a] * errArray[b])
    return covMatArray

def covMat_to_artUnc(ndata, covMatArray, is_normalized):
    epsilon = -0.0000000001
    negEValCount = 0
    psdCheck = True
    covMat = np.zeros((ndata, ndata))
    artUnc = np.zeros((ndata, ndata))
    for i in range(len(covMatArray)):
        a = i // ndata
        b = i % ndata
        covMat[a][b] = covMatArray[i]
    eigVal, eigVec = eig(covMat)
    if is_normalized == False:
        for i in range(ndata):
            for j in range(ndata):
                if eigVal[j] <= 0:
                    raise ValueError('The covariance matrix is not positive-semidefinite')
                else:
                    artUnc[i][j] = eigVec[i][j] * sqrt(eigVal[j])
    elif is_normalized == True:
        for j in range(len(eigVal)):
            if eigVal[j] < epsilon:
                psdCheck = False
            elif eigVal[j] > epsilon and eigVal[j] <= 0:
                negEValCount = negEValCount + 1
                if negEValCount == 2:
                    psdCheck = False
            elif eigVal[j] > 0:
                continue
        if psdCheck == False:
            raise ValueError('The covariance matrix is not positive-semidefinite')
        else:
            for i in range(ndata):
                for j in range(ndata):
                    if eigVal[j] < 0:
                        continue
                    else:
                        artUnc[i][j] = eigVec[i][j] * sqrt(eigVal[j]) 
    return artUnc
