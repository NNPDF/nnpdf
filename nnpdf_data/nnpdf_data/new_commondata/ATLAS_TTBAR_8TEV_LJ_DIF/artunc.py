import yaml

import numpy as np
from numpy.linalg import eig
# use #1693
from nnpdf_data.new_commondata.ATLAS_TTBAR_13TEV_HADR_DIF.utils import cormat_to_covmat as ctc
from nnpdf_data.new_commondata.ATLAS_TTBAR_13TEV_HADR_DIF.utils import covmat_to_artunc as cta
from nnpdf_data.new_commondata.ATLAS_TTBAR_13TEV_HADR_DIF.utils import percentage_to_absolute as pta
from nnpdf_data.new_commondata.ATLAS_TTBAR_13TEV_HADR_DIF.utils import concat_matrices as cm
from nnpdf_data.new_commondata.ATLAS_TTBAR_13TEV_HADR_DIF.utils import matlist_to_matrix as mtm

def artunc():
    statArr = []
    for i in [23, 29, 31, 27]:
        with open('rawdata/Table_'+str(i)+'.yaml', 'r') as file:
            input = yaml.safe_load(file)
        for j in range(len(input['dependent_variables'][0]['values'])):
            datval = input['dependent_variables'][0]['values'][j]['value']
            statperc = input['dependent_variables'][0]['values'][j]['errors'][0]['symerror']
            statArr.append(pta(statperc, datval))
    
#         mttbar(7)|  pTt (8)|  yt(5)|  yttbar(5)
# mttbar|   179       174t    170t    177t
# pTt   |   174       172     168t    173
# yt    |   170       168     167     169
# yttbar|   177       173t    169t    176
    ml179, ml174, ml170, ml177, ml172, ml168, ml167, ml173, ml169, ml176 = ([] for i in range(10))
    
    with open('rawdata/Table_179.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml179.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_174.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml174.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_170.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml170.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_177.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml177.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_172.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml172.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_168.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml168.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_167.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml167.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_173.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml173.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_169.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml169.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_176.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml176.append(input['dependent_variables'][0]['values'][i]['value'])
    
    mat179 = mtm(7, 7, ml179)
    mat174 = mtm(8, 7, ml174)
    mat174t = mat174.transpose()
    mat170 = mtm(5, 7, ml170)
    mat170t = mat170.transpose()
    mat177 = mtm(5, 7, ml177)
    mat177t = mat177.transpose()
    mat172 = mtm(8, 8, ml172)
    mat168 = mtm(5, 8, ml168)
    mat168t = mat168.transpose()
    mat167 = mtm(5, 5, ml167)
    mat173 = mtm(8, 5, ml173)
    mat173t = mat173.transpose()
    mat169 = mtm(5, 5, ml169)
    mat169t = mat169.transpose()
    mat176 = mtm(5, 5, ml176)

    cormatlist = cm(4, 4, [mat179, mat174t, mat170t, mat177t, mat174, mat172, mat168t, mat173, 
                           mat170, mat168, mat167, mat169, mat177, mat173t, mat169t, mat176])

    covmatlist = ctc(statArr, cormatlist)
    artunc = cta(25, covmatlist)
    return artunc


def artunc_norm():
    statArr = []
    for i in [24, 30, 32, 28]:
        with open('rawdata/Table_'+str(i)+'.yaml', 'r') as file:
            input = yaml.safe_load(file)
        for j in range(len(input['dependent_variables'][0]['values'])):
            datval = input['dependent_variables'][0]['values'][j]['value']
            statperc = input['dependent_variables'][0]['values'][j]['errors'][0]['symerror']
            statArr.append(pta(statperc, datval))
    
#         mttbar(7)|  pTt (8)|  yt(5)|  yttbar(5)
# mttbar|   234       229t    225t    232t
# pTt   |   229       227     223t    228
# yt    |   225       223     222     224
# yttbar|   232       228t    224t    231
    ml234, ml229, ml225, ml232, ml227, ml223, ml222, ml228, ml224, ml231 = ([] for i in range(10))
    
    with open('rawdata/Table_234.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml234.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_229.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml229.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_225.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml225.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_232.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml232.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_227.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml227.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_223.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml223.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_222.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml222.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_228.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml228.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_224.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml224.append(input['dependent_variables'][0]['values'][i]['value'])
    
    with open('rawdata/Table_231.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml231.append(input['dependent_variables'][0]['values'][i]['value'])
    
    mat234 = mtm(7, 7, ml234)
    mat229 = mtm(8, 7, ml229)
    mat229t = mat229.transpose()
    mat225 = mtm(5, 7, ml225)
    mat225t = mat225.transpose()
    mat232 = mtm(5, 7, ml232)
    mat232t = mat232.transpose()
    mat227 = mtm(8, 8, ml227)
    mat223 = mtm(5, 8, ml223)
    mat223t = mat223.transpose()
    mat222 = mtm(5, 5, ml222)
    mat228 = mtm(8, 5, ml228)
    mat228t = mat228.transpose()
    mat224 = mtm(5, 5, ml224)
    mat224t = mat224.transpose()
    mat231 = mtm(5, 5, ml231)

    cormatlist = cm(4, 4, [mat234, mat229t, mat225t, mat232t, mat229, mat227, mat223t, mat228, 
                           mat225, mat223, mat222, mat224, mat232, mat228t, mat224t, mat231])

    covmatlist = ctc(statArr, cormatlist)
    artunc = cta(25, covmatlist, 4)
    return artunc
