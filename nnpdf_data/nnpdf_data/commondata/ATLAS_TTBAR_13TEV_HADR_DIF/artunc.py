import numpy as np
import yaml

from nnpdf_data.filter_utils.utils import concat_matrices as cm
from nnpdf_data.filter_utils.utils import covmat_to_artunc as cta
from nnpdf_data.filter_utils.utils import matlist_to_matrix as mtm

covmat_mtt = []
covmat_mtt_norm = []
covmat_ytt = []
covmat_ytt_norm = []
covmat_501 = []
covmat_502 = []
covmat_503 = []
covmat_504 = []
covmat_505 = []
covmat_506 = []
covmat_492 = []
covmat_493 = []
covmat_494 = []
covmat_495 = []
covmat_496 = []
covmat_497 = []

covariance_matrix_mtt = "rawdata/Table464.yaml"
with open(covariance_matrix_mtt, 'r') as file1:
    input1 = yaml.safe_load(file1)
for i in range(81):
    covmat_mtt.append(input1['dependent_variables'][0]['values'][i]['value'])
artunc_mtt = cta(9, covmat_mtt)


covariance_matrix_mtt_norm = "rawdata/Table462.yaml"
with open(covariance_matrix_mtt_norm, 'r') as file2:
    input2 = yaml.safe_load(file2)
for i in range(81):
    covmat_mtt_norm.append(input2['dependent_variables'][0]['values'][i]['value'])
artunc_mtt_norm = cta(9, covmat_mtt_norm, 1)


covariance_matrix_ytt = "rawdata/Table476.yaml"
with open(covariance_matrix_ytt, 'r') as file3:
    input3 = yaml.safe_load(file3)
for i in range(144):
    covmat_ytt.append(input3['dependent_variables'][0]['values'][i]['value'])
artunc_ytt = cta(12, covmat_ytt)

covariance_matrix_ytt_norm = "rawdata/Table474.yaml"
with open(covariance_matrix_ytt_norm, 'r') as file4:
    input4 = yaml.safe_load(file4)
for i in range(144):
    covmat_ytt_norm.append(input4['dependent_variables'][0]['values'][i]['value'])
artunc_ytt_norm = cta(12, covmat_ytt_norm, 1)

# abs mtt-ytt
#        m1  |  m2  |  m3
#  m1 | 501    502t   504t
#  m2 | 502    503    505t
#  m3 | 504    505    506

with open("rawdata/Table501.yaml", 'r') as file5:
    input5 = yaml.safe_load(file5)
for i in range(16):
    covmat_501.append(input5['dependent_variables'][0]['values'][i]['value'])
covmat_501 = mtm(4, 4, covmat_501)
with open("rawdata/Table502.yaml", 'r') as file6:
    input6 = yaml.safe_load(file6)
for i in range(16):
    covmat_502.append(input6['dependent_variables'][0]['values'][i]['value'])
covmat_502 = mtm(4, 4, covmat_502)
covmat_502t = np.transpose(covmat_502)
with open("rawdata/Table503.yaml", 'r') as file7:
    input7 = yaml.safe_load(file7)
for i in range(16):
    covmat_503.append(input7['dependent_variables'][0]['values'][i]['value'])
covmat_503 = mtm(4, 4, covmat_503)
with open("rawdata/Table504.yaml", 'r') as file8:
    input8 = yaml.safe_load(file8)
for i in range(12):
    covmat_504.append(input8['dependent_variables'][0]['values'][i]['value'])
covmat_504 = mtm(3, 4, covmat_504)
covmat_504t = np.transpose(covmat_504)
with open("rawdata/Table505.yaml", 'r') as file9:
    input9 = yaml.safe_load(file9)
for i in range(12):
    covmat_505.append(input9['dependent_variables'][0]['values'][i]['value'])
covmat_505 = mtm(3, 4, covmat_505)
covmat_505t = np.transpose(covmat_505)
with open("rawdata/Table506.yaml", 'r') as file10:
    input10 = yaml.safe_load(file10)
for i in range(9):
    covmat_506.append(input10['dependent_variables'][0]['values'][i]['value'])
covmat_506 = mtm(3, 3, covmat_506)

covmat_mtt_ytt = cm(
    3,
    3,
    [
        covmat_501,
        covmat_502t,
        covmat_504t,
        covmat_502,
        covmat_503,
        covmat_505t,
        covmat_504,
        covmat_505,
        covmat_506,
    ],
)
artunc_mtt_ytt = cta(11, covmat_mtt_ytt)

# abs mtt-ytt-norm
#        m1  |  m2  |  m3
#  m1 | 492    493t   495t
#  m2 | 493    494    496t
#  m3 | 495    496    497

with open("rawdata/Table492.yaml", 'r') as file11:
    input11 = yaml.safe_load(file11)
for i in range(16):
    covmat_492.append(input11['dependent_variables'][0]['values'][i]['value'])
covmat_492 = mtm(4, 4, covmat_492)
with open("rawdata/Table493.yaml", 'r') as file12:
    input12 = yaml.safe_load(file12)
for i in range(16):
    covmat_493.append(input12['dependent_variables'][0]['values'][i]['value'])
covmat_493 = mtm(4, 4, covmat_493)
covmat_493t = np.transpose(covmat_493)
with open("rawdata/Table494.yaml", 'r') as file13:
    input13 = yaml.safe_load(file13)
for i in range(16):
    covmat_494.append(input13['dependent_variables'][0]['values'][i]['value'])
covmat_494 = mtm(4, 4, covmat_494)
with open("rawdata/Table495.yaml", 'r') as file14:
    input14 = yaml.safe_load(file14)
for i in range(12):
    covmat_495.append(input14['dependent_variables'][0]['values'][i]['value'])
covmat_495 = mtm(3, 4, covmat_495)
covmat_495t = np.transpose(covmat_495)
with open("rawdata/Table496.yaml", 'r') as file15:
    input15 = yaml.safe_load(file15)
for i in range(12):
    covmat_496.append(input15['dependent_variables'][0]['values'][i]['value'])
covmat_496 = mtm(3, 4, covmat_496)
covmat_496t = np.transpose(covmat_496)
with open("rawdata/Table497.yaml", 'r') as file16:
    input16 = yaml.safe_load(file16)
for i in range(9):
    covmat_497.append(input16['dependent_variables'][0]['values'][i]['value'])
covmat_497 = mtm(3, 3, covmat_497)

covmat_mtt_ytt_norm = cm(
    3,
    3,
    [
        covmat_492,
        covmat_493t,
        covmat_495t,
        covmat_493,
        covmat_494,
        covmat_496t,
        covmat_495,
        covmat_496,
        covmat_497,
    ],
)
artunc_mtt_ytt_norm = cta(11, covmat_mtt_ytt_norm, 1)

#########################################
#            lumiless covmat            #
#########################################


def lumiless_artUnc(data_array, covmat_array):

    lumi_uncert_percent = 2.1

    lumi_uncert_array = np.array([data_array * lumi_uncert_percent / 100])
    lumi_covmat = lumi_uncert_array.T @ lumi_uncert_array
    lumi_covmat_array = lumi_covmat.flatten()
    lumiless_covmat = (covmat_array - lumi_covmat_array).tolist()
    artUnc = cta(len(data_array), lumiless_covmat)

    return artUnc


mtt_data = np.array([])
with open("rawdata/Table463.yaml", 'r') as file17:
    input17 = yaml.safe_load(file17)
for i in range(len(input17['dependent_variables'][0]['values'])):
    mtt_data = np.append(mtt_data, input17['dependent_variables'][0]['values'][i]['value'])

ytt_data = np.array([])
with open("rawdata/Table475.yaml", 'r') as file18:
    input18 = yaml.safe_load(file18)
for i in range(len(input18['dependent_variables'][0]['values'])):
    ytt_data = np.append(ytt_data, input18['dependent_variables'][0]['values'][i]['value'])

mtt_ytt_data = np.array([])
with open("rawdata/Table498.yaml", 'r') as file19:
    input19 = yaml.safe_load(file19)
for i in range(len(input19['dependent_variables'][0]['values'])):
    mtt_ytt_data = np.append(mtt_ytt_data, input19['dependent_variables'][0]['values'][i]['value'])
with open("rawdata/Table499.yaml", 'r') as file20:
    input20 = yaml.safe_load(file20)
for i in range(len(input20['dependent_variables'][0]['values'])):
    mtt_ytt_data = np.append(mtt_ytt_data, input20['dependent_variables'][0]['values'][i]['value'])
with open("rawdata/Table500.yaml", 'r') as file21:
    input21 = yaml.safe_load(file21)
for i in range(len(input21['dependent_variables'][0]['values'])):
    mtt_ytt_data = np.append(mtt_ytt_data, input21['dependent_variables'][0]['values'][i]['value'])

artunc_mtt_lumiless = lumiless_artUnc(mtt_data, covmat_mtt)
artunc_ytt_lumiless = lumiless_artUnc(ytt_data, covmat_ytt)
artunc_mtt_ytt_lumiless = lumiless_artUnc(mtt_ytt_data, covmat_mtt_ytt)
