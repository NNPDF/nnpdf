import numpy as np
from numpy.linalg import eig
import yaml

from nnpdf_data.filter_utils.utils import concat_matrices as cm
from nnpdf_data.filter_utils.utils import cormat_to_covmat as ctc
from nnpdf_data.filter_utils.utils import covmat_to_artunc as cta
from nnpdf_data.filter_utils.utils import matlist_to_matrix as mtm
from nnpdf_data.filter_utils.utils import percentage_to_absolute as pta


def artunc():
    statArr = []
    for i in [618, 610, 614, 626]:
        with open('rawdata/Table' + str(i) + '.yaml', 'r') as file:
            input = yaml.safe_load(file)
        for j in range(len(input['dependent_variables'][0]['values'])):
            datval = input['dependent_variables'][0]['values'][j]['value']
            statperc = input['dependent_variables'][0]['values'][j]['errors'][0]['symerror']
            statArr.append(pta(statperc, datval))

    #         mttbar(9)|  pTt (8)|  yt(5)|  yttbar(7)
    # mttbar|   803       801t    802t    810t
    # pTt   |   801       798     799t    808t
    # yt    |   802       799     800     809t
    # yttbar|   810       808     809     812
    ml803, ml801, ml802, ml810, ml798, ml799, ml808, ml800, ml809, ml812 = ([] for i in range(10))

    with open('rawdata/Table803.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml803.append(input['dependent_variables'][0]['values'][i]['value'])

    with open('rawdata/Table801.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml801.append(input['dependent_variables'][0]['values'][i]['value'])

    with open('rawdata/Table802.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml802.append(input['dependent_variables'][0]['values'][i]['value'])

    with open('rawdata/Table810.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml810.append(input['dependent_variables'][0]['values'][i]['value'])

    with open('rawdata/Table798.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml798.append(input['dependent_variables'][0]['values'][i]['value'])

    with open('rawdata/Table799.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml799.append(input['dependent_variables'][0]['values'][i]['value'])

    with open('rawdata/Table808.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml808.append(input['dependent_variables'][0]['values'][i]['value'])

    with open('rawdata/Table800.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml800.append(input['dependent_variables'][0]['values'][i]['value'])

    with open('rawdata/Table809.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml809.append(input['dependent_variables'][0]['values'][i]['value'])

    with open('rawdata/Table812.yaml', 'r') as file:
        input = yaml.safe_load(file)
    for i in range(len(input['dependent_variables'][0]['values'])):
        ml812.append(input['dependent_variables'][0]['values'][i]['value'])

    mat803 = mtm(9, 9, ml803)
    mat801 = mtm(8, 9, ml801)
    mat801t = mat801.transpose()
    mat802 = mtm(5, 9, ml802)
    mat802t = mat802.transpose()
    mat810 = mtm(7, 9, ml810)
    mat810t = mat810.transpose()
    mat798 = mtm(8, 8, ml798)
    mat799 = mtm(5, 8, ml799)
    mat799t = mat799.transpose()
    mat808 = mtm(7, 8, ml808)
    mat808t = mat808.transpose()
    mat800 = mtm(5, 5, ml800)
    mat809 = mtm(7, 5, ml809)
    mat809t = mat809.transpose()
    mat812 = mtm(7, 7, ml812)

    cormatlist = cm(
        4,
        4,
        [
            mat803,
            mat801t,
            mat802t,
            mat810t,
            mat801,
            mat798,
            mat799t,
            mat808t,
            mat802,
            mat799,
            mat800,
            mat809t,
            mat810,
            mat808,
            mat809,
            mat812,
        ],
    )

    covmatlist = ctc(statArr, cormatlist)
    artunc = cta(29, covmatlist)
    return artunc


artunc()
