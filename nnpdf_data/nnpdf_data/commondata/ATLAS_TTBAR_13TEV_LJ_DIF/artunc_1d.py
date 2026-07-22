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
    # mttbar|   803       801    802    810
    # pTt   |   801t      798    799    808
    # yt    |   802t      799t   800    809
    # yttbar|   810t      808t   809t   812
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
    mat801 = mtm(9, 8, ml801)
    mat801t = mat801.transpose()
    mat802 = mtm(9, 5, ml802)
    mat802t = mat802.transpose()
    mat810t = mtm(7, 9, ml810)
    mat810 = mat810t.transpose()
    mat798 = mtm(8, 8, ml798)
    mat799t = mtm(5, 8, ml799)
    mat799 = mat799t.transpose()
    mat808t = mtm(7, 8, ml808)
    mat808 = mat808t.transpose()
    mat800 = mtm(5, 5, ml800)
    mat809t = mtm(7, 5, ml809)
    mat809 = mat809t.transpose()
    mat812 = mtm(7, 7, ml812)

    cormatlist = cm(
        4,
        4,
        [
            mat803,
            mat801,
            mat802,
            mat810,
            mat801t,
            mat798,
            mat799,
            mat808,
            mat802t,
            mat799t,
            mat800,
            mat809,
            mat810t,
            mat808t,
            mat809t,
            mat812,
        ],
    )

    covmatlist_stat = ctc(statArr, cormatlist)
    covmat_stat = mtm(29, 29, covmatlist_stat)
    artunc = cta(29, covmatlist_stat)
    return covmat_stat, artunc
