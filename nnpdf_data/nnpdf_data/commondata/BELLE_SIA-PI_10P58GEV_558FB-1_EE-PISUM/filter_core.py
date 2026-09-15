import yaml

from nnpdf_data.filter_utils.utils import prettify_float
from nnpdf_data.filter_utils.utils import symmetrize_errors as se

yaml.add_representer(float, prettify_float)


def magic(table, var_name):

    data_dict = {}

    with open(table, 'r') as f:
        for _ in range(2):
            next(f)
        for i, line in enumerate(f):
            data_dict['bin' + str(i)] = []
            for elements in line.split():
                data_dict['bin' + str(i)].append(float(elements))

    data_central = []
    kin = []
    error = []

    ndat = len(data_dict)

    for i in range(ndat):
        kin_min = data_dict['bin' + str(i)][0]
        kin_max = data_dict['bin' + str(i)][1]

        kin_value = {var_name: {'min': kin_min, 'mid': None, 'max': kin_max}}

        data_central_value = data_dict['bin' + str(i)][2]
        error_value = {}
        error_value['stat'] = data_dict['bin' + str(i)][3]
        data_shift1, error_value['sys_1'] = se(
            data_dict['bin' + str(i)][4], data_dict['bin' + str(i)][5]
        )
        data_shift2, error_value['sys_2'] = se(
            data_dict['bin' + str(i)][6], data_dict['bin' + str(i)][7]
        )
        data_central_value = data_central_value + data_shift1 + data_shift2

        kin.append(kin_value)
        data_central.append(data_central_value)
        error.append(error_value)

    error_definition = {}
    error_definition['stat'] = {
        'definition': 'statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    error_definition['sys_1'] = {
        'definition': 'systematic uncertainty (correlated)',
        'treatment': 'MULT',
        'type': 'CORR',
    }
    error_definition['sys_2'] = {
        'definition': 'systematic uncertainty (uncorrelated)',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }

    data_central_yaml = {'data_central': data_central}
    kin_yaml = {'bins': kin}
    uncertainties_yaml = {'definitions': error_definition, 'bins': error}

    return data_central_yaml, kin_yaml, uncertainties_yaml
