import yaml

from nnpdf_data.filter_utils.utils import percentage_to_absolute as pta
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def magic(table, var_name):
    with open(table, 'r') as f:
        input = yaml.safe_load(f)

    data_central = []
    kin = []
    error = []
    ndat = len(input['independent_variables'][0]['values'])

    values = input['dependent_variables'][0]['values']

    for i in range(ndat):
        kin_mid = input['independent_variables'][0]['values'][i]['value']
        
        kin_value = {var_name: {'min': None, 'mid': kin_mid, 'max': None}}

        data_central_value = values[i]['value']
        error_value = {}
        error_value['error'] = values[i]['errors'][0]['symerror']
        error_value['sys_norm'] = pta(values[i]['errors'][1]['symerror'], data_central_value)

        kin.append(kin_value)
        data_central.append(data_central_value)
        error.append(error_value)

    error_definition = {}
    error_definition['error'] = {
        'definition': 'total uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    error_definition['sys_norm'] = {
        'definition': 'systematic uncertainty - normalisation',
        'treatment': 'MULT',
        'type': 'CORR',
    }

    data_central_yaml = {'data_central': data_central}
    kin_yaml = {'bins': kin}
    uncertainties_yaml = {'definitions': error_definition, 'bins': error}

    return data_central_yaml, kin_yaml, uncertainties_yaml
