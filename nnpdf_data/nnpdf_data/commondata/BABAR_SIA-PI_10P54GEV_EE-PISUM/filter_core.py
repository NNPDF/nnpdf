import yaml

from nnpdf_data.filter_utils.utils import percentage_to_absolute as pta
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def magic(table, ndat, var_name, index):
    with open(table, 'r') as f:
        input = yaml.safe_load(f)

    data_central = []
    kin = []
    error = []

    values = input['dependent_variables'][index]['values']

    for i in range(ndat):
        kin_min = input['independent_variables'][0]['values'][i]['low']
        kin_max = input['independent_variables'][0]['values'][i]['high']
        
        kin_value = {var_name: {'min': kin_min, 'mid': None, 'max': kin_max}}

        data_central_value = values[i]['value']
        error_value = {}
        error_value['stat'] = values[i]['errors'][0]['symerror']
        error_value['sys'] = values[i]['errors'][1]['symerror']
        error_value['sys,norm'] = pta(values[i]['errors'][2]['symerror'], data_central_value)

        kin.append(kin_value)
        data_central.append(data_central_value)
        error.append(error_value)

    error_definition = {}
    error_definition['stat'] = {
        'definition': 'statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    error_definition['sys'] = {
        'definition': 'systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    }
    error_definition['sys,norm'] = {
        'definition': 'systematic uncertainty (normalization)',
        'treatment': 'MULT',
        'type': 'CORR',
    }

    data_central_yaml = {'data_central': data_central}
    kin_yaml = {'bins': kin}
    uncertainties_yaml = {'definitions': error_definition, 'bins': error}

    return data_central_yaml, kin_yaml, uncertainties_yaml
