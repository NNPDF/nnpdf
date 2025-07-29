import yaml

from nnpdf_data.filter_utils.utils import percentage_to_absolute as pta
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def magic1(table, ndat, var_name):
    with open(table, 'r') as f:
        input = yaml.safe_load(f)

    data_central = []
    kin = []
    error = []

    values = input['dependent_variables'][1]['values']

    for i in range(ndat):
        kin_min = input['independent_variables'][0]['values'][i]['low']
        kin_max = input['independent_variables'][0]['values'][i]['high']
        if 'value' in input['independent_variables'][0]['values'][i]:
            kin_mid = input['independent_variables'][0]['values'][i]['value']
        else:
            kin_mid = (kin_min + kin_max) / 2

        kin_value = {var_name: {'min': kin_min, 'mid': kin_mid, 'max': kin_max}}

        data_central_value = values[i]['value']
        error_value = {}
        error_value['stat'] = values[i]['errors'][0]['symerror']
        error_value['sys'] = values[i]['errors'][1]['symerror']
        error_value['sys_norm'] = pta(1, data_central_value)

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
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    error_definition['sys_norm'] = {
        'definition': 'systematic uncertainty (normalization)',
        'treatment': 'MULT',
        'type': 'CORR',
    }

    data_central_yaml = {'data_central': data_central}
    kin_yaml = {'bins': kin}
    uncertainties_yaml = {'definitions': error_definition, 'bins': error}

    return data_central_yaml, kin_yaml, uncertainties_yaml


def magic2(table, ndat, var_name):
    with open(table, 'r') as f:
        input = yaml.safe_load(f)

    data_central_uds = []
    kin_uds = []
    error_uds = []
    data_central_c = []
    kin_c = []
    error_c = []
    data_central_b = []
    kin_b = []
    error_b = []

    values_uds = input['dependent_variables'][0]['values']
    values_c = input['dependent_variables'][1]['values']
    values_b = input['dependent_variables'][2]['values']
    for i in range(ndat):
        kin_mid = input['independent_variables'][0]['values'][i]['value']

        kin_value = {var_name: {'min': None, 'mid': kin_mid, 'max': None}}

        data_central_uds_value = values_uds[i]['value']
        error_uds_value = {}
        error_uds_value['error'] = values_uds[i]['errors'][0]['symerror']
        error_uds_value['sys_norm'] = pta(1, data_central_uds_value)

        data_central_c_value = values_c[i]['value']
        error_c_value = {}
        error_c_value['error'] = values_c[i]['errors'][0]['symerror']
        error_c_value['sys_norm'] = pta(1, data_central_c_value)

        data_central_b_value = values_b[i]['value']
        error_b_value = {}
        error_b_value['error'] = values_b[i]['errors'][0]['symerror']
        error_b_value['sys_norm'] = pta(1, data_central_b_value)

        kin_uds.append(kin_value)
        data_central_uds.append(data_central_uds_value)
        error_uds.append(error_uds_value)

        kin_c.append(kin_value)
        data_central_c.append(data_central_c_value)
        error_c.append(error_c_value)

        kin_b.append(kin_value)
        data_central_b.append(data_central_b_value)
        error_b.append(error_b_value)

    error_definition = {}
    error_definition['error'] = {
        'definition': 'total uncertainty added in quadrature',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    error_definition['sys_norm'] = {
        'definition': 'systematic uncertainty (normalization)',
        'treatment': 'MULT',
        'type': 'CORR',
    }

    data_central_uds_yaml = {'data_central': data_central_uds}
    kin_uds_yaml = {'bins': kin_uds}
    uncertainties_uds_yaml = {'definitions': error_definition, 'bins': error_uds}

    data_central_c_yaml = {'data_central': data_central_c}
    kin_c_yaml = {'bins': kin_c}
    uncertainties_c_yaml = {'definitions': error_definition, 'bins': error_c}

    data_central_b_yaml = {'data_central': data_central_b}
    kin_b_yaml = {'bins': kin_b}
    uncertainties_b_yaml = {'definitions': error_definition, 'bins': error_b}

    return (
        data_central_uds_yaml,
        kin_uds_yaml,
        uncertainties_uds_yaml,
        data_central_c_yaml,
        kin_c_yaml,
        uncertainties_c_yaml,
        data_central_b_yaml,
        kin_b_yaml,
        uncertainties_b_yaml,
    )
