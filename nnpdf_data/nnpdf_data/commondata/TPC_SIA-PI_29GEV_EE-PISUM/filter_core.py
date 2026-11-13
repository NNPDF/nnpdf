import yaml

from nnpdf_data.filter_utils.utils import percentage_to_absolute as pta
from nnpdf_data.filter_utils.utils import prettify_float
from nnpdf_data.filter_utils.utils import symmetrize_errors as se

yaml.add_representer(float, prettify_float)


def magic(table, ndat, var_name, index):
    with open(table, 'r') as f:
        input = yaml.safe_load(f)

    data_central = []
    kin = []
    error = []

    values = input['dependent_variables'][index]['values']
    total_bins = len(values)

    for i in range(total_bins):
        if values[i]['value'] == '-':
            continue
        else:
            kin_min = input['independent_variables'][0]['values'][i]['low']
            kin_max = input['independent_variables'][0]['values'][i]['high']

            kin_value = {var_name: {'min': kin_min, 'mid': None, 'max': kin_max}}

            data_central_value = values[i]['value']
            error_value = {}
            error_value['error'] = values[i]['errors'][0]['symerror']

            kin.append(kin_value)
            data_central.append(data_central_value)
            error.append(error_value)

    if len(data_central) != ndat:
        raise ValueError(
            f"Number of data points {len(data_central)} does not match expected {ndat}"
        )

    error_definition = {}
    error_definition['error'] = {
        'definition': 'all uncertainties added in quadrature',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }

    data_central_yaml = {'data_central': data_central}
    kin_yaml = {'bins': kin}
    uncertainties_yaml = {'definitions': error_definition, 'bins': error}

    return data_central_yaml, kin_yaml, uncertainties_yaml
