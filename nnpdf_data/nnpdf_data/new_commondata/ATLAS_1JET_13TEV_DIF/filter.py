from math import sqrt

import yaml


def symmetrize_errors(delta_plus, delta_minus):
    r"""Compute the symmterized uncertainty and the shift in data point.
    Parameters
    ----------
    delta_plus : float
        The top/plus uncertainty with sign
    delta_minus : float
        The bottom/minus uncertainty with sign

    Returns
    -------
    se_delta : float
        The value to be added to the data point
    se_sigma : float
        The symmetrized uncertainty to be used in commondata
    """
    semi_diff = (delta_plus + delta_minus) / 2
    average = (delta_plus - delta_minus) / 2
    se_delta = semi_diff
    se_sigma = sqrt(average * average + 2 * semi_diff * semi_diff)
    return se_delta, se_sigma


def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    tables = metadata['implemented_observables'][0]['tables']
    tables_altcorr1 = metadata['implemented_observables'][1]['tables']

    data_central = []
    kin = []
    error = []
    data_central_altcorr1 = []
    kin_altcorr1 = []
    error_altcorr1 = []

    # jet data

    for i in tables:
        if i == 1:
            y_min = 0
            y_max = 0.5
        elif i == 2:
            y_min = 0.5
            y_max = 1
        elif i == 3:
            y_min = 1
            y_max = 1.5
        elif i == 4:
            y_min = 1.5
            y_max = 2
        elif i == 5:
            y_min = 2
            y_max = 2.5
        elif i == 6:
            y_min = 2.5
            y_max = 3
        y_central = None
        hepdata_tables = "rawdata/atlas_inclusive_jet2015_r04_eta" + str(i) + ".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)

        values = input['dependent_variables'][0]['values']
        sqrts = input['dependent_variables'][0]['qualifiers'][1]['value']

        for j in range(len(values)):
            pT_min = input['independent_variables'][0]['values'][j]['low']
            pT_max = input['independent_variables'][0]['values'][j]['high']
            value_delta = 0
            error_value = {}
            for k in range(len(values[j]['errors'])):
                se_delta, se_sigma = symmetrize_errors(
                    values[j]['errors'][k]['asymerror']['plus'],
                    values[j]['errors'][k]['asymerror']['minus'],
                )
                value_delta = value_delta + se_delta
                error_label = str(values[j]['errors'][k]['label'])
                error_value[error_label] = se_sigma
            data_central_value = values[j]['value'] + value_delta
            data_central.append(data_central_value)
            error.append(error_value)
            kin_value = {
                'sqrts': {'min': None, 'mid': sqrts, 'max': None},
                'pT': {'min': pT_min, 'mid': None, 'max': pT_max},
                'y': {'min': y_min, 'mid': y_central, 'max': y_max},
            }
            kin.append(kin_value)

    hepdata_tables = "rawdata/atlas_inclusive_jet2015_r04_eta1.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    error_definition = {}
    error_definition['stat'] = {
        'description': 'total statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
        error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
        error_definition[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}

    data_central_yaml = {'data_central': data_central}
    kinematics_yaml = {'bins': kin}
    uncertainties_yaml = {'definitions': error_definition, 'bins': error}

    with open('data.yaml', 'w') as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    with open('uncertainties.yaml', 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)

    # jet altcorr1 data

    for i in tables_altcorr1:
        if i == 1:
            y_min = 0
            y_max = 0.5
        elif i == 2:
            y_min = 0.5
            y_max = 1
        elif i == 3:
            y_min = 1
            y_max = 1.5
        elif i == 4:
            y_min = 1.5
            y_max = 2
        elif i == 5:
            y_min = 2
            y_max = 2.5
        elif i == 6:
            y_min = 2.5
            y_max = 3
        y_central = None
        hepdata_tables = "rawdata/atlas_inclusive_jet2015_r04_altcorr1_eta" + str(i) + ".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)

        values = input['dependent_variables'][0]['values']
        sqrts = input['dependent_variables'][0]['qualifiers'][1]['value']

        for j in range(len(values)):
            pT_min = input['independent_variables'][0]['values'][j]['low']
            pT_max = input['independent_variables'][0]['values'][j]['high']
            value_delta = 0
            error_value = {}
            for k in range(len(values[j]['errors'])):
                se_delta, se_sigma = symmetrize_errors(
                    values[j]['errors'][k]['asymerror']['plus'],
                    values[j]['errors'][k]['asymerror']['minus'],
                )
                value_delta = value_delta + se_delta
                error_label = str(values[j]['errors'][k]['label'])
                error_value[error_label] = se_sigma
            data_central_value = values[j]['value'] + value_delta
            data_central_altcorr1.append(data_central_value)
            error_altcorr1.append(error_value)
            kin_value = {
                'sqrts': {'min': None, 'mid': sqrts, 'max': None},
                'pT': {'min': pT_min, 'mid': None, 'max': pT_max},
                'y': {'min': y_min, 'mid': y_central, 'max': y_max},
            }
            kin_altcorr1.append(kin_value)

    hepdata_tables = "rawdata/atlas_inclusive_jet2015_r04_altcorr1_eta1.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    error_definition_altcorr1 = {}
    error_definition_altcorr1['stat'] = {
        'description': 'total statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
        error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
        error_definition_altcorr1[error_name] = {
            'description': '',
            'treatment': 'MULT',
            'type': 'CORR',
        }

    data_central_altcorr1_yaml = {'data_central': data_central_altcorr1}
    kinematics_altcorr1_yaml = {'bins': kin_altcorr1}
    uncertainties_altcorr1_yaml = {'definitions': error_definition_altcorr1, 'bins': error_altcorr1}

    with open('data_altcorr1.yaml', 'w') as file:
        yaml.dump(data_central_altcorr1_yaml, file, sort_keys=False)

    with open('kinematics_altcorr1.yaml', 'w') as file:
        yaml.dump(kinematics_altcorr1_yaml, file, sort_keys=False)

    with open('uncertainties_altcorr1.yaml', 'w') as file:
        yaml.dump(uncertainties_altcorr1_yaml, file, sort_keys=False)


processData()
