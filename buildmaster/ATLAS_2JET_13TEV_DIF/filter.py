# implemented by Tanishq Sharma

import yaml
from validphys.commondata_utils import symmetrize_errors as se

def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    tables = metadata['implemented_observables'][0]['tables']

    data_central = []
    kin = []
    error = []

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
        hepdata_tables="rawdata/atlas_mjj_jet2015_r04_ystar"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)

        values = input['dependent_variables'][0]['values']
        sqrt_s = input['dependent_variables'][0]['qualifiers'][1]['value']

        for j in range(len(values)):
            m_jj_min = input['independent_variables'][0]['values'][j]['low']
            m_jj_max = input['independent_variables'][0]['values'][j]['high']
            value_delta = 0
            error_value = {}
            for k in range(len(values[j]['errors'])):
                se_delta, se_sigma = se(values[j]['errors'][k]['asymerror']['plus'], values[j]['errors'][k]['asymerror']['minus'])
                value_delta = value_delta + se_delta
                error_label = str(values[j]['errors'][k]['label'])
                error_value[error_label] = se_sigma
            data_central_value = values[j]['value'] + value_delta
            data_central.append(data_central_value)
            error.append(error_value)
            kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'm_jj_sqr': {'min': m_jj_min**2, 'mid': None, 'max': m_jj_max**2}, 'y': {'min': y_min, 'mid': y_central, 'max': y_max}}
            kin.append(kin_value)

    hepdata_tables="rawdata/atlas_mjj_jet2015_r04_ystar1.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    error_definition = {}
    error_definition['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
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

processData()
