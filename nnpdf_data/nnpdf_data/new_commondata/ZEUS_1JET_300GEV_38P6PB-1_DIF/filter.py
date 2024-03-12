import yaml
from validphys.commondata_utils import symmetrize_errors as se

def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    tables_q2_et = metadata['implemented_observables'][0]['tables']

    data_central_q2_et = []
    kin_q2_et = []
    error_q2_et = []

# q2_et data

    for i in tables_q2_et:

        hepdata_tables="rawdata/Table"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)

        for j in 0, 1:
            if i == 4:
                if j == 0:
                    Q2_min = 125
                    Q2_max = 250
                elif j == 1:
                    Q2_min = 250
                    Q2_max = 500
            elif i == 5:
                if j == 0:
                    Q2_min = 500
                    Q2_max = 1000
                elif j == 1:
                    Q2_min = 1000
                    Q2_max = 2000
            elif i == 6:
                if j == 0:
                    Q2_min = 2000
                    Q2_max = 5000
                if j == 1:
                    Q2_min = 5000
                    Q2_max = 100000

            sqrts = float(input['dependent_variables'][j]['qualifiers'][5]['value'])
            values = input['dependent_variables'][j]['values']

            for k in range(len(values)):
                data_central_value = values[k]['value']
                ET_max = input['independent_variables'][0]['values'][k]['high']
                ET_min = input['independent_variables'][0]['values'][k]['low']
                kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'Q2': {'min': Q2_min, 'mid': None, 'max': Q2_max}, 'ET': {'min': ET_min, 'mid': None, 'max': ET_max}}
                kin_q2_et.append(kin_value)
                value_delta = 0
                error_value = {}
                if 'symerror' in values[k]['errors'][0]:
                    error_value['stat'] = values[k]['errors'][0]['symerror']
                else:
                    se_delta, se_sigma = se(values[k]['errors'][0]['asymerror']['plus'], values[k]['errors'][0]['asymerror']['minus'])
                    value_delta = value_delta + se_delta
                    error_value['stat'] = se_sigma
                if 'symerror' in values[k]['errors'][1]:
                    error_value['sys'] = values[k]['errors'][1]['symerror']
                else:
                    se_delta, se_sigma = se(values[k]['errors'][1]['asymerror']['plus'], values[k]['errors'][1]['asymerror']['minus'])
                    value_delta = value_delta + se_delta
                    error_value['sys'] = se_sigma
                if 'symerror' in values[k]['errors'][2]:
                    error_value['jet_es'] = values[k]['errors'][2]['symerror']
                else:
                    se_delta, se_sigma = se(values[k]['errors'][2]['asymerror']['plus'], values[k]['errors'][2]['asymerror']['minus'])
                    value_delta = value_delta + se_delta
                    error_value['jet_es'] = se_sigma
                data_central_value = data_central_value + value_delta
                data_central_q2_et.append(data_central_value)
                error_q2_et.append(error_value)

    error_definition_q2_et = {
        'stat': {'description': 'statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'sys': {'description': 'systematic uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'jet_es': {'description': 'jet energy scale uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    }

    data_central_q2_et_yaml = {'data_central': data_central_q2_et}
    kinematics_q2_et_yaml = {'bins': kin_q2_et}
    uncertainties_q2_et_yaml = {'definitions': error_definition_q2_et, 'bins': error_q2_et}

    with open('data_q2_et.yaml', 'w') as file:
        yaml.dump(data_central_q2_et_yaml, file, sort_keys=False)

    with open('kinematics_q2_et.yaml', 'w') as file:
        yaml.dump(kinematics_q2_et_yaml, file, sort_keys=False)

    with open('uncertainties_q2_et.yaml', 'w') as file:
        yaml.dump(uncertainties_q2_et_yaml, file, sort_keys=False)

processData()
