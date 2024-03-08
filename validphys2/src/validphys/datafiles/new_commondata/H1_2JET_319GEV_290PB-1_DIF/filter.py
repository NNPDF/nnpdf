import artUnc
import yaml
# use #1693
from validphys.commondata_utils import percentage_to_absolute as pta
from validphys.commondata_utils import symmetrize_errors as se
from math import sqrt

def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    tables = metadata['implemented_observables'][0]['tables']
    tables_norm = metadata['implemented_observables'][1]['tables']

    data_central = []
    kin = []
    error = []
    data_central_norm = []
    kin_norm = []
    error_norm = []

    artUncMatr = artUnc.artunc()
    artUncMatr_norm = artUnc.artunc_norm()

# dijet data

    for i in tables:
        if i == 9:
            Q2_min = 5.5
            Q2_max = 8
        elif i == 10:
            Q2_min = 8
            Q2_max = 11
        elif i == 11:
            Q2_min = 11
            Q2_max = 16
        elif i == 12:
            Q2_min = 16
            Q2_max = 22
        elif i == 13:
            Q2_min = 22
            Q2_max = 30
        elif i == 14:
            Q2_min = 30
            Q2_max = 42
        elif i == 15:
            Q2_min = 42
            Q2_max = 60
        elif i == 16:
            Q2_min = 60
            Q2_max = 80

        hepdata_tables="rawdata/data"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)

        sqrts = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
        values = input['dependent_variables'][0]['values']

        for j in range(len(values)):
            data_central_value = float(values[j]['value'])
            pT_max = input['independent_variables'][0]['values'][j]['high']
            pT_min = input['independent_variables'][0]['values'][j]['low']
            kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'Q2': {'min': Q2_min, 'mid': None, 'max': Q2_max},'pT': {'min': pT_min, 'mid': None, 'max': pT_max}}
            kin.append(kin_value)
            value_delta = 0
            error_value = {}
            for k in 0, 1, 5, 6, 7, 8, 9, 10, 11:
                if 'symerror' in values[j]['errors'][k]:
                    error_value[values[j]['errors'][k]['label']] = pta(values[j]['errors'][k]['symerror'], data_central_value)
                else:
                    se_delta, se_sigma = se(pta(values[j]['errors'][k]['asymerror']['plus'], data_central_value), pta(values[j]['errors'][k]['asymerror']['minus'], data_central_value))
                    value_delta = value_delta + se_delta
                    error_value[values[j]['errors'][k]['label']] = se_sigma
            for k in 2, 3, 4:
                if 'symerror' in values[j]['errors'][k]:
                    error_value[values[j]['errors'][k]['label']+'_1'] = pta(values[j]['errors'][k]['symerror'], data_central_value)/sqrt(2)
                    error_value[values[j]['errors'][k]['label']+'_2'] = pta(values[j]['errors'][k]['symerror'], data_central_value)/sqrt(2)
                else:
                    se_delta, se_sigma = se(pta(values[j]['errors'][k]['asymerror']['plus'], data_central_value)/sqrt(2), pta(values[j]['errors'][k]['asymerror']['minus'], data_central_value)/sqrt(2))
                    value_delta = value_delta + se_delta + se_delta
                    error_value[values[j]['errors'][k]['label']+'_1'] = se_sigma
                    error_value[values[j]['errors'][k]['label']+'_2'] = se_sigma
            data_central_value = data_central_value + value_delta
            data_central.append(data_central_value)
            for k in range(96):
                error_value['ArtUnc_'+str(k+1)] = float(artUncMatr[j + 48][k])
            error_value['stat'] = 0
            error.append(error_value)

    error_definition = {
        'stat':{'description': 'statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'Uncorr':{'description': 'systematic uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'Model_1':{'description': 'MC model uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'Model_2':{'description': 'MC model uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'ModelRW_1':{'description': 'reweighting uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'ModelRW_2':{'description': 'reweighting uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'JES_1':{'description': 'jet energy scale uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'JES_2':{'description': 'jet energy scale uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'RCES':{'description': 'remaining cluster energy scale uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'ElEn':{'description': 'electron energy uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'ElTh':{'description': 'electron theta uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'Lumi':{'description': 'luminosity', 'treatment': 'MULT', 'type': 'CORR'},
        'LArN':{'description': 'lar noise', 'treatment': 'MULT', 'type': 'CORR'},
        'StatMC':{'description': 'MC statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'RadErr':{'description': 'radiative uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    }

    for i in range(96):
        error_definition['ArtUnc_'+str(i+1)] = {'description': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'H1JETS161103421unc'+str(i+1)}

    data_central_yaml = {'data_central': data_central}
    kinematics_yaml = {'bins': kin}
    uncertainties_yaml = {'definitions': error_definition, 'bins': error}

    with open('data.yaml', 'w') as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    with open('uncertainties.yaml', 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)

# dijet_norm data

    for i in tables_norm:
        if i == 33:
            Q2_min = 5.5
            Q2_max = 8
        elif i == 34:
            Q2_min = 8
            Q2_max = 11
        elif i == 35:
            Q2_min = 11
            Q2_max = 16
        elif i == 36:
            Q2_min = 16
            Q2_max = 22
        elif i == 37:
            Q2_min = 22
            Q2_max = 30
        elif i == 38:
            Q2_min = 30
            Q2_max = 42
        elif i == 39:
            Q2_min = 42
            Q2_max = 60
        elif i == 40:
            Q2_min = 60
            Q2_max = 80

        hepdata_tables="rawdata/data"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)

        sqrts = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
        values = input['dependent_variables'][0]['values']

        for j in range(len(values)):
            data_central_value = float(values[j]['value'])
            pT_max = input['independent_variables'][0]['values'][j]['high']
            pT_min = input['independent_variables'][0]['values'][j]['low']
            kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'Q2': {'min': Q2_min, 'mid': None, 'max': Q2_max},'pT': {'min': pT_min, 'mid': None, 'max': pT_max}}
            kin_norm.append(kin_value)
            value_delta = 0
            error_value = {}
            for k in 0, 1, 5, 6, 7, 8, 9, 10, 11:
                if 'symerror' in values[j]['errors'][k]:
                    error_value[values[j]['errors'][k]['label']] = pta(values[j]['errors'][k]['symerror'], data_central_value)
                else:
                    se_delta, se_sigma = se(pta(values[j]['errors'][k]['asymerror']['plus'], data_central_value), pta(values[j]['errors'][k]['asymerror']['minus'], data_central_value))
                    value_delta = value_delta + se_delta
                    error_value[values[j]['errors'][k]['label']] = se_sigma
            for k in 2, 3, 4:
                if 'symerror' in values[j]['errors'][k]:
                    error_value[values[j]['errors'][k]['label']+'_1'] = pta(values[j]['errors'][k]['symerror'], data_central_value)/sqrt(2)
                    error_value[values[j]['errors'][k]['label']+'_2'] = pta(values[j]['errors'][k]['symerror'], data_central_value)/sqrt(2)
                else:
                    se_delta, se_sigma = se(pta(values[j]['errors'][k]['asymerror']['plus'], data_central_value)/sqrt(2), pta(values[j]['errors'][k]['asymerror']['minus'], data_central_value)/sqrt(2))
                    value_delta = value_delta + se_delta + se_delta
                    error_value[values[j]['errors'][k]['label']+'_1'] = se_sigma
                    error_value[values[j]['errors'][k]['label']+'_2'] = se_sigma
            data_central_value = data_central_value + value_delta
            data_central_norm.append(data_central_value)
            for k in range(96):
                error_value['ArtUnc_'+str(k+1)] = float(artUncMatr_norm[j + 48][k])
            error_value['stat'] = 0
            error_norm.append(error_value)

    error_definition_norm = {
        'stat':{'description': 'statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'Uncorr':{'description': 'systematic uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'Model_1':{'description': 'MC model uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'Model_2':{'description': 'MC model uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'ModelRW_1':{'description': 'reweighting uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'ModelRW_2':{'description': 'reweighting uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'JES_1':{'description': 'jet energy scale uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'JES_2':{'description': 'jet energy scale uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'RCES':{'description': 'remaining cluster energy scale uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'ElEn':{'description': 'electron energy uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'ElTh':{'description': 'electron theta uncertainty', 'treatment': 'MULT', 'type': 'CORR'},
        'Lumi':{'description': 'luminosity', 'treatment': 'MULT', 'type': 'CORR'},
        'LArN':{'description': 'lar noise', 'treatment': 'MULT', 'type': 'CORR'},
        'StatMC':{'description': 'MC statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'},
        'RadErr':{'description': 'radiative uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    }
    for i in range(96):
        error_definition_norm['ArtUnc_'+str(i+1)] = {'description': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'H1JETS161103421NORMunc'+str(i+1)}

    data_central_norm_yaml = {'data_central': data_central_norm}
    kinematics_norm_yaml = {'bins': kin_norm}
    uncertainties_norm_yaml = {'definitions': error_definition_norm, 'bins': error_norm}

    with open('data_norm.yaml', 'w') as file:
        yaml.dump(data_central_norm_yaml, file, sort_keys=False)

    with open('kinematics_norm.yaml', 'w') as file:
        yaml.dump(kinematics_norm_yaml, file, sort_keys=False)

    with open('uncertainties_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_norm_yaml, file, sort_keys=False)

processData()
