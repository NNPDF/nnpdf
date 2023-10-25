import artUnc
import yaml
from validphys.commondata_utils import percentage_to_absolute as pta
from validphys.commondata_utils import symmetrize_errors as se
from math import sqrt

def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    tables = metadata['implemented_observables'][0]['tables']
    tables_norm = metadata['implemented_observables'][1]['tables']
    tables_highQ2 = metadata['implemented_observables'][2]['tables']
    tables_highQ2_norm = metadata['implemented_observables'][3]['tables']

    data_central = []
    kin = []
    error = []
    data_central_norm = []
    kin_norm = []
    error_norm = []
    data_central_highQ2 = []
    kin_highQ2 = []
    error_highQ2 = []
    data_central_highQ2_norm = []
    kin_highQ2_norm = []
    error_highQ2_norm = []

    artUncMatr = artUnc.artunc()
    artUncMatr_norm = artUnc.artunc_norm()

# jet data

    for i in tables:
        if i == 1:
            q_sqr_min = 5.5
            q_sqr_max = 8
        elif i == 2:
            q_sqr_min = 8
            q_sqr_max = 11
        elif i == 3:
            q_sqr_min = 11
            q_sqr_max = 16
        elif i == 4:
            q_sqr_min = 16
            q_sqr_max = 22
        elif i == 5:
            q_sqr_min = 22
            q_sqr_max = 30
        elif i == 6:
            q_sqr_min = 30
            q_sqr_max = 42
        elif i == 7:
            q_sqr_min = 42
            q_sqr_max = 60
        elif i == 8:
            q_sqr_min = 60
            q_sqr_max = 80

        hepdata_tables="rawdata/data"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)

        sqrt_s = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
        values = input['dependent_variables'][0]['values']

        for j in range(len(values)):
            data_central_value = float(values[j]['value'])
            pT_max = input['independent_variables'][0]['values'][j]['high']
            pT_min = input['independent_variables'][0]['values'][j]['low']
            kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'q_sqr': {'min': q_sqr_min, 'mid': None, 'max': q_sqr_max},'pT': {'min': pT_min, 'mid': None, 'max': pT_max}}
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
                error_value['ArtUnc_'+str(k+1)] = float(artUncMatr[j][k])
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
        error_definition['ArtUnc_'+str(i+1)] = {'description': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'JET12'}

    data_central_yaml = {'data_central': data_central}
    kinematics_yaml = {'bins': kin}
    uncertainties_yaml = {'definitions': error_definition, 'bins': error}

    with open('data.yaml', 'w') as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    with open('uncertainties.yaml', 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)

# jet_norm data

    for i in tables_norm:
        if i == 25:
            q_sqr_min = 5.5
            q_sqr_max = 8
        elif i == 26:
            q_sqr_min = 8
            q_sqr_max = 11
        elif i == 27:
            q_sqr_min = 11
            q_sqr_max = 16
        elif i == 28:
            q_sqr_min = 16
            q_sqr_max = 22
        elif i == 29:
            q_sqr_min = 22
            q_sqr_max = 30
        elif i == 30:
            q_sqr_min = 30
            q_sqr_max = 42
        elif i == 31:
            q_sqr_min = 42
            q_sqr_max = 60
        elif i == 32:
            q_sqr_min = 60
            q_sqr_max = 80

        hepdata_tables="rawdata/data"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)

        sqrt_s = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
        values = input['dependent_variables'][0]['values']

        for j in range(len(values)):
            data_central_value = float(values[j]['value'])
            pT_max = input['independent_variables'][0]['values'][j]['high']
            pT_min = input['independent_variables'][0]['values'][j]['low']
            kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'q_sqr': {'min': q_sqr_min, 'mid': None, 'max': q_sqr_max},'pT': {'min': pT_min, 'mid': None, 'max': pT_max}}
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
                error_value['ArtUnc_'+str(k+1)] = float(artUncMatr_norm[j][k])
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
        error_definition_norm['ArtUnc_'+str(i+1)] = {'description': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'JET12'}

    data_central_norm_yaml = {'data_central': data_central_norm}
    kinematics_norm_yaml = {'bins': kin_norm}
    uncertainties_norm_yaml = {'definitions': error_definition_norm, 'bins': error_norm}

    with open('data_norm.yaml', 'w') as file:
        yaml.dump(data_central_norm_yaml, file, sort_keys=False)

    with open('kinematics_norm.yaml', 'w') as file:
        yaml.dump(kinematics_norm_yaml, file, sort_keys=False)

    with open('uncertainties_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_norm_yaml, file, sort_keys=False)

# jet_highQ2 data

    hepdata_tables="rawdata/data51.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    
    sqrt_s = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
    pT_min = 5
    pT_max = 7
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = float(values[i]['value'])
        q_sqr_max = input['independent_variables'][0]['values'][i]['high']
        q_sqr_min = input['independent_variables'][0]['values'][i]['low']
        kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'q_sqr': {'min': q_sqr_min, 'mid': None, 'max': q_sqr_max},'pT': {'min': pT_min, 'mid': None, 'max': pT_max}}
        kin_highQ2.append(kin_value)
        value_delta = 0
        error_value = {}
        for k in 0, 3, 4, 5, 6, 7, 8:
                if 'symerror' in values[j]['errors'][k]:
                    error_value[values[j]['errors'][k]['label']] = pta(values[j]['errors'][k]['symerror'], data_central_value)
                else:
                    se_delta, se_sigma = se(pta(values[j]['errors'][k]['asymerror']['plus'], data_central_value), pta(values[j]['errors'][k]['asymerror']['minus'], data_central_value))
                    value_delta = value_delta + se_delta
                    error_value[values[j]['errors'][k]['label']] = se_sigma
        for k in 1, 2:
            if 'symerror' in values[j]['errors'][k]:
                error_value[values[j]['errors'][k]['label']+'_1'] = pta(values[j]['errors'][k]['symerror'], data_central_value)/sqrt(2)
                error_value[values[j]['errors'][k]['label']+'_2'] = pta(values[j]['errors'][k]['symerror'], data_central_value)/sqrt(2)
            else:
                se_delta, se_sigma = se(pta(values[j]['errors'][k]['asymerror']['plus'], data_central_value)/sqrt(2), pta(values[j]['errors'][k]['asymerror']['minus'], data_central_value)/sqrt(2))
                value_delta = value_delta + se_delta + se_delta
                error_value[values[j]['errors'][k]['label']+'_1'] = se_sigma
                error_value[values[j]['errors'][k]['label']+'_2'] = se_sigma
        data_central_value = data_central_value + value_delta
        data_central_highQ2.append(data_central_value)
        error_highQ2.append(error_value)
    
    error_definition_highQ2 = {
        'stat':{'description':'statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR' },
        'Model_1':{'description':'MC model uncertainty', 'treatment': 'ADD', 'type': 'UNCORR' },
        'Model_2':{'description':'MC model uncertainty', 'treatment': 'MULT', 'type': 'CORR' },
        'JES_1':{'description':'jet energy scale uncertainty', 'treatment': 'ADD', 'type': 'UNCORR' },
        'JES_2':{'description':'jet energy scale uncertainty', 'treatment': 'MULT', 'type': 'CORR' },
        'RCES':{'description':'remaining cluster energy scale uncertainty', 'treatment': 'MULT', 'type': 'CORR' },
        '$E_{e^\prime}$':{'description':'electron energy', 'treatment': 'MULT', 'type': 'CORR' },
        '$\theta_{e^\prime}$':{'description':'electron theta', 'treatment': 'MULT', 'type': 'CORR' },
        'ID(e)':{'description': 'electron identification', 'treatment': 'MULT', 'type': 'CORR'},
        'LArNoise':{'description':'lar noice', 'treatment': 'MULT', 'type': 'CORR' },
        'Norm':{'description': 'normalization uncertainty', 'treatment': 'MULT', 'type': 'CORR'}
    }

    data_central_highQ2_yaml = {'data_central': data_central_highQ2}
    kinematics_highQ2_yaml = {'bins': kin_highQ2}
    uncertainties_highQ2_yaml = {'definitions': error_definition_highQ2, 'bins': error_highQ2}


    with open('data_highQ2.yaml', 'w') as file:
        yaml.dump(data_central_highQ2_yaml, file, sort_keys=False)

    with open('kinematics_highQ2.yaml', 'w') as file:
        yaml.dump(kinematics_highQ2_yaml, file, sort_keys=False)

    with open('uncertainties_highQ2.yaml', 'w') as file:
        yaml.dump(uncertainties_highQ2_yaml, file, sort_keys=False)

# jet_highQ2_norm data

    hepdata_tables="rawdata/data52.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    
    sqrt_s = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
    pT_min = 5
    pT_max = 7
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = float(values[i]['value'])
        q_sqr_max = input['independent_variables'][0]['values'][i]['high']
        q_sqr_min = input['independent_variables'][0]['values'][i]['low']
        kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'q_sqr': {'min': q_sqr_min, 'mid': None, 'max': q_sqr_max},'pT': {'min': pT_min, 'mid': None, 'max': pT_max}}
        kin_highQ2_norm.append(kin_value)
        value_delta = 0
        error_value = {}
        for k in 0, 3, 4, 5, 6:
                if 'symerror' in values[j]['errors'][k]:
                    error_value[values[j]['errors'][k]['label']] = pta(values[j]['errors'][k]['symerror'], data_central_value)
                else:
                    se_delta, se_sigma = se(pta(values[j]['errors'][k]['asymerror']['plus'], data_central_value), pta(values[j]['errors'][k]['asymerror']['minus'], data_central_value))
                    value_delta = value_delta + se_delta
                    error_value[values[j]['errors'][k]['label']] = se_sigma
        for k in 1, 2:
            if 'symerror' in values[j]['errors'][k]:
                error_value[values[j]['errors'][k]['label']+'_1'] = pta(values[j]['errors'][k]['symerror'], data_central_value)/sqrt(2)
                error_value[values[j]['errors'][k]['label']+'_2'] = pta(values[j]['errors'][k]['symerror'], data_central_value)/sqrt(2)
            else:
                se_delta, se_sigma = se(pta(values[j]['errors'][k]['asymerror']['plus'], data_central_value)/sqrt(2), pta(values[j]['errors'][k]['asymerror']['minus'], data_central_value)/sqrt(2))
                value_delta = value_delta + se_delta + se_delta
                error_value[values[j]['errors'][k]['label']+'_1'] = se_sigma
                error_value[values[j]['errors'][k]['label']+'_2'] = se_sigma
        data_central_value = data_central_value + value_delta
        data_central_highQ2_norm.append(data_central_value)
        error_highQ2_norm.append(error_value)
    
    error_definition_highQ2_norm = {
        'stat':{'description':'statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR' },
        'Model_1':{'description':'MC model uncertainty', 'treatment': 'ADD', 'type': 'UNCORR' },
        'Model_2':{'description':'MC model uncertainty', 'treatment': 'MULT', 'type': 'CORR' },
        'JES_1':{'description':'jet energy scale uncertainty', 'treatment': 'ADD', 'type': 'UNCORR' },
        'JES_2':{'description':'jet energy scale uncertainty', 'treatment': 'MULT', 'type': 'CORR' },
        'RCES':{'description':'remaining cluster energy scale uncertainty', 'treatment': 'MULT', 'type': 'CORR' },
        '$E_{e^\prime}$':{'description':'electron energy', 'treatment': 'MULT', 'type': 'CORR' },
        '$\theta_{e^\prime}$':{'description':'electron theta', 'treatment': 'MULT', 'type': 'CORR' },
        'LArNoise':{'description':'lar noice', 'treatment': 'MULT', 'type': 'CORR' }
    }

    data_central_highQ2_norm_yaml = {'data_central': data_central_highQ2_norm}
    kinematics_highQ2_norm_yaml = {'bins': kin_highQ2_norm}
    uncertainties_highQ2_norm_yaml = {'definitions': error_definition_highQ2_norm, 'bins': error_highQ2_norm}


    with open('data_highQ2_norm.yaml', 'w') as file:
        yaml.dump(data_central_highQ2_norm_yaml, file, sort_keys=False)

    with open('kinematics_highQ2_norm.yaml', 'w') as file:
        yaml.dump(kinematics_highQ2_norm_yaml, file, sort_keys=False)

    with open('uncertainties_highQ2_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_highQ2_norm_yaml, file, sort_keys=False)

processData()
