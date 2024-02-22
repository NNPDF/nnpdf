import artunc
import yaml

from validphys.commondata_utils import percentage_to_absolute as pta
from validphys.commondata_utils import symmetrize_errors as se

def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    data_central_dSig_dmttBar = []
    kin_dSig_dmttBar = []
    error_dSig_dmttBar = []
    data_central_dSig_dmttBar_norm = []
    kin_dSig_dmttBar_norm = []
    error_dSig_dmttBar_norm = []
    data_central_dSig_dpTt = []
    kin_dSig_dpTt = []
    error_dSig_dpTt = []
    data_central_dSig_dpTt_norm = []
    kin_dSig_dpTt_norm = []
    error_dSig_dpTt_norm = []
    data_central_dSig_dyt = []
    kin_dSig_dyt = []
    error_dSig_dyt = []
    data_central_dSig_dyt_norm = []
    kin_dSig_dyt_norm = []
    error_dSig_dyt_norm = []
    data_central_dSig_dyttBar = []
    kin_dSig_dyttBar = []
    error_dSig_dyttBar = []
    data_central_dSig_dyttBar_norm = []
    kin_dSig_dyttBar_norm = []
    error_dSig_dyttBar_norm = []

    artUnc = artunc.artunc()
    artUnc_norm = artunc.artunc_norm()
# dSig_dmttBar

    hepdata_tables='rawdata/Table_23.yaml'
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    
    sqrts = 8000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = values[i]['value']
        m_ttbar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttbar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(25):
            error_value['ArtUnc_'+str(j+1)] = artUnc[i][j]
        value_delta = 0
        for j in range(1, len(input['dependent_variables'][1]['values'][i]['errors'])-1):
            if 'symerror' in input['dependent_variables'][1]['values'][i]['errors'][j]:
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = pta(input['dependent_variables'][1]['values'][i]['errors'][j]['symerror'], data_central_value)
            else:
                se_delta, se_sigma = se(pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['plus'], data_central_value), pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['minus'], data_central_value))
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = se_sigma
                value_delta = value_delta + se_delta
        error_value['lumi'] = pta(values[i]['errors'][2]['symerror'], data_central_value)
        data_central_value = data_central_value + value_delta
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'm_ttBar': {'min': m_ttbar_min, 'mid': None, 'max': m_ttbar_max}}
        data_central_dSig_dmttBar.append(data_central_value)
        kin_dSig_dmttBar.append(kin_value)
        error_dSig_dmttBar.append(error_value)

    error_definition_dSig_dmttBar = {}
    for i in range(25):
        error_definition_dSig_dmttBar['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'ATLAS8TEVTTB151104716'}
    for i in range(1, len(input['dependent_variables'][1]['values'][0]['errors'])-1):
        error_definition_dSig_dmttBar[input['dependent_variables'][1]['values'][0]['errors'][i]['label'].replace(" ", "")] = {'definition': '', 'treatment': 'MULT', 'type': 'CORR'}
    error_definition_dSig_dmttBar['lumi'] =  {'definition': 'luminosity uncertainty', 'treatment': 'MULT', 'type': 'ATLASLUMI8'}

    data_central_dSig_dmttBar_yaml = {'data_central': data_central_dSig_dmttBar}
    kinematics_dSig_dmttBar_yaml = {'bins': kin_dSig_dmttBar}
    uncertainties_dSig_dmttBar_yaml = {'definitions': error_definition_dSig_dmttBar, 'bins': error_dSig_dmttBar}

    with open('data_dSig_dmttBar.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_yaml, file, sort_keys=False)
  
# dSig_dmttBar_norm

    hepdata_tables='rawdata/Table_24.yaml'
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    
    sqrts = 8000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = values[i]['value']
        m_ttbar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttbar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(25):
            error_value['ArtUnc_'+str(j+1)] = artUnc_norm[i][j]
        value_delta = 0
        for j in range(1, len(input['dependent_variables'][1]['values'][i]['errors'])):
            if 'symerror' in input['dependent_variables'][1]['values'][i]['errors'][j]:
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = pta(input['dependent_variables'][1]['values'][i]['errors'][j]['symerror'], data_central_value)
            else:
                se_delta, se_sigma = se(pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['plus'], data_central_value), pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['minus'], data_central_value))
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = se_sigma
                value_delta = value_delta + se_delta
        data_central_value = data_central_value + value_delta
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'm_ttBar': {'min': m_ttbar_min, 'mid': None, 'max': m_ttbar_max}}
        data_central_dSig_dmttBar_norm.append(data_central_value)
        kin_dSig_dmttBar_norm.append(kin_value)
        error_dSig_dmttBar_norm.append(error_value)

    error_definition_dSig_dmttBar_norm = {}
    for i in range(25):
        error_definition_dSig_dmttBar_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'ATLAS8TEVTTB151104716'}
    for i in range(1, len(input['dependent_variables'][1]['values'][0]['errors'])):
        error_definition_dSig_dmttBar_norm[input['dependent_variables'][1]['values'][0]['errors'][i]['label'].replace(" ", "")] = {'definition': '', 'treatment': 'MULT', 'type': 'CORR'}

    data_central_dSig_dmttBar_norm_yaml = {'data_central': data_central_dSig_dmttBar_norm}
    kinematics_dSig_dmttBar_norm_yaml = {'bins': kin_dSig_dmttBar_norm}
    uncertainties_dSig_dmttBar_norm_yaml = {'definitions': error_definition_dSig_dmttBar_norm, 'bins': error_dSig_dmttBar_norm}

    with open('data_dSig_dmttBar_norm.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar_norm.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_norm_yaml, file, sort_keys=False)

# dSig_dpTt

    hepdata_tables='rawdata/Table_29.yaml'
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    
    sqrts = 8000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = values[i]['value']
        pT_t_min = input['independent_variables'][0]['values'][i]['low']
        pT_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(25):
            error_value['ArtUnc_'+str(j+1)] = artUnc[i+7][j]
        value_delta = 0
        for j in range(1, len(input['dependent_variables'][1]['values'][i]['errors'])-1):
            if 'symerror' in input['dependent_variables'][1]['values'][i]['errors'][j]:
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = pta(input['dependent_variables'][1]['values'][i]['errors'][j]['symerror'], data_central_value)
            else:
                se_delta, se_sigma = se(pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['plus'], data_central_value), pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['minus'], data_central_value))
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = se_sigma
                value_delta = value_delta + se_delta
        error_value['lumi'] = pta(values[i]['errors'][2]['symerror'], data_central_value)
        data_central_value = data_central_value + value_delta
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max}}
        data_central_dSig_dpTt.append(data_central_value)
        kin_dSig_dpTt.append(kin_value)
        error_dSig_dpTt.append(error_value)

    error_definition_dSig_dpTt = {}
    for i in range(25):
        error_definition_dSig_dpTt['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'ATLAS8TEVTTB151104716'}
    for i in range(1, len(input['dependent_variables'][1]['values'][0]['errors'])-1):
        error_definition_dSig_dpTt[input['dependent_variables'][1]['values'][0]['errors'][i]['label'].replace(" ", "")] = {'definition': '', 'treatment': 'MULT', 'type': 'CORR'}
    error_definition_dSig_dpTt['lumi'] =  {'definition': 'luminosity uncertainty', 'treatment': 'MULT', 'type': 'ATLASLUMI8'}

    data_central_dSig_dpTt_yaml = {'data_central': data_central_dSig_dpTt}
    kinematics_dSig_dpTt_yaml = {'bins': kin_dSig_dpTt}
    uncertainties_dSig_dpTt_yaml = {'definitions': error_definition_dSig_dpTt, 'bins': error_dSig_dpTt}

    with open('data_dSig_dpTt.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dpTt_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dpTt.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dpTt_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dpTt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_yaml, file, sort_keys=False)

# dSig_dpTt_norm

    hepdata_tables='rawdata/Table_30.yaml'
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    
    sqrts = 8000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = values[i]['value']
        pT_t_min = input['independent_variables'][0]['values'][i]['low']
        pT_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(25):
            error_value['ArtUnc_'+str(j+1)] = artUnc_norm[i+7][j]
        value_delta = 0
        for j in range(1, len(input['dependent_variables'][1]['values'][i]['errors'])):
            if 'symerror' in input['dependent_variables'][1]['values'][i]['errors'][j]:
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = pta(input['dependent_variables'][1]['values'][i]['errors'][j]['symerror'], data_central_value)
            else:
                se_delta, se_sigma = se(pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['plus'], data_central_value), pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['minus'], data_central_value))
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = se_sigma
                value_delta = value_delta + se_delta
        data_central_value = data_central_value + value_delta
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max}}
        data_central_dSig_dpTt_norm.append(data_central_value)
        kin_dSig_dpTt_norm.append(kin_value)
        error_dSig_dpTt_norm.append(error_value)

    error_definition_dSig_dpTt_norm = {}
    for i in range(25):
        error_definition_dSig_dpTt_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'ATLAS8TEVTTB151104716'}
    for i in range(1, len(input['dependent_variables'][1]['values'][0]['errors'])):
        error_definition_dSig_dpTt_norm[input['dependent_variables'][1]['values'][0]['errors'][i]['label'].replace(" ", "")] = {'definition': '', 'treatment': 'MULT', 'type': 'CORR'}

    data_central_dSig_dpTt_norm_yaml = {'data_central': data_central_dSig_dpTt_norm}
    kinematics_dSig_dpTt_norm_yaml = {'bins': kin_dSig_dpTt_norm}
    uncertainties_dSig_dpTt_norm_yaml = {'definitions': error_definition_dSig_dpTt_norm, 'bins': error_dSig_dpTt_norm}

    with open('data_dSig_dpTt_norm.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dpTt_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dpTt_norm.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dpTt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dpTt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_norm_yaml, file, sort_keys=False)

# dSig_dyt

    hepdata_tables='rawdata/Table_31.yaml'
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    
    sqrts = 8000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = values[i]['value']
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(25):
            error_value['ArtUnc_'+str(j+1)] = artUnc[i+15][j]
        value_delta = 0
        for j in range(1, len(input['dependent_variables'][1]['values'][i]['errors'])-1):
            if 'symerror' in input['dependent_variables'][1]['values'][i]['errors'][j]:
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = pta(input['dependent_variables'][1]['values'][i]['errors'][j]['symerror'], data_central_value)
            else:
                se_delta, se_sigma = se(pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['plus'], data_central_value), pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['minus'], data_central_value))
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = se_sigma
                value_delta = value_delta + se_delta
        error_value['lumi'] = pta(values[i]['errors'][2]['symerror'], data_central_value)
        data_central_value = data_central_value + value_delta
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max}}
        data_central_dSig_dyt.append(data_central_value)
        kin_dSig_dyt.append(kin_value)
        error_dSig_dyt.append(error_value)

    error_definition_dSig_dyt = {}
    for i in range(25):
        error_definition_dSig_dyt['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'ATLAS8TEVTTB151104716'}
    for i in range(1, len(input['dependent_variables'][1]['values'][0]['errors'])-1):
        error_definition_dSig_dyt[input['dependent_variables'][1]['values'][0]['errors'][i]['label'].replace(" ", "")] = {'definition': '', 'treatment': 'MULT', 'type': 'CORR'}
    error_definition_dSig_dyt['lumi'] =  {'definition': 'luminosity uncertainty', 'treatment': 'MULT', 'type': 'ATLASLUMI8'}

    data_central_dSig_dyt_yaml = {'data_central': data_central_dSig_dyt}
    kinematics_dSig_dyt_yaml = {'bins': kin_dSig_dyt}
    uncertainties_dSig_dyt_yaml = {'definitions': error_definition_dSig_dyt, 'bins': error_dSig_dyt}

    with open('data_dSig_dyt.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dyt_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyt.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dyt_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_yaml, file, sort_keys=False)

# dSig_dyt_norm

    hepdata_tables='rawdata/Table_32.yaml'
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    
    sqrts = 8000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = values[i]['value']
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(25):
            error_value['ArtUnc_'+str(j+1)] = artUnc_norm[i+15][j]
        value_delta = 0
        for j in range(1, len(input['dependent_variables'][1]['values'][i]['errors'])):
            if 'symerror' in input['dependent_variables'][1]['values'][i]['errors'][j]:
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = pta(input['dependent_variables'][1]['values'][i]['errors'][j]['symerror'], data_central_value)
            else:
                se_delta, se_sigma = se(pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['plus'], data_central_value), pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['minus'], data_central_value))
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = se_sigma
                value_delta = value_delta + se_delta
        data_central_value = data_central_value + value_delta
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max}}
        data_central_dSig_dyt_norm.append(data_central_value)
        kin_dSig_dyt_norm.append(kin_value)
        error_dSig_dyt_norm.append(error_value)

    error_definition_dSig_dyt_norm = {}
    for i in range(25):
        error_definition_dSig_dyt_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'ATLAS8TEVTTB151104716'}
    for i in range(1, len(input['dependent_variables'][1]['values'][0]['errors'])):
        error_definition_dSig_dyt_norm[input['dependent_variables'][1]['values'][0]['errors'][i]['label'].replace(" ", "")] = {'definition': '', 'treatment': 'MULT', 'type': 'CORR'}

    data_central_dSig_dyt_norm_yaml = {'data_central': data_central_dSig_dyt_norm}
    kinematics_dSig_dyt_norm_yaml = {'bins': kin_dSig_dyt_norm}
    uncertainties_dSig_dyt_norm_yaml = {'definitions': error_definition_dSig_dyt_norm, 'bins': error_dSig_dyt_norm}

    with open('data_dSig_dyt_norm.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dyt_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyt_norm.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dyt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_norm_yaml, file, sort_keys=False)

# dSig_dyttBar

    hepdata_tables='rawdata/Table_27.yaml'
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    
    sqrts = 8000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = values[i]['value']
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(25):
            error_value['ArtUnc_'+str(j+1)] = artUnc[i+20][j]
        value_delta = 0
        for j in range(1, len(input['dependent_variables'][1]['values'][i]['errors'])-1):
            if 'symerror' in input['dependent_variables'][1]['values'][i]['errors'][j]:
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = pta(input['dependent_variables'][1]['values'][i]['errors'][j]['symerror'], data_central_value)
            else:
                se_delta, se_sigma = se(pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['plus'], data_central_value), pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['minus'], data_central_value))
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = se_sigma
                value_delta = value_delta + se_delta
        error_value['lumi'] = pta(values[i]['errors'][2]['symerror'], data_central_value)
        data_central_value = data_central_value + value_delta
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max}}
        data_central_dSig_dyttBar.append(data_central_value)
        kin_dSig_dyttBar.append(kin_value)
        error_dSig_dyttBar.append(error_value)

    error_definition_dSig_dyttBar = {}
    for i in range(25):
        error_definition_dSig_dyttBar['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'ATLAS8TEVTTB151104716'}
    for i in range(1, len(input['dependent_variables'][1]['values'][0]['errors'])-1):
        error_definition_dSig_dyttBar[input['dependent_variables'][1]['values'][0]['errors'][i]['label'].replace(" ", "")] = {'definition': '', 'treatment': 'MULT', 'type': 'CORR'}
    error_definition_dSig_dyttBar['lumi'] =  {'definition': 'luminosity uncertainty', 'treatment': 'MULT', 'type': 'ATLASLUMI8'}

    data_central_dSig_dyttBar_yaml = {'data_central': data_central_dSig_dyttBar}
    kinematics_dSig_dyttBar_yaml = {'bins': kin_dSig_dyttBar}
    uncertainties_dSig_dyttBar_yaml = {'definitions': error_definition_dSig_dyttBar, 'bins': error_dSig_dyttBar}

    with open('data_dSig_dyttBar.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyttBar.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_yaml, file, sort_keys=False)

# dSig_dyttBar_norm

    hepdata_tables='rawdata/Table_28.yaml'
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)
    
    sqrts = 8000.0
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        data_central_value = values[i]['value']
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(25):
            error_value['ArtUnc_'+str(j+1)] = artUnc_norm[i+20][j]
        value_delta = 0
        for j in range(1, len(input['dependent_variables'][1]['values'][i]['errors'])):
            if 'symerror' in input['dependent_variables'][1]['values'][i]['errors'][j]:
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = pta(input['dependent_variables'][1]['values'][i]['errors'][j]['symerror'], data_central_value)
            else:
                se_delta, se_sigma = se(pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['plus'], data_central_value), pta(input['dependent_variables'][1]['values'][i]['errors'][j]['asymerror']['minus'], data_central_value))
                error_value[input['dependent_variables'][1]['values'][i]['errors'][j]['label'].replace(" ", "")] = se_sigma
                value_delta = value_delta + se_delta
        data_central_value = data_central_value + value_delta
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}, 'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max}}
        data_central_dSig_dyttBar_norm.append(data_central_value)
        kin_dSig_dyttBar_norm.append(kin_value)
        error_dSig_dyttBar_norm.append(error_value)

    error_definition_dSig_dyttBar_norm = {}
    for i in range(25):
        error_definition_dSig_dyttBar_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'ATLAS8TEVTTB151104716'}
    for i in range(1, len(input['dependent_variables'][1]['values'][0]['errors'])):
        error_definition_dSig_dyttBar_norm[input['dependent_variables'][1]['values'][0]['errors'][i]['label'].replace(" ", "")] = {'definition': '', 'treatment': 'MULT', 'type': 'CORR'}

    data_central_dSig_dyttBar_norm_yaml = {'data_central': data_central_dSig_dyttBar_norm}
    kinematics_dSig_dyttBar_norm_yaml = {'bins': kin_dSig_dyttBar_norm}
    uncertainties_dSig_dyttBar_norm_yaml = {'definitions': error_definition_dSig_dyttBar_norm, 'bins': error_dSig_dyttBar_norm}

    with open('data_dSig_dyttBar_norm.yaml', 'w') as file:
         yaml.dump(data_central_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyttBar_norm.yaml', 'w') as file:
         yaml.dump(kinematics_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_norm_yaml, file, sort_keys=False)

processData()
