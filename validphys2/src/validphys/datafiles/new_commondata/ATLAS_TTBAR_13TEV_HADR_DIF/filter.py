import yaml
from artunc import artunc_mtt, artunc_mtt_norm, artunc_ytt, artunc_ytt_norm, artunc_mtt_ytt, artunc_mtt_ytt_norm

def processData():

    data_central_dSig_dmttBar = []
    kin_dSig_dmttBar = []
    error_dSig_dmttBar = []
    data_central_dSig_dmttBar_norm = []
    kin_dSig_dmttBar_norm = []
    error_dSig_dmttBar_norm = []
    data_central_dSig_dyttBar = []
    kin_dSig_dyttBar = []
    error_dSig_dyttBar = []
    data_central_dSig_dyttBar_norm = []
    kin_dSig_dyttBar_norm = []
    error_dSig_dyttBar_norm = []
    data_central_d2Sig_dmttBar_dyttBar = []
    kin_d2Sig_dmttBar_dyttBar = []
    error_d2Sig_dmttBar_dyttBar = []
    data_central_d2Sig_dmttBar_dyttBar_norm = []
    kin_d2Sig_dmttBar_dyttBar_norm = []
    error_d2Sig_dmttBar_dyttBar_norm = []

# mttbar

    hepdata_tables1="rawdata/Table463.yaml"
    with open(hepdata_tables1, 'r') as file1:
        input1 = yaml.safe_load(file1)
    
    sqrts = 13000.0
    m_t2 = 29756.25
    values = input1['dependent_variables'][0]['values']

    for i in range(len(values)):
        mttbar_min = input1['independent_variables'][0]['values'][i]['low']
        mttbar_max = input1['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(len(values)):
            error_value['ArtUnc_'+str(j+1)] = artunc_mtt[i][j]
        data_central_value = values[i]['value']
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_ttBar': {'min': mttbar_min, 'mid': None, 'max': mttbar_max}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}}
        data_central_dSig_dmttBar.append(data_central_value)
        kin_dSig_dmttBar.append(kin_value)
        error_dSig_dmttBar.append(error_value)

    error_definition_dSig_dmttBar = {}
    for i in range (len(values)):
        error_definition_dSig_dmttBar['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}
    
    data_central_dSig_dmttBar_yaml = {'data_central': data_central_dSig_dmttBar}
    kinematics_dSig_dmttBar_yaml = {'bins': kin_dSig_dmttBar}
    uncertainties_dSig_dmttBar_yaml = {'definitions': error_definition_dSig_dmttBar, 'bins': error_dSig_dmttBar}

    with open('data_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_yaml, file, sort_keys=False)


# mttbar_norm

    hepdata_tables2="rawdata/Table461.yaml"
    with open(hepdata_tables2, 'r') as file2:
        input2 = yaml.safe_load(file2)
    
    sqrts = 13000.0
    m_t2 = 29756.25
    values = input2['dependent_variables'][0]['values']

    for i in range(len(values)):
        mttbar_min = input2['independent_variables'][0]['values'][i]['low']
        mttbar_max = input2['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(len(values)):
            error_value['ArtUnc_'+str(j+1)] = artunc_mtt_norm[i][j]
        data_central_value = values[i]['value']
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_ttBar': {'min': mttbar_min, 'mid': None, 'max': mttbar_max}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}}
        data_central_dSig_dmttBar_norm.append(data_central_value)
        kin_dSig_dmttBar_norm.append(kin_value)
        error_dSig_dmttBar_norm.append(error_value)

    error_definition_dSig_dmttBar_norm = {}
    for i in range (len(values)):
        error_definition_dSig_dmttBar_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}
    
    data_central_dSig_dmttBar_norm_yaml = {'data_central': data_central_dSig_dmttBar_norm}
    kinematics_dSig_dmttBar_norm_yaml = {'bins': kin_dSig_dmttBar_norm}
    uncertainties_dSig_dmttBar_norm_yaml = {'definitions': error_definition_dSig_dmttBar_norm, 'bins': error_dSig_dmttBar_norm}

    with open('data_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_norm_yaml, file, sort_keys=False)

# yttbar

    hepdata_tables3="rawdata/Table475.yaml"
    with open(hepdata_tables3, 'r') as file3:
        input3 = yaml.safe_load(file3)
    
    sqrts = 13000.0
    m_t2 = 29756.25
    values = input3['dependent_variables'][0]['values']

    for i in range(len(values)):
        yttbar_min = input3['independent_variables'][0]['values'][i]['low']
        yttbar_max = input3['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(len(values)):
            error_value['ArtUnc_'+str(j+1)] = artunc_ytt[i][j]
        data_central_value = values[i]['value']
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'y_ttBar': {'min': yttbar_min, 'mid': None, 'max': yttbar_max}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}}
        data_central_dSig_dyttBar.append(data_central_value)
        kin_dSig_dyttBar.append(kin_value)
        error_dSig_dyttBar.append(error_value)

    error_definition_dSig_dyttBar = {}
    for i in range (len(values)):
        error_definition_dSig_dyttBar['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}
    
    data_central_dSig_dyttBar_yaml = {'data_central': data_central_dSig_dyttBar}
    kinematics_dSig_dyttBar_yaml = {'bins': kin_dSig_dyttBar}
    uncertainties_dSig_dyttBar_yaml = {'definitions': error_definition_dSig_dyttBar, 'bins': error_dSig_dyttBar}

    with open('data_dSig_dyttBar.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyttBar.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_yaml, file, sort_keys=False)

# yttbar_norm

    hepdata_tables4="rawdata/Table473.yaml"
    with open(hepdata_tables4, 'r') as file4:
        input4 = yaml.safe_load(file4)
    
    sqrts = 13000.0
    m_t2 = 29756.25
    values = input4['dependent_variables'][0]['values']

    for i in range(len(values)):
        yttbar_min = input4['independent_variables'][0]['values'][i]['low']
        yttbar_max = input4['independent_variables'][0]['values'][i]['high']
        error_value = {}
        for j in range(len(values)):
            error_value['ArtUnc_'+str(j+1)] = artunc_ytt_norm[i][j]
        data_central_value = values[i]['value']
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'y_ttBar': {'min': yttbar_min, 'mid': None, 'max': yttbar_max}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}}
        data_central_dSig_dyttBar_norm.append(data_central_value)
        kin_dSig_dyttBar_norm.append(kin_value)
        error_dSig_dyttBar_norm.append(error_value)

    error_definition_dSig_dyttBar_norm = {}
    for i in range (len(values)):
        error_definition_dSig_dyttBar_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}
    
    data_central_dSig_dyttBar_norm_yaml = {'data_central': data_central_dSig_dyttBar_norm}
    kinematics_dSig_dyttBar_norm_yaml = {'bins': kin_dSig_dyttBar_norm}
    uncertainties_dSig_dyttBar_norm_yaml = {'definitions': error_definition_dSig_dyttBar_norm, 'bins': error_dSig_dyttBar_norm}

    with open('data_dSig_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_norm_yaml, file, sort_keys=False)

# mttbar_yttbar

    hepdata_tables5="rawdata/Table498.yaml"
    with open(hepdata_tables5, 'r') as file5:
        input5 = yaml.safe_load(file5)

    sqrts = 13000.0
    m_t2 = 29756.25
    values = input5['dependent_variables'][0]['values']

    for i in range(len(values)):
        yttbar_min = input5['independent_variables'][0]['values'][i]['low']
        yttbar_max = input5['independent_variables'][0]['values'][i]['high']
        mttbar_min = 0.0
        mttbar_max = 700.0
        error_value = {}
        for j in range(11):
            error_value['ArtUnc_'+str(j+1)] = artunc_mtt_ytt[i][j]
        data_central_value = values[i]['value']
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_ttBar': {'min': mttbar_min, 'mid': None, 'max': mttbar_max}, 'y_ttBar': {'min': yttbar_min, 'mid': None, 'max': yttbar_max}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}}
        data_central_d2Sig_dmttBar_dyttBar.append(data_central_value)
        kin_d2Sig_dmttBar_dyttBar.append(kin_value)
        error_d2Sig_dmttBar_dyttBar.append(error_value)

    hepdata_tables6="rawdata/Table499.yaml"
    with open(hepdata_tables6, 'r') as file6:
        input6 = yaml.safe_load(file6)

    sqrts = 13000.0
    m_t2 = 29756.25
    values = input6['dependent_variables'][0]['values']

    for i in range(len(values)):
        yttbar_min = input6['independent_variables'][0]['values'][i]['low']
        yttbar_max = input6['independent_variables'][0]['values'][i]['high']
        mttbar_min = 700.0
        mttbar_max = 970.0
        error_value = {}
        for j in range(11):
            error_value['ArtUnc_'+str(j+1)] = artunc_mtt_ytt[i+4][j]
        data_central_value = values[i]['value']
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_ttBar': {'min': mttbar_min, 'mid': None, 'max': mttbar_max}, 'y_ttBar': {'min': yttbar_min, 'mid': None, 'max': yttbar_max}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}}
        data_central_d2Sig_dmttBar_dyttBar.append(data_central_value)
        kin_d2Sig_dmttBar_dyttBar.append(kin_value)
        error_d2Sig_dmttBar_dyttBar.append(error_value)

    hepdata_tables7="rawdata/Table500.yaml"
    with open(hepdata_tables7, 'r') as file7:
        input7 = yaml.safe_load(file7)

    sqrts = 13000.0
    m_t2 = 29756.25
    values = input7['dependent_variables'][0]['values']

    for i in range(len(values)):
        yttbar_min = input7['independent_variables'][0]['values'][i]['low']
        yttbar_max = input7['independent_variables'][0]['values'][i]['high']
        mttbar_min = 970.0
        mttbar_max = 3000.0
        error_value = {}
        for j in range(11):
            error_value['ArtUnc_'+str(j+1)] = artunc_mtt_ytt[i+8][j]
        data_central_value = values[i]['value']
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_ttBar': {'min': mttbar_min, 'mid': None, 'max': mttbar_max}, 'y_ttBar': {'min': yttbar_min, 'mid': None, 'max': yttbar_max}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}}
        data_central_d2Sig_dmttBar_dyttBar.append(data_central_value)
        kin_d2Sig_dmttBar_dyttBar.append(kin_value)
        error_d2Sig_dmttBar_dyttBar.append(error_value)

    error_definition_d2Sig_dmttBar_dyttBar = {}
    for i in range (11):
        error_definition_d2Sig_dmttBar_dyttBar['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_d2Sig_dmttBar_dyttBar_yaml = {'data_central': data_central_d2Sig_dmttBar_dyttBar}
    kinematics_d2Sig_dmttBar_dyttBar_yaml = {'bins': kin_d2Sig_dmttBar_dyttBar}
    uncertainties_d2Sig_dmttBar_dyttBar_yaml = {'definitions': error_definition_d2Sig_dmttBar_dyttBar, 'bins': error_d2Sig_dmttBar_dyttBar}

    with open('data_d2Sig_dmttBar_dyttBar.yaml', 'w') as file:
        yaml.dump(data_central_d2Sig_dmttBar_dyttBar_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dmttBar_dyttBar.yaml', 'w') as file:
        yaml.dump(kinematics_d2Sig_dmttBar_dyttBar_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dmttBar_dyttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dmttBar_dyttBar_yaml, file, sort_keys=False)

# mttbar_yttbar_norm
    
    hepdata_tables8="rawdata/Table489.yaml"
    with open(hepdata_tables8, 'r') as file8:
        input8 = yaml.safe_load(file8)

    sqrts = 13000.0
    m_t2 = 29756.25
    values = input8['dependent_variables'][0]['values']

    for i in range(len(values)):
        yttbar_min = input8['independent_variables'][0]['values'][i]['low']
        yttbar_max = input8['independent_variables'][0]['values'][i]['high']
        mttbar_min = 0.0
        mttbar_max = 700.0
        error_value = {}
        for j in range(11):
            error_value['ArtUnc_'+str(j+1)] = artunc_mtt_ytt_norm[i][j]
        data_central_value = values[i]['value']
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_ttBar': {'min': mttbar_min, 'mid': None, 'max': mttbar_max}, 'y_ttBar': {'min': yttbar_min, 'mid': None, 'max': yttbar_max}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}}
        data_central_d2Sig_dmttBar_dyttBar_norm.append(data_central_value)
        kin_d2Sig_dmttBar_dyttBar_norm.append(kin_value)
        error_d2Sig_dmttBar_dyttBar_norm.append(error_value)

    hepdata_tables9="rawdata/Table490.yaml"
    with open(hepdata_tables9, 'r') as file9:
        input9 = yaml.safe_load(file9)

    sqrts = 13000.0
    m_t2 = 29756.25
    values = input9['dependent_variables'][0]['values']

    for i in range(len(values)):
        yttbar_min = input9['independent_variables'][0]['values'][i]['low']
        yttbar_max = input9['independent_variables'][0]['values'][i]['high']
        mttbar_min = 700.0
        mttbar_max = 970.0
        error_value = {}
        for j in range(11):
            error_value['ArtUnc_'+str(j+1)] = artunc_mtt_ytt_norm[i+4][j]
        data_central_value = values[i]['value']
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_ttBar': {'min': mttbar_min, 'mid': None, 'max': mttbar_max}, 'y_ttBar': {'min': yttbar_min, 'mid': None, 'max': yttbar_max}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}}
        data_central_d2Sig_dmttBar_dyttBar_norm.append(data_central_value)
        kin_d2Sig_dmttBar_dyttBar_norm.append(kin_value)
        error_d2Sig_dmttBar_dyttBar_norm.append(error_value)

    hepdata_tables10="rawdata/Table491.yaml"
    with open(hepdata_tables10, 'r') as file10:
        input10 = yaml.safe_load(file10)

    sqrts = 13000.0
    m_t2 = 29756.25
    values = input10['dependent_variables'][0]['values']

    for i in range(len(values)):
        yttbar_min = input10['independent_variables'][0]['values'][i]['low']
        yttbar_max = input10['independent_variables'][0]['values'][i]['high']
        mttbar_min = 970.0
        mttbar_max = 3000.0
        error_value = {}
        for j in range(11):
            error_value['ArtUnc_'+str(j+1)] = artunc_mtt_ytt_norm[i+8][j]
        data_central_value = values[i]['value']
        kin_value = {'sqrts': {'min': None, 'mid': sqrts, 'max': None}, 'm_ttBar': {'min': mttbar_min, 'mid': None, 'max': mttbar_max}, 'y_ttBar': {'min': yttbar_min, 'mid': None, 'max': yttbar_max}, 'm_t2': {'min': None, 'mid': m_t2, 'max': None}}
        data_central_d2Sig_dmttBar_dyttBar_norm.append(data_central_value)
        kin_d2Sig_dmttBar_dyttBar_norm.append(kin_value)
        error_d2Sig_dmttBar_dyttBar_norm.append(error_value)

    error_definition_d2Sig_dmttBar_dyttBar_norm = {}
    for i in range (11):
        error_definition_d2Sig_dmttBar_dyttBar_norm['ArtUnc_'+str(i+1)] = {'definition': 'artificial uncertainty '+str(i+1), 'treatment': 'ADD', 'type': 'CORR'}

    data_central_d2Sig_dmttBar_dyttBar_norm_yaml = {'data_central': data_central_d2Sig_dmttBar_dyttBar_norm}
    kinematics_d2Sig_dmttBar_dyttBar_norm_yaml = {'bins': kin_d2Sig_dmttBar_dyttBar_norm}
    uncertainties_d2Sig_dmttBar_dyttBar_norm_yaml = {'definitions': error_definition_d2Sig_dmttBar_dyttBar_norm, 'bins': error_d2Sig_dmttBar_dyttBar_norm}

    with open('data_d2Sig_dmttBar_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(data_central_d2Sig_dmttBar_dyttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dmttBar_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(kinematics_d2Sig_dmttBar_dyttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dmttBar_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dmttBar_dyttBar_norm_yaml, file, sort_keys=False)

processData()
