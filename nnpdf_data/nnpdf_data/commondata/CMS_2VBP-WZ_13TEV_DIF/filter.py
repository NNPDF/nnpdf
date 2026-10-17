import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def magic():
    with open('metadata.yaml', 'r') as f:
        metadata = yaml.safe_load(f)

    ndat_ptz_comb = metadata['implemented_observables'][0]['ndata']
    ndat_ptz_eee = metadata['implemented_observables'][1]['ndata']
    ndat_ptz_eem = metadata['implemented_observables'][2]['ndata']
    ndat_ptz_emm = metadata['implemented_observables'][3]['ndata']
    ndat_ptz_mmm = metadata['implemented_observables'][4]['ndata']

    data_central_ptz_comb = []
    kin_ptz_comb = []
    error_ptz_comb = []

    data_central_ptz_eee = []
    kin_ptz_eee = []
    error_ptz_eee = []

    data_central_ptz_eem = []
    kin_ptz_eem = []
    error_ptz_eem = []

    data_central_ptz_emm = []
    kin_ptz_emm = []
    error_ptz_emm = []

    data_central_ptz_mmm = []
    kin_ptz_mmm = []
    error_ptz_mmm = []

    ptZ_bins_edges = [0, 10, 20, 30, 50, 70, 90, 110, 130, 160, 200, 300]

    # ptz_comb data

    hepdata_tab1 = "rawdata/table_9.yaml"
    with open(hepdata_tab1, 'r') as f:
        input1 = yaml.safe_load(f)

    for i in range(ndat_ptz_comb):
        data_central_value = input1['dependent_variables'][0]['values'][i]['value']
        data_central_ptz_comb.append(float(data_central_value))
        ptZ_min = ptZ_bins_edges[i]
        ptZ_max = ptZ_bins_edges[i + 1]
        kin_value = {'ptZ': {'min': ptZ_min, 'mid': None, 'max': ptZ_max}}
        kin_ptz_comb.append(kin_value)
        error_value = {}
        error_value['stat'] = input1['dependent_variables'][1]['values'][i]['errors'][0]['symerror']
        error_value['bgr'] = input1['dependent_variables'][2]['values'][i]['errors'][0]['symerror']
        error_value['sys'] = input1['dependent_variables'][3]['values'][i]['errors'][0]['symerror']
        error_ptz_comb.append(error_value)

    error_definition_ptz_comb = {}
    error_definition_ptz_comb['stat'] = {
        'definition': 'statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    error_definition_ptz_comb['bgr'] = {
        'definition': 'background uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    }
    error_definition_ptz_comb['sys'] = {
        'definition': 'systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    }

    data_central_ptz_comb_yaml = {'data_cental': data_central_ptz_comb}
    kinematics_ptz_comb_yaml = {'bins': kin_ptz_comb}
    uncertainties_ptz_comb_yaml = {'definitions': error_definition_ptz_comb, 'bins': error_ptz_comb}

    with open('data_central_ptz_comb.yaml', 'w') as f:
        yaml.dump(data_central_ptz_comb_yaml, f, sort_keys=False)

    with open('kinematics_ptz_comb.yaml', 'w') as f:
        yaml.dump(kinematics_ptz_comb_yaml, f, sort_keys=False)

    with open('uncertainties_ptz_comb.yaml', 'w') as f:
        yaml.dump(uncertainties_ptz_comb_yaml, f, sort_keys=False)

    # ptz_eee data

    hepdata_tab2 = "rawdata/table_7.yaml"
    with open(hepdata_tab2, 'r') as f:
        input2 = yaml.safe_load(f)

    for i in range(ndat_ptz_eee):
        data_central_value = input2['dependent_variables'][0]['values'][i + 1]['value']
        data_central_ptz_eee.append(float(data_central_value))
        ptZ_min = ptZ_bins_edges[i]
        ptZ_max = ptZ_bins_edges[i + 1]
        kin_value = {'ptZ': {'min': ptZ_min, 'mid': None, 'max': ptZ_max}}
        kin_ptz_eee.append(kin_value)
        error_value = {}
        error_value['stat'] = input2['dependent_variables'][1]['values'][i + 1]['errors'][0][
            'symerror'
        ]
        error_value['bgr'] = input2['dependent_variables'][2]['values'][i + 1]['errors'][0][
            'symerror'
        ]
        error_value['sys'] = input2['dependent_variables'][3]['values'][i + 1]['errors'][0][
            'symerror'
        ]
        error_ptz_eee.append(error_value)

    error_definition_ptz_eee = {}
    error_definition_ptz_eee['stat'] = {
        'definition': 'statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    error_definition_ptz_eee['bgr'] = {
        'definition': 'background uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    }
    error_definition_ptz_eee['sys'] = {
        'definition': 'systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    }

    data_central_ptz_eee_yaml = {'data_cental': data_central_ptz_eee}
    kinematics_ptz_eee_yaml = {'bins': kin_ptz_eee}
    uncertainties_ptz_eee_yaml = {'definitions': error_definition_ptz_eee, 'bins': error_ptz_eee}

    with open('data_central_ptz_eee.yaml', 'w') as f:
        yaml.dump(data_central_ptz_eee_yaml, f, sort_keys=False)

    with open('kinematics_ptz_eee.yaml', 'w') as f:
        yaml.dump(kinematics_ptz_eee_yaml, f, sort_keys=False)

    with open('uncertainties_ptz_eee.yaml', 'w') as f:
        yaml.dump(uncertainties_ptz_eee_yaml, f, sort_keys=False)

    # ptz_eem data

    for i in range(ndat_ptz_eem):
        data_central_value = input2['dependent_variables'][0]['values'][i + 2 + ndat_ptz_eee][
            'value'
        ]
        data_central_ptz_eem.append(float(data_central_value))
        ptZ_min = ptZ_bins_edges[i]
        ptZ_max = ptZ_bins_edges[i + 1]
        kin_value = {'ptZ': {'min': ptZ_min, 'mid': None, 'max': ptZ_max}}
        kin_ptz_eem.append(kin_value)
        error_value = {}
        error_value['stat'] = input2['dependent_variables'][1]['values'][i + 2 + ndat_ptz_eee][
            'errors'
        ][0]['symerror']
        error_value['bgr'] = input2['dependent_variables'][2]['values'][i + 2 + ndat_ptz_eee][
            'errors'
        ][0]['symerror']
        error_value['sys'] = input2['dependent_variables'][3]['values'][i + 2 + ndat_ptz_eee][
            'errors'
        ][0]['symerror']
        error_ptz_eem.append(error_value)

    data_central_ptz_eem_yaml = {'data_cental': data_central_ptz_eem}
    kinematics_ptz_eem_yaml = {'bins': kin_ptz_eem}
    uncertainties_ptz_eem_yaml = {'definitions': error_definition_ptz_eee, 'bins': error_ptz_eem}

    with open('data_central_ptz_eem.yaml', 'w') as f:
        yaml.dump(data_central_ptz_eem_yaml, f, sort_keys=False)

    with open('kinematics_ptz_eem.yaml', 'w') as f:
        yaml.dump(kinematics_ptz_eem_yaml, f, sort_keys=False)

    with open('uncertainties_ptz_eem.yaml', 'w') as f:
        yaml.dump(uncertainties_ptz_eem_yaml, f, sort_keys=False)

    # ptz_emm data

    hepdata_tab3 = "rawdata/table_8.yaml"
    with open(hepdata_tab3, 'r') as f:
        input3 = yaml.safe_load(f)

    for i in range(ndat_ptz_emm):
        data_central_value = input3['dependent_variables'][0]['values'][i + 1]['value']
        data_central_ptz_emm.append(float(data_central_value))
        ptZ_min = ptZ_bins_edges[i]
        ptZ_max = ptZ_bins_edges[i + 1]
        kin_value = {'ptZ': {'min': ptZ_min, 'mid': None, 'max': ptZ_max}}
        kin_ptz_emm.append(kin_value)
        error_value = {}
        error_value['stat'] = input3['dependent_variables'][1]['values'][i + 1]['errors'][0][
            'symerror'
        ]
        error_value['bgr'] = input3['dependent_variables'][2]['values'][i + 1]['errors'][0][
            'symerror'
        ]
        error_value['sys'] = input3['dependent_variables'][3]['values'][i + 1]['errors'][0][
            'symerror'
        ]
        error_ptz_emm.append(error_value)

    error_definition_ptz_emm = {}
    error_definition_ptz_emm['stat'] = {
        'definition': 'statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    error_definition_ptz_emm['bgr'] = {
        'definition': 'background uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    }
    error_definition_ptz_emm['sys'] = {
        'definition': 'systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    }

    data_central_ptz_emm_yaml = {'data_cental': data_central_ptz_emm}
    kinematics_ptz_emm_yaml = {'bins': kin_ptz_emm}
    uncertainties_ptz_emm_yaml = {'definitions': error_definition_ptz_emm, 'bins': error_ptz_emm}

    with open('data_central_ptz_emm.yaml', 'w') as f:
        yaml.dump(data_central_ptz_emm_yaml, f, sort_keys=False)

    with open('kinematics_ptz_emm.yaml', 'w') as f:
        yaml.dump(kinematics_ptz_emm_yaml, f, sort_keys=False)

    with open('uncertainties_ptz_emm.yaml', 'w') as f:
        yaml.dump(uncertainties_ptz_emm_yaml, f, sort_keys=False)

    # ptz_mmm data

    for i in range(ndat_ptz_mmm):
        data_central_value = input3['dependent_variables'][0]['values'][i + 2 + ndat_ptz_emm][
            'value'
        ]
        data_central_ptz_mmm.append(float(data_central_value))
        ptZ_min = ptZ_bins_edges[i]
        ptZ_max = ptZ_bins_edges[i + 1]
        kin_value = {'ptZ': {'min': ptZ_min, 'mid': None, 'max': ptZ_max}}
        kin_ptz_mmm.append(kin_value)
        error_value = {}
        error_value['stat'] = input3['dependent_variables'][1]['values'][i + 2 + ndat_ptz_emm][
            'errors'
        ][0]['symerror']
        error_value['bgr'] = input3['dependent_variables'][2]['values'][i + 2 + ndat_ptz_emm][
            'errors'
        ][0]['symerror']
        error_value['sys'] = input3['dependent_variables'][3]['values'][i + 2 + ndat_ptz_emm][
            'errors'
        ][0]['symerror']
        error_ptz_mmm.append(error_value)

    data_central_ptz_mmm_yaml = {'data_cental': data_central_ptz_mmm}
    kinematics_ptz_mmm_yaml = {'bins': kin_ptz_mmm}
    uncertainties_ptz_mmm_yaml = {'definitions': error_definition_ptz_emm, 'bins': error_ptz_mmm}

    with open('data_central_ptz_mmm.yaml', 'w') as f:
        yaml.dump(data_central_ptz_mmm_yaml, f, sort_keys=False)

    with open('kinematics_ptz_mmm.yaml', 'w') as f:
        yaml.dump(kinematics_ptz_mmm_yaml, f, sort_keys=False)

    with open('uncertainties_ptz_mmm.yaml', 'w') as f:
        yaml.dump(uncertainties_ptz_mmm_yaml, f, sort_keys=False)


magic()
