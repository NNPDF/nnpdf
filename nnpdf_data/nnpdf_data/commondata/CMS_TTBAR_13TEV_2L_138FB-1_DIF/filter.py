import yaml

from nnpdf_data.filter_utils.utils import prettify_float
from nnpdf_data.filter_utils.utils import symmetrize_errors as se

yaml.add_representer(float, prettify_float)


def the_function():

    data_central_dSig_dmttBar = []
    kin_dSig_dmttBar = []
    error_dSig_dmttBar = []
    data_central_dSig_dpTt = []
    kin_dSig_dpTt = []
    error_dSig_dpTt = []
    data_central_dSig_dyt = []
    kin_dSig_dyt = []
    error_dSig_dyt = []
    data_central_d2Sig_dmttBar_dyttBar = []
    kin_d2Sig_dmttBar_dyttBar = []
    error_d2Sig_dmttBar_dyttBar = []

    m_t2 = 172.5 * 172.5

    # dSig_dmttBar

    with open("rawdata/mtt_abs_parton.yaml", "r") as f:
        input1 = yaml.safe_load(f)

    values1 = input1['dependent_variables'][0]['values']

    for i in range(len(values1)):
        data_central_value = values1[i]['value']
        kin_low = input1['independent_variables'][0]['values'][i]['low']
        kin_high = input1['independent_variables'][0]['values'][i]['high']
        kin_value = {
            'm_ttBar': {'min': kin_low, 'mid': None, 'max': kin_high},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
        }
        value_delta = 0
        error_value = {}
        error_value['stat'] = values1[i]['errors'][0]['symerror']
        for j in range(1, len(values1[i]['errors'])):
            se_delta, se_sigma = se(
                values1[i]['errors'][j]['asymerror']['plus'],
                values1[i]['errors'][j]['asymerror']['minus'],
            )
            value_delta += se_delta
            error_value[values1[i]['errors'][j]['label']] = se_sigma
        data_central_value = data_central_value + value_delta
        data_central_dSig_dmttBar.append(data_central_value)
        kin_dSig_dmttBar.append(kin_value)
        error_dSig_dmttBar.append(error_value)

    error_definition_dSig_dmttBar = {}
    error_definition_dSig_dmttBar['stat'] = {
        'description': 'statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    for i in range(1, len(values1[0]['errors'])):
        error_definition_dSig_dmttBar[values1[0]['errors'][i]['label']] = {
            'description': '',
            'treatment': 'MULT',
            'type': 'CORR',
        }

    data_central_dSig_dmttBar_yaml = {'data_central': data_central_dSig_dmttBar}
    kin_dSig_dmttBar_yaml = {'bins': kin_dSig_dmttBar}
    uncertainties_dSig_dmttBar_yaml = {
        'definitions': error_definition_dSig_dmttBar,
        'bins': error_dSig_dmttBar,
    }

    with open('data_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dmttBar_yaml, file, sort_keys=False)
    with open('kinematics_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(kin_dSig_dmttBar_yaml, file, sort_keys=False)
    with open('uncertainties_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_yaml, file, sort_keys=False)

    # dSig_dpTt

    with open("rawdata/ptt_abs_parton.yaml", "r") as f:
        input2 = yaml.safe_load(f)

    values2 = input2['dependent_variables'][0]['values']

    for i in range(len(values2)):
        data_central_value = values2[i]['value']
        kin_low = input2['independent_variables'][0]['values'][i]['low']
        kin_high = input2['independent_variables'][0]['values'][i]['high']
        kin_value = {
            'pT_t': {'min': kin_low, 'mid': None, 'max': kin_high},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
        }
        value_delta = 0
        error_value = {}
        error_value['stat'] = values2[i]['errors'][0]['symerror']
        for j in range(1, len(values2[i]['errors'])):
            se_delta, se_sigma = se(
                values2[i]['errors'][j]['asymerror']['plus'],
                values2[i]['errors'][j]['asymerror']['minus'],
            )
            value_delta += se_delta
            error_value[values2[i]['errors'][j]['label']] = se_sigma
        data_central_value = data_central_value + value_delta
        data_central_dSig_dpTt.append(data_central_value)
        kin_dSig_dpTt.append(kin_value)
        error_dSig_dpTt.append(error_value)

    error_definition_dSig_dpTt = {}
    error_definition_dSig_dpTt['stat'] = {
        'description': 'statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    for i in range(1, len(values2[0]['errors'])):
        error_definition_dSig_dpTt[values2[0]['errors'][i]['label']] = {
            'description': '',
            'treatment': 'MULT',
            'type': 'CORR',
        }

    data_central_dSig_dpTt_yaml = {'data_central': data_central_dSig_dpTt}
    kin_dSig_dpTt_yaml = {'bins': kin_dSig_dpTt}
    uncertainties_dSig_dpTt_yaml = {
        'definitions': error_definition_dSig_dpTt,
        'bins': error_dSig_dpTt,
    }

    with open('data_dSig_dpTt.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dpTt_yaml, file, sort_keys=False)
    with open('kinematics_dSig_dpTt.yaml', 'w') as file:
        yaml.dump(kin_dSig_dpTt_yaml, file, sort_keys=False)
    with open('uncertainties_dSig_dpTt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_yaml, file, sort_keys=False)

    # dSig_dyt

    with open("rawdata/yt_abs_parton.yaml", "r") as f:
        input3 = yaml.safe_load(f)
    values3 = input3['dependent_variables'][0]['values']

    for i in range(len(values3)):
        data_central_value = values3[i]['value']
        kin_low = input3['independent_variables'][0]['values'][i]['low']
        kin_high = input3['independent_variables'][0]['values'][i]['high']
        kin_value = {
            'y_t': {'min': kin_low, 'mid': None, 'max': kin_high},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
        }
        value_delta = 0
        error_value = {}
        error_value['stat'] = values3[i]['errors'][0]['symerror']
        for j in range(1, len(values3[i]['errors'])):
            se_delta, se_sigma = se(
                values3[i]['errors'][j]['asymerror']['plus'],
                values3[i]['errors'][j]['asymerror']['minus'],
            )
            value_delta += se_delta
            error_value[values3[i]['errors'][j]['label']] = se_sigma
        data_central_value = data_central_value + value_delta
        data_central_dSig_dyt.append(data_central_value)
        kin_dSig_dyt.append(kin_value)
        error_dSig_dyt.append(error_value)

    error_definition_dSig_dyt = {}
    error_definition_dSig_dyt['stat'] = {
        'description': 'statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    for i in range(1, len(values3[0]['errors'])):
        error_definition_dSig_dyt[values3[0]['errors'][i]['label']] = {
            'description': '',
            'treatment': 'MULT',
            'type': 'CORR',
        }

    data_central_dSig_dyt_yaml = {'data_central': data_central_dSig_dyt}
    kin_dSig_dyt_yaml = {'bins': kin_dSig_dyt}
    uncertainties_dSig_dyt_yaml = {'definitions': error_definition_dSig_dyt, 'bins': error_dSig_dyt}

    with open('data_dSig_dyt.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dyt_yaml, file, sort_keys=False)
    with open('kinematics_dSig_dyt.yaml', 'w') as file:
        yaml.dump(kin_dSig_dyt_yaml, file, sort_keys=False)
    with open('uncertainties_dSig_dyt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_yaml, file, sort_keys=False)

    # d2Sig_dmttBar_dyttBar

    with open("rawdata/mttytt-abs_parton.yaml", "r") as f:
        input4 = yaml.safe_load(f)
    values4 = input4['dependent_variables'][0]['values']

    for i in range(len(values4)):
        kin_low_mttBar = input4['independent_variables'][1]['values'][i]['low']
        kin_high_mttBar = input4['independent_variables'][1]['values'][i]['high']
        kin_low_yttBar = input4['independent_variables'][0]['values'][i]['low']
        kin_high_yttBar = input4['independent_variables'][0]['values'][i]['high']
        binwidth_mtt = kin_high_mttBar - kin_low_mttBar

        # At the time of implementation (23/07/25), the hepdata entry was not normalized to the bin width.
        # We correct for this here by dividing the error and central value by the bin width.
        data_central_value = values4[i]['value'] / binwidth_mtt

        kin_value = {
            'm_ttBar': {'min': kin_low_mttBar, 'mid': None, 'max': kin_high_mttBar},
            'y_ttBar': {'min': kin_low_yttBar, 'mid': None, 'max': kin_high_yttBar},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
        }
        value_delta = 0
        error_value = {}
        error_value['stat'] = values4[i]['errors'][0]['symerror'] /binwidth_mtt
        for j in range(1, len(values4[i]['errors'])):
            se_delta, se_sigma = se(
                values4[i]['errors'][j]['asymerror']['plus'],
                values4[i]['errors'][j]['asymerror']['minus'],
            )
            se_delta /= binwidth_mtt
            se_sigma /= binwidth_mtt
            value_delta += se_delta
            error_value[values4[i]['errors'][j]['label']] = se_sigma
        data_central_value = data_central_value + value_delta
        data_central_d2Sig_dmttBar_dyttBar.append(data_central_value)
        kin_d2Sig_dmttBar_dyttBar.append(kin_value)
        error_d2Sig_dmttBar_dyttBar.append(error_value)

    error_definition_d2Sig_dmttBar_dyttBar = {}
    error_definition_d2Sig_dmttBar_dyttBar['stat'] = {
        'description': 'statistical uncertainty',
        'treatment': 'ADD',
        'type': 'UNCORR',
    }
    for i in range(1, len(values4[0]['errors'])):
        error_definition_d2Sig_dmttBar_dyttBar[values4[0]['errors'][i]['label']] = {
            'description': '',
            'treatment': 'MULT',
            'type': 'CORR',
        }

    data_central_d2Sig_dmttBar_dyttBar_yaml = {'data_central': data_central_d2Sig_dmttBar_dyttBar}
    kin_d2Sig_dmttBar_dyttBar_yaml = {'bins': kin_d2Sig_dmttBar_dyttBar}
    uncertainties_d2Sig_dmttBar_dyttBar_yaml = {
        'definitions': error_definition_d2Sig_dmttBar_dyttBar,
        'bins': error_d2Sig_dmttBar_dyttBar,
    }

    with open('data_d2Sig_dmttBar_dyttBar.yaml', 'w') as file:
        yaml.dump(data_central_d2Sig_dmttBar_dyttBar_yaml, file, sort_keys=False)
    with open('kinematics_d2Sig_dmttBar_dyttBar.yaml', 'w') as file:
        yaml.dump(kin_d2Sig_dmttBar_dyttBar_yaml, file, sort_keys=False)
    with open('uncertainties_d2Sig_dmttBar_dyttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dmttBar_dyttBar_yaml, file, sort_keys=False)

    # manual_decorr
    for unc in [uncertainties_dSig_dmttBar_yaml, 
                uncertainties_dSig_dpTt_yaml, 
                uncertainties_dSig_dyt_yaml, 
                uncertainties_d2Sig_dmttBar_dyttBar_yaml]:
        for i in ['16', '17', '18']:
            unc['definitions']['BTAG_STATISTIC_'+i]['treatment'] = 'ADD'
            unc['definitions']['BTAG_STATISTIC_'+i]['type'] = 'UNCORR'

    with open('uncertainties_md_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_yaml, file, sort_keys=False)
    with open('uncertainties_md_dSig_dpTt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_yaml, file, sort_keys=False)
    with open('uncertainties_md_dSig_dyt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_yaml, file, sort_keys=False)
    with open('uncertainties_md_d2Sig_dmttBar_dyttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dmttBar_dyttBar_yaml, file, sort_keys=False)


the_function()
