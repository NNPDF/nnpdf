import artunc_1d
from artunc_2d import (
    d2Sig_dmttBar_dpTt_artunc,
    d2Sig_dmttBar_dpTt_norm_artunc,
    d2Sig_dpTt_dyt_artunc,
    d2Sig_dpTt_dyt_norm_artunc,
)
from lumiless_covmat import llcm_mtt, llcm_ptt, llcm_yt, llcm_ytt
import numpy as np
import yaml

from nnpdf_data.filter_utils.utils import covmat_to_artunc as cta
from nnpdf_data.filter_utils.utils import percentage_to_absolute as pta
from nnpdf_data.filter_utils.utils import prettify_float
from nnpdf_data.filter_utils.utils import symmetrize_errors as se

yaml.add_representer(float, prettify_float)


def processData():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    ndata_dSig_dmttBar = metadata['implemented_observables'][0]['ndata']
    ndata_dSig_dmttBar_norm = metadata['implemented_observables'][1]['ndata']
    ndata_dSig_dpTt = metadata['implemented_observables'][2]['ndata']
    ndata_dSig_dpTt_norm = metadata['implemented_observables'][3]['ndata']
    ndata_dSig_dyt = metadata['implemented_observables'][4]['ndata']
    ndata_dSig_dyt_norm = metadata['implemented_observables'][5]['ndata']
    ndata_dSig_dyttBar = metadata['implemented_observables'][6]['ndata']
    ndata_dSig_dyttBar_norm = metadata['implemented_observables'][7]['ndata']
    ndata_d2Sig_dpTt_dyt = metadata['implemented_observables'][8]['ndata']
    ndata_d2Sig_dpTt_dyt_norm = metadata['implemented_observables'][9]['ndata']
    ndata_d2Sig_dmttBar_dpTt = metadata['implemented_observables'][10]['ndata']
    ndata_d2Sig_dmttBar_dpTt_norm = metadata['implemented_observables'][11]['ndata']

    data_central_dSig_dmttBar = []
    data_central_dSig_dmttBar_interspectra = []
    kin_dSig_dmttBar = []
    error_dSig_dmttBar = []
    error_dSig_dmttBar_lumiless = []
    error_dSig_dmttBar_interspectra = []
    data_central_dSig_dmttBar_norm = []
    kin_dSig_dmttBar_norm = []
    error_dSig_dmttBar_norm = []
    data_central_dSig_dpTt = []
    data_central_dSig_dpTt_interspectra = []
    kin_dSig_dpTt = []
    error_dSig_dpTt = []
    error_dSig_dpTt_lumiless = []
    error_dSig_dpTt_interspectra = []
    data_central_dSig_dpTt_norm = []
    kin_dSig_dpTt_norm = []
    error_dSig_dpTt_norm = []
    data_central_dSig_dyt = []
    data_central_dSig_dyt_interspectra = []
    kin_dSig_dyt = []
    error_dSig_dyt = []
    error_dSig_dyt_lumiless = []
    error_dSig_dyt_interspectra = []
    data_central_dSig_dyt_norm = []
    kin_dSig_dyt_norm = []
    error_dSig_dyt_norm = []
    data_central_dSig_dyttBar = []
    data_central_dSig_dyttBar_interspectra = []
    kin_dSig_dyttBar = []
    error_dSig_dyttBar = []
    error_dSig_dyttBar_lumiless = []
    error_dSig_dyttBar_interspectra = []
    data_central_dSig_dyttBar_norm = []
    kin_dSig_dyttBar_norm = []
    error_dSig_dyttBar_norm = []
    data_central_d2Sig_dpTt_dyt = []
    kin_d2Sig_dpTt_dyt = []
    error_d2Sig_dpTt_dyt = []
    data_central_d2Sig_dpTt_dyt_norm = []
    kin_d2Sig_dpTt_dyt_norm = []
    error_d2Sig_dpTt_dyt_norm = []
    data_central_d2Sig_dmttBar_dpTt = []
    kin_d2Sig_dmttBar_dpTt = []
    error_d2Sig_dmttBar_dpTt = []
    data_central_d2Sig_dmttBar_dpTt_norm = []
    kin_d2Sig_dmttBar_dpTt_norm = []
    error_d2Sig_dmttBar_dpTt_norm = []

    covMatArray_dSig_dmttBar_total = []
    covMatArray_dSig_dmttBar_total_norm = []
    covMatArray_dSig_dpTt_total = []
    covMatArray_dSig_dpTt_total_norm = []
    covMatArray_dSig_dyt_total = []
    covMatArray_dSig_dyt_total_norm = []
    covMatArray_dSig_dyttBar_total = []
    covMatArray_dSig_dyttBar_total_norm = []

    covMat_stat, artUnc_stat = artunc_1d.artunc()
    n_artUnc_stat = len(artUnc_stat)

    # dSig_dmttBar data
    hepdata_tables = "rawdata/Table618.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix_total = "rawdata/Table619.yaml"  # statistical + systematics
    with open(covariance_matrix_total, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dmttBar * ndata_dSig_dmttBar):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dmttBar_total.append(covMatEl)

    artUncMat_dSig_dmttBar = cta(ndata_dSig_dmttBar, covMatArray_dSig_dmttBar_total)
    artUncMat_dSig_dmttBar_lumiless = cta(ndata_dSig_dmttBar, llcm_mtt)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        error_value_lumiless = {}
        error_value_interspectra = {}
        # stat error is in the covariance matrix, so we skip the first error in the list, which is the total stat error
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        for j in range(1, len(values[i]['errors'])):
            plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
            minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
            se_delta, se_sigma = se(plus, minus)
            value_delta = value_delta + se_delta
            error_label = values[i]['errors'][j]['label']
            error_value_interspectra[error_label] = se_sigma
        data_central_value_interspectra = values[i]['value'] + value_delta
        data_central_value = values[i]['value']

        data_central_dSig_dmttBar.append(data_central_value)
        data_central_dSig_dmttBar_interspectra.append(data_central_value_interspectra)
        for j in range(ndata_dSig_dmttBar):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dmttBar[i][j])
            error_value_lumiless['ArtUnc_' + str(j + 1)] = float(
                artUncMat_dSig_dmttBar_lumiless[i][j]
            )
        for j in range(n_artUnc_stat):
            error_value_interspectra['ArtUnc_stat_' + str(j + 1)] = float(artUnc_stat[i][j])
        error_dSig_dmttBar.append(error_value)
        error_dSig_dmttBar_lumiless.append(error_value_lumiless)
        error_dSig_dmttBar_interspectra.append(error_value_interspectra)
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        kin_dSig_dmttBar.append(kin_value)

    error_definition_dSig_dmttBar = {}
    error_definition_dSig_dmttBar_interspectra = {}

    # error_definition_dSig_dmttBar['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
        error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
        error_definition_dSig_dmttBar_interspectra[error_name] = {
            'description': '',
            'treatment': 'MULT',
            'type': 'CORR',
        }
    for i in range(ndata_dSig_dmttBar):
        error_definition_dSig_dmttBar['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }
    for i in range(n_artUnc_stat):
        error_definition_dSig_dmttBar_interspectra['ArtUnc_stat_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'ATLAS13TEVTTBAR190807305unc' + str(i + 1),
        }

    data_central_dSig_dmttBar_yaml = {'data_central': data_central_dSig_dmttBar}
    data_central_dSig_dmttBar_interspectra_yaml = {
        'data_central': data_central_dSig_dmttBar_interspectra
    }
    kinematics_dSig_dmttBar_yaml = {'bins': kin_dSig_dmttBar}
    uncertainties_dSig_dmttBar_yaml = {
        'definitions': error_definition_dSig_dmttBar,
        'bins': error_dSig_dmttBar,
    }

    uncertainties_dSig_dmttBar_wo_lumi_yaml = {
        'definitions': error_definition_dSig_dmttBar,
        'bins': error_dSig_dmttBar_lumiless,
    }
    uncertainties_dSig_dmttBar_interspectra_yaml = {
        'definitions': error_definition_dSig_dmttBar_interspectra,
        'bins': error_dSig_dmttBar_interspectra,
    }

    with open('data_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('data_dSig_dmttBar_interspectra.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dmttBar_interspectra_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar_wo-lumi.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_wo_lumi_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar_interspectra.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_interspectra_yaml, file, sort_keys=False)

    # dSig_dmttBar_norm data
    hepdata_tables = "rawdata/Table616.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix = "rawdata/Table617.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dmttBar_norm * ndata_dSig_dmttBar_norm):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dmttBar_total_norm.append(covMatEl)
    artUncMat_dSig_dmttBar_norm = cta(ndata_dSig_dmttBar_norm, covMatArray_dSig_dmttBar_total_norm)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        m_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        m_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value']  # + value_delta
        data_central_dSig_dmttBar_norm.append(data_central_value)
        for j in range(ndata_dSig_dmttBar_norm):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dmttBar_norm[i][j])
        error_dSig_dmttBar_norm.append(error_value)
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        kin_dSig_dmttBar_norm.append(kin_value)

    error_definition_dSig_dmttBar_norm = {}
    # error_definition_dSig_dmttBar_norm['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dmttBar_norm[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dmttBar_norm):
        error_definition_dSig_dmttBar_norm['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }

    data_central_dSig_dmttBar_norm_yaml = {'data_central': data_central_dSig_dmttBar_norm}
    kinematics_dSig_dmttBar_norm_yaml = {'bins': kin_dSig_dmttBar_norm}
    uncertainties_dSig_dmttBar_norm_yaml = {
        'definitions': error_definition_dSig_dmttBar_norm,
        'bins': error_dSig_dmttBar_norm,
    }

    with open('data_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dmttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dmttBar_norm_yaml, file, sort_keys=False)

    # dSig_dpTt data
    hepdata_tables = "rawdata/Table610.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix_total = "rawdata/Table611.yaml"
    with open(covariance_matrix_total, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dpTt * ndata_dSig_dpTt):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dpTt_total.append(covMatEl)

    artUncMat_dSig_dpTt = cta(ndata_dSig_dpTt, covMatArray_dSig_dpTt_total)
    artUncMat_dSig_dpTt_lumiless = cta(ndata_dSig_dpTt, llcm_ptt)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        pT_t_min = input['independent_variables'][0]['values'][i]['low']
        pT_t_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        error_value_lumiless = {}
        error_value_interspectra = {}
        # stat error is in the covariance matrix, so we skip the first error in the list, which is the total stat error
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        for j in range(1, len(values[i]['errors'])):
            plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
            minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
            se_delta, se_sigma = se(plus, minus)
            value_delta = value_delta + se_delta
            error_label = values[i]['errors'][j]['label']
            error_value_interspectra[error_label] = se_sigma
        data_central_value_interspectra = values[i]['value'] + value_delta
        data_central_value = values[i]['value']

        data_central_dSig_dpTt.append(data_central_value)
        data_central_dSig_dpTt_interspectra.append(data_central_value_interspectra)
        for j in range(ndata_dSig_dpTt):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dpTt[i][j])
            error_value_lumiless['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dpTt_lumiless[i][j])
        for j in range(n_artUnc_stat):
            error_value_interspectra['ArtUnc_stat_' + str(j + 1)] = float(artUnc_stat[i][j])
        error_dSig_dpTt.append(error_value)
        error_dSig_dpTt_lumiless.append(error_value_lumiless)
        error_dSig_dpTt_interspectra.append(error_value_interspectra)
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
        }
        kin_dSig_dpTt.append(kin_value)

    error_definition_dSig_dpTt = {}
    error_definition_dSig_dpTt_interspectra = {}

    # error_definition_dSig_dpTt['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
        error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
        error_definition_dSig_dpTt_interspectra[error_name] = {
            'description': '',
            'treatment': 'MULT',
            'type': 'CORR',
        }
    for i in range(ndata_dSig_dpTt):
        error_definition_dSig_dpTt['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }
    for i in range(n_artUnc_stat):
        error_definition_dSig_dpTt_interspectra['ArtUnc_stat_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'ATLAS13TEVTTBAR190807305unc' + str(i + 1),
        }

    data_central_dSig_dpTt_yaml = {'data_central': data_central_dSig_dpTt}
    data_central_dSig_dpTt_interspectra_yaml = {'data_central': data_central_dSig_dpTt_interspectra}
    kinematics_dSig_dpTt_yaml = {'bins': kin_dSig_dpTt}
    uncertainties_dSig_dpTt_yaml = {
        'definitions': error_definition_dSig_dpTt,
        'bins': error_dSig_dpTt,
    }

    uncertainties_dSig_dpTt_wo_lumi_yaml = {
        'definitions': error_definition_dSig_dpTt,
        'bins': error_dSig_dpTt_lumiless,
    }
    uncertainties_dSig_dpTt_interspectra_yaml = {
        'definitions': error_definition_dSig_dpTt_interspectra,
        'bins': error_dSig_dpTt_interspectra,
    }

    with open('data_dSig_dpTt.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dpTt_yaml, file, sort_keys=False)

    with open('data_dSig_dpTt_interspectra.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dpTt_interspectra_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dpTt.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dpTt_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dpTt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dpTt_wo-lumi.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_wo_lumi_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dpTt_interspectra.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_interspectra_yaml, file, sort_keys=False)

    # dSig_dpTt_norm data
    hepdata_tables = "rawdata/Table608.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix = "rawdata/Table609.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dpTt_norm * ndata_dSig_dpTt_norm):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dpTt_total_norm.append(covMatEl)
    artUncMat_dSig_dpTt_norm = cta(ndata_dSig_dpTt_norm, covMatArray_dSig_dpTt_total_norm)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        pT_t_min = input['independent_variables'][0]['values'][i]['low']
        pT_t_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value']  # + value_delta
        data_central_dSig_dpTt_norm.append(data_central_value)
        for j in range(ndata_dSig_dpTt_norm):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dpTt_norm[i][j])
        error_dSig_dpTt_norm.append(error_value)
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
        }
        kin_dSig_dpTt_norm.append(kin_value)

    error_definition_dSig_dpTt_norm = {}
    # error_definition_dSig_dpTt_norm['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dpTt_norm[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dpTt_norm):
        error_definition_dSig_dpTt_norm['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }

    data_central_dSig_dpTt_norm_yaml = {'data_central': data_central_dSig_dpTt_norm}
    kinematics_dSig_dpTt_norm_yaml = {'bins': kin_dSig_dpTt_norm}
    uncertainties_dSig_dpTt_norm_yaml = {
        'definitions': error_definition_dSig_dpTt_norm,
        'bins': error_dSig_dpTt_norm,
    }

    with open('data_dSig_dpTt_norm.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dpTt_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dpTt_norm.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dpTt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dpTt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dpTt_norm_yaml, file, sort_keys=False)

    # dSig_dyt data
    hepdata_tables = "rawdata/Table614.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix_total = "rawdata/Table615.yaml"
    with open(covariance_matrix_total, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dyt * ndata_dSig_dyt):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyt_total.append(covMatEl)

    artUncMat_dSig_dyt = cta(ndata_dSig_dyt, covMatArray_dSig_dyt_total)
    artUncMat_dSig_dyt_lumiless = cta(ndata_dSig_dyt, llcm_yt)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        error_value_lumiless = {}
        error_value_interspectra = {}
        # stat error is in the covariance matrix, so we skip the first error in the list, which is the total stat error
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        for j in range(1, len(values[i]['errors'])):
            plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
            minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
            se_delta, se_sigma = se(plus, minus)
            value_delta = value_delta + se_delta
            error_label = values[i]['errors'][j]['label']
            error_value_interspectra[error_label] = se_sigma
        data_central_value_interspectra = values[i]['value'] + value_delta
        data_central_value = values[i]['value']

        data_central_dSig_dyt.append(data_central_value)
        data_central_dSig_dyt_interspectra.append(data_central_value_interspectra)
        for j in range(ndata_dSig_dyt):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dyt[i][j])
            error_value_lumiless['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dyt_lumiless[i][j])
        for j in range(n_artUnc_stat):
            error_value_interspectra['ArtUnc_stat_' + str(j + 1)] = float(artUnc_stat[i][j])
        error_dSig_dyt.append(error_value)
        error_dSig_dyt_lumiless.append(error_value_lumiless)
        error_dSig_dyt_interspectra.append(error_value_interspectra)
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max},
        }
        kin_dSig_dyt.append(kin_value)

    error_definition_dSig_dyt = {}
    error_definition_dSig_dyt_interspectra = {}

    # error_definition_dSig_dyt['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
        error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
        error_definition_dSig_dyt_interspectra[error_name] = {
            'description': '',
            'treatment': 'MULT',
            'type': 'CORR',
        }
    for i in range(ndata_dSig_dyt):
        error_definition_dSig_dyt['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }
    for i in range(n_artUnc_stat):
        error_definition_dSig_dyt_interspectra['ArtUnc_stat_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'ATLAS13TEVTTBAR190807305unc' + str(i + 1),
        }

    data_central_dSig_dyt_yaml = {'data_central': data_central_dSig_dyt}
    data_central_dSig_dyt_interspectra_yaml = {'data_central': data_central_dSig_dyt_interspectra}
    kinematics_dSig_dyt_yaml = {'bins': kin_dSig_dyt}
    uncertainties_dSig_dyt_yaml = {'definitions': error_definition_dSig_dyt, 'bins': error_dSig_dyt}

    uncertainties_dSig_dyt_wo_lumi_yaml = {
        'definitions': error_definition_dSig_dyt,
        'bins': error_dSig_dyt_lumiless,
    }
    uncertainties_dSig_dyt_interspectra_yaml = {
        'definitions': error_definition_dSig_dyt_interspectra,
        'bins': error_dSig_dyt_interspectra,
    }

    with open('data_dSig_dyt.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dyt_yaml, file, sort_keys=False)

    with open('data_dSig_dyt_interspectra.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dyt_interspectra_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyt.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dyt_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyt.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyt_wo-lumi.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_wo_lumi_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyt_interspectra.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_interspectra_yaml, file, sort_keys=False)

    # dSig_dyt_norm data
    hepdata_tables = "rawdata/Table612.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix = "rawdata/Table613.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dyt_norm * ndata_dSig_dyt_norm):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyt_total_norm.append(covMatEl)
    artUncMat_dSig_dyt_norm = cta(ndata_dSig_dyt_norm, covMatArray_dSig_dyt_total_norm)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_t_min = input['independent_variables'][0]['values'][i]['low']
        y_t_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value']  # + value_delta
        data_central_dSig_dyt_norm.append(data_central_value)
        for j in range(ndata_dSig_dyt_norm):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dyt_norm[i][j])
        error_dSig_dyt_norm.append(error_value)
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max},
        }
        kin_dSig_dyt_norm.append(kin_value)

    error_definition_dSig_dyt_norm = {}
    # error_definition_dSig_dyt_norm['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dyt_norm[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dyt_norm):
        error_definition_dSig_dyt_norm['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }

    data_central_dSig_dyt_norm_yaml = {'data_central': data_central_dSig_dyt_norm}
    kinematics_dSig_dyt_norm_yaml = {'bins': kin_dSig_dyt_norm}
    uncertainties_dSig_dyt_norm_yaml = {
        'definitions': error_definition_dSig_dyt_norm,
        'bins': error_dSig_dyt_norm,
    }

    with open('data_dSig_dyt_norm.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dyt_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyt_norm.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dyt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyt_norm_yaml, file, sort_keys=False)

    # dSig_dyttBar data
    hepdata_tables = "rawdata/Table626.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix_total = "rawdata/Table627.yaml"
    with open(covariance_matrix_total, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dyttBar * ndata_dSig_dyttBar):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyttBar_total.append(covMatEl)

    artUncMat_dSig_dyttBar = cta(ndata_dSig_dyttBar, covMatArray_dSig_dyttBar_total)
    artUncMat_dSig_dyttBar_lumiless = cta(ndata_dSig_dyttBar, llcm_ytt)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        error_value_lumiless = {}
        error_value_interspectra = {}
        # stat error is in the covariance matrix, so we skip the first error in the list, which is the total stat error
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        for j in range(1, len(values[i]['errors'])):
            plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
            minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
            se_delta, se_sigma = se(plus, minus)
            value_delta = value_delta + se_delta
            error_label = values[i]['errors'][j]['label']
            error_value_interspectra[error_label] = se_sigma
        data_central_value_interspectra = values[i]['value'] + value_delta
        data_central_value = values[i]['value']

        data_central_dSig_dyttBar.append(data_central_value)
        data_central_dSig_dyttBar_interspectra.append(data_central_value_interspectra)
        for j in range(ndata_dSig_dyttBar):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dyttBar[i][j])
            error_value_lumiless['ArtUnc_' + str(j + 1)] = float(
                artUncMat_dSig_dyttBar_lumiless[i][j]
            )
        for j in range(n_artUnc_stat):
            error_value_interspectra['ArtUnc_stat_' + str(j + 1)] = float(artUnc_stat[i][j])
        error_dSig_dyttBar.append(error_value)
        error_dSig_dyttBar_lumiless.append(error_value_lumiless)
        error_dSig_dyttBar_interspectra.append(error_value_interspectra)
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max},
        }
        kin_dSig_dyttBar.append(kin_value)

    error_definition_dSig_dyttBar = {}
    error_definition_dSig_dyttBar_interspectra = {}

    # error_definition_dSig_dyttBar['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
        error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
        error_definition_dSig_dyttBar_interspectra[error_name] = {
            'description': '',
            'treatment': 'MULT',
            'type': 'CORR',
        }
    for i in range(ndata_dSig_dyttBar):
        error_definition_dSig_dyttBar['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }
    for i in range(n_artUnc_stat):
        error_definition_dSig_dyttBar_interspectra['ArtUnc_stat_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'ATLAS13TEVTTBAR190807305unc' + str(i + 1),
        }

    data_central_dSig_dyttBar_yaml = {'data_central': data_central_dSig_dyttBar}
    data_central_dSig_dyttBar_interspectra_yaml = {
        'data_central': data_central_dSig_dyttBar_interspectra
    }
    kinematics_dSig_dyttBar_yaml = {'bins': kin_dSig_dyttBar}
    uncertainties_dSig_dyttBar_yaml = {
        'definitions': error_definition_dSig_dyttBar,
        'bins': error_dSig_dyttBar,
    }

    uncertainties_dSig_dyttBar_wo_lumi_yaml = {
        'definitions': error_definition_dSig_dyttBar,
        'bins': error_dSig_dyttBar_lumiless,
    }
    uncertainties_dSig_dyttBar_interspectra_yaml = {
        'definitions': error_definition_dSig_dyttBar_interspectra,
        'bins': error_dSig_dyttBar_interspectra,
    }

    with open('data_dSig_dyttBar.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('data_dSig_dyttBar_interspectra.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dyttBar_interspectra_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyttBar.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar_wo-lumi.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_wo_lumi_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar_interspectra.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_interspectra_yaml, file, sort_keys=False)

    # dSig_dyttBar_norm data
    hepdata_tables = "rawdata/Table624.yaml"
    with open(hepdata_tables, 'r') as file:
        input = yaml.safe_load(file)

    covariance_matrix = "rawdata/Table625.yaml"
    with open(covariance_matrix, 'r') as file2:
        input2 = yaml.safe_load(file2)
    for i in range(ndata_dSig_dyttBar_norm * ndata_dSig_dyttBar_norm):
        covMatEl = input2['dependent_variables'][0]['values'][i]['value']
        covMatArray_dSig_dyttBar_total_norm.append(covMatEl)
    artUncMat_dSig_dyttBar_norm = cta(ndata_dSig_dyttBar_norm, covMatArray_dSig_dyttBar_total_norm)

    sqrts = float(input['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values = input['dependent_variables'][0]['values']

    for i in range(len(values)):
        y_ttBar_min = input['independent_variables'][0]['values'][i]['low']
        y_ttBar_max = input['independent_variables'][0]['values'][i]['high']
        value_delta = 0
        error_value = {}
        # error_label = str(values[i]['errors'][0]['label'])
        # error_value[error_label] = pta(values[i]['errors'][0]['symerror'], values[i]['value'])
        # for j in range(1, len(values[i]['errors'])):
        #     plus = pta(values[i]['errors'][j]['asymerror']['plus'], values[i]['value'])
        #     minus = pta(values[i]['errors'][j]['asymerror']['minus'], values[i]['value'])
        #     se_delta, se_sigma = se(plus, minus)
        #     value_delta = value_delta + se_delta
        #     error_label = values[i]['errors'][j]['label']
        #     error_value[error_label] = se_sigma
        data_central_value = values[i]['value']  # + value_delta
        data_central_dSig_dyttBar_norm.append(data_central_value)
        for j in range(ndata_dSig_dyttBar_norm):
            error_value['ArtUnc_' + str(j + 1)] = float(artUncMat_dSig_dyttBar_norm[i][j])
        error_dSig_dyttBar_norm.append(error_value)
        kin_value = {
            'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'y_ttBar': {'min': y_ttBar_min, 'mid': None, 'max': y_ttBar_max},
        }
        kin_dSig_dyttBar_norm.append(kin_value)

    error_definition_dSig_dyttBar_norm = {}
    # error_definition_dSig_dyttBar_norm['stat'] = {'description': 'total statistical uncertainty', 'treatment': 'ADD', 'type': 'UNCORR'}
    # for i in range(1, len(input['dependent_variables'][0]['values'][0]['errors'])):
    #     error_name = input['dependent_variables'][0]['values'][0]['errors'][i]['label']
    #     error_definition_dSig_dyttBar_norm[error_name] = {'description': '', 'treatment': 'MULT', 'type': 'CORR'}
    for i in range(ndata_dSig_dyttBar_norm):
        error_definition_dSig_dyttBar_norm['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }

    data_central_dSig_dyttBar_norm_yaml = {'data_central': data_central_dSig_dyttBar_norm}
    kinematics_dSig_dyttBar_norm_yaml = {'bins': kin_dSig_dyttBar_norm}
    uncertainties_dSig_dyttBar_norm_yaml = {
        'definitions': error_definition_dSig_dyttBar_norm,
        'bins': error_dSig_dyttBar_norm,
    }

    with open('data_dSig_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(data_central_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    with open('kinematics_dSig_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(kinematics_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    with open('uncertainties_dSig_dyttBar_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_dSig_dyttBar_norm_yaml, file, sort_keys=False)

    # d2Sig_dpTt_dyt data
    hepdata_tables = "rawdata/Table649.yaml"
    with open(hepdata_tables, 'r') as file:
        input1 = yaml.safe_load(file)
    hepdata_tables = "rawdata/Table650.yaml"
    with open(hepdata_tables, 'r') as file:
        input2 = yaml.safe_load(file)
    hepdata_tables = "rawdata/Table651.yaml"
    with open(hepdata_tables, 'r') as file:
        input3 = yaml.safe_load(file)

    sqrts = float(input1['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values1 = input1['dependent_variables'][0]['values']
    values2 = input2['dependent_variables'][0]['values']
    values3 = input3['dependent_variables'][0]['values']

    for i in range(len(values1)):
        pT_t_min = input1['independent_variables'][0]['values'][i]['low']
        pT_t_max = input1['independent_variables'][0]['values'][i]['high']
        y_t_min = 0.0
        y_t_max = 0.75
        data_central_value = values1[i]['value']
        data_central_d2Sig_dpTt_dyt.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dpTt_dyt):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_artunc[i][j])
        error_d2Sig_dpTt_dyt.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max},
        }
        kin_d2Sig_dpTt_dyt.append(kin_value)
    for i in range(len(values2)):
        pT_t_min = input2['independent_variables'][0]['values'][i]['low']
        pT_t_max = input2['independent_variables'][0]['values'][i]['high']
        y_t_min = 0.75
        y_t_max = 1.5
        data_central_value = values2[i]['value']
        data_central_d2Sig_dpTt_dyt.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dpTt_dyt):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_artunc[i + 5][j])
        error_d2Sig_dpTt_dyt.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max},
        }
        kin_d2Sig_dpTt_dyt.append(kin_value)
    for i in range(len(values3)):
        pT_t_min = input3['independent_variables'][0]['values'][i]['low']
        pT_t_max = input3['independent_variables'][0]['values'][i]['high']
        y_t_min = 1.5
        y_t_max = 2.5
        data_central_value = values3[i]['value']
        data_central_d2Sig_dpTt_dyt.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dpTt_dyt):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_artunc[i + 9][j])
        error_d2Sig_dpTt_dyt.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max},
        }
        kin_d2Sig_dpTt_dyt.append(kin_value)

    error_definition_d2Sig_dpTt_dyt = {}
    for i in range(ndata_d2Sig_dpTt_dyt):
        error_definition_d2Sig_dpTt_dyt['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }

    data_central_d2Sig_dpTt_dyt_yaml = {'data_central': data_central_d2Sig_dpTt_dyt}
    kinematics_d2Sig_dpTt_dyt_yaml = {'bins': kin_d2Sig_dpTt_dyt}
    uncertainties_d2Sig_dpTt_dyt_yaml = {
        'definitions': error_definition_d2Sig_dpTt_dyt,
        'bins': error_d2Sig_dpTt_dyt,
    }

    with open('data_d2Sig_dpTt_dyt.yaml', 'w') as file:
        yaml.dump(data_central_d2Sig_dpTt_dyt_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dpTt_dyt.yaml', 'w') as file:
        yaml.dump(kinematics_d2Sig_dpTt_dyt_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dpTt_dyt.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dpTt_dyt_yaml, file, sort_keys=False)

    # d2Sig_dpTt_dyt_norm data
    hepdata_tables = "rawdata/Table640.yaml"
    with open(hepdata_tables, 'r') as file:
        input1 = yaml.safe_load(file)
    hepdata_tables = "rawdata/Table641.yaml"
    with open(hepdata_tables, 'r') as file:
        input2 = yaml.safe_load(file)
    hepdata_tables = "rawdata/Table642.yaml"
    with open(hepdata_tables, 'r') as file:
        input3 = yaml.safe_load(file)

    sqrts = float(input1['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values1 = input1['dependent_variables'][0]['values']
    values2 = input2['dependent_variables'][0]['values']
    values3 = input3['dependent_variables'][0]['values']

    for i in range(len(values1)):
        pT_t_min = input1['independent_variables'][0]['values'][i]['low']
        pT_t_max = input1['independent_variables'][0]['values'][i]['high']
        y_t_min = 0.0
        y_t_max = 0.75
        data_central_value = values1[i]['value']
        data_central_d2Sig_dpTt_dyt_norm.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dpTt_dyt_norm):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_norm_artunc[i][j])
        error_d2Sig_dpTt_dyt_norm.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max},
        }
        kin_d2Sig_dpTt_dyt_norm.append(kin_value)
    for i in range(len(values2)):
        pT_t_min = input2['independent_variables'][0]['values'][i]['low']
        pT_t_max = input2['independent_variables'][0]['values'][i]['high']
        y_t_min = 0.75
        y_t_max = 1.5
        data_central_value = values2[i]['value']
        data_central_d2Sig_dpTt_dyt_norm.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dpTt_dyt_norm):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_norm_artunc[i + 5][j])
        error_d2Sig_dpTt_dyt_norm.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max},
        }
        kin_d2Sig_dpTt_dyt_norm.append(kin_value)
    for i in range(len(values3)):
        pT_t_min = input3['independent_variables'][0]['values'][i]['low']
        pT_t_max = input3['independent_variables'][0]['values'][i]['high']
        y_t_min = 1.5
        y_t_max = 2.5
        data_central_value = values3[i]['value']
        data_central_d2Sig_dpTt_dyt_norm.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dpTt_dyt_norm):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_norm_artunc[i + 9][j])
        error_d2Sig_dpTt_dyt_norm.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'y_t': {'min': y_t_min, 'mid': None, 'max': y_t_max},
        }
        kin_d2Sig_dpTt_dyt_norm.append(kin_value)

    error_definition_d2Sig_dpTt_dyt_norm = {}
    for i in range(ndata_d2Sig_dpTt_dyt_norm):
        error_definition_d2Sig_dpTt_dyt_norm['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }

    data_central_d2Sig_dpTt_dyt_norm_yaml = {'data_central': data_central_d2Sig_dpTt_dyt_norm}
    kinematics_d2Sig_dpTt_dyt_norm_yaml = {'bins': kin_d2Sig_dpTt_dyt_norm}
    uncertainties_d2Sig_dpTt_dyt_norm_yaml = {
        'definitions': error_definition_d2Sig_dpTt_dyt_norm,
        'bins': error_d2Sig_dpTt_dyt_norm,
    }

    with open('data_d2Sig_dpTt_dyt_norm.yaml', 'w') as file:
        yaml.dump(data_central_d2Sig_dpTt_dyt_norm_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dpTt_dyt_norm.yaml', 'w') as file:
        yaml.dump(kinematics_d2Sig_dpTt_dyt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dpTt_dyt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dpTt_dyt_norm_yaml, file, sort_keys=False)

    # d2Sig_dmttBar_dpTt data
    hepdata_tables = "rawdata/Table700.yaml"
    with open(hepdata_tables, 'r') as file:
        input1 = yaml.safe_load(file)
    hepdata_tables = "rawdata/Table701.yaml"
    with open(hepdata_tables, 'r') as file:
        input2 = yaml.safe_load(file)
    hepdata_tables = "rawdata/Table702.yaml"
    with open(hepdata_tables, 'r') as file:
        input3 = yaml.safe_load(file)
    hepdata_tables = "rawdata/Table703.yaml"
    with open(hepdata_tables, 'r') as file:
        input4 = yaml.safe_load(file)

    sqrts = float(input1['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values1 = input1['dependent_variables'][0]['values']
    values2 = input2['dependent_variables'][0]['values']
    values3 = input3['dependent_variables'][0]['values']
    values4 = input4['dependent_variables'][0]['values']

    for i in range(len(values1)):
        pT_t_min = input1['independent_variables'][0]['values'][i]['low']
        pT_t_max = input1['independent_variables'][0]['values'][i]['high']
        m_ttBar_min = 325
        m_ttBar_max = 500
        data_central_value = values1[i]['value']
        data_central_d2Sig_dmttBar_dpTt.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dmttBar_dpTt):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_artunc[i][j])
        error_d2Sig_dmttBar_dpTt.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        kin_d2Sig_dmttBar_dpTt.append(kin_value)
    for i in range(len(values2)):
        pT_t_min = input2['independent_variables'][0]['values'][i]['low']
        pT_t_max = input2['independent_variables'][0]['values'][i]['high']
        m_ttBar_min = 500
        m_ttBar_max = 700
        data_central_value = values2[i]['value']
        data_central_d2Sig_dmttBar_dpTt.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dmttBar_dpTt):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_artunc[i + 3][j])
        error_d2Sig_dmttBar_dpTt.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        kin_d2Sig_dmttBar_dpTt.append(kin_value)
    for i in range(len(values3)):
        pT_t_min = input3['independent_variables'][0]['values'][i]['low']
        pT_t_max = input3['independent_variables'][0]['values'][i]['high']
        m_ttBar_min = 700
        m_ttBar_max = 1000
        data_central_value = values3[i]['value']
        data_central_d2Sig_dmttBar_dpTt.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dmttBar_dpTt):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_artunc[i + 7][j])
        error_d2Sig_dmttBar_dpTt.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        kin_d2Sig_dmttBar_dpTt.append(kin_value)
    for i in range(len(values4)):
        pT_t_min = input4['independent_variables'][0]['values'][i]['low']
        pT_t_max = input4['independent_variables'][0]['values'][i]['high']
        m_ttBar_min = 1000
        m_ttBar_max = 2000
        data_central_value = values4[i]['value']
        data_central_d2Sig_dmttBar_dpTt.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dmttBar_dpTt):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_artunc[i + 12][j])
        error_d2Sig_dmttBar_dpTt.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        kin_d2Sig_dmttBar_dpTt.append(kin_value)

    error_definition_d2Sig_dmttBar_dpTt = {}
    for i in range(ndata_d2Sig_dmttBar_dpTt):
        error_definition_d2Sig_dmttBar_dpTt['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }

    data_central_d2Sig_dmttBar_dpTt_yaml = {'data_central': data_central_d2Sig_dmttBar_dpTt}
    kinematics_d2Sig_dmttBar_dpTt_yaml = {'bins': kin_d2Sig_dmttBar_dpTt}
    uncertainties_d2Sig_dmttBar_dpTt_yaml = {
        'definitions': error_definition_d2Sig_dmttBar_dpTt,
        'bins': error_d2Sig_dmttBar_dpTt,
    }

    with open('data_d2Sig_dmttBar_dpTt.yaml', 'w') as file:
        yaml.dump(data_central_d2Sig_dmttBar_dpTt_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dmttBar_dpTt.yaml', 'w') as file:
        yaml.dump(kinematics_d2Sig_dmttBar_dpTt_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dmttBar_dpTt.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dmttBar_dpTt_yaml, file, sort_keys=False)

    # d2Sig_dmttBar_dpTt_norm data
    hepdata_tables = "rawdata/Table686.yaml"
    with open(hepdata_tables, 'r') as file:
        input1 = yaml.safe_load(file)
    hepdata_tables = "rawdata/Table687.yaml"
    with open(hepdata_tables, 'r') as file:
        input2 = yaml.safe_load(file)
    hepdata_tables = "rawdata/Table688.yaml"
    with open(hepdata_tables, 'r') as file:
        input3 = yaml.safe_load(file)
    hepdata_tables = "rawdata/Table689.yaml"
    with open(hepdata_tables, 'r') as file:
        input4 = yaml.safe_load(file)

    sqrts = float(input1['dependent_variables'][0]['qualifiers'][0]['value'])
    m_t2 = 29756.25
    values1 = input1['dependent_variables'][0]['values']
    values2 = input2['dependent_variables'][0]['values']
    values3 = input3['dependent_variables'][0]['values']
    values4 = input4['dependent_variables'][0]['values']

    for i in range(len(values1)):
        pT_t_min = input1['independent_variables'][0]['values'][i]['low']
        pT_t_max = input1['independent_variables'][0]['values'][i]['high']
        m_ttBar_min = 325
        m_ttBar_max = 500
        data_central_value = values1[i]['value']
        data_central_d2Sig_dmttBar_dpTt_norm.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dmttBar_dpTt_norm):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_norm_artunc[i][j])
        error_d2Sig_dmttBar_dpTt_norm.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        kin_d2Sig_dmttBar_dpTt_norm.append(kin_value)
    for i in range(len(values2)):
        pT_t_min = input2['independent_variables'][0]['values'][i]['low']
        pT_t_max = input2['independent_variables'][0]['values'][i]['high']
        m_ttBar_min = 500
        m_ttBar_max = 700
        data_central_value = values2[i]['value']
        data_central_d2Sig_dmttBar_dpTt_norm.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dmttBar_dpTt_norm):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_norm_artunc[i + 3][j])
        error_d2Sig_dmttBar_dpTt_norm.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        kin_d2Sig_dmttBar_dpTt_norm.append(kin_value)
    for i in range(len(values3)):
        pT_t_min = input3['independent_variables'][0]['values'][i]['low']
        pT_t_max = input3['independent_variables'][0]['values'][i]['high']
        m_ttBar_min = 700
        m_ttBar_max = 1000
        data_central_value = values3[i]['value']
        data_central_d2Sig_dmttBar_dpTt_norm.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dmttBar_dpTt_norm):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_norm_artunc[i + 7][j])
        error_d2Sig_dmttBar_dpTt_norm.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        kin_d2Sig_dmttBar_dpTt_norm.append(kin_value)
    for i in range(len(values4)):
        pT_t_min = input4['independent_variables'][0]['values'][i]['low']
        pT_t_max = input4['independent_variables'][0]['values'][i]['high']
        m_ttBar_min = 1000
        m_ttBar_max = 2000
        data_central_value = values4[i]['value']
        data_central_d2Sig_dmttBar_dpTt_norm.append(data_central_value)
        error_value = {}
        for j in range(ndata_d2Sig_dmttBar_dpTt_norm):
            error_value['ArtUnc_' + str(j + 1)] = float(d2Sig_dmttBar_dpTt_norm_artunc[i + 12][j])
        error_d2Sig_dmttBar_dpTt_norm.append(error_value)
        kin_value = {
            # 'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            'm_t2': {'min': None, 'mid': m_t2, 'max': None},
            'pT_t': {'min': pT_t_min, 'mid': None, 'max': pT_t_max},
            'm_ttBar': {'min': m_ttBar_min, 'mid': None, 'max': m_ttBar_max},
        }
        kin_d2Sig_dmttBar_dpTt_norm.append(kin_value)

    error_definition_d2Sig_dmttBar_dpTt_norm = {}
    for i in range(ndata_d2Sig_dmttBar_dpTt_norm):
        error_definition_d2Sig_dmttBar_dpTt_norm['ArtUnc_' + str(i + 1)] = {
            'definition': 'artificial uncertainty ' + str(i + 1),
            'treatment': 'ADD',
            'type': 'CORR',
        }

    data_central_d2Sig_dmttBar_dpTt_norm_yaml = {
        'data_central': data_central_d2Sig_dmttBar_dpTt_norm
    }
    kinematics_d2Sig_dmttBar_dpTt_norm_yaml = {'bins': kin_d2Sig_dmttBar_dpTt_norm}
    uncertainties_d2Sig_dmttBar_dpTt_norm_yaml = {
        'definitions': error_definition_d2Sig_dmttBar_dpTt_norm,
        'bins': error_d2Sig_dmttBar_dpTt_norm,
    }

    with open('data_d2Sig_dmttBar_dpTt_norm.yaml', 'w') as file:
        yaml.dump(data_central_d2Sig_dmttBar_dpTt_norm_yaml, file, sort_keys=False)

    with open('kinematics_d2Sig_dmttBar_dpTt_norm.yaml', 'w') as file:
        yaml.dump(kinematics_d2Sig_dmttBar_dpTt_norm_yaml, file, sort_keys=False)

    with open('uncertainties_d2Sig_dmttBar_dpTt_norm.yaml', 'w') as file:
        yaml.dump(uncertainties_d2Sig_dmttBar_dpTt_norm_yaml, file, sort_keys=False)


processData()
