import yaml

from nnpdf_data.filter_utils.utils import cormat_to_covmat as ctc
from nnpdf_data.filter_utils.utils import covmat_to_artunc as cta
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def errorinator(error_array, corr_mat, cov_mat_switch):

    if cov_mat_switch:
        total_error_array = []
        for i in range(len(error_array)):
            bin_err_arr = error_array[i]
            bin_error = 0
            for j in range(len(bin_err_arr)):
                bin_error += bin_err_arr[j][1] ** 2
            bin_error = bin_error**0.5
            total_error_array.append(bin_error)

        cov_mat = ctc(total_error_array, corr_mat)
        art_unc = cta(len(error_array), cov_mat, 0)
        error_definition = {}
        error_value = []
        for i in range(len(art_unc)):
            error_definition['ArtUnc_' + str(i + 1)] = {
                'definition': 'artificial uncertainty ' + str(i + 1),
                'treatment': 'ADD',
                'type': 'CORR',
            }
            error_value_bin = {}
            for j in range(len(art_unc[i])):
                error_value_bin['ArtUnc_' + str(j + 1)] = float(art_unc[i][j])
            error_value.append(error_value_bin)

        uncertainties = {'definitions': error_definition, 'bins': error_value}

    else:
        error_definition = {}
        error_definition['stat'] = {
            'definition': 'statistical uncertainty',
            'treatment': 'ADD',
            'type': 'UNCORR',
        }
        for i in range(1, len(error_array[0])):
            error_definition[error_array[0][i][0]] = {
                'definition': error_array[0][i][0],
                'treatment': 'MULT',
                'type': 'CORR',
            }
        error_value = []
        for i in range(len(error_array)):
            error_value_bin = {}
            for j in range(len(error_array[i])):
                error_value_bin[error_array[i][j][0]] = float(error_array[i][j][1])
            error_value.append(error_value_bin)
        uncertainties = {'definitions': error_definition, 'bins': error_value}

    return uncertainties


def potato():

    with open('metadata.yaml', 'r') as f:
        metadata = yaml.safe_load(f)

    ndata_mem = metadata['implemented_observables'][0]['ndata']
    ndata_ptem = metadata['implemented_observables'][1]['ndata']
    ndata_yem = metadata['implemented_observables'][2]['ndata']

    data_central_mem = []
    kin_mem = []
    error_mem = []

    data_central_ptem = []
    kin_ptem = []
    error_ptem = []

    data_central_yem = []
    kin_yem = []
    error_yem = []

    # mem data

    hepdata_tab1 = "rawdata/data7.yaml"
    with open(hepdata_tab1, 'r') as f:
        input1 = yaml.safe_load(f)

    hepdata_tab2 = "rawdata/data9.yaml"
    with open(hepdata_tab2, 'r') as f:
        cormat1 = yaml.safe_load(f)

    error_array = []
    cormat_list = []
    for i in range(ndata_mem):
        data_central_value = input1['dependent_variables'][0]['values'][i]['value']
        data_central_mem.append(data_central_value)
        m_em_min = input1['independent_variables'][0]['values'][i]['low']
        m_em_max = input1['independent_variables'][0]['values'][i]['high']
        kin_value = {'m_em': {'min': m_em_min, 'mid': None, 'max': m_em_max}}
        kin_mem.append(kin_value)
        error_bin = []
        for j in range(len(input1['dependent_variables'][0]['values'][i]['errors'])):
            error_bin.append(
                (
                    input1['dependent_variables'][0]['values'][i]['errors'][j]['label'],
                    input1['dependent_variables'][0]['values'][i]['errors'][j]['symerror'],
                )
            )
        error_array.append(error_bin)
    for i in range(len(cormat1['dependent_variables'][0]['values'])):
        cormat_list.append(cormat1['dependent_variables'][0]['values'][i]['value'])
    if len(cormat_list) != len(error_array) ** 2:
        raise ValueError(
            "The length of the correlation matrix does not match the number of data points."
        )

    data_central_mem_yaml = {'data_central': data_central_mem}
    kinematics_mem_yaml = {'bins': kin_mem}
    uncertainties_mem_yaml = errorinator(error_array, cormat_list, True)

    with open('data_mem.yaml', 'w') as f:
        yaml.dump(data_central_mem_yaml, f, sort_keys=False)
    with open('kinematics_mem.yaml', 'w') as f:
        yaml.dump(kinematics_mem_yaml, f, sort_keys=False)
    with open('uncertainties_mem.yaml', 'w') as f:
        yaml.dump(uncertainties_mem_yaml, f, sort_keys=False)

    # ptem data
    hepdata_tab3 = "rawdata/data10.yaml"
    with open(hepdata_tab3, 'r') as f:
        input2 = yaml.safe_load(f)

    hepdata_tab4 = "rawdata/data12.yaml"
    with open(hepdata_tab4, 'r') as f:
        cormat2 = yaml.safe_load(f)

    error_array = []
    cormat_list = []
    for i in range(ndata_ptem):
        data_central_value = input2['dependent_variables'][0]['values'][i]['value']
        data_central_ptem.append(data_central_value)
        pt_em_min = input2['independent_variables'][0]['values'][i]['low']
        pt_em_max = input2['independent_variables'][0]['values'][i]['high']
        kin_value = {'pt_em': {'min': pt_em_min, 'mid': None, 'max': pt_em_max}}
        kin_ptem.append(kin_value)
        error_bin = []
        for j in range(len(input2['dependent_variables'][0]['values'][i]['errors'])):
            error_bin.append(
                (
                    input2['dependent_variables'][0]['values'][i]['errors'][j]['label'],
                    input2['dependent_variables'][0]['values'][i]['errors'][j]['symerror'],
                )
            )
        error_array.append(error_bin)
    for i in range(len(cormat2['dependent_variables'][0]['values'])):
        cormat_list.append(cormat2['dependent_variables'][0]['values'][i]['value'])
    if len(cormat_list) != len(error_array) ** 2:
        raise ValueError(
            "The length of the correlation matrix does not match the number of data points."
        )

    data_central_ptem_yaml = {'data_central': data_central_ptem}
    kinematics_ptem_yaml = {'bins': kin_ptem}
    uncertainties_ptem_yaml = errorinator(error_array, cormat_list, True)

    with open('data_ptem.yaml', 'w') as f:
        yaml.dump(data_central_ptem_yaml, f, sort_keys=False)
    with open('kinematics_ptem.yaml', 'w') as f:
        yaml.dump(kinematics_ptem_yaml, f, sort_keys=False)
    with open('uncertainties_ptem.yaml', 'w') as f:
        yaml.dump(uncertainties_ptem_yaml, f, sort_keys=False)

    # yem data
    hepdata_tab5 = "rawdata/data13.yaml"
    with open(hepdata_tab5, 'r') as f:
        input3 = yaml.safe_load(f)
    hepdata_tab6 = "rawdata/data15.yaml"
    with open(hepdata_tab6, 'r') as f:
        cormat3 = yaml.safe_load(f)

    error_array = []
    cormat_list = []

    for i in range(ndata_yem):
        data_central_value = input3['dependent_variables'][0]['values'][i]['value']
        data_central_yem.append(data_central_value)
        y_em_min = input3['independent_variables'][0]['values'][i]['low']
        y_em_max = input3['independent_variables'][0]['values'][i]['high']
        kin_value = {'y_em': {'min': y_em_min, 'mid': None, 'max': y_em_max}}
        kin_yem.append(kin_value)
        error_bin = []
        for j in range(len(input3['dependent_variables'][0]['values'][i]['errors'])):
            error_bin.append(
                (
                    input3['dependent_variables'][0]['values'][i]['errors'][j]['label'],
                    input3['dependent_variables'][0]['values'][i]['errors'][j]['symerror'],
                )
            )
        error_array.append(error_bin)
    for i in range(len(cormat3['dependent_variables'][0]['values'])):
        cormat_list.append(cormat3['dependent_variables'][0]['values'][i]['value'])
    if len(cormat_list) != len(error_array) ** 2:
        raise ValueError(
            "The length of the correlation matrix does not match the number of data points."
        )

    data_central_yem_yaml = {'data_central': data_central_yem}
    kinematics_yem_yaml = {'bins': kin_yem}
    uncertainties_yem_yaml = errorinator(error_array, cormat_list, True)

    with open('data_yem.yaml', 'w') as f:
        yaml.dump(data_central_yem_yaml, f, sort_keys=False)
    with open('kinematics_yem.yaml', 'w') as f:
        yaml.dump(kinematics_yem_yaml, f, sort_keys=False)
    with open('uncertainties_yem.yaml', 'w') as f:
        yaml.dump(uncertainties_yem_yaml, f, sort_keys=False)


potato()
