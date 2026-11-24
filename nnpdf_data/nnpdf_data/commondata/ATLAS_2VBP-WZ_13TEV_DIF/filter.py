import yaml

from nnpdf_data.filter_utils.utils import cormat_to_covmat as ctc
from nnpdf_data.filter_utils.utils import covmat_to_artunc as cta
from nnpdf_data.filter_utils.utils import percentage_to_absolute as pta
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
        error_definition['sys,uncor'] ={
            'definition': 'systematic uncertainty, uncorrelated',
            'treatment': 'ADD',
            'type': 'UNCORR',
        }
        for i in range(2, len(error_array[0])):
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


def brrr():

    with open('metadata.yaml', 'r') as f:
        metadata = yaml.safe_load(f)

    ndata_ptz = metadata['implemented_observables'][0]['ndata']
    ndata_ptw = metadata['implemented_observables'][1]['ndata']
    ndata_mtwz = metadata['implemented_observables'][2]['ndata']

    data_central_ptz = []
    kin_ptz = []
    error_ptz = []

    data_central_ptw = []
    kin_ptw = []
    error_ptw = []

    data_central_mtwz = []
    kin_mtwz = []
    error_mtwz = []

    # ptz data

    hepdata_tab1 = "rawdata/data8.yaml"
    with open(hepdata_tab1, 'r') as f:
        input1 = yaml.safe_load(f)
    
    hepdata_tab2 = "rawdata/data9.yaml"
    with open(hepdata_tab2, 'r') as f:
        cormat1 = yaml.safe_load(f)

    error_array = []
    cormat_list = []
    for i in range(ndata_ptz):
        data_central_value = float(input1['dependent_variables'][0]['values'][i]['value'])
        data_central_ptz.append(data_central_value)
        ptZ_min = input1['independent_variables'][0]['values'][i]['low']
        ptZ_max = input1['independent_variables'][0]['values'][i]['high']
        kin_value = {
            'ptZ': {'min': ptZ_min, 'mid': None, 'max': ptZ_max},
        }
        kin_ptz.append(kin_value)
        error_bin = []
        for j in range(len(input1['dependent_variables'][0]['values'][i]['errors'])):
            error_bin.append(
                (
                    input1['dependent_variables'][0]['values'][i]['errors'][j]['label'],
                    pta(input1['dependent_variables'][0]['values'][i]['errors'][j]['symerror'], data_central_value),
                )
            )
        error_array.append(error_bin)
    for i in range(len(cormat1['dependent_variables'][0]['values'])):
        cormat_list.append(float(cormat1['dependent_variables'][0]['values'][i]['value']))
    if len(cormat_list) != len(error_array)**2:
        raise ValueError(
            "The number of elements in the correlation matrix does not match the number of errors in the data."
        )
    
    data_central_ptz_yaml = {'data_central': data_central_ptz}
    kinematics_ptz_yaml = {'bins': kin_ptz}
    uncertainties_ptz_yaml = errorinator(error_array, cormat_list, True)

    with open('data_ptz.yaml', 'w') as f:
        yaml.dump(data_central_ptz_yaml, f, sort_keys=False)
    with open('kinematics_ptz.yaml', 'w') as f:
        yaml.dump(kinematics_ptz_yaml, f, sort_keys=False)
    with open('uncertainties_ptz.yaml', 'w') as f:
        yaml.dump(uncertainties_ptz_yaml, f, sort_keys=False)

    # ptw data

    hepdata_tab3 = "rawdata/data10.yaml"
    with open(hepdata_tab3, 'r') as f:
        input2 = yaml.safe_load(f)
    hepdata_tab4 = "rawdata/data11.yaml"
    with open(hepdata_tab4, 'r') as f:
        cormat2 = yaml.safe_load(f)

    error_array = []
    cormat_list = []

    for i in range(ndata_ptw):
        data_central_value = float(input2['dependent_variables'][0]['values'][i]['value'])
        data_central_ptw.append(data_central_value)
        ptW_min = input2['independent_variables'][0]['values'][i]['low']
        ptW_max = input2['independent_variables'][0]['values'][i]['high']
        kin_value = {
            'ptW': {'min': ptW_min, 'mid': None, 'max': ptW_max},
        }
        kin_ptw.append(kin_value)
        error_bin = []
        for j in range(len(input2['dependent_variables'][0]['values'][i]['errors'])):
            error_bin.append(
                (
                    input2['dependent_variables'][0]['values'][i]['errors'][j]['label'],
                    pta(input2['dependent_variables'][0]['values'][i]['errors'][j]['symerror'], data_central_value),
                )
            )
        error_array.append(error_bin)
    for i in range(len(cormat2['dependent_variables'][0]['values'])):
        cormat_list.append(float(cormat2['dependent_variables'][0]['values'][i]['value']))
    if len(cormat_list) != len(error_array)**2:
        raise ValueError(
            "The number of elements in the correlation matrix does not match the number of errors in the data."
        )
    
    data_central_ptw_yaml = {'data_central': data_central_ptw}
    kinematics_ptw_yaml = {'bins': kin_ptw}
    uncertainties_ptw_yaml = errorinator(error_array, cormat_list, True)

    with open('data_ptw.yaml', 'w') as f:
        yaml.dump(data_central_ptw_yaml, f, sort_keys=False)
    with open('kinematics_ptw.yaml', 'w') as f:
        yaml.dump(kinematics_ptw_yaml, f, sort_keys=False)
    with open('uncertainties_ptw.yaml', 'w') as f:
        yaml.dump(uncertainties_ptw_yaml, f, sort_keys=False)

    # mtwz data

    hepdata_tab5 = "rawdata/data12.yaml"
    with open(hepdata_tab5, 'r') as f:
        input3 = yaml.safe_load(f)
    hepdata_tab6 = "rawdata/data13.yaml"
    with open(hepdata_tab6, 'r') as f:
        cormat3 = yaml.safe_load(f)

    error_array = []
    cormat_list = []

    for i in range(ndata_mtwz):
        data_central_value = float(input3['dependent_variables'][0]['values'][i]['value'])
        data_central_mtwz.append(data_central_value)
        mtwz_min = input3['independent_variables'][0]['values'][i]['low']
        mtwz_max = input3['independent_variables'][0]['values'][i]['high']
        kin_value = {
            'mtwZ': {'min': mtwz_min, 'mid': None, 'max': mtwz_max},
        }
        kin_mtwz.append(kin_value)
        error_bin = []
        for j in range(len(input3['dependent_variables'][0]['values'][i]['errors'])):
            error_bin.append(
                (
                    input3['dependent_variables'][0]['values'][i]['errors'][j]['label'],
                    pta(input3['dependent_variables'][0]['values'][i]['errors'][j]['symerror'], data_central_value),
                )
            )
        error_array.append(error_bin)
    for i in range(len(cormat3['dependent_variables'][0]['values'])):
        cormat_list.append(float(cormat3['dependent_variables'][0]['values'][i]['value']))
    if len(cormat_list) != len(error_array)**2:
        raise ValueError(
            "The number of elements in the correlation matrix does not match the number of errors in the data."
        )
    
    data_central_mtwz_yaml = {'data_central': data_central_mtwz}
    kinematics_mtwz_yaml = {'bins': kin_mtwz}
    uncertainties_mtwz_yaml = errorinator(error_array, cormat_list, True)

    with open('data_mtwz.yaml', 'w') as f:
        yaml.dump(data_central_mtwz_yaml, f, sort_keys=False)
    with open('kinematics_mtwz.yaml', 'w') as f:
        yaml.dump(kinematics_mtwz_yaml, f, sort_keys=False)
    with open('uncertainties_mtwz.yaml', 'w') as f:
        yaml.dump(uncertainties_mtwz_yaml, f, sort_keys=False)

brrr()
        

