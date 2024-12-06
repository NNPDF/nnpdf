import numpy as np
import yaml


def lumiless_covmat(data, covmat, lumi_uncert_percent):
    with open(data, 'r') as f:
        data = yaml.safe_load(f)

    with open(covmat, 'r') as f:
        covmat = yaml.safe_load(f)

    data_array = np.array([])
    for i in range(len(data['dependent_variables'][0]['values'])):
        data_array = np.append(data_array, data['dependent_variables'][0]['values'][i]['value'])

    covmat_array = np.array([])
    for i in range(len(covmat['dependent_variables'][0]['values'])):
        covmat_array = np.append(
            covmat_array, covmat['dependent_variables'][0]['values'][i]['value']
        )

    lumi_uncert_array = np.array([data_array * lumi_uncert_percent / 100])
    lumi_covmat = lumi_uncert_array.T @ lumi_uncert_array
    lumi_covmat_array = lumi_covmat.flatten()
    lumiless_covmat = (covmat_array - lumi_covmat_array).tolist()

    return lumiless_covmat


llcm_mtt = lumiless_covmat('rawdata/Table618.yaml', 'rawdata/Table619.yaml', 2.1)
llcm_ptt = lumiless_covmat('rawdata/Table610.yaml', 'rawdata/Table611.yaml', 2.1)
llcm_yt = lumiless_covmat('rawdata/Table614.yaml', 'rawdata/Table615.yaml', 2.1)
llcm_ytt = lumiless_covmat('rawdata/Table626.yaml', 'rawdata/Table627.yaml', 2.1)
