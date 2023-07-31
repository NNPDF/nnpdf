"""
Filter for CMS_1JET_8TEV

Created on Apr  2023

@author: Mark N. Costantini
"""

import yaml
from filter_utils import (
    get_data_values,
    get_kinematics,
    get_stat_uncertainties,
    block_diagonal_corr,
    correlation_to_covariance,
    uncertainties_df,
)


def filter_CMS_1JET_8TEV_data_kinetic():
    """
    writes kinetic and data central values
    to kinematics.yaml and data.yaml files
    respectively
    """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']

    # get kinematics from hepdata tables
    kin = get_kinematics(tables, version)

    # get central values from hepdata tables
    data_central = get_data_values(tables, version)

    data_central_yaml = {'data_central': data_central}
    kinematics_yaml = {'bins': kin}

    # write central values and kinematics to yaml file
    with open('data.yaml', 'w') as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_CMS_1JET_8TEV_uncertainties():
    """ """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']

    # get dataframe of uncertainties
    df_unc = uncertainties_df(tables)

    # construct block diagonal statistical covariance matrix
    bd_stat_cov = correlation_to_covariance(block_diagonal_corr(tables), get_stat_uncertainties())

    # luminosity
    lum_unc = df_unc['Luminosity']
    print(lum_unc)

    covmat = bd_stat_cov

    return covmat


if __name__ == "__main__":
    # write data central values and kinematics to file
    # filter_CMS_1JET_8TEV_data_kinetic()
    covmat = filter_CMS_1JET_8TEV_uncertainties()

    from validphys.api import API
    from validphys.commondataparser import parse_commondata
    from validphys.loader import Loader
    from validphys.covmats import dataset_inputs_covmat_from_systematics

    l = Loader()
    dat_file = (
        "/Users/markcostantini/codes/nnpdfgit/nnpdf/nnpdfcpp/data/commondata/DATA_CMS_1JET_8TEV.dat"
    )
    systype_file = "/Users/markcostantini/codes/nnpdfgit/nnpdf/nnpdfcpp/data/commondata/systypes/SYSTYPE_CMS_1JET_8TEV_DEFAULT.dat"
    setname = "CMS_1JET_8TEV"

    cd = parse_commondata(dat_file, systype_file, setname)
    dataset_input = l.check_dataset(setname, theoryid=200)

    cmat = dataset_inputs_covmat_from_systematics(
        dataset_inputs_loaded_cd_with_cuts=[cd],
        data_input=[dataset_input],
        use_weights_in_covmat=False,
        norm_threshold=None,
        _list_of_central_values=None,
        _only_additive=False,
    )

    print(covmat / cmat)

    # why does the API give a 185 x 185 shaped covmat
    # i.e. why is it ignoring the low pt stuff

    # inps = [{'dataset': "CMS_1JET_8TEV"}]
    # inp = dict(dataset_inputs=inps, theoryid=200, use_cuts="internal")
    # cmat = API.dataset_inputs_covmat_from_systematics(**inp)
