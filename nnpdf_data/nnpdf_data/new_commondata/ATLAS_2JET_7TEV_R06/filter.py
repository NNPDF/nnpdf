# Filter for ATLAS_2JET_7TEV_R06
"""
Created on Mar  2023

@author: Mark N. Costantini
"""

from filter_utils import decompose_covmat, fill_df, range_str_to_floats
import numpy as np
import pandas as pd
from scipy.linalg import block_diag
import yaml

from validphys.covmats import dataset_inputs_covmat_from_systematics
from validphys.loader import Loader


def filter_ATLAS_2JET_7TEV_R06_data_kinetic():
    """
    writes kinetic and data central values
    to yaml file
    """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']

    data_central = []
    kin = []

    for table in tables:
        hepdata_tables = f"rawdata/HEPData-ins1268975-v{version}-Table_{table}.yaml"

        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)

        # half of rapidity separation, index 8 of qualifiers list
        ystar = range_str_to_floats(input['dependent_variables'][0]['qualifiers'][8]['value'])
        # center of mass energy, index 7 of qualifiers list
        sqrts = float(input['dependent_variables'][0]['qualifiers'][7]['value'])

        # measurements of several dijet mass values
        values = input['dependent_variables'][0]['values']

        # loop over different dijet mass bins
        m12_values = input['independent_variables'][0]['values']

        for value, m_jj in zip(values, m12_values):
            # central values
            data_central_value = value['value']
            data_central.append(data_central_value)

            # kinematics
            m_jj['low'], m_jj['high'] = 1e3 * m_jj['low'], 1e3 * m_jj['high']
            m_jj['mid'] = float(f"{0.5 * (m_jj['low']+m_jj['high']):.3f}")

            kin_value = {
                'ystar': {'min': ystar['min'], 'mid': ystar['mid'], 'max': ystar['max']},
                'm_jj': {'min': m_jj['low'], 'mid': m_jj['mid'], 'max': m_jj['high']},
                'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            }

            kin.append(kin_value)

    data_central_yaml = {'data_central': data_central}
    kinematics_yaml = {'bins': kin}

    with open('data.yaml', 'w') as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_ATLAS_2JET_7TEV_R06_uncertainties(scenario='nominal'):
    """
    write uncertainties to yaml file depending on
    the correlation scenario

    Parameters
    ----------

    scenario : string
            can be either nominal, stronger orweaker
            default is nominal

    """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']

    dfs = []
    for table in tables:
        hepdata_tables = f"rawdata/HEPData-ins1268975-v{version}-Table_{table}.yaml"

        # uncertainties dataframe
        df = fill_df(hepdata_tables, scenario=scenario)
        dfs.append(df)

    # Construct Covariance matrix for Systematics
    Asys = pd.concat([df.drop(['lum'], axis=1) for df in dfs], axis=0).to_numpy()
    Csys = np.einsum('ij,kj->ik', Asys, Asys)

    # Construct Special Sys (Lum) Cov matrix
    Alum = pd.concat([df[['lum']] for df in dfs], axis=0).to_numpy()
    Clum = np.einsum('ij,kj->ik', Alum, Alum)

    # construct Block diagonal statistical Covariance matrix
    ndata = [21, 21, 19, 17, 8, 4]
    stat = pd.read_csv(f"rawdata/dijet_statcov/hepcov_R06_Eta0.txt", sep=" ", header=None)
    BD_stat = stat.loc[13:, 1 : ndata[0]].to_numpy().astype(float)

    for i in range(1, 6):
        stat = pd.read_csv(f"rawdata/dijet_statcov/hepcov_R06_Eta{i}.txt", sep=" ", header=None)
        stat = stat.loc[13:, 1 : ndata[i]].to_numpy().astype(float)
        BD_stat = block_diag(BD_stat, stat)

    # covariance matrix without the special systematics, that is, ATLASLUMI11
    covmat_no_lum = BD_stat  # Csys + BD_stat

    # generate artificial systematics
    A_art_sys = decompose_covmat(covmat=covmat_no_lum)

    # error definition
    error_definition = {
        f"art_sys_{i}": {
            "description": f"artificial systematic {i}",
            "treatment": "ADD",
            "type": "CORR",
        }
        for i in range(1, A_art_sys.shape[0] + 1)
    }

    for i in range(1, Asys.shape[1] + 1):
        error_definition[f"sys_{i}"] = {
            "description": f"sys {i}",
            "treatment": "MULT",
            "type": "CORR",
        }

    error_definition["luminosity_uncertainty"] = {
        "description": "luminosity uncertainty",
        "treatment": "MULT",
        "type": "ATLASLUMI11",
    }

    # store error in dict
    error = []
    for n1, n2 in zip(range(A_art_sys.shape[0]), range(Asys.shape[0])):
        error_value = {}
        for m in range(A_art_sys.shape[1]):
            error_value[f"art_sys_{m+1}"] = float(A_art_sys[n1, m])

        for m in range(Asys.shape[1]):
            error_value[f"sys_{m+1}"] = float(Asys[n2, m])

        error_value["luminosity_uncertainty"] = float(Alum[n2])
        error.append(error_value)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    if scenario == 'nominal':
        with open(f"uncertainties.yaml", 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)
    else:
        with open(f"uncertainties_{scenario}.yaml", 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)

    return covmat_no_lum + Clum + Csys


if __name__ == "__main__":
    # # generate kinematics and central data values
    filter_ATLAS_2JET_7TEV_R06_data_kinetic()

    # generate all uncertainties scenarios
    filter_ATLAS_2JET_7TEV_R06_uncertainties(scenario='nominal')
    filter_ATLAS_2JET_7TEV_R06_uncertainties(scenario='stronger')
    filter_ATLAS_2JET_7TEV_R06_uncertainties(scenario='weaker')

    setname = "ATLAS_2JET_7TEV_R06"
    l = Loader()
    cd = l.check_commondata(setname=setname).load_commondata_instance()
    dataset_input = l.check_dataset(setname, theoryid=200)
    from validphys.commondataparser import parse_commondata

    dat_file = '/Users/markcostantini/codes/nnpdfgit/nnpdf/buildmaster/results/DATA_ATLAS_2JET_7TEV_R06.dat'
    sys_file = '/Users/markcostantini/codes/nnpdfgit/nnpdf/buildmaster/results/systypes/SYSTYPE_ATLAS_2JET_7TEV_R06_DEFAULT.dat'

    cd = parse_commondata(dat_file, sys_file, setname)

    cmat = dataset_inputs_covmat_from_systematics(
        dataset_inputs_loaded_cd_with_cuts=[cd],
        data_input=[dataset_input],
        use_weights_in_covmat=False,
        norm_threshold=None,
        _list_of_central_values=None,
        _only_additive=False,
    )

    covmat = filter_ATLAS_2JET_7TEV_R06_uncertainties(scenario='nominal')

    ones = cmat / covmat
    print(ones[0, :])
    print(f"min covmat/cov = {np.min(ones)}")
    print(f"max covmat/cov = {np.max(ones)}")
    print(f"mean covmat / cov = {np.mean(ones)}")
    print(f"np.allclose: {np.allclose(covmat,cmat, rtol = 1e-5, atol =1e-5)}")
