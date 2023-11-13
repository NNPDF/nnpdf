"""
Utils to be used in CMS_1JET_8TEV/filter.py
"""

import yaml
import numpy as np
import logging
from scipy.linalg import block_diag
import pandas as pd


log = logging.getLogger(__name__)

TABLE_TO_RAPIDITY = {
    1: [0.0, 0.5],
    2: [0.5, 1.0],
    3: [1.0, 1.5],
    4: [1.5, 2],
    5: [2.0, 2.5],
    6: [2.5, 3.0],
}

COLUMN_NAMES = [
    'binFlag',
    'ylow',
    'yhigh',
    'ptlow',
    'pthigh',
    'Sigma',
    'NPCorr',
    'npcorerr+',
    'npcorerr-',
    'ignore',
    'Unfolding+',
    'Unfolding-',
    'AbsoluteStat+',
    'AbsoluteStat-',
    'AbsoluteScale+',
    'AbsoluteScale-',
    'AbsoluteMPFBias+',
    'AbsoluteMPFBias-',
    'Fragmentation+',
    'Fragmentation-',
    'SinglePionECAL+',
    'SinglePionECAL-',
    'SinglePionHCAL+',
    'SinglePionHCAL-',
    'FlavorQCD+',
    'FlavorQCD-',
    'RelativeJEREC1+',
    'RelativeJEREC1-',
    'RelativeJEREC2+',
    'RelativeJEREC2-',
    'RelativeJERHF+',
    'RelativeJERHF-',
    'RelativePtBB+',
    'RelativePtBB-',
    'RelativePtEC1+',
    'RelativePtEC1-',
    'RelativePtEC2+',
    'RelativePtEC2-',
    'RelativePtHF+',
    'RelativePtHF-',
    'RelativeFSR+',
    'RelativeFSR-',
    'RelativeStatEC2+',
    'RelativeStatEC2-',
    'RelativeStatHF+',
    'RelativeStatHF-',
    'PileUpDataMC+',
    'PileUpDataMC-',
    'PileUpPtRef+',
    'PileUpPtRef-',
    'PileUpPtBB+',
    'PileUpPtBB-',
    'PileUpPtEC1+',
    'PileUpPtEC1-',
    'PileUpPtEC2+',
    'PileUpPtEC2-',
    'PileUpPtHF+',
    'PileUpPtHF-',
    'RelativeStatFSR+',
    'RelativeStatFSR-',
    'Luminosity',
    'stat',
    'uncor',
]

TABLE_DATA_SHAPE = {1: 37, 2: 37, 3: 36, 4: 32, 5: 25, 6: 18}


def get_data_values(tables, version):
    """
    returns the central data

    Parameters
    ----------
    tables : list
            list that enumerates the table number

    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    Returns
    -------
    list
        list containing the central values for all
        hepdata tables

    """

    data_central = []
    for table in tables:
        hepdata_table = f"rawdata/HEPData-ins1487277-v{version}-Table_{table}.yaml"

        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)

        values = input['dependent_variables'][0]['values']
        pt_values = input['independent_variables'][0]['values']
        
        for pt, value in zip(pt_values, values):
            data_central.append(value['value'])
            

    return data_central


def get_kinematics(tables, version):
    """
    returns the relevant kinematics values

    Parameters
    ----------
    tables : list
            list that enumerates the table number

    version : int
            integer read from metadata.yaml that
            indicated the version of the hepdata
            tables

    Returns
    -------
    list
        list containing the kinematic values for all
        hepdata tables
    """
    kin = []

    for table in tables:
        hepdata_table = f"rawdata/HEPData-ins1487277-v{version}-Table_{table}.yaml"

        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)

        # rapidity
        rapidity_interval = TABLE_TO_RAPIDITY[table]
        rap = {}
        rap['min'], rap['max'] = rapidity_interval[0], rapidity_interval[1]
        rap['mid'] = 0.5 * (rap['min'] + rap['max'])

        # center of mass energy
        sqrts = float(input['dependent_variables'][0]['qualifiers'][2]['value'])

        # transverse momentum
        jet_kt_bins = input['independent_variables'][0]['values']
        KT = {}
        for kt in jet_kt_bins:
        
            KT['min'], KT['max'] = kt['low'], kt['high']
            KT['mid'] = float(f"{0.5 * (kt['low'] + kt['high']):.3f}")

            kin_value = {
                'y': {'min': rap['min'], 'mid': rap['mid'], 'max': rap['max']},
                'pT': {'min': KT['min'], 'mid': KT['mid'], 'max': KT['max']},
                'sqrt_s': {'min': None, 'mid': sqrts, 'max': None},
            }

            kin.append(kin_value)


    return kin
    

def get_stat_correlations(table):
    """
    Read a CMS_8TeV_jets_Ybin#___CMS_8TeV_jets_Ybin#.dat
    file and output the statistical correlation matrix.


    Parameters
    ----------
    table : int
            number of the table

    Returns
    -------
    np.array
            2-D array

    """
    with open(f'rawdata/CMS_8TeV_jets_Ybin{table}___CMS_8TeV_jets_Ybin{table}.dat', 'r') as file:
        card = file.readlines()
    
    # get shape of matrix
    shape_mat = TABLE_DATA_SHAPE[table] 
    
    stat_corr = np.zeros((shape_mat, shape_mat))
    
    # correlation rows always start at row 18
    for j in range(shape_mat):
        # fill rows of correlation matrix
        stat_corr[j, :] = np.array(
            [card[(17 + shape_mat * j) + k].split()[-1] for k in range(shape_mat)]
        )
    
    # add zeros for points in the pt<74 kinematic region 
    # these points should be cut (it is always 9 pt bins in the pt < 74 region)
    stat_corr = block_diag(np.zeros((9,9)), stat_corr)
    return stat_corr


def block_diagonal_corr(tables):
    """
    forms block diagonal correlation matrix
    for stat uncertainties. Each block corresponds
    to a rapidity bin.
    
    Parameters
    ----------
    tables : list
            list of integers numbering the tables

    Returns
    -------
    np.array
        block diagonal matrix of dim ndata x ndata
    """
    bd_corr = get_stat_correlations(tables[0])

    for table in tables[1:]:
        
        bd_corr = block_diag(bd_corr, get_stat_correlations(table))

    return bd_corr


def correlation_to_covariance(correlation, uncertainties):
    """
    Converts a correlation matrix into a covariance matrix
    using a list of uncertainties.

    Parameters:
    -----------
    correlation : np.ndarray
        A square matrix of correlations.
    uncertainties : np.ndarray
        A 1D array of uncertainties.

    Returns:
    --------
    np.ndarray
        The corresponding covariance matrix.
    """
    covariance = np.outer(uncertainties, uncertainties) * correlation
    return covariance


def get_stat_uncertainties():
    """
    function used to get the statistical
    uncertainty from the HEPdata tables.

    Returns
    -------
    np.array
        array with ordered stat errors for all tables

    """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables = metadata['hepdata']['tables']

    stat_err = []

    for table in tables:
        hepdata_tables = f"rawdata/HEPData-ins1487277-v{version}-Table_{table}.yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
        
        # discard pT < 74 GeV entries

        for err, pt in zip(input['dependent_variables'][0]['values'], input['independent_variables'][0]['values']):
            stat_err.append(err['errors'][0]['symerror'])
    
    return np.array(stat_err)


def get_uncertainties_df(table):
    """ 
    """

    # read dat file into dataframe by skipping the first 41 metadata rows
    df = pd.read_csv(
        f"rawdata/CMS_8TeV_jets_Ybin{table}.dat", sep="\s+", skiprows=41, names=COLUMN_NAMES
    )

    # reindex
    df = df.reset_index(drop=True)

    df = df[1:-1] # discard last one as it is repeated
    
    return df


def uncertainties_df(tables):
    """ """
    dfs = []

    for table in tables:
        dfs.append(get_uncertainties_df(table))
    df = pd.concat(dfs, axis=0)
    return df

def process_err(df):
    """
    Given the uncertainties dataframe, if the two variations in the pair 
    (of uncertainties) have the same sign, only the largest (in absolute value) 
    is retained, while the other is set to zero
    
    """
    for col_idx in np.arange(0,len(df.columns),2):

        for row_idx, (val1, val2) in enumerate(zip(df.iloc[:,col_idx], df.iloc[:, col_idx+1])):
            if np.sign(val1) ==  np.sign(val2):
    
                if np.abs(val1) > np.abs(val2):
                    df.iloc[row_idx, col_idx+1] = 0
                elif np.abs(val1) < np.abs(val2):
                    df.iloc[row_idx, col_idx] = 0
    return df

def decompose_covmat(covmat):
     """Given a covmat it return an array sys with shape (ndat,ndat)
     giving ndat correlated systematics for each of the ndat point.
     The original covmat is obtained by doing sys@sys.T"""

     lamb, mat = np.linalg.eig(covmat)
     sys = np.multiply(np.sqrt(lamb), mat)
     return sys

if __name__ == "__main__":
    # print(get_kinematics(tables=[1],version=1))

    # print(get_stat_correlations())
    # bd = block_diagonal_corr(tables=[1,2,3,4,5,6])
    # stat = get_stat_uncertainties()
    # print(stat)
    df = uncertainties_df([1, 2, 3, 4, 5, 6])
    print(df)
    print(df.columns)
    print(df.Sigma)
