"""
Utils to be used in CMS_1JET_8TEV/filter.py
"""

import yaml
import numpy as np
import logging
from scipy.linalg import block_diag
import pandas as pd


log = logging.getLogger(__name__)

table_to_rapidity = {1: [0.,0.5], 2: [0.5,1.], 3:[1.,1.5], 4:[1.5,2], 5:[2.,2.5], 6:[2.5,3.0]}

def get_data_values(tables,version):
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

        for value in values:
            data_central.append(value['value'])

    return data_central


def get_kinematics(tables,version):
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
        rapidity_interval = table_to_rapidity[table]
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
            

            kin_value = {'y' : {'min': rap['min'], 'mid': rap['mid'] , 'max': rap['max']}, 
                        'kt' : {'min': KT['min'], 'mid': KT['mid'] , 'max': KT['max']} ,
                         'sqrt_s' : {'min': None, 'mid': sqrts , 'max': None}}

            kin.append(kin_value)
    
    return kin

def get_shape(correlation_file):
    """
    returns the shape of the statistical correlation
    matrix.

    Parameters
    ----------
    correlation_file : list
                    list whose entries are the rows of the
                    correlation file
    
    Returns
    -------
    int
        integer giving the shape of the corr matrix
    
    """
    shape_list = []
    for i, line in enumerate(correlation_file):
        # 74 is the lowest pt and is common to all .dat files
        if len(line.split()) >= 5 and line.split()[-2] == '74':
            shape_list.append(i)
            if len(shape_list) == 2:
                return shape_list[1]-shape_list[0]
        

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
    with open(f'rawdata/CMS_8TeV_jets_Ybin{table}___CMS_8TeV_jets_Ybin{table}.dat','r') as file:
        card = file.readlines()
    
    # get shape of matrix 
    shape_mat = get_shape(card)

    stat_corr =  np.zeros((shape_mat,shape_mat))

    # correlation rows always start at row 18
    for j in range(shape_mat):
        # fill rows of correlation matrix
        stat_corr[j,:] = np.array([card[(17 + shape_mat * j)+ k].split()[-1] for k in range(shape_mat)])

    return stat_corr

def block_diagonal_corr(tables):
    """
    forms block diagonal correlation matrix
    for stat uncertainties. Each block corresponds
    to a rapidity bin.
    Null correlations are associated to the points 
    below the nominal pT cut (pT<74 GeV) for which the
    information on correlations is not available.
    
    Parameters
    ----------
    tables : list
            list of integers numbering the tables
        
    Returns
    -------
    np.array
        block diagonal matrix of dim ndata x ndata
    """

    # null correlations for points below nominal pT cut 
    null_correlations = np.zeros((9,9))

    bd_corr = block_diag(null_correlations,get_stat_correlations(tables[0]))
    
    for table in tables[1:]:
        tmp_corr = block_diag(null_correlations,get_stat_correlations(table))
        bd_corr = block_diag(bd_corr,tmp_corr)

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
    tables  = metadata['hepdata']['tables']
    
    stat_err = []

    for table in tables:
        
        hepdata_tables = f"rawdata/HEPData-ins1487277-v{version}-Table_{table}.yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)

        for err in input['dependent_variables'][0]['values']:
            stat_err.append(err['errors'][0]['symerror'])
        
    return np.array(stat_err)

def get_uncertainties_df(table):
    """
    
    """

    with open(f"rawdata/CMS_8TeV_jets_Ybin{table}.dat",'r') as file:
        card = file.readlines()

    for line in card:
        if "ColumnName" in line:
            columns = line.split()
            columns = [col.strip(", '") for col in columns[2:]]
                
                

            
            

    # ignore ColumnName, =, and binFlag
    
    # columns = columns[3:]
    
    # read dat file into dataframe by skipping the first 41 metadata rows
    df = pd.read_csv(f"rawdata/CMS_8TeV_jets_Ybin{table}.dat", sep="\s+", skiprows = 41, names = columns)

    # reindex
    df = df.reset_index(drop=True)

    return df

def uncertainties_df(tables):
    """
    """
    dfs = []
    
    for table in tables:
        dfs.append(get_uncertainties_df(table))
    df = pd.concat(dfs,axis = 0)
    return df


if __name__ == "__main__":
    # print(get_kinematics(tables=[1],version=1))

    # print(get_stat_correlations())
    # bd = block_diagonal_corr(tables=[1,2,3,4,5,6])
    # stat = get_stat_uncertainties()
    # print(stat)
    df = uncertainties_df([1,2,3,4,5,6])
    print(df)
    print(df.columns)
    print(df.Sigma)
    
