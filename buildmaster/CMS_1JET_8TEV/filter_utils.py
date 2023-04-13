"""
Utils to be used in CMS_1JET_8TEV/filter.py
"""

import yaml
import numpy as np
import logging


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
        

def get_stat_correlations():
    """
    
    Returns
    -------
    np.array
    
    """
    with open('rawdata/CMS_8TeV_jets_Ybin6___CMS_8TeV_jets_Ybin6.dat','r') as file:
        card = file.readlines()
    
    # get shape of matrix
    shape_mat = get_shape(card)
    print(shape_mat)

    stat_corr =  np.zeros((shape_mat,shape_mat))

    # correlation rows always start at row 18
    for j in range(shape_mat):
        # fill rows of correlation matrix
        stat_corr[j,:] = np.array([card[(17 + shape_mat * j)+ k].split()[-1] for k in range(shape_mat)])


    return stat_corr

            

if __name__ == "__main__":
    # print(get_kinematics(tables=[1],version=1))

    print(get_stat_correlations())