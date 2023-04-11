"""
Utils to be used in CMS_2JET_7TEV.filter.py
"""
import numpy as np
import yaml
import pandas as pd
import logging

log = logging.getLogger(__name__)

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
        
        hepdata_table = f"rawdata/HEPData-ins1208923-v{version}-Table_{table}.yaml"
        
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
        
        hepdata_table = f"rawdata/HEPData-ins1208923-v{version}-Table_{table}.yaml"
        
        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)
        
        rapidity_interval = input['dependent_variables'][0]['qualifiers'][0]['value']
        ydiff = {}
        if rapidity_interval == '< 0.5':
            ydiff['min'], ydiff['max'], ydiff['mid'] = 0.0, 0.5, 0.25
        else:
            ydiff = range_str_to_floats(rapidity_interval)
        
        sqrts = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
            
        dijet_mass_bins = input['independent_variables'][0]['values']
        
        for m12 in dijet_mass_bins:

            # kinematics
            m12['low'], m12['high'] = m12['low'], m12['high']
            m12['mid'] = float(f"{0.5 * (m12['low']+m12['high']):.3f}")

            kin_value = {'ydiff' : {'min': ydiff['min'], 'mid': ydiff['mid'] , 'max': ydiff['max']}, 
                        'm12' : {'min': m12['low'], 'mid': m12['mid'] , 'max': m12['high']} ,
                         'sqrt_s' : {'min': None, 'mid': sqrts , 'max': None}}

            kin.append(kin_value)
    
    return kin


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

def decompose_covmat(covmat):
     """Given a covmat it return an array sys with shape (ndat,ndat)
     giving ndat correlated systematics for each of the ndat point.
     The original covmat is obtained by doing sys@sys.T"""

     lamb, mat = np.linalg.eig(covmat)
     sys = np.multiply(np.sqrt(lamb), mat)
     return sys

def range_str_to_floats(str_range):
    """
    converts a string range to a list,
    e.g. "0.5 - 1.0" --> [0.5,1.0]
    and returns a dictionary
    """
    # Split the range string into two parts
    str_nums = str_range.split('-')
    # Convert the parts to floats
    min = float(str_nums[0])
    max = float(str_nums[1])
    mid = float(f"{0.5 * (min + max):.3f}")
    # Return a dict
    return {"min": min, "mid": mid, "max": max}


def get_corr_dat_file(filename):
    """
    read out correlation matrices from the
    dijet_corr.dat file

    Parameters
    ----------
    filename : str
        Takes path to dijet_corr.dat

    Returns
    -------
    list
        list, each element of which is a 2D np array
    """

    with open(filename) as file:
        lines = file.readlines()
    
    # store the number of the rows where the correlation matrix
    # is printed
    begin_rows = []
    end_rows = []

    for i,line in enumerate(lines):

        if "Statistical correlation" in line and begin_rows == []:
            begin_rows.append(i+2)

        elif "Statistical correlation" in line:
            begin_rows.append(i+2)
            end_rows.append(i-2)

        elif i == len(lines)-1:
            end_rows.append(i)

    correlation_matrices = []
    for begin_row, end_row in zip(begin_rows,end_rows):
        
        size_mat = end_row-begin_row+1
        stat_corr = np.zeros((size_mat,size_mat))

        i = 0
        for idx in range(begin_row,end_row+1):
            # ignore first two columns as these give the bin kin
            stat_corr[i] = np.fromstring(lines[idx], sep=' ')[2:] 
            i+=1
        
        correlation_matrices.append(stat_corr)
        
    return correlation_matrices



def get_stat_uncertainties():
    """
    function used to get the statistical
    uncertainty from the HEPdata tables.

    Returns
    -------
    dict
        dictionary with keys = number of table
        value = list of statistical uncertainties
    
    """
    
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']
    
    stat_err = {}

    for table in tables:
        stat = []
        hepdata_tables = f"rawdata/HEPData-ins1208923-v{version}-Table_{table}.yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)

        for err in input['dependent_variables'][0]['values']:
            stat.append(err['errors'][0]['symerror'])
        stat_err[table] = stat

    return stat_err


def dat_file_to_df():
    """
    from dijet_sys.dat table return a pandas
    DataFrame with index given by Ndata,
    columns by the uncertainties and 
    np.nan entries

    Returns
    -------
    list
        list of dataframes

    """
    
    with open("rawdata/dijet_sys.dat", 'r') as file:
        lines = file.readlines()

    # get rows of numeric tables in dat file
    begin_lines = []
    end_lines = []

    for i,line in enumerate(lines):
        
        if "mlo" in line:
            begin_lines.append(i+2)

        if "The NP correction" in line and len(begin_lines)>=1:
            end_lines.append(i-1)
    
    end_lines.append(len(lines))

    # define dataframe
    columns = ["JEC0-","JEC0+","JEC1-","JEC1+","JEC2-","JEC2+","JEC3-","JEC3+","JEC4-","JEC4+",
                "JEC5-","JEC5+","JEC6-","JEC6+","JEC7-","JEC7+","JEC8-","JEC8+",
                "JEC9-","JEC9+","JEC10-","JEC10+","JEC11-","JEC11+","JEC12-","JEC12+",
                "JEC13-","JEC13+","Lumi-","Lumi+", "Unfolding+","Unfolding-","Bin-by-bin+","Bin-by-bin-"]

    dfs = []
    for beg, end in zip(begin_lines,end_lines):

        df = pd.DataFrame(columns=columns,index=range(beg,end))
        j=0
        for i in range(beg,end):
            # do not consider NP uncertainty
            col_vals = np.fromstring(lines[i], sep=' ')[5:]

            df.iloc[j] = col_vals
            j+=1

        dfs.append(df)
    
    return dfs


def JEC_covmat():
    """
    Jet Energy Scale (JET): 14 Asymmetric uncertainties correlated across all 
    rapidity and mass bins (CORR). This uncertainty is not always presented as 
    (left<0 and right>0), e.g. [-delta_left,+delta_right]
    Hence the D'Agostini prescription for symmetrising errors  
    is not valid here because it works with the only case displayed above.
    Instead, we use here the experimentalists prescription, where we take every 
    subpart of the uncertainty left and right as independent source of 
    uncertainty. This is motivated by taking the average 
    of the left and right uncertainty, hence the origin of the sqrt(2) 
    that we divide by.


    Returns
    -------
    np.array
        covariance matrix for JET energy scale uncertainty

    """
    dfs = dat_file_to_df()
    JEC_err = []
    for df in dfs:
        JEC_err.append(df.filter(like = "JEC"))

    # divide by sqrt(2) since treating each unc of asymm as independent
    jec = pd.concat(JEC_err,axis=0) / np.sqrt(2)
    
    # get central value to convert mult error
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']
    cv = get_data_values(tables,version)
    cv = np.array(cv)
    
    # convert mult error to absolute
    jec = jec.multiply(cv,axis = 0)
    jec = jec.to_numpy().astype(float)
    
    return np.einsum('ij,kj->ik',jec,jec)


def lumi_covmat():
    """
    Luminosity uncertainty: this is a symmetric uncertainty of 2.2% correlated 
    accross all mass and rapidity bins and all CMS datasets at 7 TeV, hence the 
    keyword (CMSLUMI11).

    NOTE: this function is needed to test only whether the full covmat coincides
    with the old implementation and can be removed at some point

    Returns
    -------
    np.array
        covariance matrix for luminosity uncertainty
    
    """
    dfs = dat_file_to_df()
    lumi_err = []
    for df in dfs:
        lumi_err.append(df.filter(like = "Lumi+"))
    
    lumi = pd.concat(lumi_err,axis=0)

    # get central value to convert mult error
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']
    cv = get_data_values(tables,version)
    cv = np.array(cv)

    # convert mult to abs
    lumi = lumi.multiply(cv, axis = 0)
    lumi = lumi.to_numpy().astype(float)
    
    return np.einsum('ij,kj->ik',lumi,lumi)

def unfolding_covmat():
    """
    Unfolding uncertainty: this asymmetric is correlated across all rapidity 
    and mass bins (CORR).


    Returns
    -------
    np.array
        covariance matrix for unfolding uncertainty
    """
    dfs = dat_file_to_df()
    unfold_err = []
    for df in dfs:
        unfold_err.append(df.filter(like = "Unfolding"))
    
    unfold = pd.concat(unfold_err,axis=0) / np.sqrt(2)

    # get central value to convert mult error
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']
    cv = get_data_values(tables,version)
    cv = np.array(cv)

    # convert mult to abs
    unfold = unfold.multiply(cv, axis = 0)
    unfold = unfold.to_numpy().astype(float)
    
    return np.einsum('ij,kj->ik',unfold,unfold)

def bin_by_bin_covmat():
    """
    Bin-by-Bin uncertainty: this is a symmetric uncertainty fully uncorrelated 
    accross bins of mass and rapidity (UNCORR)

    NOTE: this function is needed to test only whether the full covmat coincides
    with the old implementation and can be removed at some point

    Returns
    -------
    np.array
        covariance matrix for bin by bin uncertainty
    """
    dfs = dat_file_to_df()
    bin_err = []
    for df in dfs:
        bin_err.append(df.filter(like = "Bin-by-bin-")) # symm so choose only one
    
    bin = pd.concat(bin_err,axis=0)
    
    # get central value to convert mult error
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']
    cv = get_data_values(tables,version)
    cv = np.array(cv)

    # convert mult to abs
    bin = bin.multiply(cv, axis = 0)
    
    bin = bin.to_numpy().astype(float)

    # fully uncorrelated
    bin_cov = np.diag(bin.reshape(bin.shape[0])**2)
    return bin_cov
