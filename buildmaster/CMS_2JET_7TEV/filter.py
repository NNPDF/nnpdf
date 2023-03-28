import yaml
import numpy as np
from scipy.linalg import block_diag

def filter_CMS_2JET_7TRV_data_kinetic():
    """
    writes kinetic and data central values
    to yaml file
    """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']

    data_central = []
    kin = []

    for table in tables:
        
        hepdata_tables = f"rawdata/HEPData-ins1208923-v{version}-Table_{table}.yaml"
        
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
        
        rapidity_interval = input['dependent_variables'][0]['qualifiers'][0]['value']
        ydiff = {}
        if rapidity_interval == '< 0.5':
            ydiff['min'], ydiff['max'], ydiff['mid'] = 0.0, 0.5, 0.25
        else:
            ydiff = range_str_to_floats(rapidity_interval)
        
        sqrts = float(input['dependent_variables'][0]['qualifiers'][2]['value'])
            
        dijet_mass_bins = input['independent_variables'][0]['values']

        values = input['dependent_variables'][0]['values']
        
        for value, m12 in zip(values, dijet_mass_bins):

            # central values
            data_central_value = value['value']
            data_central.append(data_central_value)

            # kinematics
            m12['low'], m12['high'] = m12['low'], m12['high']
            m12['mid'] = float(f"{0.5 * (m12['low']+m12['high']):.3f}")

            kin_value = {'ydiff' : {'min': ydiff['min'], 'mid': ydiff['mid'] , 'max': ydiff['max']}, 
                        'm12' : {'min': m12['low'], 'mid': m12['mid'] , 'max': m12['high']} ,
                         'sqrt_s' : {'min': None, 'mid': sqrts , 'max': None}}

            kin.append(kin_value)
    
    data_central_yaml  = { 'data_central' : data_central }
    kinematics_yaml    = { 'bins' : kin }

    with open('data.yaml', 'w') as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)




def filterCMS_2JET_7TEV_uncertainties():
    """
    
    """

    # generate block diagonal statistical covariance matrix
    # Statistical uncertainty correlated between the mass bins in 
    # the same rapidity range

    # get correlation matrix for statistical uncertainties
    corr_matrices = get_corr_dat_file('rawdata/dijet_corr.dat')
    # get statistical uncertainties from each HEPData table
    stat_uncertainties = get_stat_uncertainties()


    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']

    stat_cov_mats = []
    
    for i,table in enumerate(tables):
        if corr_matrices[i].shape[0] != np.array(stat_uncertainties[table]).shape[0]:
            raise("Shapes of correlation matrix and uncertainties array are not compatible")
        # convert correlation matrices to covariance matrices
        stat_cov_mats.append(correlation_to_covariance(corr_matrices[i],stat_uncertainties[table]))

    BD_stat = stat_cov_mats[0]

    for i in range(1,len(stat_cov_mats)):
        stat = stat_cov_mats[i]
        BD_stat = block_diag(BD_stat,stat)



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
    


if __name__ == "__main__":

    # save kinematics and central values
    # filter_CMS_2JET_7TRV_data_kinetic()
    
    filterCMS_2JET_7TEV_uncertainties()

    # print(get_stat_uncertainties())
    

    # C = filterCMS_2JET_7TEV_uncertainties()

    # for c in C:
    #     print()
    #     print(f"shape = {c.shape}")
    #     print()
    #     print(c)
        