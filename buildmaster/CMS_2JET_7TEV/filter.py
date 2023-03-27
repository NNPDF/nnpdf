import yaml
import numpy as np

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

    # read correlation matrix
    corr_matrices = read_dat_file('rawdata/dijet_corr.dat')
    print(len(corr_matrices))

    # construct covariance matrix from correlation matrix
    # TODOO

    return corr_matrices

    

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

def read_dat_file(filename):
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
            stat_corr[i] = np.fromstring(lines[idx], sep=' ')[2:]
            i+=1
        
        correlation_matrices.append(stat_corr)
        
    return correlation_matrices
    


if __name__ == "__main__":

    # save kinematics and central values
    # filter_CMS_2JET_7TRV_data_kinetic()
    
    filterCMS_2JET_7TEV_uncertainties()
    # C = filterCMS_2JET_7TEV_uncertainties()

    # for c in C:
    #     print()
    #     print(f"shape = {c.shape}")
    #     print()
    #     print(c)
        