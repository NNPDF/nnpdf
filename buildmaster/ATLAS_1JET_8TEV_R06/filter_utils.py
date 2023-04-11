import yaml





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
        hepdata_table = f"rawdata/HEPData-ins1604271-v{version}-Table_{table}.yaml"
        
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
        
        hepdata_table = f"rawdata/HEPData-ins1604271-v{version}-Table_{table}.yaml"
        
        with open(hepdata_table, 'r') as file:
            input = yaml.safe_load(file)
        
        # rapidity
        rapidity_interval = table_to_rapidity[table]
        rap = {}
        rap['min'], rap['max'] = rapidity_interval[0], rapidity_interval[1]
        rap['mid'] = 0.5 * (rap['min'] + rap['max'])
        
        # center of mass energy
        sqrts = float(input['dependent_variables'][0]['qualifiers'][1]['value'])
        
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


if __name__ == "__main__":
    cv = get_data_values(tables = [1],version= 1)
    kin = get_kinematics(tables=[1],version=1)
    print(kin)