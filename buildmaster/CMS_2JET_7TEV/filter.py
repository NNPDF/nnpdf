import yaml


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



if __name__ == "__main__":

    # save kinematics and central values
    filter_CMS_2JET_7TRV_data_kinetic()