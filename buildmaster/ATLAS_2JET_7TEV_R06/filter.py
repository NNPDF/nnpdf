# Filter for ATLAS_2JET_7TEV_R06
"""
Created on Mar  2023

@author: Mark N. Costantini
"""

import yaml

def range_str_to_floats(str_range):
    """
    converts a string range to a list,
    e.g. "0.5 - 1.0" --> [0.5,1.0]
    """
    # Split the range string into two parts
    str_nums = str_range.split('-')
    # Convert the parts to floats
    min = float(str_nums[0])
    max = float(str_nums[1])
    mid = float(f"{0.5 * (min + max):.3f}")
    # Return a dict
    return {"min": min, "mid": mid, "max": max}

from_TeV_to_GeV = lambda x: 1e3 * x


def filter_ATLAS_2JET_7TEV_R06():
    """
    
    """

    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']

    data_central = []
    kin = []
    
    for table in tables:
        hepdata_tables = f"rawdata/HEPData-ins1268975-v{version}-Table_{table}.yaml"
        
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
        
        ystar = range_str_to_floats(input['dependent_variables'][0]['qualifiers'][8]['value'])
        sqrts = float(input['dependent_variables'][0]['qualifiers'][7]['value'])

        # measurements of several dijet mass values
        values = input['dependent_variables'][0]['values']

        # loop over different diijet mass bins
        for i,value in enumerate(values):

            data_central_value = value['value']
            data_central.append(data_central_value)
              
            m12 = input['independent_variables'][0]['values'][i]
            m12['low'], m12['high'] = 1e3 * m12['low'], 1e3 * m12['high']
            m12['mid'] = float(f"{0.5 * (m12['low']+m12['high']):.3f}")

            kin_value = {'ystar' : {'min': ystar['min'], 'mid': ystar['mid'] , 'max': ystar['max']}, 
                        'm12' : {'min': m12['low'], 'mid': m12['mid'] , 'max': m12['high']} ,
                         'sqrt_s' : {'min': None, 'mid': sqrts , 'max': None}}

            kin.append(kin_value)
            
    #         error_value = {'stat_1': input['dependent_variables'][0]['values'][j]['errors'][0]['symerror'], 'syst_1': input['dependent_variables'][0]['values'][j]['errors'][1]['symerror']}
    #         error.append(error_value)

    #         error_value_nuc = {'nuclear': input['dependent_variables'][0]['values'][j]['errors'][0]['symerror']}
    #         error_nuc.append(error_value_nuc)

    # error_definition = {'stat_1': {'description': "total statistical uncertainty", 'treatment': "ADD", 'type': "UNCORR"},
    #                     'syst_1': {'description': "total systematic uncertainty",  'treatment': "ADD", 'type': "CORR"}}

    # error_definition_dw = {'nuclear': {'description': "nuclear uncertainty (deweighted)", 'treatment': "ADD", 'type': "NUC_DW"}}
    # error_definition_sh = {'nuclear': {'description': "nuclear uncertainty (shifted)", 'treatment': "ADD", 'type': "NUC_SH"}}
    
    data_central_yaml  = { 'data_central' : data_central }
    kinematics_yaml    = { 'bins' : kin }
    # uncertainties_yaml = { 'definition': error_definition, 'bins' : error }
    # uncertainties_dw_yaml = { 'definition': error_definition_dw, 'bins' : error_nuc }
    # uncertainties_sh_yaml = { 'definition': error_definition_sh, 'bins' : error_nuc }
    
    with open('data.yaml', 'w') as file:
         yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
         yaml.dump(kinematics_yaml, file, sort_keys=False)

    # with open('uncertainties.yaml', 'w') as file:
    #     yaml.dump(uncertainties_yaml, file, sort_keys=False)

    # with open('uncertainties_dw.yaml', 'w') as file:
    #     yaml.dump(uncertainties_dw_yaml, file, sort_keys=False)

    # with open('uncertainties_sh.yaml', 'w') as file:
    #     yaml.dump(uncertainties_sh_yaml, file, sort_keys=False)    

if __name__ == "__main__":
    filter_ATLAS_2JET_7TEV_R06()