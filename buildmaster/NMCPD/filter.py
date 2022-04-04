#Filter for NMCPD

import sys
import yaml

def filter_NMCPD():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']

    data_central = []
    kin = []
    error = []
    
    for i in tables:
        hepdata_tables="rawdata/HEPData-ins426595-"+version+"-Table_"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
            
        values = input['dependent_variables'][0]['values']
        x = float(input['dependent_variables'][0]['qualifiers'][4]['value'])
        sqrts = float(input['dependent_variables'][0]['qualifiers'][3]['value'])
        
        for j in range(len(values)):
            
            data_central_value = input['dependent_variables'][0]['values'][j]['value']    
            data_central.append(data_central_value)
            Q2 = input['independent_variables'][0]['values'][j]['value']
            y = Q2 / ( sqrts * sqrts * x )
            kin_value = {'x' : {'min': None, 'mid': x , 'max': None}, 'q2' : {'min': None, 'mid': Q2 , 'max': None} , 'y' : {'min': None, 'mid': y , 'max': None}}
            kin.append(kin_value)
            
            error_value = {'stat_1': input['dependent_variables'][0]['values'][j]['errors'][0]['symerror'], 'syst_1': input['dependent_variables'][0]['values'][j]['errors'][1]['symerror']}
            error.append(error_value)

    error_definition = {'stat_1': {'description': "total statistical uncertainty", 'treatment': "ADD", 'type': "UNCORR"},
                        'syst_1': {'description': "total systematic uncertainty",  'treatment': "ADD", 'type': "CORR"}}
    
    data_central_yaml  = { 'data_central' : data_central }
    kinematics_yaml    = { 'bins' : kin }
    uncertainties_yaml = { 'definition': error_definition, 'bins' : error }
    
    with open('data.yaml', 'w') as file:
         yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
         yaml.dump(kinematics_yaml, file, sort_keys=False)

    with open('uncertainties.yaml', 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)
         
filter_NMCPD()
