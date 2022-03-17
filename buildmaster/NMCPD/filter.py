#Filter for NMCPD

import yaml

def filter_NMCPD():
    with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

    version = metadata['hepdata']['version']
    tables  = metadata['hepdata']['tables']

    data_central = []
    kin_avg = []
    kin_min = []
    kin_max = []
    stat_error = []
    sys_error = []
    index = 1
    
    for i in tables:
        hepdata_tables="rawdata/HEPData-ins426595-"+version+"-Table_"+str(i)+".yaml"
        with open(hepdata_tables, 'r') as file:
            input = yaml.safe_load(file)
            
        values = input['dependent_variables'][0]['values']
        x = float(input['dependent_variables'][0]['qualifiers'][4]['value'])
        sqrts = float(input['dependent_variables'][0]['qualifiers'][3]['value'])
        
        for j in range(len(values)):
            data_central_value = {'index' : index, 'value' : input['dependent_variables'][0]['values'][j]['value']}      
            data_central.append(data_central_value)
            Q2 = input['independent_variables'][0]['values'][j]['value']
            y = Q2 / ( sqrts * sqrts * x )
            kin_avg_value = {'index' : index , 'x' : x , 'q2' : Q2 , 'y' : y}
            kin_avg.append(kin_avg_value)
            kin_min_value = {'index' : index , 'x' : x , 'q2' : Q2 , 'y' : y}
            kin_min.append(kin_min_value)
            kin_max_value = {'index' : index , 'x' : x , 'q2' : Q2 , 'y' : y}
            kin_max.append(kin_max_value)
            
            stat_error_value = {'index' : index, 'value' : input['dependent_variables'][0]['values'][j]['errors'][0]['symerror']}
            stat_error.append(stat_error_value)
            sys_error_value  = {'index' : index, 'value' : input['dependent_variables'][0]['values'][j]['errors'][1]['symerror']}
            sys_error.append(sys_error_value)
            
            index = index + 1

    data_central_yaml = {'data_central' : data_central}
    kinematics_yaml = {'kin_avg' : kin_avg, 'kin_min' : kin_min, 'kin_max' : kin_max}
    uncertainties_yaml = {'stat' : stat_error, 'sys_1' : {'mode' : ["ADD", "CORR"] , 'errors' : sys_error} }
    
    with open('data.yaml', 'w') as file:
         yaml.dump(data_central_yaml, file)

    with open('kinematics.yaml', 'w') as file:
         yaml.dump(kinematics_yaml, file)

    with open('uncertainties.yaml', 'w') as file:
        yaml.dump(uncertainties_yaml, file)

    with open('uncertainties_dw.yaml', 'w') as file:
        yaml.dump(uncertainties_yaml, file)
         
filter_NMCPD()
