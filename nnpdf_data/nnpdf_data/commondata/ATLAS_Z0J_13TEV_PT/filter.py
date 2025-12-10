"""
This file contains the piece of code needed to implement the ATLAS ZpT 
measurement at 13 TeV. We consider the combined electron-muon measurement, for
which a total uncorrelated and a total correlated ucnertainties are given. The 
measurement is normalised to the fiducial cross section, therefore there is no
luminosity ucnertainty. The first three bins in pT are cropped, because of
the impossiblity of producing theoretical predictions.
"""

import yaml

def get_tables():
    """
    Get the Hepdata table
    """
    hepdata_tables = ["rawdata/HEPData-ins1768911-v3-Table_4a.yaml"]

    return hepdata_tables

def get_all():
    """
    Returns data, kinematics and uncertainties for dumping in the .yaml files
    """
    data_central = []
    kinematics = []
    uncertainties = []

    hepdata_tables = get_tables()
    for table in hepdata_tables:
        with open(table, 'r') as f:
            input = yaml.safe_load(f)
            
        # Central values
        data_values = input["dependent_variables"][0]["values"]
        for data_value in data_values:
            data_central.append(data_value["value"])
        # Kinematic bins
        kin_values = input["independent_variables"][0]["values"]
        for kin_value in kin_values:
            kin = {
                'pT': {'min': kin_value['low'],
                       'mid': 0.5 * (kin_value['low'] + kin_value['high']),
                       'max': kin_value['high']},
                'm_Z2': {'min': None, 'mid': 8317.44, 'max': None},
                'sqrts':  {'min': None, 'mid': 13000, 'max': None}}
               
            kinematics.append(kin)
        # Uncertainties
        i = 0
        for data_value in data_values:
            errors = data_value["errors"]
            uncertainty = {}
            for error in errors:
                uncertainty[error["label"]] = float(error["symerror"].replace('%',''))*data_central[i]/100.
                uncertainty.update(uncertainty)
                
            uncertainties.append(uncertainty)
            i = i+1

    n=3
    return (data_central[n:], kinematics[n:], uncertainties[n:])

def filter_ATLAS_Z0J_13TEV_PT():
    """
    Dumps data, kinematics, and uncertainties on .yaml files
    """
    central_values, kinematics, uncertainties = get_all()
    # Central values
    data_central_yaml = {"data_central": central_values}
    # Kinematics
    kinematics_yaml = {"bins": kinematics}
    # Uncertainties
    treatment = {"correlated uncertainty": "ADD",
                 "uncorrelated uncertainty": "ADD",}
    correlation = {"correlated uncertainty": "CORR",
                 "uncorrelated uncertainty": "UNCORR",}
    definitions = {}
    for key,value in uncertainties[0].items():
        definition = {key :
                      {"description": key,
                       "treatment": treatment[key],
                       "type": correlation[key]}}
        definitions.update(definition)
    uncertainties_yaml = {"definitions": definitions,"bins": uncertainties}

    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)
    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)
    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)
        
if __name__ == "__main__":
    filter_ATLAS_Z0J_13TEV_PT()
