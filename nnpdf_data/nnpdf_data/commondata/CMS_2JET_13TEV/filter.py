"""

"""

import yaml
import numpy as np

def get_tables(obs=None):
    """
    Get the Hepdata tables, given the tables and version specified in metadata 
    """
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    prefix = "rawdata/"
    if obs==None:
        print("Observable not available. Please choose between on of the two:")
        print("2d_mass_ak8")
        print("3d_madd_ak8")
        exit(1)
    else:
        table_cv = prefix + "crosssection_" + obs + ".yaml"
        table_er = prefix + "statcorr_" + obs + ".yaml"

    return [table_cv, table_er]

def get_all(obs=None):
    """
    Returns data, kinematics and uncertainties for dumping in the .yaml files
    """
    data_central = []
    kinematics = []
    uncertainties = []
    
    table_cv, table_er = get_tables(obs)
    
    with open(table_cv, 'r') as f:
        input = yaml.safe_load(f)
        
    # Central values
    data_values = input["dependent_variables"][0]["values"]
    for data_value in data_values:
        data_central.append(data_value["value"])
        
    # Kinematic bins
    kin_1_values = input["independent_variables"][0]["values"]
    kin_2_values = input["independent_variables"][1]["values"]
    if obs == "3d_mass_ak8":
        kin_3_values = input["independent_variables"][2]["values"]
        
    for i in range(len(kin_1_values)):
        if obs == "2d_mass_ak8":
            kin = {
                'ymax': {'min': kin_2_values[i]['low'],
                         'mid': 0.5 * (kin_2_values[i]['low'] + kin_2_values[i]['high']),
                         'max': kin_2_values[i]['high']},
                'm_jj': {'min': kin_1_values[i]['low'],
                         'mid': 0.5 * (kin_1_values[i]['low'] + kin_1_values[i]['high']),
                         'max': kin_1_values[i]['high']},
                'sqrts': {'min': None, 'mid': 13000.0, 'max': None},}
        elif obs == "3d_mass_ak8":
            kin = {
                'yb': {'min': kin_2_values[i]['low'],
                       'mid': 0.5 * (kin_2_values[i]['low'] + kin_2_values[i]['high']),
                       'max': kin_2_values[i]['high']},
                'y_star': {'min': kin_3_values[i]['low'],
                           'mid': 0.5 * (kin_3_values[i]['low'] + kin_3_values[i]['high']),
                           'max': kin_3_values[i]['high']},
                'm_jj': {'min': kin_1_values[i]['low'],
                         'mid': 0.5 * (kin_1_values[i]['low'] + kin_1_values[i]['high']),
                         'max': kin_1_values[i]['high']},}
        kinematics.append(kin)

    # Uncertainties
    shifts = []
    for data_value in data_values:
        errors = data_value["errors"]
        uncertainty = {}
        shift = 0
        for error in errors:
            unc = {}
            if "asymerror" in error.keys():
                delta_plus = error["asymerror"]["plus"]
                delta_minus = error["asymerror"]["minus"]
                if delta_plus > 0 and delta_minus > 0:
                    delta_minus = 0
                elif delta_plus < 0 and delta_minus < 0:
                    delta_plus = 0
                delta = ( delta_plus + delta_minus ) / 2.
                Deltabar = ( delta_plus - delta_minus ) / 2.
                shift = shift + delta
                # Symmetrise asymmetric uncertainties
                unc[error["label"]] = Deltabar * (1. + delta*delta / (Deltabar * Deltabar))

            elif "symerror" in error.keys():
                unc[error["label"]] = error["symerror"]

            uncertainty.update(unc)
        
        uncertainties.append(uncertainty)
        shifts.append(shift)

    # Shift the data due to asymmetric uncertainties
    data_central_shifted = [sum(x) for x in zip(data_central, shifts)]
    return (data_central_shifted, kinematics, uncertainties)

def filter_CMS_2JET_13TEV():
    """
    Dumps data, kinematics, and uncertainties on .yaml files
    """
    observables = ["2d_mass_ak8","3d_mass_ak8"]
    for observable in observables:
        data_file = "data_" + observable + ".yaml"
        kin_file = "kinematics_" + observable + ".yaml"
        unc_file = "uncertainties_" + observable + ".yaml"
        central_values, kinematics, uncertainties = get_all(observable)
        # Central values
        data_central_yaml = {"data_central": central_values}
        # Kinematics
        kinematics_yaml = {"bins": kinematics}
        # Uncertainties
        definitions = {}
        correlation = ""
        treatment = ""
        for key in uncertainties[0].keys():
            if key=="Stat":
                correlation = "UNCORR"
                treatment = "ADD"
            elif key=="Lumi":
                correlation = "CMSLUMI16"
                treatment = "MULT"
            else:
                correlation = "CORR"
                treatment = "MULT"   
            
            definition = {key:
                          {
                              "description": key + " unc. from Hepdata",
                              "treatment": treatment,
                              "type": correlation,}}
            definitions.update(definition)

        uncertainties_yaml = {"definitions": definitions, "bins": uncertainties}
        with open(data_file, "w") as file:
            yaml.dump(data_central_yaml, file, sort_keys=False)
        with open(kin_file, "w") as file:
            yaml.dump(kinematics_yaml, file, sort_keys=False)
        with open(unc_file, "w") as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)

if __name__ == "__main__":
    filter_CMS_2JET_13TEV()
