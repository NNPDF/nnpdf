"""
This file contains the piece of code needed to implement the CMS AFB 
measurement at 13 TeV. Uncertainties are obtained by combining the correlation
matrix and the total uncertainty provided on Hepdata, after which a covariance
matrix is constructed, which is finally decomposed into Ndat artificial
uncertainties
"""

import yaml
from nnpdf_data.filter_utils.utils import cormat_to_covmat,covmat_to_artunc

def get_tables():
    """
    Get the Hepdata tables, given the tables and version specified in metadata 
    """
    prefix = "rawdata/HEPData-ins2818125"
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)
        
    version = metadata["hepdata"]["version"]
    tables = metadata["implemented_observables"][0]["tables"]
    hepdata_tables = []
    
    for table in tables:
        hepdata_tables.append(f"{prefix}-v{version}-{table}.yaml")
    return hepdata_tables

def get_all():
    """
    Returns data, kinematics and uncertainties for dumping in the .yaml files
    """
    data_central = []
    kinematics = []
    uncertainties = []
    art_uncert = []
    art_uncertainties = []
    correlations = []
    hepdata_tables = get_tables()

    table=hepdata_tables[0]
    with open(table, 'r') as f:
        input = yaml.safe_load(f)
    
    # Central values
    data_values = input["dependent_variables"][4]["values"]
    for data_value in data_values:
        data_central.append(data_value["value"])

    # Kinematic bins
    kin_values_yll_min = input["dependent_variables"][0]["values"]
    kin_values_yll_max = input["dependent_variables"][1]["values"]
    kin_values_mll_min = input["dependent_variables"][2]["values"]
    kin_values_mll_max = input["dependent_variables"][3]["values"]
    
    for i in range(len(kin_values_yll_min)):
        kin = {
            'y': {'min': kin_values_yll_min[i]["value"],
                  'mid': 0.5 * (kin_values_yll_min[i]["value"] + kin_values_yll_max[i]["value"]),
                  'max': kin_values_yll_max[i]["value"]},
            'mll': {'min': kin_values_mll_min[i]["value"],
                    'mid': 0.5 * (kin_values_mll_min[i]["value"] + kin_values_mll_max[i]["value"]),
                    'max': kin_values_mll_max[i]["value"]}}
        kinematics.append(kin)
        
    # Artificial uncertainties (from correlation matrix)
    # Errors
    for data_value in data_values:
        errors = data_value["errors"]
        uncertainty = {}
        for error in errors:
            uncertainties.append(error["symerror"])

    # Correlation coefficients
    table_corr = hepdata_tables[1]
    with open(table_corr, 'r') as f:
        input = yaml.safe_load(f)
    corr_coeffs = input["dependent_variables"][0]["values"]
    for corr_coeff in corr_coeffs:
        correlations.append(corr_coeff["value"])

    covmat = cormat_to_covmat(uncertainties,correlations)
    art_uncertainties = covmat_to_artunc(len(uncertainties),covmat,1)
    
    # Remap uncertainties
    for uncertainties in art_uncertainties:
        j = 1
        artificial_unc={}
        for uncertainty in uncertainties:
            key = "artificial uncertainty " + str(j)
            art_unc = {key: uncertainty}
            j = j + 1
            
            artificial_unc.update(art_unc)

        art_uncert.append(artificial_unc)
            
    return(data_central, kinematics, art_uncert)
 
def filter_CMS_Z0_13TEV_PT():
    """
    Dumps data, kinematics, and uncertainties on .yaml files
    """
    central_values, kinematics, uncertainties = get_all()
    # Central values
    data_central_yaml = {"data_central": central_values}
    # Kinematics
    kinematics_yaml = {"bins": kinematics}
    # Uncertainties
    treatment = {"artificial uncertainty": "ADD"}
    correlation = {"artificial uncertainty": "CORR"}
    definitions = {}
    for key,value in uncertainties[0].items():
        definition = {key :
                      {"description": key,
                       "treatment": "ADD",
                       "type": "CORR"}}
        definitions.update(definition)
    uncertainties_yaml = {"definitions": definitions,"bins": uncertainties}

    with open("data_AFB.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)  
    with open("kinematics_AFB.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)
    with open("uncertainties_AFB.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)

if __name__ == "__main__":
    filter_CMS_Z0_13TEV_PT()
