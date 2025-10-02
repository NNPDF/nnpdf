"""
This file contains the piece of code needed to implement the CMS ZpT measurement
at 13 TeV. Systematic uncertainties are implemented starting from the breakdown
available on HepData. The correlation treatment follows the approach mentioned
in the paper (see the lines immediately before Sect. 6): "The systematic
and statistical uncertainties are obtained using the linear combination method 
described in Ref. [77], considering as fully correlated the uncertainties in 
the jet energy scale and resolution, the pileup, the background subtraction, 
b tagging, and the integrated luminosity. Other uncertainties are considered 
as uncorrelated." Note that correlations are kept not only across different 
pT(ll) bins, but also across different m(ll) bins. Covariance matrices, albeit 
only for covariances across pT(ll) bins, are available on HepData. These are
disregarded, because they do not include any correlations in m(ll). Their 
inspection confirms that indeed the aforementioned uncertainties are fully
correlated. The other uncertainties are weakly correlated or anti-correlated. 
"""

import yaml, re

def get_tables():
    """
    Get the Hepdata tables, given the tables and version specified in metadata 
    """
    prefix = "rawdata/HEPData-ins2079374"
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)
        
    version = metadata["hepdata"]["version"]
    tables = metadata["implemented_observables"][0]["tables"]
    hepdata_tables = []
    
    for table in tables:
        hepdata_tables.append(f"{prefix}-v{version}-pT_ll_mass_{table}.yaml")

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
        ranges = list(map(int, re.search(r"mass_(\d+)-(\d+)", table).groups()))
        
        min_mll2 = ranges[0] ** 2
        mid_mll2 = (0.5 * (ranges[0] + ranges[1])) ** 2
        max_mll2 = ranges[1] ** 2
            
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
                'm_ll2': {'min': min_mll2, 'mid': mid_mll2, 'max': max_mll2},}

            kinematics.append(kin)
        # Uncertainties
        for data_value in data_values:
            errors = data_value["errors"]
            uncertainty = {}
            for error in errors:
                uncertainty[error["label"]] = error["symerror"]
                uncertainty.update(uncertainty)

            uncertainties.append(uncertainty)

    return (data_central, kinematics, uncertainties)

def filter_CMS_Z0J_13TEV_PT():
    """
    Dumps data, kinematics, and uncertainties on .yaml files
    """
    central_values, kinematics, uncertainties = get_all()
    # Central values
    data_central_yaml = {"data_central": central_values}
    # Kinematics
    kinematics_yaml = {"bins": kinematics}
    # Uncertainties
    treatment = {"Data stat.": "ADD",
                 "Unfolding stat.": "ADD",
                 "Unfolding model": "ADD",
                 "Int. luminosity": "MULT",
                 "Lepton energy": "ADD",
                 "Efficiency": "ADD",
                 "Backgrounds": "MULT",
                 "Jet energy": "MULT",
                 "Others": "MULT"}
    correlation = {"Data stat.": "UNCORR",
                   "Unfolding stat.": "UNCORR",
                   "Unfolding model": "UNCORR",
                   "Int. luminosity": "CMSLUMI16",
                   "Lepton energy": "UNCORR",
                   "Efficiency": "UNCORR",
                   "Backgrounds": "CORR",
                   "Jet energy": "CORR",
                   "Others": "CORR"}
    definitions = {}
    for key,value in uncertainties[0].items():
        definition = {key :
                      {"description": key + " unc. from HepData",
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
    filter_CMS_Z0J_13TEV_PT()
