"""

"""

import yaml
import numpy as np
from nnpdf_data.filter_utils.utils import covmat_to_artunc

def get_tables(observable=None):
    """
    Get the Hepdata tables, given the tables and version specified in metadata 
    """
    prefix = "rawdata/HEPData-ins2628732"
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]

    if observable == "WMWP-D":
        tables = metadata["implemented_observables"][0]["tables"]
    elif observable == "WMWP-Dstar":
        tables = metadata["implemented_observables"][1]["tables"]
    else:
        print("Observable not implemented.")
        print("Choose one of the following observables:")
        print("- WMWP-D")
        print("- WMWP-Dstar")
        
    hepdata_tables = []

    for table in tables:
        hepdata_tables.append(f"{prefix}-v{version}-Table_{table}.yaml")

    return hepdata_tables

def get_all(observable=None):
    """
    Returns data, kinematics and uncertainties for dumping in the .yaml files
    """
    data_central = []
    kinematics = []
    uncertainties = []

    hepdata_tables = get_tables(observable)
    data_tables = hepdata_tables[:-1]
    for table in data_tables:
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
                'abs_eta': {'min': kin_value['low'],
                            'mid': 0.5 * (kin_value['low'] + kin_value['high']),
                            'max': kin_value['high']},
                'm_W2': {'min': None, 'mid': 6.46174823e+03, 'max': None},}
            kinematics.append(kin)

    ndata = len(data_central)
    
    # Uncertainties
    # Construct luminosity covariance matrix
    lumi_unc = 0.83 #%
    lumi_uncs = []
    lumi_cov = []
    tot_cov = []
    for data in data_central:
        lumi_uncs.append(data * lumi_unc / 100.)
    for lumi_i in lumi_uncs:
        for lumi_j in lumi_uncs:
            lumi_cov.append(lumi_i * lumi_j)

    # Read total covariance matrix
    with open(hepdata_tables[2], 'r') as f:
            input = yaml.safe_load(f)
    cov_values = input["dependent_variables"][0]["values"]
    for cov_value in cov_values:
        tot_cov.append(cov_value["value"])

    # Compute covariance matrix without luminosity uncertainty    
    partial_cov = np.subtract(tot_cov,lumi_cov)

    # Generate artifical systematic uncertainties form partial_cov
    art_unc = covmat_to_artunc(ndata, partial_cov, 0)

    for i in range(len(art_unc)):
        errors = art_unc[i]
        uncertainty = {}
        for j in range(len(errors)):
            unc = {"art. sys. " + f"{j+1}" : errors[j]}
            uncertainty.update(unc)

        lumi_unc = {"luminosity": lumi_uncs[i] }
        uncertainty.update(lumi_unc)
        uncertainties.append(uncertainty)
        
    return (data_central, kinematics, uncertainties)

def filter_ATLAS_WCHARM_13TEV(observable=None):
    """
    Dumps data, kinematics, and uncertainties on .yaml files
    """
    central_values, kinematics, uncertainties = get_all(observable)
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
        if key == "luminosity":
            definition = {key :
                          {"description": key + " unc. from HepData",
                           "treatment": "MULT",
                           "type": "ATLASLUMI16"}}
        else:             
            definition = {key :
                          {"description": key + " unc. from HepData",
                           "treatment": "ADD",
                           "type": "CORR"}}                        
        definitions.update(definition)
    uncertainties_yaml = {"definitions": definitions,"bins": uncertainties}
    
    with open("data_" + observable + ".yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)
    with open("kinematics_" + observable + ".yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)
    with open("uncertainties_" + observable + ".yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)

if __name__ == "__main__":
    filter_ATLAS_WCHARM_13TEV("WMWP-D")
    filter_ATLAS_WCHARM_13TEV("WMWP-Dstar")


















"""
from filter_utils import (
    get_data_values,
    get_kinematics,
    get_artificial_uncertainties,
    get_uncertainties,
)
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def filter_ATLAS_WCHARM_13TEV_data_kinematic():
    

    central_values = get_data_values()

    kin = get_kinematics()

    data_central_yaml_WMDP = {"data_central": central_values[:5]}
    data_central_yaml_WPDM = {"data_central": central_values[5:10]}
    data_central_yaml_WMDP_star = {"data_central": central_values[10:15]}
    data_central_yaml_WPDM_star = {"data_central": central_values[15:20]}

    kinematics_yaml_WMDP = {"bins": kin[:5]}
    kinematics_yaml_WPDM = {"bins": kin[5:10]}
    kinematics_yaml_WMDP_star = {"bins": kin[10:15]}
    kinematics_yaml_WPDM_star = {"bins": kin[15:20]}

    # write central values and kinematics to yaml file
    for channel in ["WMDP", "WPDM", "WMDP_star", "WPDM_star"]:
        with open(f"data_{channel}.yaml", "w") as file:
            data_central_yaml = eval(f"data_central_yaml_{channel}")
            yaml.dump(data_central_yaml, file, sort_keys=False)

    for channel in ["WMDP", "WPDM", "WMDP_star", "WPDM_star"]:
        kinematics_yaml = eval(f"kinematics_yaml_{channel}")
        with open(f"kinematics_{channel}.yaml", "w") as file:
            yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_get_artificial_uncertainties():
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    systematics_full = get_artificial_uncertainties()

    for channel in ["WMDP", "WPDM", "WMDP_star", "WPDM_star"]:

        systematics = systematics_full[channel]

        error_definitions = {}
        errors = []

        for sys in systematics:
            if sys['name'] == 'stat':
                error_definitions[sys['name']] = {
                    "description": "Uncorrelated statistical uncertainties",
                    "treatment": "ADD",
                    "type": "UNCORR",
                }
            else:
                error_definitions[sys['name']] = {
                    "description": "Systematic uncertainty",
                    "treatment": "MULT",  # Not sure if this is correct
                    "type": "CORR",
                }

        for i in range(metadata['implemented_observables'][0]['ndata']):
            error_value = {}

            for sys in systematics:
                error_value[sys['name']] = float(sys['values'][i])

            errors.append(error_value)

        uncertainties_yaml = {"definitions": error_definitions, "bins": errors}

        # write uncertainties
        with open(f"uncertainties_covariances_{channel}.yaml", 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)


def filter_get_systematics():
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    systematics_full, deltas = get_uncertainties()
    for channel in ["WMDP", "WPDM", "WMDP_star", "WPDM_star"]:
        systematics = systematics_full[channel]
        error_definitions = {}
        errors = []
        # print("Systematics:", systematics)
        for sys in systematics:
            if sys['name'] == 'stat':
                error_definitions[sys['name']] = {
                    "description": "Uncorrelated statistical uncertainties",
                    "treatment": "ADD",
                    "type": "UNCORR",
                }
            else:
                error_definitions[sys['name']] = {
                    "description": "Systematic uncertainty",
                    "treatment": "MULT",  # Not sure if this is correct
                    "type": "CORR",
                }

        for i in range(metadata['implemented_observables'][0]['ndata']):
            error_value = {}

            for sys in systematics:

                error_value[sys['name']] = float(sys['values'][i])

            errors.append(error_value)

        uncertainties_yaml = {"definitions": error_definitions, "bins": errors}

        # write uncertainties
        with open(f"uncertainties_{channel}.yaml", 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_WCHARM_13TEV_data_kinematic()
    filter_get_artificial_uncertainties()
    filter_get_systematics()
"""
