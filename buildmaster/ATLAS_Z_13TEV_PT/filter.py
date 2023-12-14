import yaml

from filter_utils import get_kinematics, get_data_values, get_systematics

UNCORRELATED_SYS = ["Stat (Data)", "Stat (MC)", "Efficiencies (Uncorellated)"]

def filter_ATLAS_Z_13TEV_PT_data_kinetic():
    """
    writes data central values and kinematics
    to respective .yaml file
    """
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)
    
    version = metadata["hepdata"]["version"]
    # tables for Z->e+e- observable
    tables_E = metadata["implemented_observables"][0]["tables"]
    # tables for Z->mu+mu- observable
    tables_MU = metadata["implemented_observables"][1]["tables"]
    
    kin_E = get_kinematics(tables_E, version)
    central_values_E = get_data_values(tables_E, version)

    kin_MU = get_kinematics(tables_MU, version)
    central_values_MU = get_data_values(tables_MU, version)


    data_central_E_yaml = {"data_central": central_values_E}
    kinematics_E_yaml = {"bins": kin_E}

    data_central_MU_yaml = {"data_central": central_values_MU}
    kinematics_MU_yaml = {"bins": kin_MU}

    # write central values and kinematics to yaml file
    with open("data_E.yaml", "w") as file:
        yaml.dump(data_central_E_yaml, file, sort_keys=False)

    with open("kinematics_E.yaml", "w") as file:
        yaml.dump(kinematics_E_yaml, file, sort_keys=False)

    with open("data_MU.yaml", "w") as file:
        yaml.dump(data_central_MU_yaml, file, sort_keys=False)

    with open("kinematics_MU.yaml", "w") as file:
        yaml.dump(kinematics_MU_yaml, file, sort_keys=False)

def filter_ATLAS_Z_13TEV_PT_uncertainties():
    """
    writes uncertainties to respective .yaml file
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    # tables for Z->e+e- observable
    tables_E = metadata["implemented_observables"][0]["tables"]
    # tables for Z->mu+mu- observable
    tables_MU = metadata["implemented_observables"][1]["tables"]
    tables = {"E": tables_E, "MU": tables_MU}

    systematics_E = get_systematics(tables_E, version)
    systematics_MU = get_systematics(tables_MU, version)
    systematics = {"E": systematics_E, "MU": systematics_MU}

    # error definition
    error_definitions = {}
    errors = {} 

    for obs in ["E","MU"]:

        error_definitions[obs] = {}

        for sys in systematics[obs]:

            if sys[0]['name'] in UNCORRELATED_SYS:
                error_definitions[obs][sys[0]['name']] = {
                    "description": f"{sys[0]['name']} from HEPDATA",
                    "treatment": "ADD",
                    "type": "UNCORR",
                }

            else:
                error_definitions[obs][sys[0]['name']] = {
                    "description": f"{sys[0]['name']} from HEPDATA",
                    "treatment": "ADD",
                    "type": "CORR",
                }

        # TODO:
        # store error in dict
        errors[obs] = []
        central_values = get_data_values(tables[obs], version)

        for i in range(len(central_values)):
            error_value = {}
            
            for sys in systematics[obs]:
                error_value[sys[0]['name']] = float(sys[0]['values'][i])

            errors[obs].append(error_value)
        
        uncertainties_yaml = {"definitions": error_definitions[obs], "bins": errors[obs]}

        # write uncertainties
        with open(f"uncertainties_{obs}.yaml", 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)
        
    




if __name__ == "__main__":
    filter_ATLAS_Z_13TEV_PT_data_kinetic()
    filter_ATLAS_Z_13TEV_PT_uncertainties()
    
