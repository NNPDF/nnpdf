import yaml

from filter_utils import get_kinematics, get_data_values, get_systematics


def filter_CMS_WP_13TEV_data_kinetic():
    """
    writes data central values and kinematics
    to respective .yaml file
    """
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)
    
    version = metadata["hepdata"]["version"]
    tables = metadata["implemented_observables"][0]["tables"]
    
    kin = get_kinematics(version)
    central_values = get_data_values(version)


    data_central_yaml = {"data_central": central_values}
    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)



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
    filter_CMS_WP_13TEV_data_kinetic()
    
    