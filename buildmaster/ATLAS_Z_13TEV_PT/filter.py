"""

"""

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

    # tables for Z->l+l- observable
    tables = metadata["implemented_observables"][0]["tables"]

    kin = get_kinematics(tables, version)
    central_values = get_data_values(tables, version)

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
    # tables for Z->l+l- observable
    tables = metadata["implemented_observables"][0]["tables"]

    systematics_LL = get_systematics(tables, version)

    systematics = {"LL": systematics_LL}

    # error definition
    error_definitions = {}
    errors = {}

    for obs in ["LL"]:

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

        central_values = get_data_values(tables, version)

        for i in range(len(central_values)):
            error_value = {}

            for sys in systematics[obs]:
                error_value[sys[0]['name']] = float(sys[0]['values'][i])

            errors[obs].append(error_value)

        uncertainties_yaml = {"definitions": error_definitions[obs], "bins": errors[obs]}

        # write uncertainties
        with open(f"uncertainties.yaml", 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_Z_13TEV_PT_data_kinetic()
    filter_ATLAS_Z_13TEV_PT_uncertainties()
