import yaml

from filter_utils import get_kinematics, get_data_values, get_systematics


def filter_CMS_W_13TEV_data_kinetic(figure):
    """
    writes data central values and kinematics
    to respective .yaml file
    """
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    tables = metadata["implemented_observables"][0]["tables"]

    kin = get_kinematics(version, figure)
    central_values = get_data_values(version, figure)

    data_central_yaml = {"data_central": central_values}
    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    if figure == "17a":
        with open("data_WP.yaml", "w") as file:
            yaml.dump(data_central_yaml, file, sort_keys=False)

        with open("kinematics_WP.yaml", "w") as file:
            yaml.dump(kinematics_yaml, file, sort_keys=False)

    elif figure == "17b":
        with open("data_WM.yaml", "w") as file:
            yaml.dump(data_central_yaml, file, sort_keys=False)

        with open("kinematics_WM.yaml", "w") as file:
            yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_CMS_W_13TEV_uncertainties(observable):
    """
    writes uncertainties to respective .yaml file

    Parameters
    ----------
    observable : str
        eg. "W+" or "W-"
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    systematics = get_systematics(observable)

    # error definition
    error_definitions = {}
    errors = []

    for sys in systematics:
        
        error_definitions[sys[0]['name']] = {
            "description": f"{sys[0]['name']}",
            "treatment": "ADD",
            "type": "CORR",
        }

    #
    for i in range(metadata['implemented_observables'][0]['ndata']):
        error_value = {}

        for sys in systematics:
            error_value[sys[0]['name']] = float(sys[0]['values'][i])

        errors.append(error_value)

    uncertainties_yaml = {"definitions": error_definitions, "bins": errors}

    # write uncertainties
    if observable == "W+":
        with open(f"uncertainties_WP.yaml", 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)
    elif observable == "W-":
        with open(f"uncertainties_WM.yaml", 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    # WP data
    filter_CMS_W_13TEV_data_kinetic(figure="17a")
    filter_CMS_W_13TEV_uncertainties(observable="W+")

    # WM data
    filter_CMS_W_13TEV_data_kinetic(figure="17b")
    filter_CMS_W_13TEV_uncertainties(observable="W-")
