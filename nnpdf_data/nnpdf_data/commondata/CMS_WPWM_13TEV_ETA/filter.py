from filter_utils import get_data_values, get_kinematics, get_systematics
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

import numpy as np

yaml.add_representer(float, prettify_float)


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


def filter_CMS_W_13TEV_uncertainties(observable, figure):
    """
    writes uncertainties to respective .yaml file

    Parameters
    ----------
    observable : str
        eg. "W+" or "W-"
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)
    version = metadata["hepdata"]["version"]

    systematics = get_systematics(observable, version, figure)

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


def get_data_CMS_W_13TEV_ASY():

    data_central = []
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)
    version = metadata["hepdata"]["version"]
    figure = "18"
    kinematics = get_kinematics(version, figure)
    hepdata_table = f"rawdata/HEPData-ins1810913-v{version}-Figure_{figure}.yaml"

    with open(hepdata_table, "r") as f:
        input = yaml.safe_load(f)

    data_values = input["dependent_variables"][0]["values"]

    for data_value in data_values:
        data_central.append(data_value["value"])

    ndata = len(data_central)
    syst_dict = {}

    # Luminosity uncertainty?
    value_id = 0

    for point in data_values:
        for err in point["errors"]:
            label = err["label"]
            symerr = err['symerror']

            if label not in syst_dict:
                syst_dict[label] = np.zeros(ndata)

            syst_dict[label][value_id] = symerr

        value_id += 1

    sys_list = []
    for label, values in syst_dict.items():
        sys_list.append({"name": label, "values": values.tolist()})

    return data_central, kinematics, sys_list


def filter_CMS_W_13TEV_ASY():
    central_values, kinematics, uncertainties = get_data_CMS_W_13TEV_ASY()
    data_central_yaml = {"data_central": central_values}

    kinematics_yaml = {"bins": kinematics}
    definitions = {
        uncertainties[0]['name']: {
            "description": uncertainties[0]['name'],
            "treatment": "ADD",
            "type": "CORR",
        }
    }
    errors_yaml = []
    unc_name = uncertainties[0]['name']
    for bin in range(len(central_values)):
        errors_yaml.append({unc_name: uncertainties[0]['values'][bin]})

    uncertainties_yaml = {"definitions": definitions, "bins": errors_yaml}
    with open("data_ASY.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)
    with open("kinematics_ASY.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)
    with open("uncertainties_ASY.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)
    return


if __name__ == "__main__":
    # WP data
    filter_CMS_W_13TEV_data_kinetic(figure="17a")
    filter_CMS_W_13TEV_uncertainties(observable="W+", figure="17a")

    # WM data
    filter_CMS_W_13TEV_data_kinetic(figure="17b")
    filter_CMS_W_13TEV_uncertainties(observable="W-", figure="17b")

    # ASY data
    filter_CMS_W_13TEV_ASY()
