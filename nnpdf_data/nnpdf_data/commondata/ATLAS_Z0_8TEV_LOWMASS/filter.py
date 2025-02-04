"""
filter.py module for ATLAS_Z0_8TEV_LOWMASS dataset
When running `python filter.py` the relevant uncertainties , data and kinematics yaml
file will be created in the `nnpdf_data/commondata/ATLAS_Z0_8TEV_LOWMASS` directory.
"""

import yaml
from filter_utils import get_kinematics, get_data_values, get_systematics, get_systematics_light
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)


def filter_ATLAS_Z0_8TEV_LOWMASS_data_kinetic():
    """
    This function writes the central values and kinematics to yaml files.
    """

    kin = get_kinematics()
    central_values = list(get_data_values())

    data_central_yaml = {"data_central": central_values}

    kinematics_yaml = {"bins": kin}

    # write central values and kinematics to yaml file
    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_ATLAS_Z0_8TEV_LOWMASS_systematics(version=3):
    """
    This function writes the systematics to a yaml file.
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    systematics = get_systematics(version=version)

    # error definition
    error_definitions = {}
    errors = []

    for sys in systematics:
        if (sys[0]['name'] == 'stat') or (sys[0]['name'] == 'sys,uncor'):
            error_definitions[sys[0]['name']] = {
                "description": f"{sys[0]['name']}",
                "treatment": "ADD",
                "type": "UNCORR",
            }

        elif (sys[0]['name'] == 'ATLAS_LUMI') or (sys[0]['name'] == 'Lumi:M'):
            if version == 3 or version == 1:
                type_lumi = "ATLASLUMI12"
            elif version == "3_uncorr" or version == "1_uncorr":
                type_lumi = "UNCORR"
            elif version == "3_corr" or version == "1_corr":
                type_lumi = "CORR"
            error_definitions[sys[0]['name']] = {
                "description": f"{sys[0]['name']}",
                "treatment": "MULT",
                "type": type_lumi
            }

        else:
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
    if version == 1:
        
        filename = "uncertainties_new_v1_corr_lumi_all.yaml"
    elif version == 3:
        
        filename = "uncertainties.yaml"
    elif version == "3_uncorr":
        
        filename = "uncertainties_new_v3_uncorr_lumi.yaml"
    elif version == "1_uncorr":
        
        filename = "uncertainties_new_v1_uncorr_lumi.yaml"
    elif version == "3_corr":
        
        filename = "uncertainties_new_v3_corr_lumi.yaml"
    elif version == "1_corr":
        
        filename = "uncertainties_new_v1_corr_lumi.yaml"
    else:    
        raise ValueError(f"Version {version} is not supported.")
    
    with open(filename, 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)

def filter_ATLAS_Z0_8TEV_LOWMASS_systematics_light(lumi_corr):
    """
    This function writes the light version of systematics 
    to a yaml file.
    """

    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    systematics = get_systematics_light()

    # error definition
    error_definitions = {}
    errors = []

    for sys in systematics:
        # Uncorrelated uncertainties
        if (sys['name'] == 'stat') or (sys['name'] == 'syst_unc'):
            error_definitions[sys['name']] = {
                "description": f"{sys['name']}",
                "treatment": "ADD",
                "type": "UNCORR",
            }
            
        elif sys['name'] == 'ATLAS_LUMI':
            error_definitions[sys['name']] = {
                "description": f"{sys['name']}",
                "treatment": "MULT",
                "type": lumi_corr
            }

        else:
            error_definitions[sys['name']] = {
                "description": f"{sys['name']}",
                "treatment": "ADD",
                "type": "CORR",
            }

    for i in range(metadata['implemented_observables'][0]['ndata']):
        error_value = {}

        for sys in systematics:
            error_value[sys['name']] = float(sys['values'][i])

        errors.append(error_value)

    uncertainties_yaml = {"definitions": error_definitions, "bins": errors}

    # write uncertainties
    if lumi_corr == "ATLASLUMI12":
      filename = "uncertainties_sys_light_corr_all.yaml"
    elif lumi_corr == "CORR":
      filename = "uncertainties_sys_light_corr.yaml"
    elif lumi_corr == "UNCORR":
      filename = "uncertainties_sys_light_uncorr.yaml"
    else:
        raise ValueError(f"{lumi_corr} is not supported.")
    
    with open(filename, 'w') as file:
          yaml.dump(uncertainties_yaml, file, sort_keys=False)

if __name__ == "__main__":
    filter_ATLAS_Z0_8TEV_LOWMASS_data_kinetic()
    filter_ATLAS_Z0_8TEV_LOWMASS_systematics(version=3)
    filter_ATLAS_Z0_8TEV_LOWMASS_systematics(version="3_uncorr")
    filter_ATLAS_Z0_8TEV_LOWMASS_systematics(version="3_corr")
    filter_ATLAS_Z0_8TEV_LOWMASS_systematics(version=1)
    filter_ATLAS_Z0_8TEV_LOWMASS_systematics(version="1_uncorr")
    filter_ATLAS_Z0_8TEV_LOWMASS_systematics(version="1_corr")
    filter_ATLAS_Z0_8TEV_LOWMASS_systematics_light(lumi_corr="ATLASLUMI12")
    filter_ATLAS_Z0_8TEV_LOWMASS_systematics_light(lumi_corr="CORR")
    filter_ATLAS_Z0_8TEV_LOWMASS_systematics_light(lumi_corr="UNCORR")
