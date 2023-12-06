from validphys.loader import Loader
import yaml


SETNAME = "ATLAS_2JET_7TEV_R06"

def filter_ATLAS_2JET_7TEV_R06_uncertainties_bugged():
    """
    - read systematics from old CommonData format
    - write to uncertainties_bugged.yaml the old bugged systematics
    
    This reproduces the same covariance matrix as the old CommonData.
    """

    l = Loader()
    cd = l.check_commondata(setname=SETNAME).load_commondata_instance()
    
    additive_sys = cd.commondata_table.drop(['process','kin1','kin2','kin3','data', 'stat'],axis=1)['ADD'].to_numpy()
    systype = cd.systype_table

    # error definition
    error_definition = {}
    for index, row in systype.iterrows():
        error_definition[f"sys_{index}"] = {
            "description": f"sys_{index}",
            "treatment": row['type'],
            "type": row['name'],

        }

    # store error in dict
    error = []
    for n in range(additive_sys.shape[0]):
        error_value = {}
        for m in range(additive_sys.shape[1]):
            error_value[f"sys_{m+1}"] = float(additive_sys[n, m])

        # error_value["luminosity_uncertainty"] = float(atlas_lumi[n])
        error.append(error_value)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}

    # save uncertainties to .yaml file
    with open("uncertainties_bugged.yaml", 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_2JET_7TEV_R06_uncertainties_bugged()
