"""


"""

import yaml

def get_tables(observable):
    """
    get the Hepdata tables, given the tables and version specified in metadata
    """
    prefix = "rawdata/Table"
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    if observable == "PT-Y":
        tables = metadata["implemented_observables"][0]["tables"]
    elif observable == "PT-M":
        tables = metadata["implemented_observables"][1]["tables"]
    else:
        print("Data set not implemented")
        print("Available data sets are:")
        print("- ATLAS_Z0J_8TEV_PT-Y")
        print("- ATLAS_Z0J_8TEV_PT-M")
        exit()
            
    hepdata_tables = []

    for table in tables:
        hepdata_tables.append(f"{prefix}{table}.yaml")

    return hepdata_tables

def get_uncertainties(observable):
    """
    Returns uncertainties for dumping in the -yaml file
    """
    data_central = []
    uncertainties = []

    hepdata_tables = get_tables(observable)
    for table in hepdata_tables:
        with open(table, 'r') as f:
            input = yaml.safe_load(f)
        # Central values
        data_values = input["dependent_variables"][5]["values"]
        for data_value in data_values:
            data_central.append(data_value["value"])

        # Uncertainties
        for data_value in data_values:
            errors = data_value["errors"]
            uncertainty = {}
            for error in errors:
                uncertainty[error["label"]] = error["symerror"]
                uncertainty.update(uncertainty)

            uncertainties.append(uncertainty)

    return (data_central, uncertainties)

def filter_unc_ATLAS_Z0J_8TEV(MC=False):
    """
    Dumps uncertainties on .yaml files
    """
    lumi_unc = 2.8 # %
    mc_unc = 1.0 # %
    observables = ["PT-Y", "PT-M"]
    for observable in observables:
        if MC==False:
            unc_file = "uncertainties_decorr_" + observable + ".yaml"
        else:
            unc_file= "uncertainties_decorr_sys_10_" + observable + ".yaml"
        central_values, uncertainties = get_uncertainties(observable)

        for i in range(len(central_values)):
            for k in uncertainties[i]:
                uncertainties[i][k] = float(uncertainties[i][k].replace("%",""))/100. * central_values[i] * 1000.
            uncertainties[i].update({"sys_lumi_corr": lumi_unc/100 * central_values[i] * 1000.})
            if(MC==True):
                uncertainties[i].update({"sys_mc_uncorr": mc_unc/100 * central_values[i] * 1000.})  
                                 
        # Uncertainties
        treatment = {"stat": "ADD",
                     "sys,Uncorrelated": "ADD",
                     "sys,Correlated": "MULT",
                     "sys_lumi_corr": "MULT",}
        correlation = {"stat": "UNCORR",
                       "sys,Uncorrelated": "UNCORR",
                       "sys,Correlated": "CORR",
                       "sys_lumi_corr": "ATLASLUMI12",}
        if MC == True:
            treatment.update({"sys_mc_uncorr": "ADD"})
            correlation.update({"sys_mc_uncorr": "UNCORR"})

        definitions = {}
        for key,value in uncertainties[0].items():
            definition = {key :
                          {"description": key + " unc. from HepData",
                           "treatment": treatment[key],
                           "type": correlation[key]}}
            definitions.update(definition)
            uncertainties_yaml = {"definitions": definitions,"bins": uncertainties}
            
        with open(unc_file, "w") as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)
        

if __name__ == "__main__":
    filter_unc_ATLAS_Z0J_8TEV(MC=False)
    filter_unc_ATLAS_Z0J_8TEV(MC=True)
