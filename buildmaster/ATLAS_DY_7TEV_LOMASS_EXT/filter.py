import yaml

def filter_ATLAS_DY_7TEV_LOMASS_EXT():
    with open("metadata.yaml", "r") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    tables = metadata["implemented_observables"][0]["tables"]
    npoints = metadata["implemented_observables"][0]["npoints"]

    data_central = []
    kin = []
    error = []
    input = []

    hepdata_tables = [
        "rawdata/HEPData-ins1288706-v" + str(version) + "-Table_6.yaml",
        "rawdata/HEPData-ins1288706-v" + str(version) + "-Table_7.yaml",
    ]

    for table in hepdata_tables:
        with open(table, "r") as file:
            input.append(yaml.safe_load(file))

    for j in range(npoints[0]):
        data_central_value = input[0]["dependent_variables"][0]["values"][j]["value"]
        data_central.append(data_central_value)
        mll_low = float(input[0]["independent_variables"][0]["values"][j]["low"])
        mll_hig = float(input[0]["independent_variables"][0]["values"][j]["high"])
        kin_value = {
            "mll": {
                "min": mll_low,
                "mid": 0.5 * (mll_low + mll_hig),
                "max": mll_hig,
            },
            "mll2": {
                "min": mll_low * mll_low,
                "mid": (0.5 * (mll_low + mll_hig))**2,
                "max": mll_hig * mll_hig,
            },
            "sqrts": {"min": 7000.0, "mid": 7000.0, "max": 7000.0},
        }
        kin.append(kin_value)

        error_value = {
            "stat_1": float(
                input[0]["dependent_variables"][0]["values"][j]["errors"][0][
                    "symerror"
                ].rstrip("%")
            ),
            "sys_reco": float(input[1]["dependent_variables"][0]["values"][j]["value"]),
            "sys_trig": float(input[1]["dependent_variables"][1]["values"][j]["value"]),
            "sys_iso": float(input[1]["dependent_variables"][2]["values"][j]["value"]),
            "sys_mjet": float(input[1]["dependent_variables"][3]["values"][j]["value"]),
            "sys_pTsc": float(input[1]["dependent_variables"][4]["values"][j]["value"]),
            "sys_res": float(input[1]["dependent_variables"][5]["values"][j]["value"]),
            "sys_MC": float(input[1]["dependent_variables"][6]["values"][j]["value"]),
            "sys_lumi": 3.5,
        }
        error.append(error_value)

    error_definition = {
        "stat_1": {
            "description": "statistical uncertainty",
            "treatment": "MULT",
            "type": "UNCORR",
        },
        "sys_reco": {
            "description": "TODO",
            "treatment": "MULT",
            "type": "CORR",
        },
        "sys_trig": {
            "description": "TODO",
            "treatment": "MULT",
            "type": "CORR",
        },
        "sys_iso": {
            "description": "TODO",
            "treatment": "MULT",
            "type": "CORR",
        },
        "sys_mjet": {
            "description": "TODO",
            "treatment": "MULT",
            "type": "CORR",
        },
        "sys_pTsc": {
            "description": "TODO",
            "treatment": "MULT",
            "type": "CORR",
        },
        "sys_res": {
            "description": "TODO",
            "treatment": "MULT",
            "type": "CORR",
        },
        "sys_MC": {
            "description": "TODO",
            "treatment": "MULT",
            "type": "CORR",
        },
        "sys_lumi": {
            "description": "TODO",
            "treatment": "MULT",
            "type": "ATLASLUMI11",
        },
    }

    data_central_yaml = {"data_central": data_central}
    kinematics_yaml = {"bins": kin}
    uncertainties_yaml = {"definition": error_definition, "bins": error}

    with open("data.yaml", "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)

    with open("kinematics.yaml", "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)

    with open("uncertainties.yaml", "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)


filter_ATLAS_DY_7TEV_LOMASS_EXT()
