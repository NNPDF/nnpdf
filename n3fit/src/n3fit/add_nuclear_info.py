"""
This module is a really dirty way to propagate the A-dependence into
the fit. By A-dependence on means the atomic number of the nucleus
involved in the fitting. There are two main reasons why this is done
in the current way:
    1. The information on the target(s) are not stored in the current
    FK tables. Ideally, one should have all the information concerning
    what are being collided as metadata in the tables. Currently, this
    is not the case, but this will change with the advent of `pineko`.
    2. There are still ongoing works happening towards (a) stripping
    off the C++ dependence of the code, and (b) modifying various parts
    of the codes to accommodate for the new `pineko` pipeline. In this
    sense, it is not worth spending huge amount of time implementing
    something that will be changed later.

This module will become completely obsolete once the new pipeline is
in place and working.
"""


NUCLEAR_DATASPEC = {
    "NMC_AL_C": {"A": [27, 12], "Z": [13, 6]},
    "NMC_CA_D": {"A": [40, 2], "Z": [20, 1]},
    "NMC_BE_C": {"A": [9, 12], "Z": [4, 6]},
    "NMC_CA_C": {"A": [40, 12], "Z": [20, 6]},
    "NMC_CA_LI": {"A": [40, 6], "Z": [20, 3]},
    "NMC_C_LI": {"A": [12, 6], "Z": [6, 3]},
    "NMC_C_D": {"A": [12, 2], "Z": [6, 1]},
    "NMC_FE_C": {"A": [56, 12], "Z": [26, 6]},
    "NMC_PB_C": {"A": [208, 12], "Z": [82, 6]},
    "NMC_SN_C": {"A": [119, 12], "Z": [50, 6]},
    "NMC_LI_D": {"A": [6, 2], "Z": [3, 1]},
    "NMC_HE_D": {"A": [4, 2], "Z": [2, 1]},
    "NMC_p_D": {"A": [2, 1], "Z": [1, 1]},
    "SLAC_AG_D": {"A": [108, 2], "Z": [47, 1]},
    "SLAC_AL_D": {"A": [27, 2], "Z": [13, 1]},
    "SLAC_AU_D": {"A": [197, 2], "Z": [79, 1]},
    "SLAC_BE_D": {"A": [9, 2], "Z": [4, 1]},
    "SLAC_CA_D": {"A": [40, 2], "Z": [20, 1]},
    "SLAC_C_D": {"A": [12, 2], "Z": [6, 1]},
    "SLAC_FE_D": {"A": [56, 2], "Z": [26, 1]},
    "SLAC_HE_D": {"A": [4, 2], "Z": [2, 1]},
    "EMC_CU_D": {"A": [64, 2], "Z": [29, 1]},
    "EMC_FE_D": {"A": [56, 2], "Z": [26, 1]},
    "EMC_SN_D": {"A": [119, 2], "Z": [50, 1]},
    "EMC_C_D": {"A": [12, 2], "Z": [6, 1]},
    "BCDMS_FE_D": {"A": [56, 2], "Z": [26, 1]},
    "BCDMS_N_D": {"A": [14, 2], "Z": [7, 1]},
    "FNAL_C_D": {"A": [12, 2], "Z": [6, 1]},
    "FNAL_C_D": {"A": [12, 2], "Z": [6, 1]},
    "FNAL_PB_D": {"A": [208, 2], "Z": [82, 1]},
    "FNAL_XE_D": {"A": [131, 2], "Z": [54, 1]},
    "SLAC_D": {"A": [2], "Z": [1]},
    "BCDMS_D": {"A": [2], "Z": [1]},
}


def add_nuclear_dependence(datasetinfo: list) -> list:
    """Takes the usual inputs for n3fit (raw datasets from validphys) and
    add the A-values to the dictionaries.

    Parameters:
    -----------
    datasetinfo: list
        list of dictionaries containing information of a given dataset

    Returns:
    --------
    list:
        the same list as the input but with extra-information on the
        nuclear datasets
    """
    for dataset_group in datasetinfo:
        for dataset in dataset_group["datasets"]:
            if dataset["name"] in NUCLEAR_DATASPEC.keys():
                dataset['A'] = NUCLEAR_DATASPEC[dataset["name"]]["A"]
                dataset['Z'] = NUCLEAR_DATASPEC[dataset["name"]]["Z"]
    return datasetinfo


def list_active_nuclei(datasetinfo: list) -> list:
    """Take the new list from `add_nuclear_dependence` in which the info
    on the atomic mass number A and the atomic number Z is included and 
    returns an ordered list of the A included in the fit.

    Parameters:
    -----------
    datasetinfo: list
        list of dictionaries containing information of a given dataset

    Returns:
    --------
        ordered list of the active A involved in the fit
    """
    nuclear_lists = []
    for dataset_group in datasetinfo:
        for dataset in dataset_group["datasets"]:
            nuclear_lists.append(dataset.get("A", [1]))
    merged = [item for sublist in nuclear_lists for item in sublist]
    return sorted(list(set(merged)))
