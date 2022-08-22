"""
commondata.py

Module containing actions which return loaded commondata, leverages utils
found in :py:mod:`validphys.commondataparser`, and returns objects from
:py:mod:`validphys.coredata`

"""
from reportengine import collect

from validphys.commondataparser import load_commondata

def loaded_commondata_with_cuts(commondata, cuts):
    """Load the commondata and apply cuts.

    Parameters
    ----------
    commondata: validphys.core.CommonDataSpec
        commondata to load and cut.
    cuts: validphys.core.cuts, None
        valid cuts, used to cut loaded commondata.

    Returns
    -------
    loaded_cut_commondata: validphys.coredata.CommonData

    """
    lcd = load_commondata(commondata)
    return lcd.with_cuts(cuts)

dataset_inputs_loaded_cd_with_cuts = collect(
    "loaded_commondata_with_cuts", ("data_input",)
)

groups_dataset_inputs_loaded_cd_with_cuts = collect(
    "loaded_commondata_with_cuts", ("group_dataset_inputs_by_fitting_group", "data_input")
)

