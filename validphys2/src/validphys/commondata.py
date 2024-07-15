"""
commondata.py

Module containing actions which return loaded commondata, leverages utils
found in :py:mod:`nnpdf_data.commondataparser`, and returns objects from
:py:mod:`nnpdf_data.coredata`

"""

import functools

from nnpdf_data.commondataparser import load_commondata
from reportengine import collect


@functools.cache
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
    loaded_cut_commondata: nnpdf_data.coredata.CommonData

    """
    lcd = load_commondata(commondata.metadata)
    return lcd.with_cuts(cuts)


dataset_inputs_loaded_cd_with_cuts = collect("loaded_commondata_with_cuts", ("data_input",))

groups_dataset_inputs_loaded_cd_with_cuts = collect(
    "loaded_commondata_with_cuts", ("group_dataset_inputs_by_metadata", "data_input")
)

groups_dataset_inputs_loaded_cd_with_cuts_byprocess = collect(
     "loaded_commondata_with_cuts", ("group_dataset_inputs_by_process", "data")
 )