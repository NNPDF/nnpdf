"""
theorycovarianceutils.py

Low level utilities for theorycovariance module
"""

import logging

from reportengine.checks import check, make_argcheck
from validphys.loader import Loader
from validphys.plotoptions.core import get_info

log = logging.getLogger(__name__)


def check_correct_theory_combination_internal(
    theoryids, point_prescription: (str, type(None)) = None
):
    """Checks that a valid theory combination corresponding to an existing
    prescription has been inputted"""
    l = len(theoryids)
    check(
        l in {3, 5, 7, 9, 62, 64, 66, 70, 19, 23},
        f"Expecting exactly 3, 5, 7, 9, 62, 64, 66, 23, 19 or 70 theories, but got {l}.",
    )
    opts = {"bar", "nobar"}
    xifs = [theoryid.get_description()["XIF"] for theoryid in theoryids]
    xirs = [theoryid.get_description()["XIR"] for theoryid in theoryids]
    if l == 3:
        if point_prescription == "3f point":
            correct_xifs = [1.0, 2.0, 0.5]
            correct_xirs = [1.0, 1.0, 1.0]
        elif point_prescription == "3r point":
            correct_xifs = [1.0, 1.0, 1.0]
            correct_xirs = [1.0, 2.0, 0.5]
        elif "n3lo" in point_prescription or "missing" in point_prescription:
            correct_xifs = [1.0, 1.0, 1.0]
            correct_xirs = [1.0, 0.5, 2.0]
        else:
            correct_xifs = [1.0, 2.0, 0.5]
            correct_xirs = [1.0, 2.0, 0.5]
    elif l == 5:
        if point_prescription == "5 point":
            correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0]
            correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5]
        elif point_prescription == "5bar point":
            correct_xifs = [1.0, 2.0, 0.5, 2.0, 0.5]
            correct_xirs = [1.0, 2.0, 0.5, 0.5, 2.0]
    elif l == 7:
        correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0, 2.0, 0.5]
        correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5]
    elif l in [62, 64, 66, 70]:
        # check Anomalous dimensions variations
        id_max = None
        if l == 70:
            id_max = -8
        elif l == 66:
            id_max = -4
        elif l == 64:
            id_max = -2
        n3lo_vars_list = []
        for theoryid in theoryids[:id_max]:
            n3lo_ad_variation = theoryid.get_description()["n3lo_ad_variation"]
            # Only take the first 4, the last three are NS and only relevant for FHMRUVV
            if any(n3lo_ad_variation[4:]):
                raise ValueError(
                    f"Theory {theoryid.id} has non-zero entries in the last three (NS) elements "
                    "of n3lo_ad_variation, the covmat construction does not support this!"
                )
            n3lo_vars_list.append(theoryid.get_description()["n3lo_ad_variation"][:4])
        full_var_list = [[0, 0, 0, 0]]
        n3lo_vars_dict = {"gg": 19, "gq": 21, "qg": 15, "qq": 6}
        for entry, max_var in enumerate(n3lo_vars_dict.values()):
            for idx in range(1, max_var + 1):
                base_var = [0, 0, 0, 0]
                base_var[entry] = idx
                full_var_list.append(base_var)
        check(
            n3lo_vars_list == full_var_list,
            f"Theories do not include the full list of N3LO variation but {n3lo_vars_list}",
        )
        if l == 70:
            # check Scale variations
            varied_xifs = [xifs[0]]
            varied_xirs = [xirs[0]]
            varied_xifs.extend(xifs[-8:-2])
            varied_xirs.extend(xirs[-8:-2])
            correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0, 2.0, 0.5]
            correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5]
            check(
                varied_xifs == correct_xifs and varied_xirs == correct_xirs,
                "Choice of input theories does not correspond to a valid "
                "prescription for theory covariance matrix calculation",
            )
        elif l == 66:
            # check Scale variations
            varied_xifs = [xifs[0]]
            varied_xirs = [xirs[0]]
            varied_xifs.extend(xifs[-4:-2])
            varied_xirs.extend(xirs[-4:-2])
            correct_xifs = [1.0, 1.0, 1.0]
            correct_xirs = [1.0, 0.5, 2.0]
            check(
                varied_xifs == correct_xifs and varied_xirs == correct_xirs,
                "Choice of input theories does not correspond to a valid "
                "prescription for theory covariance matrix calculation",
            )
        return
    elif l in [19, 23]:
        if l == 23:
            # check Scale variations
            varied_xifs = [xifs[0]]
            varied_xirs = [xirs[0]]
            varied_xifs.extend(xifs[-8:-2])
            varied_xirs.extend(xirs[-8:-2])
            correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0, 2.0, 0.5]
            correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5]
            check(
                varied_xifs == correct_xifs and varied_xirs == correct_xirs,
                "Choice of input theories does not correspond to a valid "
                "prescription for theory covariance matrix calculation",
            )
        elif l == 19:
            # check Scale variations
            varied_xifs = [xifs[0]]
            varied_xirs = [xirs[0]]
            varied_xifs.extend(xifs[-4:-2])
            varied_xirs.extend(xirs[-4:-2])
            correct_xifs = [1.0, 1.0, 1.0]
            correct_xirs = [1.0, 0.5, 2.0]
            check(
                varied_xifs == correct_xifs and varied_xirs == correct_xirs,
                "Choice of input theories does not correspond to a valid "
                "prescription for theory covariance matrix calculation",
            )
        return
    else:
        correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5]
        correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5, 0.5, 2.0]
    check(
        xifs == correct_xifs and xirs == correct_xirs,
        "Choice of input theories does not correspond to a valid "
        "prescription for theory covariance matrix calculation",
    )


check_correct_theory_combination = make_argcheck(check_correct_theory_combination_internal)


@make_argcheck
def check_fit_dataset_order_matches_grouped(
    group_dataset_inputs_by_metadata, data_input, processed_metadata_group
):
    """
    Check for use with theory covmat generation.

    Makes sure that the order of datasets listed in the fit runcard is the same
    as that specified by the metadata grouping. Otherwise there can be a
    misalignment between the experiment covmat and theory covmat.
    """
    data_input_iter = iter(data_input)
    for group in group_dataset_inputs_by_metadata:
        for dsinput in group["data_input"]:
            grouped_ds = dsinput.name
            input_ds = next(data_input_iter).name
            check(
                grouped_ds == input_ds,
                "Dataset ordering is changed by grouping, this will cause "
                "errors when running fits with theory covmat. Datasets should "
                f"be ordered by {processed_metadata_group} in the runcard.",
            )


def process_lookup(name):
    """
    Returns the `nnpdf31_process` of the corresponding dataset.
    """
    cd = Loader().check_commondata(setname=name)
    proc = get_info(cd).nnpdf31_process
    return proc
