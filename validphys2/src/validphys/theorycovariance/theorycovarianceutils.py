"""
theorycovarianceutils.py

Low level utilities for theorycovariance module
"""
import logging

from reportengine.checks import make_argcheck, check
from validphys.loader import Loader
from validphys.plotoptions import get_info

log = logging.getLogger(__name__)


def check_correct_theory_combination_internal(
    theoryids, fivetheories, point_prescription: (str, type(None)) = None
):
    """Checks that a valid theory combination corresponding to an existing
    prescription has been inputted"""
    l = len(theoryids)
    check(l in {3, 5, 7, 9, 71, 81}, f"Expecting exactly 3, 5, 7 or 9 or 81 theories, but got {l}.")
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
        else:
            correct_xifs = [1.0, 2.0, 0.5]
            correct_xirs = [1.0, 2.0, 0.5]
    elif l == 5:
        check(
            fivetheories is not None,
            "For five input theories a prescription bar or nobar"
            "for the flag fivetheories must be specified.",
        )
        check(
            fivetheories in opts,
            "Invalid choice of prescription for 5 points",
            fivetheories,
            opts,
        )
        if fivetheories == "nobar":
            correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0]
            correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5]
        elif fivetheories == "bar":
            correct_xifs = [1.0, 2.0, 0.5, 2.0, 0.5]
            correct_xirs = [1.0, 2.0, 0.5, 0.5, 2.0]
    elif l == 7:
        correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0, 2.0, 0.5]
        correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5]
    elif l == 71 or l == 81:
        # check Anomalous dimensions variations
        n3lo_vars_dict = {
            "P_gg": 18,
            "P_gq": 24,
            "P_qg": 20,
            "P_qq": 8,
        }
        # TODO: for the moment fish the n3lo_ad_variation from the comments
        n3lo_vars_list = []
        id_max = -10 if l == 81 else None
        for theoryid in theoryids[:id_max]:
            n3lo_vars_list.append([int(val) for val in theoryid.get_description()["Comments"][28:-1].split(",")])
        full_var_list = [[0,0,0,0]]
        for entry, max_var in enumerate(n3lo_vars_dict.values()):
            for idx in range(1, max_var + 1):
                base_var = [0,0,0,0]
                base_var[entry] = idx
                full_var_list.append(base_var)
        check(
            n3lo_vars_list == full_var_list, 
            f"Theories do not include the full list of N3LO variation but {n3lo_vars_list}"
        )
        if l == 81:
            # check Scale variations
            varied_xifs = [xifs[0]]
            varied_xirs = [xirs[0]]
            varied_xifs.extend(xifs[-10:-2])
            varied_xirs.extend(xirs[-10:-2])
            correct_xifs = [1.0, 0.5, 2.0, 0.5, 1.0, 2.0, 0.5, 1.0, 2.0]
            correct_xirs = [1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 2.0, 2.0, 2.0]
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


check_correct_theory_combination = make_argcheck(
    check_correct_theory_combination_internal
)


@make_argcheck
def check_correct_theory_combination_theoryconfig(collected_theoryids, fivetheories):
    check_correct_theory_combination_internal(collected_theoryids[0], fivetheories)


@make_argcheck
def check_correct_theory_combination_dataspecs(dataspecs_theoryids, fivetheories):
    """Like check_correct_theory_combination but for matched dataspecs."""
    return check_correct_theory_combination.__wrapped__(
        dataspecs_theoryids, fivetheories
    )


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
