"""
construction.py
Tools for constructing theory covariance matrices and computing their chi2s.
"""
from __future__ import generator_stop

from collections import defaultdict, namedtuple
import logging

import numpy as np
import pandas as pd
import scipy.linalg as la

from reportengine import collect
from reportengine.checks import check
from reportengine.table import table
from validphys.calcutils import all_chi2_theory, calc_chi2, central_chi2_theory
from validphys.checks import check_using_theory_covmat
from validphys.results import Chi2Data, procs_central_values, procs_central_values_no_table, results
from validphys.theorycovariance.theorycovarianceutils import (
    check_correct_theory_combination,
    check_fit_dataset_order_matches_grouped,
    process_lookup,
)

log = logging.getLogger(__name__)

theoryids_procs_central_values = collect(procs_central_values, ("theoryids",))

theoryids_procs_central_values_no_table = collect(procs_central_values_no_table, ("theoryids",))

collected_theoryids = collect("theoryids", ["theoryconfig"])


def make_scale_var_covmat(predictions):
    raise Exception(
        "If you are seeing this error, please open an issue. This function should not be used"
    )


@check_correct_theory_combination
def theory_covmat_singleprocess_no_table(
    theoryids_procs_central_values_no_table, procs_index, theoryids, fivetheories
):
    """Calculates the theory covariance matrix for scale variations.
    The matrix is a dataframe indexed by procs_index."""
    s = make_scale_var_covmat(theoryids_procs_central_values_no_table)
    df = pd.DataFrame(s, index=procs_index, columns=procs_index)
    return df


@table
@check_correct_theory_combination
def theory_covmat_singleprocess(theory_covmat_singleprocess_no_table, fivetheories):
    """Duplicate of theory_covmat_singleprocess_no_table but with a table decorator."""
    return theory_covmat_singleprocess_no_table


results_bytheoryids = collect(results, ("theoryids",))
each_dataset_results_bytheory = collect(
    "results_bytheoryids", ("group_dataset_inputs_by_process", "data")
)


@check_using_theory_covmat
def theory_covmat_dataset(
    results,
    results_bytheoryids,
    use_theorycovmat,  # for the check
    point_prescription,
    fivetheories=None,
    seventheories=None,
):
    """
    Compute the theory covmat for a collection of theoryids for a single dataset.

    In general this will come from some point prescription and it could be guessed from the input
    however, it takes as input all relevant variables for generality
    """
    _, theory_results = zip(*results_bytheoryids)
    _, central_th_result = results
    l = len(results_bytheoryids)

    # Remove the central theory from the list if it was included
    theory_results = [i for i in theory_results if i._theoryid != central_th_result._theoryid]
    cv = central_th_result.central_value

    # Compute the theory contribution to the covmats
    deltas = list(t.central_value - cv for t in theory_results)
    thcovmat = compute_covs_pt_prescrip(
        point_prescription, l, "A", deltas, fivetheories=fivetheories, seventheories=seventheories
    )

    return thcovmat


@check_correct_theory_combination
def theory_covmat_datasets(each_dataset_results_bytheory, fivetheories):
    """Produces an array of theory covariance matrices. Each matrix corresponds
    to a different dataset, which must be specified in the runcard."""
    dataset_covmats = []
    for dataset in each_dataset_results_bytheory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        dataset_covmats.append(s)
    return dataset_covmats


@check_correct_theory_combination
def total_covmat_datasets(each_dataset_results_bytheory, fivetheories):
    """Produces an array of total covariance matrices; the sum of experimental
    and scale-varied theory covariance matrices. Each matrix corresponds
    to a different dataset, which must be specified in the runcard.
    These are needed for calculation of chi2 per dataset."""
    dataset_covmats = []
    for dataset in each_dataset_results_bytheory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        sigma = dataset[0][0].covmat
        cov = s + sigma
        dataset_covmats.append(cov)
    return dataset_covmats


@check_correct_theory_combination
def total_covmat_diagtheory_datasets(each_dataset_results_bytheory, fivetheories):
    """Same as total_covmat_theory_datasets but for diagonal theory only"""
    dataset_covmats = []
    for dataset in each_dataset_results_bytheory:
        theory_centrals = [x[1].central_value for x in dataset]
        s = make_scale_var_covmat(theory_centrals)
        # Initialise array of zeros and set precision to same as FK tables
        s_diag = np.zeros((len(s), len(s)), dtype=np.float32)
        np.fill_diagonal(s_diag, np.diag(s))
        sigma = dataset[0][0].covmat
        cov = s_diag + sigma
        dataset_covmats.append(cov)
    return dataset_covmats


@table
def theory_block_diag_covmat(theory_covmat_datasets, procs_index):
    """Takes the theory covariance matrices for individual datasets and
    returns a data frame with a block diagonal theory covariance matrix
    by dataset"""
    s = la.block_diag(*theory_covmat_datasets)
    df = pd.DataFrame(s, index=procs_index, columns=procs_index)
    return df


@table
def theory_diagonal_covmat(theory_covmat_singleprocess_no_table, procs_index):
    """Returns theory covmat with only diagonal values"""
    s = theory_covmat_singleprocess_no_table.values
    # Initialise array of zeros and set precision to same as FK tables
    s_diag = np.zeros((len(s), len(s)), dtype=np.float32)
    np.fill_diagonal(s_diag, np.diag(s))
    df = pd.DataFrame(s_diag, index=procs_index, columns=procs_index)
    return df


procs_results_theory = collect("procs_results", ("theoryids",))


@check_correct_theory_combination
def total_covmat_procs(procs_results_theory, fivetheories):
    """Same as total_covmat_datasets but per experiment rather than
    per dataset. Needed for calculation of chi2 per experiment."""
    proc_result_covmats = []
    for proc_result in zip(*procs_results_theory):
        theory_centrals = [x[1].central_value for x in proc_result]
        s = make_scale_var_covmat(theory_centrals)
        sigma = proc_result[0][0].covmat
        cov = s + sigma
        proc_result_covmats.append(cov)
    return proc_result_covmats


def dataset_names(data_input):
    """Returns a list of the names of the datasets, in the same order as
    they are inputted in the runcard"""
    return [el.name for el in data_input]


ProcessInfo = namedtuple("ProcessInfo", ("theory", "namelist", "sizes", "data"))


def combine_by_type(each_dataset_results_bytheory, dataset_names):
    """Groups the datasets according to processes and returns three objects:
    theories_by_process: the relevant theories grouped by process type
    ordered_names: dictionary with keys of process type and values being the
                   corresponding list of names of datasets, in the order they
                   are appended to theories_by_process
    dataset_size:  dictionary with keys of dataset name and values being the
                   number of points in that dataset"""
    dataset_size = defaultdict(list)
    theories_by_process = defaultdict(list)
    ordered_names = defaultdict(list)
    for dataset, name in zip(each_dataset_results_bytheory, dataset_names):
        theory_centrals = [x[1].central_value for x in dataset]
        dataset_size[name] = len(theory_centrals[0])
        proc_type = process_lookup(name)
        ordered_names[proc_type].append(name)
        theories_by_process[proc_type].append(theory_centrals)
    for key, item in theories_by_process.items():
        theories_by_process[key] = np.concatenate(item, axis=1)
    process_info = ProcessInfo(
        theory=theories_by_process, namelist=ordered_names, sizes=dataset_size, data=None
    )
    return process_info


def combine_by_type_ht(
    each_dataset_results, groups_dataset_inputs_loaded_cd_with_cuts_byprocess
):
    """same as combine_by_type but now for a single theory and including commondata info"""
    dataset_size = defaultdict(list)
    theories_by_process = defaultdict(list)
    cd_by_process = defaultdict(list)
    ordered_names = defaultdict(list)
    for dataset, cd in zip(
        each_dataset_results, groups_dataset_inputs_loaded_cd_with_cuts_byprocess
    ):
        name = cd.setname
        theory_centrals = [x.central_value for x in dataset]
        dataset_size[name] = len(theory_centrals[0])
        proc_type = process_lookup(name)
        ordered_names[proc_type].append(name)
        cd_by_process[proc_type].append(cd.kinematics.values)
        theories_by_process[proc_type].append(theory_centrals)

    for key in theories_by_process.keys():
        theories_by_process[key] = np.concatenate(theories_by_process[key], axis=1)
        cd_by_process[key] = np.concatenate(cd_by_process[key], axis=0)
    process_info = ProcessInfo(
        theory=theories_by_process, namelist=ordered_names, sizes=dataset_size, data=cd_by_process
    )
    return process_info


def process_starting_points(combine_custom):
    """Returns a dictionary of indices in the covariance matrix corresponding
    to the starting point of each process."""
    process_info = combine_custom
    running_index = 0
    start_proc = defaultdict(list)
    for name in process_info.theory:
        size = len(process_info.theory[name][0])
        start_proc[name] = running_index
        running_index += size
    return start_proc


def covmap(combine_custom, dataset_names):
    """Creates a map between the covmat indices from matrices ordered by
    process to matrices ordered by experiment as listed in the runcard"""
    mapping = defaultdict(list)
    start_exp = defaultdict(list)
    process_info = combine_custom
    running_index = 0
    for dataset in dataset_names:
        size = process_info.sizes[dataset]
        start_exp[dataset] = running_index
        running_index += size
    start = 0
    names_by_proc_list = [item for sublist in process_info.namelist.values() for item in sublist]
    for dataset in names_by_proc_list:
        for i in range(process_info.sizes[dataset]):
            mapping[start + i] = start_exp[dataset] + i
        start += process_info.sizes[dataset]
    return mapping


def covmat_3fpt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 3pt factorisation
    scale variation *only*, given two dataset names and collections
    of scale variation shifts"""
    s = 0.5 * (np.outer(deltas1[0], deltas2[0]) + np.outer(deltas1[1], deltas2[1]))
    return s


def covmat_3rpt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 3pt renormalisation
    scale variation *only*, given two dataset names and collections
    of scale variation shifts"""
    if name1 == name2:
        s = 0.5 * (np.outer(deltas1[0], deltas2[0]) + np.outer(deltas1[1], deltas2[1]))
    else:
        s = 0.25 * (np.outer((deltas1[0] + deltas1[1]), (deltas2[0] + deltas2[1])))
    return s


def covmat_3pt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 3pt prescription,
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = 0.5 * sum(np.outer(d, d) for d in deltas1)
    else:
        s = 0.25 * (np.outer((deltas1[0] + deltas1[1]), (deltas2[0] + deltas2[1])))
    return s


def covmat_5pt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 5pt prescription,
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = 0.5 * sum(np.outer(d, d) for d in deltas1)
    else:
        s = 0.5 * (np.outer(deltas1[0], deltas2[0]) + np.outer(deltas1[1], deltas2[1])) + 0.25 * (
            np.outer((deltas1[2] + deltas1[3]), (deltas2[2] + deltas2[3]))
        )
    return s


def covmat_5barpt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 5barpt prescription,
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = 0.5 * sum(np.outer(d, d) for d in deltas1)
    else:
        s = 0.25 * (
            np.outer((deltas1[0] + deltas1[2]), (deltas2[0] + deltas2[2]))
            + np.outer((deltas1[1] + deltas1[3]), (deltas2[1] + deltas2[3]))
        )
    return s


def covmat_7pt_orig(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for original 7pt prescription,
    now deprecated but kept for posterity,
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = (1 / 3) * sum(np.outer(d, d) for d in deltas1)
    else:
        s = (1 / 6) * (
            np.outer((deltas1[0] + deltas1[4]), (deltas2[0] + deltas2[4]))
            + np.outer((deltas1[1] + deltas1[5]), (deltas2[1] + deltas2[5]))
            + np.outer((deltas1[2] + deltas1[3]), (deltas2[2] + deltas2[3]))
        )
    return s


def covmat_7pt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 7pt prescription (Gavin),
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = (1 / 3) * sum(np.outer(d, d) for d in deltas1)
    else:
        s = (1 / 6) * (
            2 * (np.outer(deltas1[0], deltas2[0]) + np.outer(deltas1[1], deltas2[1]))
            + (
                np.outer((deltas1[2] + deltas1[3]), (deltas2[2] + deltas2[3]))
                + np.outer((deltas1[4] + deltas1[5]), (deltas2[4] + deltas2[5]))
            )
        )
    return s


def covmat_9pt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 9pt prescription,
    given two dataset names and collections of scale variation shifts"""
    if name1 == name2:
        s = 0.25 * sum(np.outer(d, d) for d in deltas1)
    else:
        s = (1 / 12) * (
            np.outer((deltas1[0] + deltas1[4] + deltas1[6]), (deltas2[0] + deltas2[4] + deltas2[6]))
            + np.outer(
                (deltas1[1] + deltas1[5] + deltas1[7]), (deltas2[1] + deltas2[5] + deltas2[7])
            )
        ) + (1 / 8) * (np.outer((deltas1[2] + deltas1[3]), (deltas2[2] + deltas2[3])))
    return s


def thcov_HT(combine_by_type_ht, process_starting_points, ht_coeff):
    start_proc = process_starting_points
    process_info = combine_by_type_ht
    covmats = defaultdict(list)
    for name1 in process_info.theory:
        for name2 in process_info.theory:
            central1 = process_info.theory[name1]
            central1 = central1[1]
            kin1_1 = process_info.data[name1][:, 0]
            kin2_1 = process_info.data[name1][:, 1]
            central2 = process_info.theory[name2]
            central2 = central2[1]
            kin1_2 = process_info.data[name2][:, 0]
            kin2_2 = process_info.data[name2][:, 1]
            deltas1 = central1 * ht_coeff / kin2_1**2 / (1 - kin1_1)
            deltas2 = central2 * ht_coeff / kin2_2**2 / (1 - kin1_2)
            s = np.outer(deltas1, deltas2)
            start_locs = (start_proc[name1], start_proc[name2])
            covmats[start_locs] = s
    return covmats


def covmat_n3lo_singlet(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for all the
    singlet splitting function variations.
    """
    n3lo_vars_dict = {"gg": 19, "gq": 21, "qg": 15, "qq": 6}
    s_singlet_ad = 0
    cnt = 0
    for n_var in n3lo_vars_dict.values():
        s_singlet_ad += covmat_n3lo_ad(
            name1, name2, deltas1[cnt : cnt + n_var], deltas2[cnt : cnt + n_var]
        )
        cnt += n_var
    return s_singlet_ad

def covmat_n3lo_fhmv(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for all the
    FHMV splitting function variations.
    """
    s_ad = 0
    n_var = 2
    # loop on the 7 splitting functions variations
    for cnt in range(0,14,2):
        s_ad += covmat_n3lo_ad(
            name1, name2, deltas1[cnt : cnt + n_var], deltas2[cnt : cnt + n_var]
        )
    return s_ad

def covmat_n3lo_ad(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for each of the
    singlet splitting function variations.

    Normalization is given by:
        (n_pt - 1)

    where:
        * n_pt = number of point presctiption
    """
    norm = len(deltas1)
    if name1 == name2:
        s = sum(np.outer(d, d) for d in deltas1)
    else:
        s = 0
        for i, d1 in enumerate(deltas1):
            for j, d2 in enumerate(deltas2):
                if i == j:
                    s += np.outer(d1, d2)
    return 1 / norm * s


def compute_covs_pt_prescrip(
    point_prescription,
    l,
    name1,
    deltas1,
    name2=None,
    deltas2=None,
    fivetheories=None,
    seventheories=None,
):
    """Utility to compute the covariance matrix by prescription given the
    shifts with respect to the central value for a pair of processes.

    The processes are defined by the variables ``name1`` and ``name2`` with
    ``deltas1`` and ``deltas2`` the associated shifts wrt the central prediction.

    This utility also allows for the computation of the theory covmat for a single process
    or dataset if ``names2=deltas2=None``.

    Parameters
    ---------
        point_prescription: str
            defines the point prescription to be utilized
        l: int
            Number of theory variations (counting the central theory)
        name1: str
            Process name of the first set of shifts
        deltas1: list(np.ndarray)
            list of shifts for each of the non-central theories
        name2: str
            Process name of the second set of shifts
        deltas2: list(np.ndarray)
            list of shifts for each of the non-central theories
        fivetheories: str
            5-point prescription variation
        seventheories: str
            7-point prescription variation
    """
    if name2 is None and deltas2 is not None:
        raise ValueError(
            f"Error building theory covmat: predictions have been given with no associated process/dataset name"
        )
    elif deltas2 is None and name2 is not None:
        raise ValueError(
            f"Error building theory covmat: a process/dataset name has been given {name2} with no predictions"
        )

    if name2 is None:
        name2 = name1
        deltas2 = deltas1

    if l == 3:
        if point_prescription == "3f point":
            s = covmat_3fpt(name1, name2, deltas1, deltas2)
        elif point_prescription == "3r point":
            s = covmat_3rpt(name1, name2, deltas1, deltas2)
        else:
            s = covmat_3pt(name1, name2, deltas1, deltas2)
    elif l == 5:
        # 5 point --------------------------------------------------------------
        if fivetheories == "nobar":
            s = covmat_5pt(name1, name2, deltas1, deltas2)
        # 5bar-point -----------------------------------------------------------
        else:
            s = covmat_5barpt(name1, name2, deltas1, deltas2)
    #  ---------------------------------------------------------------------
    elif l == 7:
        # Outdated 7pts implementation: left for posterity ---------------------
        if seventheories == "original":
            s = covmat_7pt_orig(name1, name2, deltas1, deltas2)
        # 7pt (Gavin) ----------------------------------------------------------
        else:
            s = covmat_7pt(name1, name2, deltas1, deltas2)
    elif l == 9:
        s = covmat_9pt(name1, name2, deltas1, deltas2)
    # n3lo ad variation prescriprion
    elif l == 62:
        s = covmat_n3lo_singlet(name1, name2, deltas1, deltas2)
    # n3lo ihou prescriprion
    elif l == 64:
        s_ad = covmat_n3lo_singlet(name1, name2, deltas1[:-2], deltas2[:-2])
        s_cf = covmat_3pt(name1, name2, deltas1[-2:], deltas2[-2:])
        s = s_ad + s_cf
    # n3lo 3 pt MHOU see also
    # see https://github.com/NNPDF/papers/blob/e2ac1832cf4a36dab83a696564eaa75a4e55f5d2/minutes/minutes-2023-08-18.txt#L148-L157
    elif l == 66:
        s_ad = covmat_n3lo_singlet(name1, name2, deltas1[:-4], deltas2[:-4])
        s_mhou = covmat_3pt(name1, name2, deltas1[-4:-2], deltas2[-4:-2])
        s_cf = covmat_3pt(name1, name2, deltas1[-2:], deltas2[-2:])
        s = s_ad + s_cf + s_mhou
    # n3lo full covmat prescriprion
    elif l == 70:
        # spit deltas and compose thcovmat
        # splittin functions variatons
        s_ad = covmat_n3lo_singlet(name1, name2, deltas1[:-8], deltas2[:-8])
        # scale variations
        s_mhou = covmat_7pt(name1, name2, deltas1[-8:-2], deltas2[-8:-2])
        # massive coefficient function variations
        s_cf = covmat_3pt(name1, name2, deltas1[-2:], deltas2[-2:])
        s = s_ad + s_cf + s_mhou
    elif l == 19:
        s_ad = covmat_n3lo_fhmv(name1, name2, deltas1[:-4], deltas2[:-4])
        s_mhou = covmat_3pt(name1, name2, deltas1[-4:-2], deltas2[-4:-2])
        s_cf = covmat_3pt(name1, name2, deltas1[-2:], deltas2[-2:])
        s = s_ad + s_cf + s_mhou
    # n3lo full covmat prescriprion
    elif l == 23:
        # spit deltas and compose thcovmat
        # splitting functions variatons
        s_ad = covmat_n3lo_fhmv(name1, name2, deltas1[:-8], deltas2[:-8])
        # scale variations
        s_mhou = covmat_7pt(name1, name2, deltas1[-8:-2], deltas2[-8:-2])
        # massive coefficient function variations
        s_cf = covmat_3pt(name1, name2, deltas1[-2:], deltas2[-2:])
        s = s_ad + s_cf + s_mhou
    return s


@check_correct_theory_combination
def covs_pt_prescrip(
    combine_by_type,
    process_starting_points,
    theoryids,
    point_prescription,
    fivetheories,
    seventheories,
):
    """Produces the sub-matrices of the theory covariance matrix according
    to a point prescription which matches the number of input theories.
    If 5 theories are provided, a scheme 'bar' or 'nobar' must be
    chosen in the runcard in order to specify the prescription. Sub-matrices
    correspond to applying the scale variation prescription to each pair of
    processes in turn, using a different procedure for the case where the
    processes are the same relative to when they are different."""
    l = len(theoryids)
    start_proc = process_starting_points
    process_info = combine_by_type
    covmats = defaultdict(list)
    for name1 in process_info.theory:
        for name2 in process_info.theory:
            central1, *others1 = process_info.theory[name1]
            deltas1 = list(other - central1 for other in others1)
            central2, *others2 = process_info.theory[name2]
            deltas2 = list(other - central2 for other in others2)
            s = compute_covs_pt_prescrip(
                point_prescription, l, name1, deltas1, name2, deltas2, fivetheories, seventheories
            )
            start_locs = (start_proc[name1], start_proc[name2])
            covmats[start_locs] = s
    return covmats


@table
def theory_covmat_custom(covmat_custom, covmap, procs_index):
    """Takes the individual sub-covmats between each two processes and assembles
    them into a full covmat. Then reshuffles the order from ordering by process
    to ordering by experiment as listed in the runcard"""
    matlength = int(
        sum([len(covmat) for covmat in covmat_custom.values()]) / int(np.sqrt(len(covmat_custom)))
    )
    # Initialise arrays of zeros and set precision to same as FK tables
    mat = np.zeros((matlength, matlength), dtype=np.float32)
    cov_by_exp = np.zeros((matlength, matlength), dtype=np.float32)
    for locs in covmat_custom:
        cov = covmat_custom[locs]
        mat[locs[0] : (len(cov) + locs[0]), locs[1] : (len(cov.T) + locs[1])] = cov
    for i in range(matlength):
        for j in range(matlength):
            cov_by_exp[covmap[i]][covmap[j]] = mat[i][j]
    df = pd.DataFrame(cov_by_exp, index=procs_index, columns=procs_index)
    return df


@table
def fromfile_covmat(covmatpath, procs_data, procs_index):
    """Reads a general theory covariance matrix from file. Then
    1: Applies cuts to match experiment covariance matrix
    2: Expands dimensions to match experiment covariance matrix
       by filling additional entries with 0."""
    # Load covmat as pandas DataFrame
    filecovmat = pd.read_csv(
        covmatpath, index_col=[0, 1, 2], header=[0, 1, 2], sep="\t|,", engine="python"
    )
    # Remove string in column id
    filecovmat.columns = filecovmat.index
    # Reordering covmat to match exp order in runcard
    # Datasets in exp covmat
    dslist = []
    for group in procs_data:
        for ds in group.datasets:
            dslist.append(ds.name)
    # Datasets in filecovmat in exp covmat order
    shortlist = []
    for ds in dslist:
        if ds in filecovmat.index.get_level_values(level="dataset"):
            shortlist.append(ds)
    filecovmat = filecovmat.reindex(shortlist, level="dataset")
    filecovmat = filecovmat.reindex(shortlist, level="dataset", axis=1)
    # ------------- #
    # 1: Apply cuts #
    # ------------- #
    # Loading cuts to apply to covariance matrix
    indextuples = []
    for group in procs_data:
        for ds in group.datasets:
            # Load cuts for each dataset in the covmat
            if ds.name in filecovmat.index.get_level_values(1):
                cuts = ds.cuts
                # Creating new index for post cuts
                lcuts = cuts.load()
                if lcuts is not None:
                    for keeploc in lcuts:
                        indextuples.append((group.name, ds.name, keeploc))
                # If loaded cuts are None, keep all index points for that dataset
                else:
                    for ind in filecovmat.index():
                        if ind[1] == ds:
                            indextuples.append(ind)
    newindex = pd.MultiIndex.from_tuples(
        indextuples, names=["group", "dataset", "index"], sortorder=0
    )
    # Reindex covmat with the new cut index
    cut_df = filecovmat.reindex(newindex).T
    cut_df = cut_df.reindex(newindex).T
    # Elements where cuts are applied will become NaN - remove these rows and columns
    cut_df = cut_df.dropna(axis=0).dropna(axis=1)
    # -------------------- #
    # 2: Expand dimensions #
    # -------------------- #
    # First make empty df of exp covmat dimensions
    empty_df = pd.DataFrame(0, index=procs_index, columns=procs_index)
    covmats = []
    # Make a piece of the covmat for each combination of two datasetes
    for ds1 in dslist:
        for ds2 in dslist:
            if (ds1 in shortlist) and (ds2 in shortlist):
                # If both datasets in the fromfile covmat, use the piece of the fromfile covmat
                covmat = (
                    cut_df.xs(ds1, level=1, drop_level=False).T.xs(ds2, level=1, drop_level=False).T
                )
            else:
                # Otherwise use a covmat of 0s
                covmat = (
                    empty_df.xs(ds1, level=1, drop_level=False)
                    .T.xs(ds2, level=1, drop_level=False)
                    .T
                )
            covmats.append(covmat)
    chunks = []
    # Arrange into chunks, each chunk is a list of pieces of covmat which are associated with
    # one dataset in particular
    for x in range(0, len(covmats), len(dslist)):
        chunk = covmats[x : x + len(dslist)]
        chunks.append(chunk)
    strips = []
    # Concatenate each chunk into a strip of the covariance matrix
    for chunk in chunks:
        strip = pd.concat(chunk, axis=1)
        strips.append(strip.T)
    # strips.reverse()
    # Concatenate the strips to make the full matrix
    full_df = pd.concat(strips, axis=1)
    # Reindex to align with experiment covmat index
    full_df = full_df.reindex(procs_index)
    full_df = full_df.reindex(procs_index, axis=1)
    return full_df


@table
def user_covmat(procs_data, procs_index, loaded_user_covmat_path):
    """
    General theory covariance matrix provided by the user.
    Useful for testing the impact of externally produced
    covariance matrices. Matrices must be produced as a
    csv of pandas DataFrame, and uploaded to the validphys
    server. The server path is then provided via
    ``user_covmat_path`` in ``theorycovmatconfig`` in the
    runcard. For more information see documentation.
    """
    return fromfile_covmat(loaded_user_covmat_path, procs_data, procs_index)


@table
@check_fit_dataset_order_matches_grouped
def total_theory_covmat(theory_covmat_custom, user_covmat):
    """
    Sum of scale variation and user covmat, where both are used.
    """
    return theory_covmat_custom + user_covmat


def theory_covmat_custom_fitting(theory_covmat_custom, procs_index_matched):
    """theory_covmat_custom but reindexed so the order of the datasets matches
    those in the experiment covmat so they are aligned when fitting."""
    df = theory_covmat_custom.reindex(procs_index_matched).T.reindex(procs_index_matched)
    return df


def total_theory_covmat_fitting(total_theory_covmat, procs_index_matched):
    """total_theory_covmat but reindexed so the order of the datasets matches
    those in the experiment covmat so they are aligned when fitting."""
    return theory_covmat_custom_fitting(total_theory_covmat, procs_index_matched)


def user_covmat_fitting(user_covmat, procs_index_matched):
    """user_covmat but reindexed so the order of the datasets matches
    those in the experiment covmat so they are aligned when fitting."""
    return theory_covmat_custom_fitting(user_covmat, procs_index_matched)


def procs_index_matched(groups_index, procs_index):
    """procs_index but matched to the dataset order given
    by groups_index."""
    # Making list with exps ordered like in groups_index
    groups_ds_order = groups_index.get_level_values(level=1).unique().tolist()
    # Tuples to make multiindex, ordered like in groups_index
    tups = []
    for ds in groups_ds_order:
        for orig in procs_index:
            if orig[1] == ds:
                tups.append(orig)

    return pd.MultiIndex.from_tuples(tups, names=("process", "dataset", "id"))


@check_correct_theory_combination
def total_covmat_diagtheory_procs(procs_results_theory, fivetheories):
    """Same as total_covmat_datasets but per proc rather than
    per dataset. Needed for calculation of chi2 per proc."""
    exp_result_covmats = []
    for exp_result in zip(*procs_results_theory):
        theory_centrals = [x[1].central_value for x in exp_result]
        s = make_scale_var_covmat(theory_centrals)
        # Initialise array of zeros and set precision to same as FK tables
        s_diag = np.zeros((len(s), len(s)), dtype=np.float32)
        np.fill_diagonal(s_diag, np.diag(s))
        sigma = exp_result[0][0].covmat
        cov = s_diag + sigma
        exp_result_covmats.append(cov)
    return exp_result_covmats


@table
def theory_corrmat_singleprocess(theory_covmat_singleprocess):
    """Calculates the theory correlation matrix for scale variations."""
    df = theory_covmat_singleprocess
    covmat = df.values
    diag_minus_half = (np.diagonal(covmat)) ** (-0.5)
    mat = diag_minus_half[:, np.newaxis] * df * diag_minus_half
    return mat


@table
def theory_blockcorrmat(theory_block_diag_covmat):
    """Calculates the theory correlation matrix for scale variations
    with block diagonal entries by dataset only"""
    mat = theory_corrmat_singleprocess(theory_block_diag_covmat)
    return mat


@table
def theory_corrmat_custom(theory_covmat_custom):
    """Calculates the theory correlation matrix for scale variations
    with variations by process type"""
    mat = theory_corrmat_singleprocess(theory_covmat_custom)
    return mat


@table
def theory_normcovmat_singleprocess(theory_covmat_singleprocess, procs_data_values):
    """Calculates the theory covariance matrix for scale variations normalised
    to data."""
    df = theory_covmat_singleprocess
    mat = df / np.outer(procs_data_values, procs_data_values)
    return mat


@table
def theory_normblockcovmat(theory_block_diag_covmat, procs_data_values):
    """Calculates the theory covariance matrix for scale variations
    normalised to data, block diagonal by dataset."""
    df = theory_block_diag_covmat
    mat = df / np.outer(procs_data_values, procs_data_values)
    return mat


@table
def theory_normcovmat_custom(theory_covmat_custom, procs_data_values):
    """Calculates the theory covariance matrix for scale variations normalised
    to data, with variations according to the relevant prescription."""
    df = theory_covmat_custom
    mat = df / np.outer(procs_data_values, procs_data_values)
    return mat


@table
def experimentplustheory_covmat_singleprocess(
    procs_covmat_no_table, theory_covmat_singleprocess_no_table
):
    """Calculates the experiment + theory covariance matrix for
    scale variations."""
    df = procs_covmat_no_table + theory_covmat_singleprocess_no_table
    return df


@table
def experimentplusblocktheory_covmat(procs_covmat, theory_block_diag_covmat):
    """Calculates the experiment + theory covariance
    matrix for scale variations."""
    df = procs_covmat + theory_block_diag_covmat
    return df


@table
def experimentplustheory_covmat_custom(procs_covmat, theory_covmat_custom):
    """Calculates the experiment + theory covariance matrix for
    scale variations correlated according to the relevant prescription."""
    df = procs_covmat + theory_covmat_custom
    return df


@table
def experimentplustheory_normcovmat_singleprocess(
    procs_covmat, theory_covmat_singleprocess, procs_data
):
    """Calculates the experiment + theory covariance matrix for scale
    variations normalised to data."""
    df = procs_covmat + theory_covmat_singleprocess
    procs_data_array = np.array(procs_data)
    mat = df / np.outer(procs_data_array, procs_data_array)
    return mat


@table
def experimentplusblocktheory_normcovmat(
    procs_covmat, theory_block_diag_covmat, procs_data_values, experimentplustheory_normcovmat
):
    """Calculates the experiment + theory covariance matrix for scale
    variations normalised to data, block diagonal by data set."""
    mat = experimentplustheory_normcovmat(procs_covmat, theory_block_diag_covmat, procs_data_values)
    return mat


@table
def experimentplustheory_normcovmat_custom(
    procs_covmat, theory_covmat_custom, procs_data_values, experimentplustheory_normcovmat
):
    """Calculates the experiment + theory covariance matrix for scale
    variations normalised to data, correlations by process type."""
    mat = experimentplustheory_normcovmat(procs_covmat, theory_covmat_custom, procs_data_values)

    return mat


@table
def experimentplustheory_corrmat_singleprocess(procs_covmat, theory_covmat_singleprocess):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices."""
    total_df = procs_covmat + theory_covmat_singleprocess
    total_cov = (procs_covmat + theory_covmat_singleprocess).values
    diag_minus_half = (np.diagonal(total_cov)) ** (-0.5)
    corrmat = diag_minus_half[:, np.newaxis] * total_df * diag_minus_half
    return corrmat


@table
def experimentplusblocktheory_corrmat(procs_covmat, theory_block_diag_covmat):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices, block diagonal by dataset."""
    corrmat = experimentplustheory_corrmat_singleprocess(procs_covmat, theory_block_diag_covmat)
    return corrmat


@table
def experimentplustheory_corrmat_custom(procs_covmat, theory_covmat_custom):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices, correlations by prescription."""
    corrmat = experimentplustheory_corrmat_singleprocess(procs_covmat, theory_covmat_custom)
    return corrmat


def chi2_impact(theory_covmat_singleprocess, procs_covmat, procs_results):
    """Returns total chi2 including theory cov mat"""
    dataresults, theoryresults = zip(*procs_results)
    dat_central_list = [x.central_value for x in dataresults]
    th_central_list = [x.central_value for x in theoryresults]
    dat_central = np.concatenate(dat_central_list)
    th_central = np.concatenate([x for x in th_central_list])
    central_diff = dat_central - th_central
    cov = theory_covmat_singleprocess.values + procs_covmat.values
    return calc_chi2(la.cholesky(cov, lower=True), central_diff) / len(central_diff)


def data_theory_diff(procs_results):
    """Returns (D-T) for central theory, for use in chi2 calculations"""
    dataresults, theoryresults = zip(*procs_results)
    dat_central_list = [x.central_value for x in dataresults]
    th_central_list = [x.central_value for x in theoryresults]
    dat_central = np.concatenate(dat_central_list)
    th_central = np.concatenate(th_central_list)
    central_diff = dat_central - th_central
    return central_diff


def chi2_block_impact(theory_block_diag_covmat, procs_covmat, procs_results):
    """Returns total chi2 including theory cov mat"""
    chi2 = chi2_impact(theory_block_diag_covmat, procs_covmat, procs_results)
    return chi2


def chi2_impact_custom(theory_covmat_custom, procs_covmat, procs_results):
    """Returns total chi2 including theory cov mat"""
    chi2 = chi2_impact(theory_covmat_custom, procs_covmat, procs_results)
    return chi2


def theory_diagcovmat(theory_covmat_singleprocess):
    """Returns theory covmat with only diagonal values"""
    s = theory_covmat_singleprocess.values
    # Initialise array of zeros and set precision to same as FK tables
    s_diag = np.zeros((len(s), len(s)), dtype=np.float32)
    np.fill_diagonal(s_diag, np.diag(s))
    return s_diag


def chi2_diag_only(theory_diagcovmat, procs_covmat, data_theory_diff):
    """Returns total chi2 including only diags of theory cov mat"""
    cov = theory_diagcovmat + procs_covmat.values
    elements = np.dot(data_theory_diff.T, np.dot(la.inv(cov), data_theory_diff))
    chi2 = (1 / len(data_theory_diff)) * np.sum(elements)
    return chi2


each_dataset_results = collect(results, ("group_dataset_inputs_by_process", "data"))


def abs_chi2_data_theory_dataset(each_dataset_results, total_covmat_datasets):
    """Returns an array of tuples (member_chi², central_chi², numpoints)
    corresponding to each data set, where theory errors are included"""
    chi2data_array = []
    for datresults, covmat in zip(each_dataset_results, total_covmat_datasets):
        data_result, th_result = datresults
        chi2s = all_chi2_theory(datresults, covmat)
        central_result = central_chi2_theory(datresults, covmat)
        chi2data_array.append(
            Chi2Data(th_result.stats_class(chi2s[:, np.newaxis]), central_result, len(data_result))
        )
    return chi2data_array


def abs_chi2_data_theory_proc(procs_results, total_covmat_procs):
    """Like abs_chi2_data_theory_dataset but for procs not datasets"""
    chi2data_array = []
    for expresults, covmat in zip(procs_results, total_covmat_procs):
        data_result, th_result = expresults
        chi2s = all_chi2_theory(expresults, covmat)
        central_result = central_chi2_theory(expresults, covmat)
        chi2data_array.append(
            Chi2Data(th_result.stats_class(chi2s[:, np.newaxis]), central_result, len(data_result))
        )
    return chi2data_array


def abs_chi2_data_diagtheory_proc(procs_results, total_covmat_diagtheory_procs):
    """For a diagonal theory covmat"""
    return abs_chi2_data_theory_proc(procs_results, total_covmat_diagtheory_procs)


def abs_chi2_data_diagtheory_dataset(each_dataset_results, total_covmat_diagtheory_datasets):
    """For a diagonal theory covmat"""
    return abs_chi2_data_theory_dataset(each_dataset_results, total_covmat_diagtheory_datasets)
