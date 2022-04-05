# -*- coding: utf-8 -*-
"""
construction.py
Tools for constructing theory covariance matrices and computing their chi2s.
"""
from __future__ import generator_stop

import logging

from collections import defaultdict, namedtuple
import numpy as np
import scipy.linalg as la
import pandas as pd

from reportengine.table import table
from reportengine import collect

from validphys.results import (
    procs_central_values,
    procs_central_values_no_table,
)
from validphys.results import Chi2Data, results
from validphys.calcutils import calc_chi2, all_chi2_theory, central_chi2_theory
from validphys.theorycovariance.theorycovarianceutils import (
    process_lookup,
    check_correct_theory_combination,
    check_fit_dataset_order_matches_grouped,
)


log = logging.getLogger(__name__)

theoryids_procs_central_values = collect(procs_central_values, ("theoryids",))

theoryids_procs_central_values_no_table = collect(
    procs_central_values_no_table, ("theoryids",)
)

collected_theoryids = collect("theoryids", ["theoryconfig",])


def make_scale_var_covmat(predictions):
    """Takes N theory predictions at different scales and applies N-pt scale
    variations to produce a covariance matrix."""
    l = len(predictions)
    central, *others = predictions
    deltas = (other - central for other in others)
    if l == 3:
        norm = 0.5
    elif l == 5:
        norm = 0.5
    elif l == 7:
        norm = 1 / 3
    elif l == 9:
        norm = 0.25
    s = norm * sum(np.outer(d, d) for d in deltas)
    return s


@check_correct_theory_combination
def theory_covmat_singleprocess_no_table(
    theoryids_procs_central_values_no_table, procs_index, theoryids, fivetheories,
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


ProcessInfo = namedtuple("ProcessInfo", ("theory", "namelist", "sizes"))


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
        theory=theories_by_process, namelist=ordered_names, sizes=dataset_size
    )
    return process_info


def process_starting_points(combine_by_type):
    """Returns a dictionary of indices in the covariance matrix corresponding
    to the starting point of each process."""
    process_info = combine_by_type
    running_index = 0
    start_proc = defaultdict(list)
    for name in process_info.theory:
        size = len(process_info.theory[name][0])
        start_proc[name] = running_index
        running_index += size
    return start_proc


def covmap(combine_by_type, dataset_names):
    """Creates a map between the covmat indices from matrices ordered by
    process to matrices ordered by experiment as listed in the runcard"""
    mapping = defaultdict(list)
    start_exp = defaultdict(list)
    process_info = combine_by_type
    running_index = 0
    for dataset in dataset_names:
        size = process_info.sizes[dataset]
        start_exp[dataset] = running_index
        running_index += size
    start = 0
    names_by_proc_list = [
        item for sublist in process_info.namelist.values() for item in sublist
    ]
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
        s = 0.5 * (
            np.outer(deltas1[0], deltas2[0]) + np.outer(deltas1[1], deltas2[1])
        ) + 0.25 * (np.outer((deltas1[2] + deltas1[3]), (deltas2[2] + deltas2[3])))
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
            np.outer(
                (deltas1[0] + deltas1[4] + deltas1[6]),
                (deltas2[0] + deltas2[4] + deltas2[6]),
            )
            + np.outer(
                (deltas1[1] + deltas1[5] + deltas1[7]),
                (deltas2[1] + deltas2[5] + deltas2[7]),
            )
        ) + (1 / 8) * (np.outer((deltas1[2] + deltas1[3]), (deltas2[2] + deltas2[3])))
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
    batches_list = []
    for name1 in process_info.theory:
        central, *others = process_info.theory[name1]
        deltas = np.array(list((other - central for other in others)))
        zeros = np.array(list(cen - cen for cen in central))
        fact05 = np.array([deltas[0], deltas[1], deltas[2]])
        fact1 = np.array([deltas[3], zeros, deltas[4]])
        fact2 = np.array([deltas[5], deltas[6], deltas[7]])
        shifts_fin = np.array([fact05,fact1,fact2])
        batches_list.append(np.transpose(shifts_fin,axes=[2,1,0]))
    shifts = shifts_vec(batches_list)
    thcmat = thcovmat(shifts) 
    return thcmat

def shifts_vec(raw: list[np.ndarray]) -> list[np.ndarray]:
    """Pump more dimensions into shifts.

    Since renormalization scale for each process has a different meaning, the
    most intrinsic and intuitive way to represent this difference is to put them
    on different dimensions:

    - a renormalization dimension is non-trivial only for the process that scale
      is relative to
    - a trivial dimension has the meaning of not affecting that datum
    - trivial dimension can be exploited by `numpy` with the broadcasting
      semantic

    Parameters
    ----------
    raw : list[np.ndarray]
        sequence of raw shifts, of dimension ``(n,3,3)``, like those generated
        by :func:`raw_shifts`

    Returns
    -------
    list[np.ndarray]
        a list of batches, blown up with new dimensions to separate on different
        dimensions the renormalization scales related to different processes

    """
    upgraded = []
    n = len(raw)
    dims = list(np.arange(n) + 2)
    for i, batch in enumerate(raw):
        newdims = dims.copy()
        newdims.remove(i + 2)
        upgraded.append(np.expand_dims(batch, newdims))

    return upgraded


def thcovmat(shifts: list[np.ndarray]) -> np.ndarray:
    """Generate theory covariance matrix from upgraded vector of shifts.

    Exploit `numpy` broadcasting to apply Eq.(4.2) of arXiv:1906.10698
    literally.

    The prescription used is always the 9-point one, so the other prescriptions
    have to be implemented on the input and the output of this funcions, i.e.:

    - the input has to be masked, to replace elements that should be missing
      with 0s (it is trivial that the effect is the same)
    - an overall normalization

    Note
    ----
    In order to apply exactly Eq.(4.2), all the renormalization scale dimensions
    should be non-trivial, and the broadcasting semantics should be implemented
    for all dimensions and on the vector itself, such that a single rectangular
    `np.ndarray` is obtained.

    This of course is very inefficient for memory and operations, since
    :math:`n-2` renormalization scales are always not involved while generating
    a block for 2 given processes.

    So for the :math:`n-2` scales it would correspond to an overall degeneracy
    factor that would simplify with normalization, so we can simply skip (in
    `numpy` this is done by the contraction on a dimension of size 1 for both
    arrays, resulting in the removal of that dimension).

    The remaining scales are at most 2, but not always 2, changing from
    on-diagonal or off-diagonal blocks. The difference has to be compensated
    with as a further degeneracy factor for the on-diagonal blocsk.

    A minimal example is the case of only to process: there would be two shift
    arrays, one each, of shapes ``[n1, nf, nr, 1]`` and ``[n2, nf, 1, nr]``. The
    internal dimensions on data is irrelevant, so we can drop, as well as the
    factorization scale one (that is always the same).
    For definiteness let's consider ``nr = 3``, so we'll have two shifts arrays
    of shapes ``[3,1]`` and ``[1, 3]``.

    - off-diagonal: the contraction is done with ``[3,1]`` and ``[1,3]``, and
      both the 1 dimensions are broadcasted to 3 before contraction, leading to
      a 9 elements contraction
    - on-diagonal: a contraction is ``[3, 1]`` with ``[3, 1]``, so the second
      dimension is not broadcasted, leading to only 3 elements, but it should be
      according to the meaning of Eq.(4.2); for this reason, a degeneracy
      factor of 3 is multiplied, in order to make it homogeneous with the
      off-diagonal blocks


    Parameters
    ----------
    raw : list[np.ndarray]
        sequence of upgraded shifts, like those generated by :func:`shifts_vec`

    Returns
    -------
    np.ndarray
        matrix of "covariances" generated out of shifts, whose shape is ``(N,
        N)``, where ``N`` is the number of all data points (i.e. ``N =
        sum(len(s) for s in shifts)``)

    Raises
    ------
    ValueError
      if not all the processes have the same number of points for their own
      renormalization scales

    """
    blockmat = []
    sumdims = list(np.arange(len(shifts) + 1) + 2)

    # for each process there is only one non-trivial renormalization scale, so
    # the other dimensions have to be 1
    murs = np.array(shifts[0].shape[2:]).prod()
    if not all(np.array(proc_shift.shape[2:]).prod() == murs for proc_shift in shifts):
        raise ValueError(
            "All the different renormalization scales should have the"
            " same number of points"
        )

    for bi in shifts:
        blockrow = []
        for bj in shifts:
            degeneracy = murs if bi is bj else 1.0
            blockrow.append(
                np.einsum(bi, [0, *sumdims], bj, [1, *sumdims], [0, 1]) * degeneracy
            )
        blockmat.append(blockrow)

    return np.block(blockmat)
           


@table
def theory_covmat_custom(covs_pt_prescrip, covmap, procs_index):
    """Takes the individual sub-covmats between each two processes and assembles
    them into a full covmat. Then reshuffles the order from ordering by process
    to ordering by experiment as listed in the runcard"""
   # matlength = int(
   #     sum([len(covmat) for covmat in covs_pt_prescrip.values()])
   #     / int(np.sqrt(len(covs_pt_prescrip)))
   # )
    # Initialise arrays of zeros and set precision to same as FK tables
    #mat = np.zeros((matlength, matlength), dtype=np.float32)
    matlength = covs_pt_prescrip.shape[0]
    cov_by_exp = np.zeros((matlength, matlength), dtype=np.float64)
    for i in range(matlength):
        for j in range(matlength):
            cov_by_exp[covmap[i]][covmap[j]] = covs_pt_prescrip[i][j]
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
    cut_df = cut_df.dropna(0).dropna(1)
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
                    cut_df.xs(ds1, level=1, drop_level=False)
                    .T.xs(ds2, level=1, drop_level=False)
                    .T
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
    df = theory_covmat_custom.reindex(procs_index_matched).T.reindex(
        procs_index_matched
    )
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
    by groups_index. """
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
    procs_covmat,
    theory_block_diag_covmat,
    procs_data_values,
    experimentplustheory_normcovmat,
):
    """Calculates the experiment + theory covariance matrix for scale
       variations normalised to data, block diagonal by data set."""
    mat = experimentplustheory_normcovmat(
        procs_covmat, theory_block_diag_covmat, procs_data_values
    )
    return mat


@table
def experimentplustheory_normcovmat_custom(
    procs_covmat,
    theory_covmat_custom,
    procs_data_values,
    experimentplustheory_normcovmat,
):
    """Calculates the experiment + theory covariance matrix for scale
       variations normalised to data, correlations by process type."""
    mat = experimentplustheory_normcovmat(
        procs_covmat, theory_covmat_custom, procs_data_values
    )

    return mat


@table
def experimentplustheory_corrmat_singleprocess(
    procs_covmat, theory_covmat_singleprocess
):
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
    corrmat = experimentplustheory_corrmat_singleprocess(
        procs_covmat, theory_block_diag_covmat
    )
    return corrmat


@table
def experimentplustheory_corrmat_custom(procs_covmat, theory_covmat_custom):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices, correlations by prescription."""
    corrmat = experimentplustheory_corrmat_singleprocess(
        procs_covmat, theory_covmat_custom
    )
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
            Chi2Data(
                th_result.stats_class(chi2s[:, np.newaxis]),
                central_result,
                len(data_result),
            )
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
            Chi2Data(
                th_result.stats_class(chi2s[:, np.newaxis]),
                central_result,
                len(data_result),
            )
        )
    return chi2data_array


def abs_chi2_data_diagtheory_proc(procs_results, total_covmat_diagtheory_procs):
    """For a diagonal theory covmat"""
    return abs_chi2_data_theory_proc(procs_results, total_covmat_diagtheory_procs)


def abs_chi2_data_diagtheory_dataset(
    each_dataset_results, total_covmat_diagtheory_datasets
):
    """For a diagonal theory covmat"""
    return abs_chi2_data_theory_dataset(
        each_dataset_results, total_covmat_diagtheory_datasets
    )
