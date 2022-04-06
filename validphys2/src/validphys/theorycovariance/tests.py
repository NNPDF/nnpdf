# -*- coding: utf-8 -*-
"""
tests.py
Tools for testing theory covariance matrices and their properties.
"""
from __future__ import generator_stop

import logging

from collections import namedtuple
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

from reportengine.figure import figure
from reportengine.table import table
from reportengine import collect
from reportengine import floatformatting

from validphys.checks import check_two_dataspecs

from validphys.theorycovariance.construction import (
    combine_by_type,
    process_starting_points,
)
from validphys.theorycovariance.construction import theory_corrmat_singleprocess
from validphys.theorycovariance.construction import (
    covmap,
    covs_pt_prescrip,
    theory_covmat_custom,
)

from validphys.theorycovariance.output import matrix_plot_labels, _get_key
from validphys.theorycovariance.theorycovarianceutils import (
    process_lookup,
    check_correct_theory_combination_theoryconfig,
)
from validphys.theorycovariance.theorycovarianceutils import (
    check_correct_theory_combination_dataspecs,
)

log = logging.getLogger(__name__)

matched_dataspecs_results = collect("results", ["dataspecs"])

LabeledShifts = namedtuple("LabeledShifts", ("process", "dataset_name", "shifts"))


@check_two_dataspecs
def dataspecs_dataset_prediction_shift(
    matched_dataspecs_results, process, dataset_name
):
    """Compute the difference in theory predictions between two dataspecs.
    This can be used in combination with `matched_datasets_from_dataspecs`
    It returns a ``LabeledShifts`` containing ``dataset_name``,
    ``process`` and ``shifts``.
    """
    r1, r2 = matched_dataspecs_results
    res = r1[1].central_value - r2[1].central_value
    return LabeledShifts(dataset_name=dataset_name, process=process, shifts=res)


matched_dataspecs_dataset_prediction_shift = collect(
    "dataspecs_dataset_prediction_shift", ["dataspecs"]
)


def shift_vector(
    matched_dataspecs_dataset_prediction_shift, matched_dataspecs_dataset_theory
):
    """Returns a DataFrame of normalised shift vectors for matched dataspecs."""
    all_shifts = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_prediction_shift]
    )
    all_theory = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_theory]
    )
    norm_shifts = all_shifts / all_theory
    dsnames = np.concatenate(
        [
            np.full(len(val.shifts), val.dataset_name, dtype=object)
            for val in matched_dataspecs_dataset_prediction_shift
        ]
    )
    point_indexes = np.concatenate(
        [
            np.arange(len(val.shifts))
            for val in matched_dataspecs_dataset_prediction_shift
        ]
    )
    index = pd.MultiIndex.from_arrays(
        [dsnames, point_indexes], names=["Dataset name", "Point"]
    )
    return pd.DataFrame(norm_shifts, index=index)


def dataspecs_dataset_theory(matched_dataspecs_results, process, dataset_name):
    """Returns a tuple of shifts processed by data set and experiment
    for matched dataspecs."""
    central = matched_dataspecs_results[0]
    res = central[1].central_value
    return LabeledShifts(dataset_name=dataset_name, process=process, shifts=res)


matched_dataspecs_dataset_theory = collect("dataspecs_dataset_theory", ["dataspecs"])


def theory_vector(matched_dataspecs_dataset_theory):
    """Returns a DataFrame of the central theory vector for
    matched dataspecs."""
    all_theory = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_theory]
    )
    dsnames = np.concatenate(
        [
            np.full(len(val.shifts), val.dataset_name, dtype=object)
            for val in matched_dataspecs_dataset_theory
        ]
    )
    point_indexes = np.concatenate(
        [np.arange(len(val.shifts)) for val in matched_dataspecs_dataset_theory]
    )
    index = pd.MultiIndex.from_arrays(
        [dsnames, point_indexes], names=["Dataset name", "Point"]
    )
    return pd.DataFrame(all_theory, index=index)


def dataspecs_dataset_alltheory(matched_dataspecs_results, process, dataset_name):
    """Returns a LabeledShifts tuple corresponding to the theory
    vectors for all the scale varied theories (not the central one),
    processed by data set and experiment for matched dataspecs."""
    others = matched_dataspecs_results[1:]
    res = [other[1].central_value for other in others]
    return LabeledShifts(dataset_name=dataset_name, process=process, shifts=res)


matched_dataspecs_dataset_alltheory = collect(
    "dataspecs_dataset_alltheory", ["dataspecs"]
)


def alltheory_vector(
    matched_dataspecs_dataset_alltheory, matched_dataspecs_dataset_theory
):
    """Returns a DataFrame with the theory vectors for matched
    dataspecs for the scale-varied theories (not the central one)."""
    all_theory = np.concatenate(
        [val.shifts for val in matched_dataspecs_dataset_alltheory], axis=1
    )
    dsnames = np.concatenate(
        [
            np.full(len(val.shifts), val.dataset_name, dtype=object)
            for val in matched_dataspecs_dataset_theory
        ]
    )
    point_indexes = np.concatenate(
        [np.arange(len(val.shifts)) for val in matched_dataspecs_dataset_theory]
    )
    index = pd.MultiIndex.from_arrays(
        [dsnames, point_indexes], names=["Dataset name", "Point"]
    )
    theory_vectors = []
    for theoryvector in all_theory:
        theory_vectors.append(pd.DataFrame(theoryvector, index=index))
    return theory_vectors


all_matched_results = collect("matched_dataspecs_results", ["dataspecs"])


def combine_by_type_dataspecs(all_matched_results, matched_dataspecs_dataset_name):
    """Like combine_by_type but for matched dataspecs"""
    return combine_by_type(all_matched_results, matched_dataspecs_dataset_name)


dataspecs_theoryids = collect("theoryid", ["theoryconfig", "original", "dataspecs"])


def process_starting_points_dataspecs(combine_by_type_dataspecs):
    """Like process_starting_points but for matched dataspecs."""
    return process_starting_points(combine_by_type_dataspecs)


@check_correct_theory_combination_dataspecs
def covs_pt_prescrip_dataspecs(
    combine_by_type_dataspecs,
    process_starting_points_dataspecs,
    dataspecs_theoryids,
    point_prescription,
    fivetheories,
    seventheories,
):
    """Like covs_pt_prescrip but for matched dataspecs."""
    return covs_pt_prescrip(
        combine_by_type_dataspecs,
        process_starting_points_dataspecs,
        dataspecs_theoryids,
        point_prescription,
        fivetheories,
        seventheories,
    )


def covmap_dataspecs(combine_by_type_dataspecs, matched_dataspecs_dataset_name):
    """Like covmap but for matched dataspecs."""
    return covmap(combine_by_type_dataspecs, matched_dataspecs_dataset_name)


matched_dataspecs_process = collect("process", ["dataspecs"])
matched_dataspecs_dataset_name = collect("dataset_name", ["dataspecs"])
matched_cuts_datasets = collect("dataset", ["dataspecs"])
all_matched_datasets = collect("matched_cuts_datasets", ["dataspecs"])


def all_matched_data_lengths(all_matched_datasets):
    """Returns a list of the data sets lengths."""
    lens = []
    for rlist in all_matched_datasets:
        lens.append(rlist[0].load().GetNData())
    return lens


def matched_experiments_index(matched_dataspecs_dataset_name, all_matched_data_lengths):
    """Returns MultiIndex composed of data set name and
    starting point of data set."""
    dsnames = matched_dataspecs_dataset_name
    lens = all_matched_data_lengths
    dsnames = np.concatenate(
        [np.full(l, dsname, dtype=object) for (l, dsname) in zip(lens, dsnames)]
    )
    point_indexes = np.concatenate([np.arange(l) for l in lens])
    index = pd.MultiIndex.from_arrays(
        [dsnames, point_indexes], names=["Dataset name", "Point"]
    )
    return index


@table
def theory_covmat_custom_dataspecs(
    covs_pt_prescrip_dataspecs, covmap_dataspecs, matched_experiments_index
):
    """Like theory_covmat_custom but for matched dataspecs."""
    return theory_covmat_custom(
        covs_pt_prescrip_dataspecs, covmap_dataspecs, matched_experiments_index
    )


thx_corrmat = collect(
    "theory_corrmat_custom_dataspecs",
    ["combined_shift_and_theory_dataspecs", "theoryconfig"],
)

shx_corrmat = collect(
    "matched_datasets_shift_matrix_correlations",
    ["combined_shift_and_theory_dataspecs", "shiftconfig"],
)

thx_covmat = collect(
    "theory_covmat_custom_dataspecs",
    ["combined_shift_and_theory_dataspecs", "theoryconfig"],
)

combined_dataspecs_results = collect(
    "all_matched_results", ["combined_shift_and_theory_dataspecs", "theoryconfig"]
)

shx_vector = collect(
    "shift_vector", ["combined_shift_and_theory_dataspecs", "shiftconfig"]
)

thx_vector = collect(
    "theory_vector", ["combined_shift_and_theory_dataspecs", "theoryconfig"]
)

allthx_vector = collect(
    "alltheory_vector", ["combined_shift_and_theory_dataspecs", "theoryconfig"]
)


def theory_matrix_threshold(theory_threshold: (int, float) = 0):
    """Returns the threshold below which theory correlation elements are set to
    zero when comparing to shift correlation matrix"""
    return theory_threshold


@table
def theory_corrmat_custom_dataspecs(theory_covmat_custom_dataspecs):
    """Calculates the theory correlation matrix for scale variations
    with variations by process type"""
    mat = theory_corrmat_singleprocess(theory_covmat_custom_dataspecs)
    return mat


def _shuffle_list(l, shift):
    """Function that moves list elements left by 'shift' entries"""
    i = 0
    newlist = l.copy()
    while i <= (shift - 1):
        newlist.append(newlist.pop(0))
        i = i + 1
    return newlist


def vectors_3pt(splitdiffs):
    """Returns the linearly independent vectors for 3pt prescription"""
    # N.B. mu_0 correlated with mu_i
    xs = []
    num_pts = 3
    pps = splitdiffs[:: (num_pts - 1)]
    mms = _shuffle_list(splitdiffs, 1)[:: (num_pts - 1)]
    # Constructing (+, +, +, ...)
    xs.append(sum(pps))
    # Constructing the p vectors with one minus
    # (-, +, + ...) + cyclic
    # for each process i:
    for procloc, mm in enumerate(mms):
        # Start with a copy of one (+,+) vector for dataframe structure
        newvec = pps[0].copy()
        # Set all elements to 0 so the dataframe is empty
        newvec.loc[:] = 0
        # Copy the set of vectors [(+,+)_1, (+,+)_2, ...],
        # i.e. one vector for each process type
        subpps = pps.copy()
        # Delete (+,+)_i from the list
        del subpps[procloc]
        # Now add the corresponding (-,-)_i vector to all the other
        # (+,+) vectors
        # This gives one "(-, +, + ...) & cyclic" vector
        newvec = newvec + sum(subpps) + mm
        # Append this vector to the list of vectors
        xs.append(newvec)
    return xs


def vectors_5pt(splitdiffs):
    """Returns the linearly independent vectors for 5pt prescription"""
    num_pts = 5
    pzs = splitdiffs[:: (num_pts - 1)]
    mzs = _shuffle_list(splitdiffs, 1)[:: (num_pts - 1)]
    zps = _shuffle_list(splitdiffs, 2)[:: (num_pts - 1)]
    zms = _shuffle_list(splitdiffs, 3)[:: (num_pts - 1)]
    xs = []
    # Constructing (+; 0, 0, 0 ...)
    #              (-; 0, 0, 0 ...)
    #              (0; +, +, + ...)
    xs.append(sum(pzs))
    xs.append(sum(mzs))
    xs.append(sum(zps))
    # Constructing the p vectors with one minus
    # (0; -, +, + ...) + cyclic
    for procloc, zm in enumerate(zms):
        newvec = zps[0].copy()
        newvec.loc[:] = 0
        subzps = zps.copy()
        del subzps[procloc]
        newvec = newvec + sum(subzps) + zm
        xs.append(newvec)
    return xs


def vectors_5barpt(splitdiffs):
    """Returns the linearly independent vectors for 5barpt prescription"""
    num_pts = 5
    pps = splitdiffs[:: (num_pts - 1)]
    mms = _shuffle_list(splitdiffs, 1)[:: (num_pts - 1)]
    pms = _shuffle_list(splitdiffs, 2)[:: (num_pts - 1)]
    mps = _shuffle_list(splitdiffs, 3)[:: (num_pts - 1)]
    xs = []
    # Constructing (+/-; +, + ...)
    xs.append(sum(pps))
    xs.append(sum(mps))
    # Constructing the 2p vectors with one minus
    # (+; -, +, + ...) + cyclic
    for procloc, pm in enumerate(pms):
        newvec = pms[0].copy()
        newvec.loc[:] = 0
        subpps = pps.copy()
        del subpps[procloc]
        newvec = newvec + sum(subpps) + pm
        xs.append(newvec)
    # (-; -, +, + ...) + cyclic
    for procloc, mm in enumerate(mms):
        newvec = mms[0].copy()
        newvec.loc[:] = 0
        submps = mps.copy()
        del submps[procloc]
        newvec = newvec + sum(submps) + mm
        xs.append(newvec)
    return xs


def vectors_7pt(splitdiffs):
    """Returns the linearly independent vectors for 7pt prescription"""
    num_pts = 7
    pzs = splitdiffs[:: (num_pts - 1)]
    mzs = _shuffle_list(splitdiffs, 1)[:: (num_pts - 1)]
    zps = _shuffle_list(splitdiffs, 2)[:: (num_pts - 1)]
    zms = _shuffle_list(splitdiffs, 3)[:: (num_pts - 1)]
    pps = _shuffle_list(splitdiffs, 4)[:: (num_pts - 1)]
    mms = _shuffle_list(splitdiffs, 5)[:: (num_pts - 1)]
    xs = []
    # 7pt is the sum of 3pts and 5pts
    # 3pt-like part:
    xs.append(sum(pps))
    for procloc, mm in enumerate(mms):
        newvec = pps[0].copy()
        newvec.loc[:] = 0
        subpps = pps.copy()
        del subpps[procloc]
        newvec = newvec + sum(subpps) + mm
        xs.append(newvec)
    # 5pt-like part:
    xs.append(sum(pzs))
    xs.append(sum(mzs))
    xs.append(sum(zps))
    for procloc, zm in enumerate(zms):
        newvec = zps[0].copy()
        newvec.loc[:] = 0
        subzps = zps.copy()
        del subzps[procloc]
        newvec = newvec + sum(subzps) + zm
        xs.append(newvec)
    return xs


def vectors_9pt(splitdiffs):
    """Returns the linearly independent vectors for 9pt prescription"""
    num_pts = 9
    pzs = splitdiffs[:: (num_pts - 1)]
    mzs = _shuffle_list(splitdiffs, 1)[:: (num_pts - 1)]
    zps = _shuffle_list(splitdiffs, 2)[:: (num_pts - 1)]
    zms = _shuffle_list(splitdiffs, 3)[:: (num_pts - 1)]
    pps = _shuffle_list(splitdiffs, 4)[:: (num_pts - 1)]
    mms = _shuffle_list(splitdiffs, 5)[:: (num_pts - 1)]
    pms = _shuffle_list(splitdiffs, 6)[:: (num_pts - 1)]
    mps = _shuffle_list(splitdiffs, 7)[:: (num_pts - 1)]
    xs = []
    # Constructing (+/-/0; +, +, ...)
    xs.append(sum(pps))
    xs.append(sum(mps))
    xs.append(sum(zps))
    # Constructing (+/-/0; -, +, + ...) + cyclic
    # -- Constructing (0; -, +, + ...) + cyclic
    for procloc, zm in enumerate(zms):
        newvec = zps[0].copy()
        newvec.loc[:] = 0
        subzps = zps.copy()
        del subzps[procloc]
        newvec = newvec + sum(subzps) + zm
        xs.append(newvec)
    # -- Constructing (+; -, +, + ...) + cyclic
    for procloc, pm in enumerate(pms):
        newvec = pps[0].copy()
        newvec.loc[:] = 0
        subpps = pps.copy()
        del subpps[procloc]
        newvec = newvec + sum(subpps) + pm
        xs.append(newvec)
    # -- Constructing (-; -, +, + ...) + cyclic
    for procloc, mm in enumerate(mms):
        newvec = mps[0].copy()
        newvec.loc[:] = 0
        submps = mps.copy()
        del submps[procloc]
        newvec = newvec + sum(submps) + mm
        xs.append(newvec)
    # Constructing (+/-; 0, +, +, ...) + cyclic
    # -- Constructing (+; 0, +, + ...) + cyclic
    for procloc, pz in enumerate(pzs):
        newvec = pps[0].copy()
        newvec.loc[:] = 0
        subpps = pps.copy()
        del subpps[procloc]
        newvec = newvec + sum(subpps) + pz
        xs.append(newvec)
    # -- Constructing (-; 0, +, + ...) + cyclic
    for procloc, mz in enumerate(mzs):
        newvec = mps[0].copy()
        newvec.loc[:] = 0
        submps = mps.copy()
        del submps[procloc]
        newvec = newvec + sum(submps) + mz
        xs.append(newvec)
    return xs


@check_correct_theory_combination_theoryconfig
def evals_nonzero_basis(
    allthx_vector,
    thx_covmat,
    thx_vector,
    collected_theoryids,
    fivetheories,
    seventheories: (str, type(None)) = None,
    orthonormalisation: (str, type(None)) = None,
):
    """Projects the theory covariance matrix from the data space into
    the basis of non-zero eigenvalues, dependent on point-prescription.
    Then returns the eigenvalues (w) and eigenvectors (v)
    in the data space. There are three methods to linearly
    orthonormalise the basis vectors for the covariance matrix,
    and the choice must be specified using the "orthonormalisation"
    flag in the runcard. The choices are: gs, the Gram-Schmidt
    method; qr, QR decomposition; svd, singular value decomposition.
    QR is the method which should be used as standard; the others
    exist for testing purposes."""

    covmat = thx_covmat[0] / (np.outer(thx_vector[0], thx_vector[0]))
    # constructing vectors of shifts due to scale variation
    diffs = [
        ((thx_vector[0] - scalevarvector) / thx_vector[0])
        for scalevarvector in allthx_vector[0]
    ]
    # number of points in point prescription
    num_pts = len(diffs) + 1
    # constructing dictionary of datasets in each process type
    indexlist = list(diffs[0].index.values)
    procdict = {}
    for index in indexlist:
        name = index[0]
        proc = process_lookup(name)
        if proc not in list(procdict.keys()):
            procdict[proc] = [name]
        elif name not in procdict[proc]:
            procdict[proc].append(name)
    # splitting up the scale-varied shift vectors into different spaces per process
    splitdiffs = []
    for process in procdict.keys():
        alldatasets = [y for x in list(procdict.values()) for y in x]
        otherdatasets = [x for x in alldatasets if x not in procdict[process]]
        for diff in diffs:
            splitdiff = diff.copy()
            for ds in otherdatasets:
                splitdiff.loc[ds] = 0
            splitdiffs.append(splitdiff)
    # --------------------------------------------------
    # CONSTRUCTING THE LINEARLY INDEPENDENT VECTORS
    # treating each prescription on a case-by-case basis
    # Notation:
    # e.g. pp => (mu_0; mu_i) = (+;+)
    #      mz => (mu_0; mu_i) = (-;0)
    #      zp => (mu_0; mu_i) = (0;+) ...
    # for a process i,
    # and total vectors are notated like
    # (mu_0; mu_1, mu_2, ..., mu_p)
    if num_pts == 3:
        xs = vectors_3pt(splitdiffs)
    elif (num_pts == 5) and (fivetheories == "nobar"):
        xs = vectors_5pt(splitdiffs)
    elif (num_pts == 5) and (fivetheories == "bar"):
        xs = vectors_5barpt(splitdiffs)
    elif (num_pts == 7) and (seventheories != "original"):
        xs = vectors_7pt(splitdiffs)
    elif num_pts == 9:
        xs = vectors_9pt(splitdiffs)
    # ------------------------------------------------
    # Orthonormalising vectors according to Gram-Schmidt
    if orthonormalisation == "gs":
        ys = [x / np.linalg.norm(x) for x in xs]
        for i in range(1, len(xs)):
            for j in range(0, i):
                ys[i] = ys[i] - (ys[i].T.dot(ys[j]))[0][0] * ys[j] / np.linalg.norm(
                    ys[j]
                )
                ys[i] = ys[i] / np.linalg.norm(ys[i])
        p = pd.concat(ys, axis=1)
    # Orthonormalising vectors according to singular value decomposition
    elif orthonormalisation == "svd":
        xsmatrix = pd.concat(xs, axis=1)
        p = la.orth(xsmatrix)
    # Orthonormalising vectors according to QR decomposition
    elif orthonormalisation == "qr":
        xsmatrix = pd.concat(xs, axis=1)
        p = np.linalg.qr(xsmatrix)[0]
    # Projecting covariance matrix onto subspace of non-zero eigenvalues
    projected_matrix = (p.T).dot(covmat.dot(p))
    cond_num = np.linalg.cond(projected_matrix)
    w, v_projected = la.eigh(projected_matrix)
    # Finding eigenvectors in data space
    v = p.dot(v_projected)
    return w, v, cond_num


def projected_condition_num(evals_nonzero_basis):
    cond_num = evals_nonzero_basis[2]
    return cond_num


def theory_shift_test(shx_vector, evals_nonzero_basis):
    """Compares the NNLO-NLO shift, f, with the eigenvectors and eigenvalues of the
    theory covariance matrix, and returns the component of the NNLO-NLO shift
    space which is missed by the covariance matrix space: fmiss, as well as the
    projections of the shift vector onto each of the eigenvectors: projectors."""
    w, v = evals_nonzero_basis[:2]
    v = np.real(v)
    # NNLO-NLO shift vector
    f = -shx_vector[0].values.T[0]
    # Projecting the shift vector onto each of the eigenvectors
    projectors = np.sum(f * v.T, axis=1)
    # Initialise array of zeros and set precision to same as FK tables
    projected_evectors = np.zeros((len(projectors), (len(f))), dtype=np.float32)
    for i, projector in enumerate(projectors):
        projected_evectors[i] = projector * v[:, i]
    fmiss = f - np.sum(projected_evectors, axis=0)
    return w, v, projectors, f, fmiss


@table
def theory_covmat_eigenvalues(theory_shift_test):
    """Returns a table of s = sqrt(eigenvalue), the projector and
    the ratio of the two, ordered by largest eigenvalue."""
    w = theory_shift_test[0]
    projectors = theory_shift_test[2]
    s = np.sqrt(np.abs(w))
    projectors = np.ndarray.tolist(projectors)
    ratio = projectors / s
    table = pd.DataFrame(
        [s[::-1], projectors[::-1], ratio[::-1]],
        index=[r"$s_a$", r"$\delta_a$", r"$\delta_a/s_a$"],
        columns=np.arange(1, len(s) + 1, 1),
    )
    return table


def efficiency(theory_shift_test):
    """Returns (efficiency = 1 - fmiss/f) with which the theory
    covariance matrix encapsulates the NNLO-NLO shift."""
    f = theory_shift_test[3]
    fmiss = theory_shift_test[4]
    fs = f - fmiss
    fmod = np.sqrt(np.sum(f ** 2))
    fs_mod = np.sqrt(np.sum(fs ** 2))
    efficiency = fs_mod / fmod
    print(f"efficiency = {efficiency}")
    return efficiency


def validation_theory_chi2(theory_shift_test):
    """Returns the theory chi2 for comparing NNLO-NLO shift
    with theory covariance matrix."""
    projectors = theory_shift_test[2]
    evals = theory_shift_test[0]
    ratio = projectors / np.sqrt(np.abs(evals))
    th_chi2 = 1 / len(evals) * np.sum(ratio ** 2)
    print(f"Theory chi2 = {th_chi2}")
    return th_chi2


def theta(theory_shift_test):
    """Returns the angle between the NNLO-NLO
    shift vector and the component of this which is captured
    by the theory covariance matrix"""
    f = theory_shift_test[3]
    fmiss = theory_shift_test[4]
    fs = f - fmiss
    fmod = np.sqrt(np.sum(f ** 2))
    fs_mod = np.sqrt(np.sum(fs ** 2))
    costheta = f @ fs / (fmod * fs_mod)
    th = np.arccos(costheta)
    return th


@figure
def projector_eigenvalue_ratio(theory_shift_test):
    """Produces a plot of the ratio between the projectors and the square roots
    of the corresponding eigenvalues."""
    evals = theory_shift_test[0][::-1]
    projectors = theory_shift_test[2][::-1]
    fmiss = theory_shift_test[4]
    fmiss_mod = np.sqrt(np.sum(fmiss ** 2))
    ratio = np.abs(projectors) / np.sqrt(np.abs(evals))
    # Initialise array of zeros and set precision to same as FK tables
    # Ordering according to shift vector
    mask = np.argsort(np.abs(projectors))[::-1]
    evals = np.asarray(evals)[mask]
    projectors = projectors[mask]
    ratio = ratio[mask]
    xvals = np.arange(1, len(evals) + 1, 1)
    # Plotting
    fig, (ax1, ax2) = plt.subplots(2, figsize=(5, 5))
    ax1.plot(xvals, np.abs(projectors), "s", label=r"|$\delta_a$|")
    ax1.plot(xvals, np.sqrt(np.abs(evals)), "o", label=r"$|s_a|$")
    ax1.plot(0, fmiss_mod, "*", label=r"$|\delta_{miss}|$", color="b")
    ax2.plot(xvals, ratio, "D", color="red")
    ax2.plot(0, 0, ".", color="w")
    ax1.set_title(f"Number of eigenvalues = {len(evals)}", fontsize=10)
    ax1.set_yscale("log")
    ax2.set_yscale("log")
    ax1.legend()
    labels = [item.get_text() for item in ax1.get_xticklabels()]
    ax1.set_xticklabels(labels)
    ax2.set_xticklabels(labels)
    ax2.axhline(y=1, color="k", label=r"|$\delta_a$/$s_a$| = 1")
    ax2.legend()
    ax2.set_ylabel(r"|$\delta_a$/$s_a$|")
    print(f"Subspace dimension = {len(evals)}")
    return fig


@figure
def eigenvector_plot(evals_nonzero_basis, shx_vector):
    """Produces a plot of the eigenvectors for the
    projected matrix, transformed back to the data space."""
    evals = evals_nonzero_basis[0][::-1]
    evecs = evals_nonzero_basis[1].T[::-1]
    f = shx_vector[0]
    indexlist = list(f.index.values)
    # adding process index for plotting, and reindexing matrices and vectors
    dsnames = []
    processnames = []
    ids = []
    for index in indexlist:
        name = index[0]
        i = index[1]
        dsnames.append(name)
        ids.append(i)
        proc = process_lookup(name)
        processnames.append(proc)
    tripleindex = pd.MultiIndex.from_arrays(
        [processnames, dsnames, ids], names=("process", "dataset", "id")
    )
    f = pd.DataFrame(f.values, index=tripleindex)
    f.sort_index(0, inplace=True)
    oldindex = f.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    f = f.reindex(newindex)
    fig, axes = plt.subplots(nrows=len(evecs), figsize=(10, 2 * len(evecs)))
    fig.subplots_adjust(hspace=0.8)
    for ax, evec, eval in zip(axes.flatten(), evecs, evals):
        eval_3sf = floatformatting.significant_digits(eval.item(), 3)
        evec = pd.DataFrame(evec, index=tripleindex)
        evec = evec.reindex(newindex)
        ax.plot(-f.values, color="k", label="NNLO-NLO shift")
        ax.plot(evec.values, label="Eigenvector")
        ticklocs, ticklabels, startlocs = matrix_plot_labels(evec)
        # Shift startlocs elements 0.5 to left so lines are between indexes
        startlocs_lines = [x - 0.5 for x in startlocs]
        ax.vlines(
            startlocs_lines, ax.get_ylim()[0], ax.get_ylim()[1], linestyles="dashed"
        )
        ax.margins(x=0, y=0)
        # Adding eigenvalue to legend
        extraString = f"Eigenvalue = {eval_3sf}"
        handles, labels = ax.get_legend_handles_labels()
        handles.append(mpatches.Patch(color="none", label=extraString))
        ax.legend(handles=handles)
        ax.set_xticks(ticklocs)
        ax.set_xticklabels(ticklabels, rotation=45, fontsize=10)
    return fig


@figure
def deltamiss_plot(theory_shift_test, allthx_vector, evals_nonzero_basis, shx_vector):
    """Produces a plot of the missing component of the
    shift vector, transformed back to the data space."""
    # Define l, which is the number of points in the point prescription being used
    l = len(allthx_vector[0]) + 1
    # Minus sign changes it from NLO-NNLO shift to NNLO-NLO shift (convention)
    f = -shx_vector[0]
    fmiss = theory_shift_test[4]
    indexlist = list(f.index.values)
    # adding process index for plotting, and reindexing matrices and vectors
    dsnames = []
    processnames = []
    ids = []
    for index in indexlist:
        name = index[0]
        i = index[1]
        dsnames.append(name)
        ids.append(i)
        proc = process_lookup(name)
        processnames.append(proc)
    # Index and reindex f and fmiss
    tripleindex = pd.MultiIndex.from_arrays(
        [processnames, dsnames, ids], names=("process", "dataset", "id")
    )
    f = pd.DataFrame(f.values, index=tripleindex)
    f.sort_index(0, inplace=True)
    oldindex = f.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    f = f.reindex(newindex)
    fmiss = pd.DataFrame(fmiss, index=tripleindex)
    fmiss.sort_index(0, inplace=True)
    fmiss = fmiss.reindex(newindex)
    # Plotting
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.plot(f.values * 100, ".-", label="NNLO-NLO Shift", color="black")
    ax.plot(
        fmiss.values * 100, ".-", label=r"$\delta_{miss}$" + f" ({l} pt)", color="blue"
    )
    ticklocs, ticklabels, startlocs = matrix_plot_labels(f)
    plt.xticks(ticklocs, ticklabels, rotation=45, fontsize=20)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ax.vlines(startlocs_lines, -70, 70, linestyles="dashed")
    ax.margins(x=0, y=0)
    ax.set_ylabel(r"% wrt central theory $T_i^{(0)}$", fontsize=20)
    ax.set_ylim(-20, 20)
    ax.legend(fontsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    return fig


@figure
def shift_diag_cov_comparison(allthx_vector, shx_vector, thx_covmat, thx_vector):
    """Produces a plot of a comparison between the NNLO-NLO shift and the
    envelope given by the diagonal elements of the theory covariance matrix."""
    l = len(allthx_vector[0]) + 1
    matrix = thx_covmat[0] / (np.outer(thx_vector[0], thx_vector[0]))
    fnorm = -shx_vector[0]
    indexlist = list(matrix.index.values)
    # adding process index for plotting, and reindexing matrices and vectors
    dsnames = []
    processnames = []
    ids = []
    for index in indexlist:
        name = index[0]
        i = index[1]
        dsnames.append(name)
        ids.append(i)
        proc = process_lookup(name)
        processnames.append(proc)
    tripleindex = pd.MultiIndex.from_arrays(
        [processnames, dsnames, ids], names=("process", "dataset", "id")
    )
    matrix = pd.DataFrame(matrix.values, index=tripleindex, columns=tripleindex)
    matrix.sort_index(0, inplace=True)
    matrix.sort_index(1, inplace=True)
    oldindex = matrix.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    matrix = matrix.reindex(newindex)
    matrix = (matrix.T.reindex(newindex)).T
    sqrtdiags = np.sqrt(np.diag(matrix))
    fnorm = pd.DataFrame(fnorm.values, index=tripleindex)
    fnorm.sort_index(0, inplace=True)
    fnorm = fnorm.reindex(newindex)
    # Plotting
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.plot(sqrtdiags * 100, ".-", label=f"MHOU ({l} pt)", color="red")
    ax.plot(-sqrtdiags * 100, ".-", color="red")
    ax.plot(fnorm.values * 100, ".-", label="NNLO-NLO Shift", color="black")
    ticklocs, ticklabels, startlocs = matrix_plot_labels(matrix)
    plt.xticks(ticklocs, ticklabels, rotation=45, fontsize=20)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ax.vlines(startlocs_lines, -70, 70, linestyles="dashed")
    ax.margins(x=0, y=0)
    ax.set_ylabel(r"% wrt central theory $T_i^{(0)}$", fontsize=20)
    ax.set_ylim(-35, 35)
    ax.legend(fontsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    return fig
