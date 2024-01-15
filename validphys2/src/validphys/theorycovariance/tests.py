# -*- coding: utf-8 -*-
"""
tests.py
Tools for testing theory covariance matrices and their properties.
"""
from __future__ import generator_stop

from collections import namedtuple
import logging

from matplotlib.figure import Figure
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scipy.linalg as la

from reportengine import collect, floatformatting
from reportengine.figure import figure
from reportengine.table import table
from validphys import plotutils
from validphys.theorycovariance.output import _get_key, matrix_plot_labels
from validphys.theorycovariance.theorycovarianceutils import process_lookup

log = logging.getLogger(__name__)

matched_dataspecs_results = collect("results", ["dataspecs"])

LabeledShifts = namedtuple("LabeledShifts", ("process", "dataset_name", "shifts"))


def shift_vector(matched_dataspecs_results, process, dataset_name):
    """Returns a DataFrame of normalised shift vectors and normalization vector for matched dataspecs."""
    r1, r2 = matched_dataspecs_results
    shifts = r1[1].central_value - r2[1].central_value
    norm = r1[1].central_value
    norm_shifts = shifts / norm
    dsnames = np.full(len(shifts), dataset_name, dtype=object)
    processnames = np.full(len(shifts), process, dtype=object)
    index = pd.MultiIndex.from_arrays(
        [processnames, dsnames], names=["group", "dataset"]
    )
    return pd.DataFrame({'shifts': norm_shifts, 'norm': norm}, index=index)


def dataset_alltheory(each_dataset_results_central_bytheory):
    """Returns a LabeledShifts tuple corresponding to the theory
    vectors for all the scale varied theories, including the central."""
    labeled_shifts_list = []
    for dataset in each_dataset_results_central_bytheory:
        res = [theoryid[1].central_value for theoryid in dataset]
        labeled_shifts_list.append(
            LabeledShifts(
                dataset_name=dataset[0][0].name,
                process=process_lookup(dataset[0][0].name),
                shifts=res,
            )
        )
    return labeled_shifts_list


def alltheory_vector(dataset_alltheory):
    """Returns a DataFrame with the theory vectors for matched
    dataspecs for the scale-varied theories, including the central"""
    all_theory = np.concatenate([val.shifts for val in dataset_alltheory], axis=1)
    dsnames = np.concatenate(
        [np.full(len(val.shifts[0]), val.dataset_name, dtype=object) for val in dataset_alltheory]
    )
    groupnames = np.concatenate(
        [np.full(len(val.shifts[0]), val.process, dtype=object) for val in dataset_alltheory]
    )
    index = pd.MultiIndex.from_arrays([groupnames, dsnames], names=["group", "dataset"])
    theory_vectors = []
    for theoryvector in all_theory:
        theory_vectors.append(pd.DataFrame(theoryvector, index=index))
    return theory_vectors


shx_vector = collect("shift_vector", ["matched_datasets_from_dataspecs"])


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


def tripleindex_thcovmat_complete(theory_covmat_custom):
    """The complete tripleindex of the theory covmat."""
    return theory_covmat_custom.index


def doubleindex_thcovmat(theory_covmat_custom):
    """The unique set of (group, dataset) index of the theory covmat."""
    tripleindex = theory_covmat_custom.index
    return list(dict.fromkeys([(ind[0], ind[1]) for ind in tripleindex]))


def ordered_alltheory_vector(alltheory_vector, doubleindex_thcovmat):
    """Order the group and the dataset of all the theory vectors of the scale-varied predictions in the same way the theory covmat is ordered."""
    ord_alltheory_vector = []
    for theoryid_vector in alltheory_vector:
        list_df = []
        for index in doubleindex_thcovmat:
            list_df.append(theoryid_vector.loc[index[0]].loc[index[1]])
        ord_alltheory_vector.append(pd.concat(list_df).values.T[0])
    return ord_alltheory_vector


def evals_nonzero_basis(
    ordered_alltheory_vector,
    theory_covmat_custom,
    fivetheories,
    orthonormalisation: (str, type(None)) = None,
    seventheories: (str, type(None)) = None,
):
    covmat = theory_covmat_custom / (
        np.outer(ordered_alltheory_vector[0], ordered_alltheory_vector[0])
    )
    # constructing vectors of shifts due to scale variation
    diffs = [
        pd.DataFrame(
            ((ordered_alltheory_vector[0] - scalevarvector) / ordered_alltheory_vector[0]),
            index=theory_covmat_custom.index.droplevel(0).droplevel(1),
        )
        for scalevarvector in ordered_alltheory_vector[1:]
    ]
    # number of points in point prescription
    num_pts = len(ordered_alltheory_vector)
    # constructing dictionary of datasets in each process type
    doubleindex = doubleindex_thcovmat(theory_covmat_custom)
    procdict = {}
    for proc, name in doubleindex:
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
                ys[i] = ys[i] - (ys[i].T.dot(ys[j]))[0][0] * ys[j] / np.linalg.norm(ys[j])
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

def fnorm_shifts_ordered(concatenated_shx_vector, doubleindex_thcovmat):
    """Shifts vectors ordered as the theory covmat."""
    fnorm_vector = fnorm_shifts_byprocess(concatenated_shx_vector, doubleindex_thcovmat)
    return np.array(sum(fnorm_vector, []))

def theory_shift_test(fnorm_shifts_ordered, evals_nonzero_basis):
    """Compares the NNLO-NLO shift, f, with the eigenvectors and eigenvalues of the
    theory covariance matrix, and returns the component of the NNLO-NLO shift
    space which is missed by the covariance matrix space: fmiss, as well as the
    projections of the shift vector onto each of the eigenvectors: projectors."""
    w, v = evals_nonzero_basis[:2]
    v = np.real(v)
    # NNLO-NLO shift vector
    f = -fnorm_shifts_ordered.T
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
    fmod = np.sqrt(np.sum(f**2))
    fs_mod = np.sqrt(np.sum(fs**2))
    efficiency = fs_mod / fmod
    print(f"efficiency = {efficiency}")
    return efficiency


def validation_theory_chi2(theory_shift_test):
    """Returns the theory chi2 for comparing NNLO-NLO shift
    with theory covariance matrix."""
    projectors = theory_shift_test[2]
    evals = theory_shift_test[0]
    ratio = projectors / np.sqrt(np.abs(evals))
    th_chi2 = 1 / len(evals) * np.sum(ratio**2)
    print(f"Theory chi2 = {th_chi2}")
    return th_chi2


def theta(theory_shift_test):
    """Returns the angle between the NNLO-NLO
    shift vector and the component of this which is captured
    by the theory covariance matrix"""
    f = theory_shift_test[3]
    fmiss = theory_shift_test[4]
    fs = f - fmiss
    fmod = np.sqrt(np.sum(f**2))
    fs_mod = np.sqrt(np.sum(fs**2))
    costheta = f @ fs / (fmod * fs_mod)
    th = np.arccos(costheta)
    return np.degrees(th)


@figure
def projector_eigenvalue_ratio(theory_shift_test):
    """Produces a plot of the ratio between the projectors and the square roots
    of the corresponding eigenvalues."""
    evals = theory_shift_test[0][::-1]
    projectors = theory_shift_test[2][::-1]
    fmiss = theory_shift_test[4]
    fmiss_mod = np.sqrt(np.sum(fmiss**2))
    ratio = np.abs(projectors) / np.sqrt(np.abs(evals))
    # Initialise array of zeros and set precision to same as FK tables
    # Ordering according to shift vector
    mask = np.argsort(np.abs(projectors))[::-1]
    evals = np.asarray(evals)[mask]
    projectors = projectors[mask]
    ratio = ratio[mask]
    xvals = np.arange(1, len(evals) + 1, 1)
    # Plotting
    fig = Figure(figsize=(5, 5))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)

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
def eigenvector_plot(evals_nonzero_basis, fnorm_shifts_ordered, tripleindex_thcovmat_complete, dataspecs):
    """Produces a plot of the eigenvectors for the
    projected matrix, transformed back to the data space."""
    evals = evals_nonzero_basis[0][::-1]
    evecs = evals_nonzero_basis[1].T[::-1]
    f = pd.DataFrame(fnorm_shifts_ordered, index=tripleindex_thcovmat_complete)
    # Sort index
    f.sort_index(axis=0, inplace=True)
    oldindex = f.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    f = f.reindex(newindex)
    fig, axes = plotutils.subplots(figsize=(10, 2 * len(evecs)), nrows=len(evecs))
    fig.subplots_adjust(hspace=0.8)
    for ax, evec, eval in zip(axes.flatten(), evecs, evals):
        eval_3sf = floatformatting.significant_digits(eval.item(), 3)
        evec = pd.DataFrame(evec, index=tripleindex_thcovmat_complete)
        evec = evec.reindex(newindex)
        ax.plot(-1.0 * (f.values), color="k", label=f"{dataspecs[1]['speclabel']}-{dataspecs[0]['speclabel']} shift")
        ax.plot(evec.values, label="Eigenvector")
        ticklocs, ticklabels, startlocs = matrix_plot_labels(evec)
        # Shift startlocs elements 0.5 to left so lines are between indexes
        startlocs_lines = [x - 0.5 for x in startlocs]
        ax.vlines(startlocs_lines, ax.get_ylim()[0], ax.get_ylim()[1], linestyles="dashed")
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
def deltamiss_plot(
    theory_shift_test,
    ordered_alltheory_vector,
    fnorm_shifts_ordered,
    tripleindex_thcovmat_complete,
    dataspecs,
):
    """Produces a plot of the missing component of the
    shift vector, transformed back to the data space."""
    # Define l, which is the number of points in the point prescription being used
    l = len(ordered_alltheory_vector) + 1
    # Minus sign changes it from NLO-NNLO shift to NNLO-NLO shift (convention)
    f = -1.0 * pd.DataFrame(fnorm_shifts_ordered, index=tripleindex_thcovmat_complete)
    # Sort index
    f.sort_index(axis=0, inplace=True)
    oldindex = f.index.tolist()
    newindex = sorted(oldindex, key=_get_key)
    f = f.reindex(newindex)
    fmiss = pd.DataFrame(theory_shift_test[4], index=tripleindex_thcovmat_complete)
    fmiss_reordered = fmiss.reindex(f.index)
    # Plotting
    fig, ax = plotutils.subplots(figsize=(20, 10))
    ax.plot(f.values * 100, ".-", label=f"{dataspecs[1]['speclabel']}-{dataspecs[0]['speclabel']} Shift", color="black")
    ax.plot(
        fmiss_reordered.values * 100, ".-", label=r"$\delta_{miss}$" + f" ({l} pt)", color="blue"
    )
    ticklocs, ticklabels, startlocs = matrix_plot_labels(f)
    ax.set_xticks(ticklocs)
    ax.set_xticklabels(ticklabels, rotation=45, fontsize=20)

    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ax.vlines(startlocs_lines, -70, 70, linestyles="dashed")
    ax.margins(x=0, y=0)
    ax.set_ylabel(r"% wrt central theory $T_i^{(0)}$", fontsize=20)
    ax.set_ylim(-20, 20)
    ax.legend(fontsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    return fig


def diagdf_theory_covmat(theory_covmat_custom):
    """Return a Dataframe indexed with groups and dataset of the diagonal entry of the theory covmat."""
    return pd.DataFrame(data=np.diag(theory_covmat_custom.values), index=theory_covmat_custom.index.droplevel(2))


def doubleindex_set_byprocess(group_dataset_inputs_by_process):
    """The (group, dataset) index ordered by process."""
    doubleindex = []
    for process in group_dataset_inputs_by_process:
        for dataset in process['data_input']:
            doubleindex.append((process['group_name'], dataset.name))
    return doubleindex


def concatenated_shx_vector(shx_vector):
    """Single DataFrame for all the datasets of the shift vectors."""
    return pd.concat(shx_vector)


def sqrtdiags_thcovmat_byprocess(doubleindex_set_byprocess, diagdf_theory_covmat, concatenated_shx_vector):
    """Ratio of the sqrts of the diagonal entries of the theory covmat to the normalization entries of the shift vectors, ordered by process."""
    sqrtdiags = []
    for index in doubleindex_set_byprocess:
        if isinstance(concatenated_shx_vector.loc[index[0]].loc[index[1]].norm, np.float64):
            sqrtdiags.append(
            list(
                np.sqrt(diagdf_theory_covmat.loc[index[0]].loc[index[1]].values)
                / concatenated_shx_vector.loc[index[0]].loc[index[1]].norm
            )
        )
        else:
            sqrtdiags.append(
                list(
                    np.sqrt(diagdf_theory_covmat.loc[index[0]].loc[index[1]].values.transpose()[0])
                    / concatenated_shx_vector.loc[index[0]].loc[index[1]].norm.values
                )
            )
    return sqrtdiags


def fnorm_shifts_byprocess(concatenated_shx_vector, doubleindex_set_byprocess):
    """Shift vector ordered by process."""
    fnorm = []
    for index in doubleindex_set_byprocess:
        if isinstance(concatenated_shx_vector.loc[index[0]].loc[index[1]].shifts, np.float64):
            fnorm.append([concatenated_shx_vector.loc[index[0]].loc[index[1]].shifts])
        else:
            fnorm.append(list(concatenated_shx_vector.loc[index[0]].loc[index[1]].shifts.values))
    return fnorm


def ticklocs_thcovmat(theory_covmat_custom):
    return matrix_plot_labels(theory_covmat_custom)


@figure
def shift_diag_cov_comparison(
    sqrtdiags_thcovmat_byprocess, fnorm_shifts_byprocess, point_prescription, ticklocs_thcovmat, dataspecs
):
    """Plot of the comparison of a shift between two pertubative order and the diagonal entries of the theory covmat, both normalized to the first of the two perturbative orders."""
    # Concatenation of the arrays
    fnorm_concat = sum(fnorm_shifts_byprocess, [])
    sqrtdiags_concat = sum(sqrtdiags_thcovmat_byprocess, [])
    fig, ax = plotutils.subplots(figsize=(20, 10))
    ax.plot(
        np.array(sqrtdiags_concat), ".-", label=f"MHOU ({point_prescription})", color="red"
    )
    ax.plot(-np.array(sqrtdiags_concat), ".-", color="red")
    ax.plot(
        -np.array(fnorm_concat),
        ".-",
        label=f"{dataspecs[1]['speclabel']}-{dataspecs[0]['speclabel']} Shift",
        color="black",
    )
    ticklocs, ticklabels, startlocs = ticklocs_thcovmat
    ax.set_xticks(ticklocs)
    ax.set_xticklabels(ticklabels, rotation=45, fontsize=20)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ax.vlines(startlocs_lines, -70, 70, linestyles="dashed")
    ax.margins(x=0, y=0)
    ax.set_ylabel(r"$\pm \sqrt{\hat{S_{ii}}}$, $\delta_{i}$", fontsize=20)
    ax.set_ylim(-0.35, 0.35)
    ax.legend(fontsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    return fig
