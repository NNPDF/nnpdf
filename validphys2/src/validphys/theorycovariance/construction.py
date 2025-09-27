"""
construction.py
Tools for constructing theory covariance matrices and computing their chi2s.
"""

from collections import defaultdict, namedtuple
import logging
import operator

import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.interpolate as scint

from reportengine import collect
from reportengine.table import table

pass
from validphys.results import results, results_central
from validphys.convolution import central_fk_predictions
from validphys.core import PDF
from validphys.theorycovariance.theorycovarianceutils import (
    check_correct_theory_combination,
    check_fit_dataset_order_matches_grouped,
    process_lookup,
)

log = logging.getLogger(__name__)

results_central_bytheoryids = collect(results_central, ("theoryids",))
each_dataset_results_central_bytheory = collect("results_central_bytheoryids", ("data",))


def theory_covmat_dataset(results, results_central_bytheoryids, point_prescription):
    """
    Compute the theory covmat for a collection of theoryids for a single dataset.

    In general this will come from some point prescription and it could be guessed from the input
    however, it takes as input all relevant variables for generality
    """
    _, theory_results = zip(*results_central_bytheoryids)
    _, central_th_result = results

    # Remove the central theory from the list if it was included
    theory_results = [i for i in theory_results if i._theoryid != central_th_result._theoryid]
    cv = central_th_result.central_value

    # Compute the theory contribution to the covmats
    deltas = list(t.central_value - cv for t in theory_results)
    thcovmat = compute_covs_pt_prescrip(point_prescription, "A", deltas)

    return thcovmat


ProcessInfo = namedtuple("ProcessInfo", ("preds", "namelist", "sizes", "data"))


def combine_by_type(each_dataset_results_central_bytheory):
    """Groups the datasets by process and returns an instance of the ProcessInfo class

    Parameters
    ----------
    each_dataset_results_central_bytheory: list[list[(DataResult,ThPredictionsResult)]]
        Tuples of DataResult and ThPredictionsResult (where only the second is used for the
        construction of the theory covariance matrix), wrapped in a list such that there is a tuple
        per theoryid, wrapped in another list per dataset.

    Returns
    -------
    :ProcesInfo :py:class:`validphys.theorycovariance.construction.ProcessInfo`
        Class with info needed to construct the theory covmat.
    """
    dataset_size = defaultdict(list)
    theories_by_process = defaultdict(list)
    ordered_names = defaultdict(list)
    for dataset in each_dataset_results_central_bytheory:
        name = dataset[0][0].name
        theory_centrals = [x[1].central_value for x in dataset]
        dataset_size[name] = len(theory_centrals[0])
        proc_type = process_lookup(name)
        ordered_names[proc_type].append(name)
        theories_by_process[proc_type].append(theory_centrals)
    for key, item in theories_by_process.items():
        theories_by_process[key] = np.concatenate(item, axis=1)
    process_info = ProcessInfo(
        preds=theories_by_process, namelist=ordered_names, sizes=dataset_size, data=None
    )
    return process_info

def combine_by_type_ht(each_dataset_results, groups_dataset_inputs_loaded_cd_with_cuts_byprocess):
    """same as combine_by_type but now for a single theory and including commondata info"""
    dataset_size = defaultdict(list)
    theories_by_process = defaultdict(list)
    cd_by_process = defaultdict(list)
    ordered_names = defaultdict(list)
    for dataset, cd in zip(
        each_dataset_results, groups_dataset_inputs_loaded_cd_with_cuts_byprocess
    ):
        name = cd.setname
        if name != dataset[0].name:
            raise ValueError("The underlying datasets do not match!")
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
        preds=theories_by_process, namelist=ordered_names, sizes=dataset_size, data=cd_by_process
    )
    return process_info


def thcov_ht(H2_list,
             HL_list,
             groups_data_by_process,
             pdf,
             ht_knots = list()):
      """
          Same as `thcov_HT` but implementing theory covariance method for each node of the spline.
          Note that 'groups_data_by_process' contains the same info as 'combine_by_type_ht'. At some
          point we should use only one of them.
      """
      groups_data = groups_data_by_process
      x_knots = list()
      start_proc_by_exp = defaultdict(list)
      ndata_by_exp = defaultdict(list)
      running_index_tot = 0
      included_proc = ["DIS NC"]
      excluded_exp = {"DIS NC" : []}
      deltas = defaultdict(list)

      for i in range(len(x_knots)):
        deltas[f"p({i+1}+,0)"] = np.array([])
        deltas[f"p(0,{i+1}+)"] = np.array([])
        deltas[f"d({i+1}+,0)"] = np.array([])
        deltas[f"d(0,{i+1}+)"] = np.array([])

      if len(ht_knots) == 0:
        # ABMP parametrisation
        x_knots = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
      else:
        x_knots = ht_knots

      # Check that H2_list and HL_list have the same size as x
      if (len(H2_list) != len(x_knots)) or (len(HL_list) != len(x_knots)):
          raise ValueError(f"The size of HT parameters does not match the number of nodes in the spline.")

      def ht_parametrisation(
        index: int,
        nodes: list,
        x: list,
        Q2: list,
        h_prior: list,
        reverse: bool = False
      ):
        if not reverse:
            shifted_H_list = [0 for _ in range(len(nodes))]
            shifted_H_list[index] = h_prior[index]
        else:
            shifted_H_list = h_prior.copy()
            shifted_H_list[index] = 0
            log.warning("The implemenation of the reverse 5-pt prescription must be checked.")

        H = scint.CubicSpline(nodes, shifted_H_list)
        H = np.vectorize(H)

        PC = H(x) / Q2
        return PC

      for idx_proc, group_proc in enumerate(groups_data):
        for idx_exp, exp_set in enumerate(group_proc.datasets):
              # For covmat construction
              exp_ndata = exp_set.load_commondata().ndata
              start_proc_by_exp[exp_set.name] = running_index_tot
              running_index_tot += exp_ndata
              ndata_by_exp[exp_set.name] = exp_ndata

              if group_proc.name in included_proc and exp_set.name not in excluded_exp[group_proc.name]:
                  cd_table = exp_set.load_commondata().commondata_table
                  process_type = cd_table['process'].iloc[0]
                  x = cd_table['kin1'].to_numpy()
                  q2 = cd_table['kin2'].to_numpy()
                  y = cd_table['kin3'].to_numpy()
                  if x.size != exp_ndata:
                      raise ValueError("Problem with the number of data.")

                  for i in range(len(x_knots)):
                    if process_type == "DIS_F2R":
                      cuts = exp_set.cuts
                      fkspec_F2D, fkspec_F2P = exp_set.fkspecs
                      fk_F2D = fkspec_F2D.load_with_cuts(cuts)
                      fk_F2P = fkspec_F2P.load_with_cuts(cuts)
                      F2D = central_fk_predictions(fk_F2D, pdf)
                      F2P = central_fk_predictions(fk_F2P, pdf)

                      F2D = np.concatenate(F2D.values)
                      F2P = np.concatenate(F2P.values)
                      F2_ratio = operator.truediv(F2D, F2P)
                      PC = ht_parametrisation(i, x_knots, x, q2, H2_list)

                      # NOTE
                      # Find a better way to store deltas
                      PC_2_p = np.array(operator.truediv(F2D, np.sum([F2P, PC],axis=0)) - F2_ratio)
                      PC_2_d = np.array(operator.truediv(F2D, np.sum([F2D, PC],axis=0)) - F2_ratio)
                      PC_L_p = np.zeros(F2_ratio.size)
                      PC_L_d = np.zeros(F2_ratio.size)


                    elif process_type == "DIS_F2P":
                      PC_2_p = ht_parametrisation(i, x_knots, x, q2, H2_list)
                      PC_2_d = np.zeros(exp_ndata)
                      PC_L_p = np.zeros(exp_ndata)
                      PC_L_d = np.zeros(exp_ndata)

                    elif process_type == "DIS_F2D":
                      PC_2_p = np.zeros(exp_ndata)
                      PC_2_d = ht_parametrisation(i, x_knots, x, q2, H2_list)
                      PC_L_p = np.zeros(exp_ndata)
                      PC_L_d = np.zeros(exp_ndata)

                    elif process_type == 'DIS_F2C':
                      # Iron target
                      Z = 23.403
                      A = 49.618
                      PC_2_p = ht_parametrisation(i, x_knots, x, q2, H2_list)
                      PC_2_d = 2 * (Z - A) / A * ht_parametrisation(i, x_knots, x, q2, H2_list)
                      PC_L_p = np.zeros(exp_ndata)
                      PC_L_d = np.zeros(exp_ndata)

                    elif process_type == "DIS_NCE" or "DIS_NCP":
                        yp = 1 + np.power(1 - y, 2)
                        yL = np.power(y, 2)
                        N_L = - yL / yp
                        PC_2_p = ht_parametrisation(i, x_knots, x, q2, H2_list)
                        PC_2_d = np.zeros(exp_ndata)
                        PC_L_p = N_L * ht_parametrisation(i, x_knots, x, q2, HL_list)
                        PC_L_d = np.zeros(exp_ndata)

                    else:
                      raise Exception(f"The process type `{process_type}` has not been implemented.")

                    deltas[f"p({i+1}+,0)"] = np.append(deltas[f"p({i+1}+,0)"], PC_2_p)
                    deltas[f"p(0,{i+1}+)"] = np.append(deltas[f"p(0,{i+1}+)"], PC_2_d)
                    deltas[f"d({i+1}+,0)"] = np.append(deltas[f"d({i+1}+,0)"], PC_L_p)
                    deltas[f"d(0,{i+1}+)"] = np.append(deltas[f"d(0,{i+1}+)"], PC_L_d)

              else:
                for i in range(len(x_knots)):
                  deltas[f"p({i+1}+,0)"] = np.append(deltas[f"p({i+1}+,0)"], np.zeros(exp_ndata))
                  deltas[f"p(0,{i+1}+)"] = np.append(deltas[f"p(0,{i+1}+)"], np.zeros(exp_ndata))
                  deltas[f"d({i+1}+,0)"] = np.append(deltas[f"d({i+1}+,0)"], np.zeros(exp_ndata))
                  deltas[f"d(0,{i+1}+)"] = np.append(deltas[f"d(0,{i+1}+)"], np.zeros(exp_ndata))

      # Construct the covariance matrix
      covmats = defaultdict(list)
      for exp_name_1, exp_idx_1 in start_proc_by_exp.items():
          for exp_name_2, exp_idx_2 in start_proc_by_exp.items():
              s = np.zeros(shape=(ndata_by_exp[exp_name_1], ndata_by_exp[exp_name_2]))
              for shifts in deltas.keys():
                s += np.outer(deltas[shifts][exp_idx_1: exp_idx_1 + ndata_by_exp[exp_name_1]],
                              deltas[shifts][exp_idx_2: exp_idx_2 + ndata_by_exp[exp_name_2]])

              start_locs = (exp_idx_1, exp_idx_2)
              covmats[start_locs] = s
      return covmats


def covmat_3fpt(name1, name2, deltas1, deltas2):
    """Returns theory covariance sub-matrix for 3pt factorisation
    scale variation *only*, given two dataset names and collections
    of scale variation shifts"""
    s = 0.5 * (np.outer(deltas1[0], deltas2[0]) + np.outer(deltas1[1], deltas2[1]))
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
    for cnt in range(0, 14, 2):
        s_ad += covmat_n3lo_ad(name1, name2, deltas1[cnt : cnt + n_var], deltas2[cnt : cnt + n_var])
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


def compute_covs_pt_prescrip(point_prescription, name1, deltas1, name2=None, deltas2=None):
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
        name1: str
            Process name of the first set of shifts
        deltas1: list(np.ndarray)
            list of shifts for each of the non-central theories
        name2: str
            Process name of the second set of shifts
        deltas2: list(np.ndarray)
            list of shifts for each of the non-central theories
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

    if point_prescription == "3f point":
        s = covmat_3fpt(deltas1, deltas2)
    elif point_prescription == "3r point":
        s = covmat_3pt(name1, name2, deltas1, deltas2)
    elif point_prescription == "3 point":
        s = covmat_3pt(name1, name2, deltas1, deltas2)
    elif point_prescription == "5 point":
        s = covmat_5pt(name1, name2, deltas1, deltas2)
    elif point_prescription == "5bar point":
        s = covmat_5barpt(name1, name2, deltas1, deltas2)
    elif point_prescription == "7 point":
        # 7pt (Gavin)
        s = covmat_7pt(name1, name2, deltas1, deltas2)
    elif point_prescription == "9 point":
        s = covmat_9pt(name1, name2, deltas1, deltas2)
    elif point_prescription == "ad ihou":
        # n3lo ad variation prescriprion
        s = covmat_n3lo_singlet(name1, name2, deltas1, deltas2)
    elif point_prescription == "dis ihou":
        # n3lo ihou prescriprion
        s = covmat_3pt(name1, name2, deltas1, deltas2)
    elif point_prescription == "3pt missing":
        # 3 point renormalization scale variations for subset of data
        s = covmat_3pt(name1, name2, deltas1, deltas2)
    elif point_prescription == "3pt hadronic":
        # N3LO 3 point scale variations for hadronic datasets
        s = covmat_3pt(name1, name2, deltas1, deltas2)
    elif point_prescription == "fhmv ihou":
        # n3lo full covmat prescriprion
        s = covmat_n3lo_fhmv(name1, name2, deltas1, deltas2)
    elif point_prescription.startswith("alphas"):
        # alphas is correlated for all datapoints and the covmat construction is
        # therefore equivalent to that of the factorization scale variations
        s = covmat_3fpt(deltas1, deltas2)
    return s


@check_correct_theory_combination
def covs_pt_prescrip(combine_by_type, point_prescription):
    """Produces the sub-matrices of the theory covariance matrix according
    to a point prescription which matches the number of input theories.
    If 5 theories are provided, a scheme 'bar' or 'nobar' must be
    chosen in the runcard in order to specify the prescription. Sub-matrices
    correspond to applying the scale variation prescription to each pair of
    processes in turn, using a different procedure for the case where the
    processes are the same relative to when they are different."""

    process_info = combine_by_type
    running_index = 0
    start_proc = defaultdict(list)
    for name in process_info.preds:
        size = len(process_info.preds[name][0])
        start_proc[name] = running_index
        running_index += size

    covmats = defaultdict(list)
    for name1 in process_info.preds:
        for name2 in process_info.preds:
            central1, *others1 = process_info.preds[name1]
            deltas1 = list(other - central1 for other in others1)
            central2, *others2 = process_info.preds[name2]
            deltas2 = list(other - central2 for other in others2)
            s = compute_covs_pt_prescrip(point_prescription, name1, deltas1, name2, deltas2)
            start_locs = (start_proc[name1], start_proc[name2])
            covmats[start_locs] = s
    return covmats


@table
def theory_covmat_custom(covmat_custom, procs_index, combine_by_type_custom):
    """Takes the individual sub-covmats between each two processes and assembles
    them into a full covmat. Then reshuffles the order from ordering by process
    to ordering by experiment as listed in the runcard"""
    process_info = combine_by_type_custom

    # Construct a covmat_index based on the order of experiments as they are in combine_by_type
    # NOTE: maybe the ordering of covmat_index is always the same as that of procs_index?
    # Regardless, we don't want to open ourselves up to the risk of the ordering of procs_index
    # changing and breaking this function
    indexlist = []
    for procname in process_info.preds:
        for datasetname in process_info.namelist[procname]:
            slicer = procs_index.get_locs((procname, datasetname))
            indexlist += procs_index[slicer].to_list()
    covmat_index = pd.MultiIndex.from_tuples(indexlist, names=procs_index.names)

    # Put the covariance matrices between two process into a single covariance matrix
    total_datapoints = sum(process_info.sizes.values())
    mat = np.zeros((total_datapoints, total_datapoints), dtype=np.float32)
    for locs, cov in covmat_custom.items():
        xsize, ysize = cov.shape
        mat[locs[0] : locs[0] + xsize, locs[1] : locs[1] + ysize] = cov
    df = pd.DataFrame(mat, index=covmat_index, columns=covmat_index)
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


def theory_covmat_custom_fitting(theory_covmat_custom_per_prescription, procs_index_matched):
    """theory_covmat_custom_per_prescription but reindexed so the order of the datasets matches
    those in the experiment covmat so they are aligned when fitting."""
    df = theory_covmat_custom_per_prescription.reindex(procs_index_matched).T.reindex(
        procs_index_matched
    )
    return df


theory_covmats_fitting = collect(theory_covmat_custom_per_prescription, ("point_prescriptions",))


@table
def theory_covmat_custom(theory_covmats_fitting):
    """Sum all the theory covmat listed in `point_prescriptions`."""
    return sum(theory_covmats_fitting)


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


@table
def theory_corrmat_custom(theory_covmat_custom):
    """Calculates the theory correlation matrix for scale variations
    with variations by process type"""
    df = theory_covmat_custom
    covmat = df.values
    diag_minus_half = (np.diagonal(covmat)) ** (-0.5)
    mat = diag_minus_half[:, np.newaxis] * df * diag_minus_half
    return mat


@table
def theory_normcovmat_custom(theory_covmat_custom, procs_data_values):
    """Calculates the theory covariance matrix for scale variations normalised
    to data, with variations according to the relevant prescription."""
    df = theory_covmat_custom
    vals = procs_data_values.reindex(df.index)
    mat = df / np.outer(vals, vals)
    return mat


@table
def experimentplustheory_corrmat_custom(procs_covmat, theory_covmat_custom):
    """Calculates the correlation matrix for the experimental
    plus theory covariance matrices, correlations by prescription."""
    total_df = procs_covmat + theory_covmat_custom
    diag_minus_half = (np.diagonal(total_df.values)) ** (-0.5)
    corrmat = diag_minus_half[:, np.newaxis] * total_df * diag_minus_half
    return corrmat


each_dataset_results = collect(results, ("group_dataset_inputs_by_process", "data"))
groups_data_by_process = collect("data", ("group_dataset_inputs_by_process",))