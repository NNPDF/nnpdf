"""
This module contains the utilities for the computation of the shifts in
the theoretical predictions due to the effects of power corrections. Contrary
to scale variations, in the case of power corrections the shifts are not
computed using theories. The shifts are computed at "runtime" during vp-setupfit.

The aim is that, after shifts being computed, the covmat can be constructed using
the same functions implemented for scale variations (e.g. `covs_pt_prescrip`).
The way shifts are computed depends also on the point prescription. In the case of
scale variations, the point prescription specifies the theories whereby shifts are
computed. In the case of power corrections, shifts and covmat constructions are
computed using a 5-point prescription extended to every parameter used to define
the power correction.

This module comprehends a bunch of ``factory`` functions such as `mult_dis_pc`. Each
of these functions returns another function that computes the shifts taking as arguments
the values of the parameters used to parametrise the power corrections. In
other words, these factory functions hard-code the dependence on the kinematic
and leave the dependence on the parameters free. In this way, the shifts
can be computed using different combinations of parameters (i.e. different prescriptions)
if needed.
"""

from collections import defaultdict
import operator

import numpy as np
import numpy.typing as npt

from validphys.convolution import central_fk_predictions
from validphys.core import PDF, DataSetSpec

F2P_exps = ['SLAC_NC_NOTFIXED_P_EM-F2', 'BCDMS_NC_NOTFIXED_P_EM-F2']
F2D_exps = ['SLAC_NC_NOTFIXED_D_EM-F2', 'BCDMS_NC_NOTFIXED_D_EM-F2']
NC_SIGMARED_P_EM = ['NMC_NC_NOTFIXED_P_EM-SIGMARED', 'HERA_NC_318GEV_EM-SIGMARED']
NC_SIGMARED_P_EP = [
    'HERA_NC_225GEV_EP-SIGMARED',
    'HERA_NC_251GEV_EP-SIGMARED',
    'HERA_NC_300GEV_EP-SIGMARED',
    'HERA_NC_318GEV_EP-SIGMARED',
]
NC_SIGMARED_P_EAVG = ['HERA_NC_318GEV_EAVG_CHARM-SIGMARED', 'HERA_NC_318GEV_EAVG_BOTTOM-SIGMARED']


# TODO This function will be deleted in the future
def step_function(a: npt.ArrayLike, y_shift: npt.ArrayLike, bin_edges: npt.ArrayLike) -> np.ndarray:
    """
    This function defines the step function used to construct the prior. The bins of the step
    function are constructed using pairs of consecutive points. For instance, given the set of
    points [0.0, 0.1, 0.3, 0.5], there will be three bins with edges [[0.0, 0.1], [0.1, 0.3],
    0.3, 0.5]]. Each bin is coupled with a shift, which correspond to the y-value of the bin.

    Parameters
    ----------
    a:  ArrayLike of float
      A one-dimensional array of points at which the function is evaluated.
    y_shift:  ArrayLike of float
      A one-dimensional array whose elements represent the y-value of each bin
    bin_edges: ArrayLike of float
      A one-dimensional array containing the edges of the bins. The bins are
      constructed using pairs of consecutive points.

    Return
    ------
    A one-dimensional array containing the function values evaluated at the points
    specified in `a`.
    """
    res = np.zeros_like(a)
    for shift_pos, shift in enumerate(y_shift):
        bin_low = bin_edges[shift_pos]
        bin_high = bin_edges[shift_pos + 1]
        condition = np.multiply(
            a >= bin_low, a < bin_high if shift_pos != len(y_shift) - 1 else a <= bin_high
        )
        res = np.add(res, [shift if cond else 0.0 for cond in condition])
    return res


# TODO This function will be deleted in the future
def cubic_spline_function(
    a: npt.ArrayLike, y_shift: npt.ArrayLike, nodes: npt.ArrayLike
) -> np.ndarray:
    """
    This function defines the cubic spline function used to construct the prior. The spline
    is constructed using the nodes specified in `nodes` and the y-values in `y_shift`. The
    spline is evaluated at the points specified in `a`.

    Parameters
    ----------
    a:  ArrayLike of float
      A one-dimensional array of points at which the function is evaluated.
    y_shift:  ArrayLike of float
      A one-dimensional array whose elements represent the y-value of each bin
    nodes: ArrayLike of float
      A one-dimensional array containing the nodes used to construct the spline.

    Return
    ------
    A one-dimensional array containing the function values evaluated at the points
    specified in `a`.
    """
    from scipy.interpolate import CubicSpline

    cs = CubicSpline(nodes, y_shift)
    return cs(a)


def linear_bin_function(
    a: npt.ArrayLike, y_shift: npt.ArrayLike, bin_edges: npt.ArrayLike
) -> np.ndarray:
    """
    This function defines the linear bin function used to construct the prior. Specifically,
    the prior is constructed using a triangular function whose value at the peak of the node
    is linked to the right and left nodes using a straight line.

    Parameters
    ----------
    a:  ArrayLike of float
      A one-dimensional array of points at which the function is evaluated.
    y_shift:  ArrayLike of float
      A one-dimensional array whose elements represent the y-value of each bin
    bin_nodes: ArrayLike of float
      A one-dimensional array containing the edges of the bins. The bins are
      constructed using pairs of consecutive points.

    Return
    ------
    A one-dimensional array containing the function values evaluated at the points
    specified in `a`.
    """
    res = np.zeros_like(a)
    for shift_pos, shift in enumerate(y_shift):
        if shift_pos > 0 and shift_pos < len(y_shift) - 1:
            bin_low = bin_edges[shift_pos - 1]
            bin_high = bin_edges[shift_pos + 1]
            bin_mid = bin_edges[shift_pos]
            m1 = shift / (bin_mid - bin_low)
            m2 = shift / (bin_high - bin_mid)
        elif shift_pos == 0:  # Left-most bin
            bin_high = bin_edges[shift_pos + 1]
            bin_mid = bin_edges[shift_pos]
            bin_low = bin_mid
            m1 = 0.0
            m2 = shift / (bin_high - bin_mid)
        else:  # Right-most bin
            bin_low = bin_edges[shift_pos - 1]
            bin_mid = bin_edges[shift_pos]
            bin_high = bin_mid
            m1 = shift / (bin_mid - bin_low)
            m2 = 0.0
        cond_low = np.multiply(
            a >= bin_low, a < bin_mid if shift_pos != len(y_shift) - 1 else a <= bin_mid
        )
        cond_high = np.multiply(
            a >= bin_mid, a < bin_high if shift_pos != len(y_shift) - 1 else a <= bin_high
        )
        res = np.add(res, [m1 * (val - bin_low) if cond else 0.0 for val, cond in zip(a, cond_low)])
        res = np.add(
            res, [-m2 * (val - bin_high) if cond else 0.0 for val, cond in zip(a, cond_high)]
        )
    return res


def dis_pc_func(
    delta_h: npt.ArrayLike,
    nodes: npt.ArrayLike,
    x: npt.ArrayLike,
    Q2: npt.ArrayLike,
    pc_func_type: str = "step",
) -> npt.ArrayLike:
    """
    This function builds the functional form of the power corrections for DIS-like processes.
    Power corrections are modelled using a step-function. The edges of the bins used in the
    step-function are specified by the list of nodes. The y-values for each bin are given
    by the array `delta_h`. The power corrections will be computed for the pairs (xb, Q2),
    where `xb` is the Bjorken x. The power correction for DIS processes are rescaled by Q2.

    Parameters
    ----------
    delta_h: ArrayLike
      One-dimensional array containing the shifts for each bin.
    nodes: ArrayLike
      One-dimensional array containing the edges of the bins in x-Bjorken.
    x: ArrayLike
      List of x-Bjorken points at which the power correction function is evaluated.
    Q2: ArrayLike
      List of scales where the power correction function is evaluated.

    Returns
    -------
    A one-dimensional array of power corrections for DIS-like processes where each point is
    evaluated at the kinematic pair (x,Q2).
    """
    if pc_func_type == "step":
        PC = step_function(x, delta_h, nodes) / Q2
    elif pc_func_type == "linear":
        PC = linear_bin_function(x, delta_h, nodes) / Q2
    elif pc_func_type == "cubic":
        PC = cubic_spline_function(x, delta_h, nodes) / Q2
    else:
        raise ValueError(f"Invalid function type: {pc_func_type} is not supported.")

    return PC


def jets_pc_func(
    delta_h: npt.ArrayLike,
    nodes: npt.ArrayLike,
    pT: npt.ArrayLike,
    rap: npt.ArrayLike,
    pc_func_type: str = "step",
) -> npt.ArrayLike:
    """
    Same as `dis_pc_func`, but for jet data. Here, the kinematic pair consists of the rapidity
    `rap` and the transverse momentum `pT`.

    Parameters
    ----------
    delta_h: ArrayLike
      One-dimensional array containing the shifts for each bin.
    nodes: ArrayLike
     One-dimensional array containing the edges of the bins in rapidity.
    rap: ArrayLike
      List of rapidity points at which the power correction is evaluated.
    pT: ArrayLike
      List of pT points at which the power correction is evaluated.

    Returns
    -------
    A one-dimensional array of power corrections for jet processes where each point is
    evaluated at the kinematic pair (y, pT).
    """
    if pc_func_type == "step":
        PC = step_function(rap, delta_h, nodes) / pT
    elif pc_func_type == "linear":
        PC = linear_bin_function(rap, delta_h, nodes) / pT
    elif pc_func_type == "cubic":
        PC = cubic_spline_function(rap, delta_h, nodes) / pT
    else:
        raise ValueError(f"Invalid function type: {pc_func_type} is not supported.")
    return PC


# def jet_single_par(delta_h: float, pT: npt.ArrayLike, rap: npt.ArrayLike) -> npt.ArrayLike:
#     ret = [delta_h for _ in range(rap.size)]
#     return np.array(ret) / pT

# def mult_jet_pc_single_par(dataset_sp, pdf, pc_nodes, pT, rap, pc_func_type: str = "step"):
#     """
#     As mult_jet_pc, but with one single shift for all rapidity bins.

#     This function is meant to be for development purposes only. It will either substitute
#     mult_jet_pc or be deleted in the future."""
#     cuts = dataset_sp.cuts
#     (fkspec,) = dataset_sp.fkspecs
#     fk = fkspec.load_with_cuts(cuts)
#     xsec = central_fk_predictions(fk, pdf)

#     def func(y_values):
#         assert y_values.size == 1
#         ret = [y_values[0] for _ in range(rap.size)]
#         ret = np.array(ret) / pT
#         return np.multiply(ret, xsec.to_numpy()[:, 0])

#     return func


def mult_dis_pc(nodes, x, q2, dataset_sp, pdf, pc_func_type: str = "step"):
    """
    Returns the function that computes the shift to observables due to
    power corrections. Power corrections are treated as multiplicative
    shifts. Hence if `O` is the observable, the prediction is shifted as

      O -> O * (1 + PC),

    and the shift is defined as

      Delta(O) = O * (1 + PC) - O = O * PC.

    This function returns a function that computes the shift
    given the y-values of the nodes used to define the power corrections.
    The interpolation between the nodes is specified by `pc_func_type`.
    """
    cuts = dataset_sp.cuts
    (fkspec,) = dataset_sp.fkspecs
    fk = fkspec.load_with_cuts(cuts)
    th_preds = central_fk_predictions(fk, pdf)

    def func(y_values):
        result = dis_pc_func(y_values, nodes, x, q2, pc_func_type)
        return np.multiply(result, th_preds.to_numpy()[:, 0])

    return func


def mult_dis_ratio_pc(p_nodes, d_nodes, x, q2, dataset_sp, pdf, pc_func_type: str = "step"):
    """
    Returns the function that computes the shift for the ratio of structure
    functions F2_d / F2_p. For this observable, power corrections are defined
    such that

      F2_d / F2_p -> F2_d * (1 + PC2_d) / F2_p * (1 + PC2_p) ,

    and the shift is the defined as

      Delta(F2 ratio) = F2_d / F2_p * (PC2_d - PC2_p) / (1 + PC2_d).

    As for `mult_dis_pc`, this function returns a function that computes the shift
    for the ratio of structure functions F2_d / F2_p given a set of y-values.

    Parameters
    ----------
    p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for the proton (see `dis_pc_func`).
    d_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for the deuteron (see `dis_pc_func`).
    x: list[float]
      Set of points in x-Bjorken where the power corrections will be evaluated.
    q2: list[float]
      Set of points in Q2 where the power corrections will be evaluated.
    dataset_sp: DataSetSpec
      An instance of DataSetSpec used to extract information such as cuts
      and fk tables.
    pdf: PDF
      An instance of the class PDF. This specifies the PDF to bo convoluted
      with the FK tables.

    Returns
    -------
    The function the computes the shift for this observable. It depends on the
    y-values for the parameterization of P2_d and P2_p.
    """
    cuts = dataset_sp.cuts
    fkspec_F2D, fkspec_F2P = dataset_sp.fkspecs
    fk_F2D = fkspec_F2D.load_with_cuts(cuts)
    fk_F2P = fkspec_F2P.load_with_cuts(cuts)
    F2D = central_fk_predictions(fk_F2D, pdf)
    F2P = central_fk_predictions(fk_F2P, pdf)

    F2D = np.concatenate(F2D.values)
    F2P = np.concatenate(F2P.values)
    F2_ratio = operator.truediv(F2D, F2P)

    def func(y_values_p, y_values_d):
        h2d = dis_pc_func(y_values_d, d_nodes, x, q2, pc_func_type)
        h2p = dis_pc_func(y_values_p, p_nodes, x, q2, pc_func_type)
        num = np.sum([h2d, -h2p], axis=0)
        denom = np.sum([np.ones_like(h2p), h2p], axis=0)
        result = np.multiply(operator.truediv(num, denom), F2_ratio)
        return result

    return func

def mult_jet_pc(nodes, pT, rap, dataset_sp, pdf, pc_func_type: str = "step"):
    """
    As `mult_dis_pc`, but for jet data. The power corrections are defined as

      xsec -> xsec * ( 1 + PC ),

    and the shift is defined as

      Delta(xsec) = (xsec + xsec) - xsec = PC.

    The power correction is a function of the rapidity.
    """
    cuts = dataset_sp.cuts
    (fkspec,) = dataset_sp.fkspecs
    fk = fkspec.load_with_cuts(cuts)
    xsec = central_fk_predictions(fk, pdf)

    def func(y_values):
        result = jets_pc_func(y_values, nodes, pT, rap, pc_func_type)
        return np.multiply(result, xsec.to_numpy()[:, 0])

    return func

def construct_pars_combs(parameters_dict):
    """Construct the combination of parameters (the ones that parametrize the power
    corrections) used to compute the shifts.

    Example
    -------
    Given the following dictionary that specifies that power corrections
    ```
    pc_dict = {
      'H1': {'list': [1,2],   'nodes': [0,1]} }
      'H2': {'list': [3,4,5], 'nodes': [0,1,2]} }
    }
    ```
    then this functions constructs a list as follows
    ```
    pars_combs = [
     {'label': 'H1(1)', 'comb': {'H1': [1,0], 'H2': [0,0,0]},
     {'label': 'H1(2)', 'comb': {'H1': [0,1], 'H2': [0,0,0]},
     {'label': 'H2(1)', 'comb': {'H1': [0,0], 'H2': [3,0,0]},
     {'label': 'H2(2)', 'comb': {'H1': [0,0], 'H2': [0,4,0]},
     {'label': 'H2(3)', 'comb': {'H1': [0,0], 'H2': [0,0,5]},
    ]
    ```

    Parameters
    ----------
    """
    combinations = []
    for key, values in parameters_dict.items():
        for i in range(len(values['yshift'])):
            # Create a new dictionary with all keys and zeroed-out values
            new_dict = {k: np.zeros_like(v['yshift']) for k, v in parameters_dict.items()}
            # Set the specific value for the current index
            label = key + f'({i})'
            new_dict[key][i] = values['yshift'][i]
            new_dict = {'label': label, 'comb': new_dict}
            combinations.append(new_dict)

    return combinations

# TODO `pc_func_type` will be removed
def compute_deltas_pc(dataset_sp: DataSetSpec, pdf: PDF, pc_dict: dict, pc_func_type: str):
    """
    Computes the shifts due to power corrections for a single dataset given
    the set of parameters that model the power corrections. The result is
    a dictionary containing as many arrays of shifts as the number of
    combinations of the parameters. For instance, the final dictionary
    may look like:
    ```
    deltas1 = {comb1_label: array_of_shifts1,
               comb2_label: array_of_shifts2,
               comb3_label: array_of_shifts3,
               ...}
    ```
    Note that, as of now, we don't need to specify different prescriptions.
    For that reason, the prescription adopted to construct the shifts is hard
    coded in the function `construct_pars_combs`, and the prescription used to
    compute the sub-matrix is hard-coded in `covmat_power_corrections`.
    """

    exp_name = dataset_sp.name
    process_type = dataset_sp.commondata.metadata.process_type.name
    cuts = dataset_sp.cuts.load()

    pars_combs = construct_pars_combs(pc_dict)
    deltas = defaultdict(list)

    pc_func = None
    if process_type.startswith('DIS'):
        f2_p_nodes = pc_dict["f2p"]['nodes']
        f2_d_nodes = pc_dict["f2d"]['nodes']
        dis_cc_nodes = pc_dict["dis_cc"]['nodes']

        x = dataset_sp.commondata.metadata.load_kinematics()['x'].to_numpy().reshape(-1)[cuts]
        q2 = dataset_sp.commondata.metadata.load_kinematics()['Q2'].to_numpy().reshape(-1)[cuts]

        # F2 ratio
        if exp_name == "NMC_NC_NOTFIXED_EM-F2":
            pc_func_ratio = mult_dis_ratio_pc(f2_p_nodes, f2_d_nodes, x, q2, dataset_sp, pdf, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func_ratio(pars_pc['comb']['f2p'], pars_pc['comb']['f2d'])

        # F2 proton traget
        elif exp_name in F2P_exps:
            pc_func = mult_dis_pc(f2_p_nodes, x, q2, dataset_sp, pdf, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['f2p'])

        # F2 deuteron traget
        elif exp_name in F2D_exps:
            pc_func = mult_dis_pc(f2_d_nodes, x, q2, dataset_sp, pdf, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['f2d'])

        # EMC
        elif exp_name.startswith('EMC_NC_250GEV'):
           raise NotImplementedError(
                f"The {process_type} observable for {exp_name} "
                "has not been implemented."
            )

        # HERA NC xsec
        elif exp_name in np.concatenate([NC_SIGMARED_P_EM, NC_SIGMARED_P_EP, NC_SIGMARED_P_EAVG]):
            pc_func = mult_dis_pc(f2_p_nodes, x, q2, dataset_sp, pdf, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['f2p'])

        # CHORUS
        elif exp_name.startswith('CHORUS_CC'):
            pc_func = mult_dis_pc(dis_cc_nodes, x, q2, dataset_sp, pdf, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['dis_cc'])

        # NuTeV
        elif exp_name.startswith('NUTEV_CC'):
            pc_func = mult_dis_pc(dis_cc_nodes  , x, q2, dataset_sp, pdf, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['dis_cc'])

        # HERA_CC
        elif exp_name.startswith('HERA_CC'):
            pc_func = mult_dis_pc(dis_cc_nodes, x, q2, dataset_sp, pdf, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['dis_cc'])

        else:
            raise ValueError(
                f"The {process_type} observable for {exp_name} " "has not been implemented."
            )

    elif process_type == 'JET':
        pc_jet_nodes = pc_dict["Hj"]['nodes']
        eta = dataset_sp.commondata.metadata.load_kinematics()['y'].to_numpy().reshape(-1)[cuts]
        pT = dataset_sp.commondata.metadata.load_kinematics()['pT'].to_numpy().reshape(-1)[cuts]

        pc_func = mult_jet_pc(pc_jet_nodes, pT, eta, dataset_sp, pdf, pc_func_type)
        for pars_pc in pars_combs:
            deltas[pars_pc['label']] = pc_func(pars_pc['comb']['Hj'])

    elif process_type == 'DIJET':

        if dataset_sp.commondata.metadata.experiment == 'ATLAS':
            pc_jet_nodes = pc_dict["H2j_ATLAS"]['nodes']
            eta_star = (
                dataset_sp.commondata.metadata.load_kinematics()['ystar']
                .to_numpy()
                .reshape(-1)[cuts]
            )
            m_jj = (
                dataset_sp.commondata.metadata.load_kinematics()['m_jj']
                .to_numpy()
                .reshape(-1)[cuts]
            )
            pc_func = mult_jet_pc(pc_jet_nodes, m_jj, eta_star, dataset_sp, pdf, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['H2j_ATLAS'])

        elif dataset_sp.commondata.metadata.experiment == 'CMS':
            pc_jet_nodes = pc_dict["H2j_CMS"]['nodes']
            eta_diff = (
                dataset_sp.commondata.metadata.load_kinematics()['ydiff']
                .to_numpy()
                .reshape(-1)[cuts]
            )
            m_jj = (
                dataset_sp.commondata.metadata.load_kinematics()['m_jj']
                .to_numpy()
                .reshape(-1)[cuts]
            )
            pc_func = mult_jet_pc(pc_jet_nodes, m_jj, eta_diff, dataset_sp, pdf, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['H2j_CMS'])

        else:
            raise ValueError(
                f"{dataset_sp.commondata.metadata.experiment} is not implemented for DIJET."
            )

    else:
        raise RuntimeError(f"{process_type} has not been implemented.")

    return deltas
