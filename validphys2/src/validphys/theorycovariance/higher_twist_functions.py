"""
Utilities for the computation of shifts in theoretical predictions due to power
corrections. Contrary to scale variations, the multiplicative shifts due to
power corrections are computed at runtime during vp-setupfit. Once the shifts
are computed using the function ``compute_deltas_pc``, the covmat can be
constructed using the same functions implemented for scale variations (e.g.
`covs_pt_prescrip`).

In the case of power corrections, shifts and covmat constructions are computed
using a 5-point prescription extended to every parameter used to define the
power correction.

This module extracts the dataset-to-power-correction-type routing into a
reusable function ``get_pc_type``. It also comprehends a bunch of
``factory`` functions such as `mult_dis_pc`. Each of these functions returns
another function that computes the shifts taking as arguments the values of the
parameters used to parametrise the power corrections. In other words, these
factory functions hard-code the dependence on the kinematic and leave the
dependence on the parameters free. In this way, the shifts can be computed using
different combinations of parameters (i.e. different prescriptions) if needed.

"""

from collections import defaultdict
import operator
from typing import Optional, Tuple, Union

import numpy as np
import numpy.typing as npt

from validphys.convolution import central_fk_predictions
from validphys.core import PDF, DataSetSpec

# ---------------------------------------------------------------------------
# Dataset name constants
# ---------------------------------------------------------------------------

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

# Combined list for convenience (avoids repeated np.concatenate)
NC_SIGMARED_P_ALL = NC_SIGMARED_P_EM + NC_SIGMARED_P_EP + NC_SIGMARED_P_EAVG

# CC experiment prefixes (merged from three identical branches)
DIS_CC_PREFIXES = ('CHORUS_CC', 'NUTEV_CC', 'HERA_CC')

# Dijet rapidity variable name per experiment
DIJET_RAPIDITY_VAR = {'ATLAS': 'ystar', 'CMS': 'ydiff'}


# ---------------------------------------------------------------------------
# Dataset -> PC type routing
# ---------------------------------------------------------------------------

PCTypeResult = Union[str, Tuple[str, str]]


def get_pc_type(
    exp_name: str,
    process_type: str,
    experiment: Optional[str] = None,
    pc_dict: Optional[dict] = None,
) -> PCTypeResult:
    """
    Determine which power correction type(s) apply to a given dataset.

    This is the single source of truth for mapping dataset names and process
    types to power correction parameter keys.

    Parameters
    ----------
    exp_name : str
        Dataset name (e.g. ``'BCDMS_NC_NOTFIXED_P_EM-F2'``).
    process_type : str
        Process type string (e.g. ``'DIS_NCE'``, ``'JET'``, ``'DIJET'``).
    experiment : str, optional
        Experiment name, required for DIJET (e.g. ``'ATLAS'``, ``'CMS'``).
    pc_dict : dict, optional
        Power correction parameter dictionary. When provided, enables
        fallback logic for DIJET (``H2j_ATLAS`` -> ``H2j`` if the
        experiment-specific key is absent).

    Returns
    -------
    str or tuple of (str, str)
        The PC type key(s). For the NMC ratio dataset
        (``NMC_NC_NOTFIXED_EM-F2``), returns ``("f2p", "f2d")``.
        For all other datasets, returns a single string key.

    Raises
    ------
    NotImplementedError
        For EMC datasets (not yet implemented).
    ValueError
        For unknown DIS experiments or unknown DIJET experiments.
    RuntimeError
        For unsupported process types.
    """
    if process_type.startswith('DIS'):
        return _get_dis_pc_type(exp_name)
    elif process_type == 'JET':
        return 'Hj'
    elif process_type == 'DIJET':
        if experiment is None:
            raise ValueError("The 'experiment' argument is required for DIJET process type.")
        return _get_dijet_pc_type(experiment, pc_dict)
    else:
        raise RuntimeError(f"{process_type} has not been implemented.")


def _get_dis_pc_type(exp_name: str) -> PCTypeResult:
    """Resolve DIS dataset name to PC type key(s)."""
    if exp_name == "NMC_NC_NOTFIXED_EM-F2":
        return ("f2p", "f2d")
    if exp_name in F2P_exps:
        return "f2p"
    if exp_name in F2D_exps:
        return "f2d"
    if exp_name.startswith('EMC_NC_250GEV'):
        raise NotImplementedError(f"The DIS observable for {exp_name} has not been implemented.")
    if exp_name in NC_SIGMARED_P_ALL:
        return "f2p"
    if any(exp_name.startswith(prefix) for prefix in DIS_CC_PREFIXES):
        return "dis_cc"
    raise ValueError(f"The DIS observable for {exp_name} has not been implemented.")


def _get_dijet_pc_type(experiment: str, pc_dict: Optional[dict] = None) -> str:
    """Resolve DIJET experiment to PC type key with optional fallback."""
    if experiment == 'ATLAS':
        specific_key = "H2j_ATLAS"
    elif experiment == 'CMS':
        specific_key = "H2j_CMS"
    else:
        raise ValueError(f"{experiment} is not implemented for DIJET.")

    if pc_dict is not None and not pc_dict.get(specific_key):
        return "H2j"
    return specific_key


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
    delta_h: npt.ArrayLike, nodes: npt.ArrayLike, x: npt.ArrayLike, Q2: npt.ArrayLike
) -> npt.ArrayLike:
    """
    This function builds the functional form of the power corrections for DIS-like processes.
    Power corrections are modelled using a linear function, which interpolates between the nodes
    of the parametrisation. The y-values for each node are given
    by the array `delta_h`. The power corrections will be computed for the pairs (xb, Q2),
    where `xb` is the Bjorken x. The power correction for DIS processes is divided by Q2.

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
    PC = linear_bin_function(x, delta_h, nodes) / Q2
    return PC


def jets_pc_func(
    delta_h: npt.ArrayLike, nodes: npt.ArrayLike, pT: npt.ArrayLike, rap: npt.ArrayLike
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
    PC = linear_bin_function(rap, delta_h, nodes) / pT
    return PC


def mult_dis_pc(nodes, x, q2, dataset_sp, pdf):
    """
    Returns the function that computes the shift to observables due to
    power corrections. Power corrections are treated as multiplicative
    shifts. Hence if `O` is the observable, the prediction is shifted as

      O -> O * (1 + PC),

    and the shift is defined as

      Delta(O) = O * (1 + PC) - O = O * PC.

    This function returns a function that computes the shift
    given the y-values of the nodes used to define the power corrections.
    """
    cuts = dataset_sp.cuts
    (fkspec,) = dataset_sp.fkspecs
    fk = fkspec.load_with_cuts(cuts)
    th_preds = central_fk_predictions(fk, pdf)

    def func(y_values):
        result = dis_pc_func(y_values, nodes, x, q2)
        return np.multiply(result, th_preds.to_numpy()[:, 0])

    return func


def mult_dis_ratio_pc(p_nodes, d_nodes, x, q2, dataset_sp, pdf):
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
        h2d = dis_pc_func(y_values_d, d_nodes, x, q2)
        h2p = dis_pc_func(y_values_p, p_nodes, x, q2)
        num = np.sum([h2d, -h2p], axis=0)
        denom = np.sum([np.ones_like(h2p), h2p], axis=0)
        result = np.multiply(operator.truediv(num, denom), F2_ratio)
        return result

    return func


def mult_jet_pc(nodes, pT, rap, dataset_sp, pdf):
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
        result = jets_pc_func(y_values, nodes, pT, rap)
        return np.multiply(result, xsec.to_numpy()[:, 0])

    return func


def construct_pars_combs(parameters_dict: dict) -> list[dict]:
    """Construct the combination of parameters (the ones that parametrize the power
    corrections) used to compute the shifts.

    Example
    -------
    Given the following dictionary that specifies the power corrections::

        pc_dict = {
          'H1': {'yshift': [1,2],   'nodes': [1, 2]},
          'H2': {'yshift': [3,4,5], 'nodes': [1, 2, 3]},
        }

    then this function constructs a list as follows::

        pars_combs = [
         {'label': 'H1(0)', 'comb': {'H1': [1,0], 'H2': [0,0,0]}},
         {'label': 'H1(1)', 'comb': {'H1': [0,2], 'H2': [0,0,0]}},
         {'label': 'H2(0)', 'comb': {'H1': [0,0], 'H2': [3,0,0]}},
         {'label': 'H2(1)', 'comb': {'H1': [0,0], 'H2': [0,4,0]}},
         {'label': 'H2(2)', 'comb': {'H1': [0,0], 'H2': [0,0,5]}},
        ]
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


def _apply_pars_combs(pars_combs, pc_key, pc_func):
    """Apply parameter combinations to a PC function, returning deltas dict.

    This eliminates the repeated pattern::

        for pars_pc in pars_combs:
            deltas[pars_pc['label']] = pc_func(pars_pc['comb'][key])

    Parameters
    ----------
    pars_combs : list of dict
        Output of ``construct_pars_combs``.
    pc_key : str
        The key into ``pars_pc['comb']`` to extract the parameter values.
    pc_func : callable
        A function that takes a single array of y-values and returns shifts.

    Returns
    -------
    dict
        ``{label: array_of_shifts}``
    """
    deltas = {}
    for pars_pc in pars_combs:
        deltas[pars_pc['label']] = pc_func(pars_pc['comb'][pc_key])
    return deltas


def compute_deltas_pc(dataset_sp: DataSetSpec, pdf: PDF, pc_dict: dict):
    """
    Computes the shifts due to power corrections for a single dataset given
    the set of parameters that model the power corrections. The result is
    a dictionary containing as many arrays of shifts as the number of
    combinations of the parameters. For instance, the final dictionary
    may look like::

        deltas1 = {comb1_label: array_of_shifts1,
                   comb2_label: array_of_shifts2,
                   comb3_label: array_of_shifts3,
                   ...}

    Note that, as of now, we don't need to specify different prescriptions.
    Hence, the prescription adopted to construct the shifts is hard
    coded in the function ``construct_pars_combs``, and the prescription used to
    compute the sub-matrix is hard-coded in ``covmat_power_corrections``.
    """
    exp_name = dataset_sp.name
    process_type = dataset_sp.commondata.metadata.process_type.name
    experiment = getattr(dataset_sp.commondata.metadata, 'experiment', None)
    cuts = dataset_sp.cuts.load()

    pars_combs = construct_pars_combs(pc_dict)
    deltas = defaultdict(list)

    pc_type = get_pc_type(exp_name, process_type, experiment=experiment, pc_dict=pc_dict)

    if process_type.startswith('DIS'):
        x = dataset_sp.commondata.metadata.load_kinematics()['x'].to_numpy().reshape(-1)[cuts]
        q2 = dataset_sp.commondata.metadata.load_kinematics()['Q2'].to_numpy().reshape(-1)[cuts]

        # NMC ratio: special case with two PC types (proton and deuteron)
        if isinstance(pc_type, tuple):
            p_key, d_key = pc_type
            f2_p_nodes = pc_dict[p_key]['nodes']
            f2_d_nodes = pc_dict[d_key]['nodes']
            pc_func_ratio = mult_dis_ratio_pc(f2_p_nodes, f2_d_nodes, x, q2, dataset_sp, pdf)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func_ratio(
                    pars_pc['comb'][p_key], pars_pc['comb'][d_key]
                )
        else:
            nodes = pc_dict[pc_type]['nodes']
            pc_func = mult_dis_pc(nodes, x, q2, dataset_sp, pdf)
            deltas.update(_apply_pars_combs(pars_combs, pc_type, pc_func))

    elif process_type == 'JET':
        eta = dataset_sp.commondata.metadata.load_kinematics()['y'].to_numpy().reshape(-1)[cuts]
        pT = dataset_sp.commondata.metadata.load_kinematics()['pT'].to_numpy().reshape(-1)[cuts]
        nodes = pc_dict[pc_type]['nodes']
        pc_func = mult_jet_pc(nodes, pT, eta, dataset_sp, pdf)
        deltas.update(_apply_pars_combs(pars_combs, pc_type, pc_func))

    elif process_type == 'DIJET':
        rap_var = DIJET_RAPIDITY_VAR.get(experiment)
        if rap_var is None:
            raise ValueError(f"{experiment} is not implemented for DIJET.")

        kinematics = dataset_sp.commondata.metadata.load_kinematics()
        rap = kinematics[rap_var].to_numpy().reshape(-1)[cuts]
        m_jj = kinematics['m_jj'].to_numpy().reshape(-1)[cuts]

        nodes = pc_dict[pc_type]['nodes']
        pc_func = mult_jet_pc(nodes, m_jj, rap, dataset_sp, pdf)
        deltas.update(_apply_pars_combs(pars_combs, pc_type, pc_func))

    else:
        raise RuntimeError(f"{process_type} has not been implemented.")

    return deltas
