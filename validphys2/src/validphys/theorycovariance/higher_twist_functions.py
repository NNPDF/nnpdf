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

This module comprehends a bunch of ``factory`` functions such as `DIS_F2_pc`. Each
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

GEV_CM2_CONV = 3.893793e10
GF = 1.1663787e-05  # Fermi's constant [GeV^-2]
Mh = 0.938  # Proton's mass in GeV/c^2
MW = 80.398  # W boson mass in GeV/c^2

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


def jet_single_par(delta_h: float, pT: npt.ArrayLike, rap: npt.ArrayLike) -> npt.ArrayLike:
    ret = [delta_h for _ in range(rap.size)]
    return np.array(ret) / pT


def JET_pc(dataset_sp, pdf, pc_nodes, pT, rap, pc_func_type: str = "step"):
    """
    Returns the function that computes the shift for the ratio for single
    jet cross sections. In particular, the shift is computed such that

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
        result = jets_pc_func(y_values, pc_nodes, pT, rap, pc_func_type)
        return np.multiply(result, xsec.to_numpy()[:, 0])

    return func


def JET_pc_single_par(dataset_sp, pdf, pc_nodes, pT, rap, pc_func_type: str = "step"):
    """
    As JET_pc, but with one single shift for all rapidity bins.

    This function is meant to be for development purposes only. It will either substitute
    JET_pc or be deleted in the future."""
    cuts = dataset_sp.cuts
    (fkspec,) = dataset_sp.fkspecs
    fk = fkspec.load_with_cuts(cuts)
    xsec = central_fk_predictions(fk, pdf)

    def func(y_values):
        assert y_values.size == 1
        ret = [y_values[0] for _ in range(rap.size)]
        ret = np.array(ret) / pT
        return np.multiply(ret, xsec.to_numpy()[:, 0])

    return func


# TODO Maybe we want to treat the function that parametrizes the PC
# as argument?
def DIS_F2_pc(pc2_nodes, x, q2, pc_func_type: str = "step"):
    """
    Returns the function that computes the shift for the ratio of structure
    functions F2_d / F2_p. For this observable, power corrections are defined
    such that

      F2 -> F2 + PC2,

    and the shift is defined as

      Delta(F2) = (F2 + PC2) - F2 = PC2.

    Note that, as in the case of `DIS_F2R_ht`, the shift still depends on the set
    of parameters needed to define the parameterization of PC2. Also, this function
    can be used to compute the shift for both proton and deuteron, provided that the
    correct list of parameters is passed to the **curried** function.

    The function used to parametrize the the power correction is `dis_pc_func` and
    it is hard coded.

    Parameters
    ----------

    """

    def PC_2(y_values):
        result = dis_pc_func(y_values, pc2_nodes, x, q2, pc_func_type)
        return result

    return PC_2


def DIS_F2R_pc(experiment, pdf, pc_2_p_nodes, pc_2_d_nodes, x, q2, pc_func_type: str = "step"):
    """
    Returns the function that computes the shift for the ratio of structure
    functions F2_d / F2_p. For this observable, power corrections are defined
    such that

      F2_d / F2_p -> (F2_d + PC2_d) / (F2_p + PC2_p) ,

    and the shift is the defined as

      Delta(F2 ratio) = (F2_d + PC2_d) / (F2_p + PC2_p) - F2_d / F2_p .

    The shift is computed for a given set of kinematic variables specified
    by the paris (x,Q2), but it still depends on the set of parameters need by
    the power correction terms PC2_d and PC2_p.

    Note that this function does **not** return the power corrections for the
    given kinematic points, but rather it returns another function where the
    kinematic dependence has been **curried** (i.e. hard coded). This new function
    takes as arguments the y-values of the nodes used to compute PC2_d and PC2_p
    (see `delta_h` in `dis_pc_func`). Note that these y-values are not necessarily
    the values listed in the runcard, as we can apply different point prescription.
    For instance, we may want to pass in a set of y-values where the nodes are shifted
    one at the time, leaving the others zero. The prescription is thus handled separately.
    The returning function allows thus to compute  Delta(F2 ratio)({...}_d, {...}_p), where
    `{...}_d` and `{...}_p` are the sets of y-values for the parametrisation for the proton
    and deuteron terms respectively.

    The function used to parametrize the the power correction is `dis_pc_func` and
    it is hard coded.

    Parameters
    ----------
    experiment: DataSetSpec
      An instance of DataSetSpec used to extract information such as cuts
      and fk tables.
    pdf: PDF
      An instance of the class PDF. This specifies the PDF to bo convoluted
      with the FK tables.
    pc_2_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for the proton (see `dis_pc_func`).
    pc_2_d_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for the deuteron (see `dis_pc_func`).
    x: list[float]
      Set of points in x-Bjorken where the power corrections will be evaluated.
    q2: list[float]
      Set of points in Q2 where the power corrections will be evaluated.

    Returns
    -------
    The function the computes the shift for this observable. It depends on the
    y-values for the parameterization of P2_d and P2_p.
    """
    cuts = experiment.cuts
    fkspec_F2D, fkspec_F2P = experiment.fkspecs
    fk_F2D = fkspec_F2D.load_with_cuts(cuts)
    fk_F2P = fkspec_F2P.load_with_cuts(cuts)
    F2D = central_fk_predictions(fk_F2D, pdf)
    F2P = central_fk_predictions(fk_F2P, pdf)

    F2D = np.concatenate(F2D.values)
    F2P = np.concatenate(F2P.values)
    F2_ratio = operator.truediv(F2D, F2P)

    def func(y_values_d, y_values_p):
        PC_d = dis_pc_func(y_values_d, pc_2_d_nodes, x, q2, pc_func_type)
        PC_p = dis_pc_func(y_values_p, pc_2_p_nodes, x, q2, pc_func_type)
        num = np.sum([F2D, PC_d], axis=0)
        denom = np.sum([F2P, PC_p], axis=0)
        result = np.array(operator.truediv(num, denom) - F2_ratio)
        return result

    return func


def DIS_F2C_pc(pc2_p_nodes, pc2_d_nodes, x, q2, pc_func_type: str = "step"):
    """
    Builds the function used to compute the shifts for the charm structure
    function measured by EMC. The process involved is

        mu^+ + Fe -> mu+^ + c cbar + X .

    This function works exactly as the previous functions used to compute
    nuisance shifts. In this case, the constructed function (`func` below)
    requires two lists of parameters for the proton and the deuteron
    contribution. The reason being that in this process the muon scatters off an
    iron target, and the power correction contribution is a mixture of proton
    and deuteron nucleons. Hence, proton and deuteron contribution are weighted
    by the appropriate atomic factor.

    Note that we are parametrising power corrections as proton and deuteron
    targets. If we were to parametrize such contributions using, say, proton and
    nucleon, than the weights would change.


    Nuclear target
    --------------
    The power corrections for nuclear observables, like in this case, are
    affected by the pc contribution of the protons and that of the neutrons. If
    we allow for the non-isoscalarity of the target, and combine the two
    contributions in accordance with the atomic and mass number (A and Z
    respectively), the power correction for the nuclear target can be written as
    (see  eq.(4.2.5) in
    https://nnpdf.mi.infn.it/wp-content/uploads/2021/09/thesis_master_RP.pdf)

      PC_N = 1/A (Z * PC_p + (A-Z) * PC_n) .

    The deuteron is obtained using the isoscalarity, namely

      PC_d = 1/2 (PC_p + PC_n) .

    Since we parametrise the power corrections of the proton and the deuteron,
    we can combine the above equations and write

      PC_N = 1/A * ( PC_p * (2Z - A) + 2 * PC_d * (A - Z) )

    Parameters
    ----------
    pc2_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for the proton (see `dis_pc_func`).
    pc2_d_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for the deuteron (see `dis_pc_func`).
    x: list[float]
      Set of points in x-Bjorken where the power corrections will be evaluated.
    q2: list[float]
      Set of points in Q2 where the power corrections will be evaluated.

    Returns
    -------
    The function the computes the shift for this observable. It depends on the
    y-values for the parameterization of P2_d and P2_p.
    """
    # Iron target
    Z = 23.403
    A = 49.618

    def func(y_values_d, y_values_p):
        PC2_d = dis_pc_func(y_values_d, pc2_d_nodes, x, q2, pc_func_type)
        PC2_p = dis_pc_func(y_values_p, pc2_p_nodes, x, q2, pc_func_type)
        result = (2 * Z - A) / A * PC2_p + 2 * (A - Z) / A * PC2_d
        return result

    return func


def DIS_NC_XSEC_pc(pc2_nodes, pcL_nodes, pc3_nodes, lepton, x, q2, y, pc_func_type: str = "step"):
    """
    Builds the function used to compute the shifts for the DIS NC x-secs
    delivered by HERA and NMC. The x-sec is reconstructed as calculated
    in Yadism (see https://yadism.readthedocs.io/en/latest/theory/intro.html).
    In particular, the x-sec is a linear combination of the structure functions
    F_2, F_L, and F_3. The coefficients are also computed appropriately (see
    link). The contribution of the power corrections is then

      Delta(x-sec) = x-sec_w_pc - x-sec_wo_pc =  PC_2 + N_L * PC_L + N_3 * PC_3

    where PC_i are the power corrections relative to the respective structure
    functions and the N_i the respective coefficients (as defined in Yadism).

    This function works exactly as the previous functions used to
    compute nuisance shifts. In addition, it requires the kinematic
    invariant `y` to build the shift-function.

    Note that this function can be used for both proton and deuteron targets,
    provided that the appropriate lists of nodes is given.

    Parameters
    ----------
    pc2_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_2.
    pcL_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_L.
    pc3_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_3.
    lepton: int
      Whether the scattering particle is a lepton (0) or an anti-lepton(1).
    x: list[float]
      Set of points in x-Bjorken where the power corrections will be evaluated.
    q2: list[float]
      Set of points in Q2 where the power corrections will be evaluated.
    y: list[float]
      Set of points in y where the power corrections will be evaluated.

    Returns
    -------
    The function the computes the shift for this observable. It depends on the
    y-values for the parameterization of P2 and PL and P3.
    """
    yp = 1 + np.power(1 - y, 2)
    ym = 1 - np.power(1 - y, 2)
    yL = np.power(y, 2)
    N_L = -yL / yp  # Coefficient for F_L
    N_3 = np.power(-1, lepton) * ym / yp  # Coefficient for F_3

    def func(y_values_pc2, y_values_pcL, y_values_pc3):
        PC_2 = dis_pc_func(y_values_pc2, pc2_nodes, x, q2, pc_func_type)
        PC_L = dis_pc_func(y_values_pcL, pcL_nodes, x, q2, pc_func_type)
        PC_3 = dis_pc_func(y_values_pc3, pc3_nodes, x, q2, pc_func_type)
        result = PC_2 + N_L * PC_L + N_3 * PC_3
        return result

    return func


def DIS_CC_HERA_XSEC_pc(
    pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, lepton, x, q2, y, pc_func_type: str = "step"
):
    """
    Builds the function used to compute the shifts for the DIS CC x-secs
    delivered by HERA. The x-sec is reconstructed as calculated
    in Yadism (see https://yadism.readthedocs.io/en/latest/theory/intro.html).
    In particular, the x-sec is a linear combination of the structure functions
    F_2, F_L, and F_3. The coefficients are also computed appropriately (see
    link). The contribution of the power corrections is then

      Delta(x-sec) = x-sec_w_pc - x-sec_wo_pc =  N * (PC_2 + N_L * PC_L + N_3 * PC_3)

    where PC_i are the power corrections relative to the respective structure
    functions and the N_i the respective coefficients (as defined in Yadism).
    N is the overall normalization factor.

    For the HERA_CC_318GEV dataset, the target is always a proton. However, the
    lepton may be either the electron (0) or the positron (1). This information
    is needed in order to compute the coefficient N_3.

    Parameters
    ----------
    pc2_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_2.
    pcL_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_L.
    pc3_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_3.
    lepton: int
      Whether the scattering particle is a lepton (0) or an anti-lepton(1).
    x: list[float]
      Set of points in x-Bjorken where the power corrections will be evaluated.
    q2: list[float]
      Set of points in Q2 where the power corrections will be evaluated.
    y: list[float]
      Set of points in y where the power corrections will be evaluated.

    Returns
    -------
    The function the computes the shift for this observable. It depends on the
    y-values for the parameterization of P2 and PL and P3.
    """
    yp = 1 + np.power(1 - y, 2)
    ym = 1 - np.power(1 - y, 2)
    yL = np.power(y, 2)
    N = 1 / 4 * yp  # Overall normalization
    N_L = -yL / yp  # Coefficient for F_L
    N_3 = np.power(-1, lepton) * ym / yp  # Coefficient for F_3

    def func(y_values_pc2_p, y_values_pcL_p, y_values_pc3_p):
        # Initialize power corrections for each structure function
        PC2_p = dis_pc_func(y_values_pc2_p, pc2_p_nodes, x, q2, pc_func_type)
        PCL_p = dis_pc_func(y_values_pcL_p, pcL_p_nodes, x, q2, pc_func_type)
        PC3_p = dis_pc_func(y_values_pc3_p, pc3_p_nodes, x, q2, pc_func_type)

        # Build the contribution to the x-sec of the power corrections
        result = N * (PC2_p + N_L * PCL_p + N_3 * PC3_p)
        return result

    return func


def DIS_CC_NUTEV_pc(
    pc2_p_nodes,
    pcL_p_nodes,
    pc3_p_nodes,
    pc2_d_nodes,
    pcL_d_nodes,
    pc3_d_nodes,
    lepton,
    x,
    q2,
    y,
    pc_func_type: str = "step",
):
    """
    Builds the function used to compute the shifts for the DIS CC x-secs
    delivered by NuTeV. The x-sec is reconstructed as calculated
    in Yadism (see https://yadism.readthedocs.io/en/latest/theory/intro.html).
    In particular, the x-sec is a linear combination of the structure functions
    F_2, F_L, and F_3. The coefficients are also computed appropriately (see
    link). Note that this experiment uses iron targets, and thus the coefficients
    must take into account the nuclear mixture of porton and deuteron. The contribution
    of the power corrections is then

      Delta(x-sec) = x-sec_w_pc - x-sec_wo_pc =  N * (PC_2 + N_L * PC_L + N_3 * PC_3)

    where PC_i are the power corrections relative to the respective structure
    functions (nuclear mixture implicit) and the N_i the respective coefficients (as defined in Yadism).
    N is the overall normalization factor.

    For the NuTeV CC dataset, the target is always iron. However, the
    lepton may be either the electron (0) or the positron (1). This
    information is needed in order to compute the coefficient N_3.

    Nuclear target
    --------------
    See `DIS_F2C_pc`.

    Parameters
    ----------
    pc2_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_2 of the proton.
    pcL_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_L of the proton.
    pc3_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_3 of the proton.
    pc2_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_2 of the deuteron.
    pcL_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_L of the deuteron.
    pc3_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_3 of the deuteron.
    lepton: int
      Whether the scattering particle is a lepton (0) or an anti-lepton(1).
    x: list[float]
      Set of points in x-Bjorken where the power corrections will be evaluated.
    q2: list[float]
      Set of points in Q2 where the power corrections will be evaluated.
    y: list[float]
      Set of points in y where the power corrections will be evaluated.

    Returns
    -------
    The function the computes the shift for this observable. It depends on the
    y-values for the parameterization of P2 and PL and P3 for proton and deuteron.
    """
    # Iron target
    Z = 23.403
    A = 49.618
    yp = 1 + np.power(1 - y, 2) - 2 * np.power(x * y * Mh, 2) / q2
    ym = 1 - np.power(1 - y, 2)
    yL = np.power(y, 2)
    N_L = -yL / yp  # Coefficient for F_L
    N_3 = np.power(-1, lepton) * ym / yp  # Coefficient for F_3

    MW2 = np.power(MW, 2)
    # Overall coefficient
    # TODO: cross-check
    N = 100 / 2 / np.power(1 + q2 / MW2, 2) * yp

    def func(
        y_values_pc2_p,
        y_values_pcL_p,
        y_values_pc3_p,
        y_values_pc2_d,
        y_values_pcL_d,
        y_values_pc3_d,
    ):
        PC2_p = dis_pc_func(y_values_pc2_p, pc2_p_nodes, x, q2, pc_func_type)
        PCL_p = dis_pc_func(y_values_pcL_p, pcL_p_nodes, x, q2, pc_func_type)
        PC3_p = dis_pc_func(y_values_pc3_p, pc3_p_nodes, x, q2, pc_func_type)
        PC2_d = dis_pc_func(y_values_pc2_d, pc2_d_nodes, x, q2, pc_func_type)
        PCL_d = dis_pc_func(y_values_pcL_d, pcL_d_nodes, x, q2, pc_func_type)
        PC3_d = dis_pc_func(y_values_pc3_d, pc3_d_nodes, x, q2, pc_func_type)
        tmp_2 = (2 * Z - A) / A * PC2_p + 2 * (A - Z) / A * PC2_d
        tmp_L = (2 * Z - A) / A * PCL_p + 2 * (A - Z) / A * PCL_d
        tmp_3 = (2 * Z - A) / A * PC3_p + 2 * (A - Z) / A * PC3_d
        result = N * (tmp_2 + N_L * tmp_L + N_3 * tmp_3)
        return result

    return func


# TODO This is function is really similar to the one
# defined for NUTEV CC. Can we reduce code repetitions?
def DIS_CC_CHORUS_pc(
    pc2_p_nodes,
    pcL_p_nodes,
    pc3_p_nodes,
    pc2_d_nodes,
    pcL_d_nodes,
    pc3_d_nodes,
    lepton,
    x,
    q2,
    y,
    pc_func_type: str = "step",
):
    """
    Same as DIS_CC_NUTEV_pc, but for CHORUS CC.

    Note that the difference here is in the definition of the overall
    normalization N.

    Nuclear target
    --------------
    See `DIS_F2C_pc`.

    Parameters
    ----------
    pc2_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_2 of the proton.
    pcL_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_L of the proton.
    pc3_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_3 of the proton.
    pc2_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_2 of the deuteron.
    pcL_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_L of the deuteron.
    pc3_p_nodes: list[float]
      The list of nodes in x-Bjorken used to define the parametrization of the
      power correction for F_3 of the deuteron.
    lepton: int
      Whether the scattering particle is a lepton (0) or an anti-lepton(1).
    x: list[float]
      Set of points in x-Bjorken where the power corrections will be evaluated.
    q2: list[float]
      Set of points in Q2 where the power corrections will be evaluated.
    y: list[float]
      Set of points in y where the power corrections will be evaluated.

    Returns
    -------
    The function the computes the shift for this observable. It depends on the
    y-values for the parameterization of P2 and PL and P3 for proton and deuteron.
    """
    # Lead target
    A = 208.0
    Z = 82
    yp = 1 + np.power(1 - y, 2) - 2 * np.power(x * y * Mh, 2) / q2
    ym = 1 - np.power(1 - y, 2)
    yL = np.power(y, 2)
    N_L = -yL / yp  # Coefficient for F_L
    N_3 = np.power(-1, lepton) * ym / yp  # Coefficient for F_3

    MW2 = np.power(MW, 2)
    # Overall coefficient
    # TODO: cross-check
    N = GEV_CM2_CONV * (GF**2) * Mh / (2 * np.pi * np.power(1 + q2 / MW2, 2)) * yp

    def func(
        y_values_pc2_p,
        y_values_pcL_p,
        y_values_pc3_p,
        y_values_pc2_d,
        y_values_pcL_d,
        y_values_pc3_d,
    ):
        PC2_p = dis_pc_func(y_values_pc2_p, pc2_p_nodes, x, q2, pc_func_type)
        PCL_p = dis_pc_func(y_values_pcL_p, pcL_p_nodes, x, q2, pc_func_type)
        PC3_p = dis_pc_func(y_values_pc3_p, pc3_p_nodes, x, q2, pc_func_type)
        PC2_d = dis_pc_func(y_values_pc2_d, pc2_d_nodes, x, q2, pc_func_type)
        PCL_d = dis_pc_func(y_values_pcL_d, pcL_d_nodes, x, q2, pc_func_type)
        PC3_d = dis_pc_func(y_values_pc3_d, pc3_d_nodes, x, q2, pc_func_type)
        tmp_2 = (2 * Z - A) / A * PC2_p + 2 * (A - Z) / A * PC2_d
        tmp_L = (2 * Z - A) / A * PCL_p + 2 * (A - Z) / A * PCL_d
        tmp_3 = (2 * Z - A) / A * PC3_p + 2 * (A - Z) / A * PC3_d
        result = N * (tmp_2 + N_L * tmp_L + N_3 * tmp_3)
        return result

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
        pc2_p_nodes = pc_dict["H2p"]['nodes']
        pcL_p_nodes = pc_dict["HLp"]['nodes']
        pc3_p_nodes = pc_dict["H3p"]['nodes']
        pc2_d_nodes = pc_dict["H2d"]['nodes']
        pcL_d_nodes = pc_dict["HLd"]['nodes']
        pc3_d_nodes = pc_dict["H3d"]['nodes']
        x = dataset_sp.commondata.metadata.load_kinematics()['x'].to_numpy().reshape(-1)[cuts]
        q2 = dataset_sp.commondata.metadata.load_kinematics()['Q2'].to_numpy().reshape(-1)[cuts]
        y = dataset_sp.commondata.metadata.load_kinematics()['y'].to_numpy().reshape(-1)[cuts]

        # F2 ratio
        if exp_name == "NMC_NC_NOTFIXED_EM-F2":
            pc_func = DIS_F2R_pc(dataset_sp, pdf, pc2_p_nodes, pc2_d_nodes, x, q2, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['H2p'], pars_pc['comb']['H2d'])

        # F2 proton traget
        elif exp_name in F2P_exps:
            pc_func = DIS_F2_pc(pc2_p_nodes, x, q2, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['H2p'])

        # F2 deuteron traget
        elif exp_name in F2D_exps:
            pc_func = DIS_F2_pc(pc2_d_nodes, x, q2, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['H2d'])

        # EMC
        elif exp_name.startswith('EMC_NC_250GEV'):
            pc_func = DIS_F2C_pc(pc2_p_nodes, pc2_d_nodes, x, q2, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['H2p'], pars_pc['comb']['H2d'])

        # HERA and NMC SIGMARED NC
        elif exp_name in np.concatenate([NC_SIGMARED_P_EM, NC_SIGMARED_P_EP, NC_SIGMARED_P_EAVG]):
            # Electron
            if exp_name in NC_SIGMARED_P_EM:
                pc_func = DIS_NC_XSEC_pc(
                    pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 0, x, q2, y, pc_func_type
                )
            # Positron
            elif exp_name in NC_SIGMARED_P_EP:
                pc_func = DIS_NC_XSEC_pc(
                    pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 1, x, q2, y, pc_func_type
                )
            # Average positron and electron
            # TODO
            # Check if this is correct (ach)
            elif NC_SIGMARED_P_EAVG:

                def average(y_values_pc2_p, y_values_pcL_p, y_values_pc3_p):
                    electron = DIS_NC_XSEC_pc(
                        pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 0, x, q2, y, pc_func_type
                    )
                    positron = DIS_NC_XSEC_pc(
                        pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 1, x, q2, y, pc_func_type
                    )
                    result = electron(y_values_pc2_p, y_values_pcL_p, y_values_pc3_p) + positron(
                        y_values_pc2_p, y_values_pcL_p, y_values_pc3_p
                    )
                    return result / 2

                pc_func = average
            else:
                raise ValueError(f"{exp_name} not implemented.")

            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(
                    pars_pc['comb']['H2p'], pars_pc['comb']['HLp'], pars_pc['comb']['H3p']
                )

        # CHORUS
        elif exp_name.startswith('CHORUS_CC'):
            # Nu
            if exp_name == 'CHORUS_CC_NOTFIXED_PB_NU-SIGMARED':
                pc_func = DIS_CC_CHORUS_pc(
                    pc2_p_nodes,
                    pcL_p_nodes,
                    pc3_p_nodes,
                    pc2_d_nodes,
                    pcL_d_nodes,
                    pc3_d_nodes,
                    0,
                    x,
                    q2,
                    y,
                    pc_func_type,
                )
            # Nu bar
            elif exp_name == 'CHORUS_CC_NOTFIXED_PB_NB-SIGMARED':
                pc_func = DIS_CC_CHORUS_pc(
                    pc2_p_nodes,
                    pcL_p_nodes,
                    pc3_p_nodes,
                    pc2_d_nodes,
                    pcL_d_nodes,
                    pc3_d_nodes,
                    1,
                    x,
                    q2,
                    y,
                    pc_func_type,
                )
            else:
                raise ValueError(f"{exp_name} not implemented.")
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(
                    pars_pc['comb']['H2p'],
                    pars_pc['comb']['HLp'],
                    pars_pc['comb']['H3p'],
                    pars_pc['comb']['H2d'],
                    pars_pc['comb']['HLd'],
                    pars_pc['comb']['H3d'],
                )

        # NuTeV
        elif exp_name.startswith('NUTEV_CC'):
            # Nu
            if exp_name == 'NUTEV_CC_NOTFIXED_FE_NU-SIGMARED':
                pc_func = DIS_CC_NUTEV_pc(
                    pc2_p_nodes,
                    pcL_p_nodes,
                    pc3_p_nodes,
                    pc2_d_nodes,
                    pcL_d_nodes,
                    pc3_d_nodes,
                    0,
                    x,
                    q2,
                    y,
                    pc_func_type,
                )
            # Nu bar
            elif exp_name == 'NUTEV_CC_NOTFIXED_FE_NB-SIGMARED':
                pc_func = DIS_CC_NUTEV_pc(
                    pc2_p_nodes,
                    pcL_p_nodes,
                    pc3_p_nodes,
                    pc2_d_nodes,
                    pcL_d_nodes,
                    pc3_d_nodes,
                    1,
                    x,
                    q2,
                    y,
                    pc_func_type,
                )
            else:
                raise ValueError(f"{exp_name} not implemented.")
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(
                    pars_pc['comb']['H2p'],
                    pars_pc['comb']['HLp'],
                    pars_pc['comb']['H3p'],
                    pars_pc['comb']['H2d'],
                    pars_pc['comb']['HLd'],
                    pars_pc['comb']['H3d'],
                )

        # HERA_CC
        elif exp_name.startswith('HERA_CC'):
            # electron
            if exp_name == 'HERA_CC_318GEV_EM-SIGMARED':
                pc_func = DIS_CC_HERA_XSEC_pc(
                    pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 0, x, q2, y, pc_func_type
                )
            # positron
            elif exp_name == 'HERA_CC_318GEV_EP-SIGMARED':
                pc_func = DIS_CC_HERA_XSEC_pc(
                    pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 1, x, q2, y, pc_func_type
                )
            else:
                raise ValueError(f"{exp_name} not implemented.")

            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(
                    pars_pc['comb']['H2p'], pars_pc['comb']['HLp'], pars_pc['comb']['H3p']
                )

        else:
            raise ValueError(
                f"The {process_type} observable for {exp_name} " "has not been implemented."
            )

    elif process_type == 'JET':
        pc_jet_nodes = pc_dict["Hj"]['nodes']
        eta = dataset_sp.commondata.metadata.load_kinematics()['y'].to_numpy().reshape(-1)[cuts]
        pT = dataset_sp.commondata.metadata.load_kinematics()['pT'].to_numpy().reshape(-1)[cuts]

        pc_func = JET_pc(dataset_sp, pdf, pc_jet_nodes, pT, eta, pc_func_type)
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
            pc_func = JET_pc(dataset_sp, pdf, pc_jet_nodes, m_jj, eta_star, pc_func_type)
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
            pc_func = JET_pc(dataset_sp, pdf, pc_jet_nodes, m_jj, eta_diff, pc_func_type)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['H2j_CMS'])

        else:
            raise ValueError(
                f"{dataset_sp.commondata.metadata.experiment} is not implemented for DIJET."
            )

    else:
        raise RuntimeError(f"{process_type} has not been implemented.")

    return deltas
