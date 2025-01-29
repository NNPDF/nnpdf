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
from validphys.process_options import _Process

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
    res = []
    for shift_pos, shift in enumerate(y_shift):
        bin_low = bin_edges[shift_pos]
        bin_high = bin_edges[shift_pos + 1]
        condition = np.multiply(a >= bin_low, a < bin_high)
        res.append([shift for cond in condition if cond])
    res = np.concatenate(res)
    return res


def dis_pc_func(
    delta_h: npt.ArrayLike, nodes: npt.ArrayLike, x: npt.ArrayLike, Q2: npt.ArrayLike
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
    PC = step_function(x, delta_h, nodes) / Q2
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
    PC = step_function(rap, delta_h, nodes) / pT
    return PC


def JET_pc(pc_nodes, pT, rap):
    """
    Returns the function that computes the shift for the ratio for single
    jet cross sections. In particular, the shift is computed such that

      xsec -> xsec + PC,

    and the shift is defined as

      Delta(xsec) = (xsec + xsec) - xsec = PC.

    The power correction is a function of the transverse momentum of the jet.
    """

    def func(y_values):
        result = jets_pc_func(y_values, pc_nodes, pT, rap)
        return result

    return func


# TODO Maybe we want to treat the function that parametrizes the PC
# as argument?
def DIS_F2_pc(pc2_nodes, x, q2):
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
        result = dis_pc_func(y_values, pc2_nodes, x, q2)
        return result

    return PC_2


def DIS_F2R_pc(experiment, pdf, pc_2_p_nodes, pc_2_d_nodes, x, q2):
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
        PC_d = dis_pc_func(y_values_d, pc_2_d_nodes, x, q2)
        PC_p = dis_pc_func(y_values_p, pc_2_p_nodes, x, q2)
        num = np.sum([F2D, PC_d], axis=0)
        denom = np.sum([F2P, PC_p], axis=0)
        result = np.array(operator.truediv(num, denom) - F2_ratio)
        return result

    return func


def DIS_F2C_pc(pc2_p_nodes, pc2_d_nodes, x, q2):
    """
    Builds the function used to compute the shifts for the charm
    structure function measured by EMC. The process involved is

        mu^+ + Fe -> mu+^ + c cbar + X .

    This function works exactly as the previous functions used to
    compute nuisance shifts. In this case, the constructed function
    (`func` below) requires two lists of parameters for the proton
    and the deuteron contribution. The reason being that in this process
    the muon scatters off an iron target, and the power correction
    contribution is a mixture of proton and deuteron nucleons. Hence, proton
    and deuteron contribution are weighted by the appropriate atomic factor.

    Note that we are parametrising power corrections as proton and deuteron
    targets. If we were to parametrize such contributions using, say, proton
    and nucleon, than the weights would change.


    Nuclear target
    --------------
    The power corrections for nuclear observables, like in this case, are affected
    by the pc contribution of the protons and that of the neutrons.
    If we allow for the non-iscoscalarity of the target, and combining the two
    contributions in accordance with the atomic and mass number (A and Z), the
    power correction for the nuclear target can be written as (see  eq.(4.2.5)
    in https://nnpdf.mi.infn.it/wp-content/uploads/2021/09/thesis_master_RP.pdf)

      PC_N = 1/A (Z * PC_p + (A-Z) * PC_n) .

    The deuteron is obtained using the isoscalarity, namely

      PC_c = 1/2 (PC_p + PC_n) .

    Since we parametrise the power corrections of the proton and the deuteron,
    we can combined the above equations and write

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
        PC2_d = dis_pc_func(y_values_d, pc2_d_nodes, x, q2)
        PC2_p = dis_pc_func(y_values_p, pc2_p_nodes, x, q2)
        result = (2 * Z - A) / A * PC2_p + 2 * (A - Z) / A * PC2_d
        return result

    return func


def DIS_NC_XSEC_pc(pc2_nodes, pcL_nodes, pc3_nodes, lepton, x, q2, y):
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
        PC_2 = dis_pc_func(y_values_pc2, pc2_nodes, x, q2)
        PC_L = dis_pc_func(y_values_pcL, pcL_nodes, x, q2)
        PC_3 = dis_pc_func(y_values_pc3, pc3_nodes, x, q2)
        result = PC_2 + N_L * PC_L + N_3 * PC_3
        return result

    return func


def DIS_CC_HERA_XSEC_pc(pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, lepton, x, q2, y):
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
        PC2_p = dis_pc_func(y_values_pc2_p, pc2_p_nodes, x, q2)
        PCL_p = dis_pc_func(y_values_pcL_p, pcL_p_nodes, x, q2)
        PC3_p = dis_pc_func(y_values_pc3_p, pc3_p_nodes, x, q2)

        # Build the contribution to the x-sec of the power corrections
        result = N * (PC2_p + N_L * PCL_p + N_3 * PC3_p)
        return result

    return func


def DIS_CC_NUTEV_pc(
    pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, pc2_d_nodes, pcL_d_nodes, pc3_d_nodes, lepton, x, q2, y
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
        PC2_p = dis_pc_func(y_values_pc2_p, pc2_p_nodes, x, q2)
        PCL_p = dis_pc_func(y_values_pcL_p, pcL_p_nodes, x, q2)
        PC3_p = dis_pc_func(y_values_pc3_p, pc3_p_nodes, x, q2)
        PC2_d = dis_pc_func(y_values_pc2_d, pc2_d_nodes, x, q2)
        PCL_d = dis_pc_func(y_values_pcL_d, pcL_d_nodes, x, q2)
        PC3_d = dis_pc_func(y_values_pc3_d, pc3_d_nodes, x, q2)
        tmp_2 = (2 * Z - A) / A * PC2_p + 2 * (A - Z) / A * PC2_d
        tmp_L = (2 * Z - A) / A * PCL_p + 2 * (A - Z) / A * PCL_d
        tmp_3 = (2 * Z - A) / A * PC3_p + 2 * (A - Z) / A * PC3_d
        result = N * (tmp_2 + N_L * tmp_L + N_3 * tmp_3)
        return result

    return func


# TODO This is function is really similar to the one
# defined for NUTEV CC. Can we reduce code repetitions?
def DIS_CC_CHORUS_pc(
    pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, pc2_d_nodes, pcL_d_nodes, pc3_d_nodes, lepton, x, q2, y
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
        PC2_p = dis_pc_func(y_values_pc2_p, pc2_p_nodes, x, q2)
        PCL_p = dis_pc_func(y_values_pcL_p, pcL_p_nodes, x, q2)
        PC3_p = dis_pc_func(y_values_pc3_p, pc3_p_nodes, x, q2)
        PC2_d = dis_pc_func(y_values_pc2_d, pc2_d_nodes, x, q2)
        PCL_d = dis_pc_func(y_values_pcL_d, pcL_d_nodes, x, q2)
        PC3_d = dis_pc_func(y_values_pc3_d, pc3_d_nodes, x, q2)
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


def compute_deltas_pc(dataset_sp: DataSetSpec, pdf: PDF, power_corr_dict: dict):
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
    cd_table = dataset_sp.load_commondata().commondata_table
    process_type = cd_table['process'].iloc[0]
    if isinstance(process_type, _Process):
        process_type = process_type.name

    pars_combs = construct_pars_combs(power_corr_dict)
    deltas = defaultdict(list)

    pc_func = None
    if process_type.startswith('DIS'):
        pc2_p_nodes = power_corr_dict["H2p"]['nodes']
        pcL_p_nodes = power_corr_dict["HLp"]['nodes']
        pc3_p_nodes = power_corr_dict["H3p"]['nodes']
        pc2_d_nodes = power_corr_dict["H2d"]['nodes']
        pcL_d_nodes = power_corr_dict["HLd"]['nodes']
        pc3_d_nodes = power_corr_dict["H3d"]['nodes']

        # TODO
        # AFter the data re-implementation the name of the variables
        # in the commondata table will change as indicated in the metadata.
        # When this happens, this part must be updated.
        x = cd_table['kin1'].to_numpy()
        q2 = cd_table['kin2'].to_numpy()
        y = cd_table['kin3'].to_numpy()

        # F2 ratio
        if exp_name == "NMC_NC_NOTFIXED_EM-F2":
            pc_func = DIS_F2R_pc(dataset_sp, pdf, pc2_p_nodes, pc2_d_nodes, x, q2)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['H2p'], pars_pc['comb']['H2d'])

        # F2 proton traget
        elif exp_name in F2P_exps:
            pc_func = DIS_F2_pc(pc2_p_nodes, x, q2)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['H2p'])

        # F2 deuteron traget
        elif exp_name in F2D_exps:
            pc_func = DIS_F2_pc(pc2_d_nodes, x, q2)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['H2d'])

        # EMC
        elif exp_name.startswith('EMC_NC_250GEV'):
            pc_func = DIS_F2C_pc(pc2_p_nodes, pc2_d_nodes, x, q2)
            for pars_pc in pars_combs:
                deltas[pars_pc['label']] = pc_func(pars_pc['comb']['H2p'], pars_pc['comb']['H2d'])

        # HERA and NMC SIGMARED NC
        elif exp_name in np.concatenate([NC_SIGMARED_P_EM, NC_SIGMARED_P_EP, NC_SIGMARED_P_EAVG]):
            # Electron
            if exp_name in NC_SIGMARED_P_EM:
                pc_func = DIS_NC_XSEC_pc(pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 0, x, q2, y)
            # Positron
            elif exp_name in NC_SIGMARED_P_EP:
                pc_func = DIS_NC_XSEC_pc(pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 1, x, q2, y)
            # Average positron and electron
            # TODO
            # Check if this is correct (ach)
            elif NC_SIGMARED_P_EAVG:

                def average(y_values_pc2_p, y_values_pcL_p, y_values_pc3_p):
                    electron = DIS_NC_XSEC_pc(pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 0, x, q2, y)
                    positron = DIS_NC_XSEC_pc(pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 1, x, q2, y)
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
                pc_func = DIS_CC_HERA_XSEC_pc(pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 0, x, q2, y)
            # positron
            elif exp_name == 'HERA_CC_318GEV_EP-SIGMARED':
                pc_func = DIS_CC_HERA_XSEC_pc(pc2_p_nodes, pcL_p_nodes, pc3_p_nodes, 1, x, q2, y)
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
        pc_jet_nodes = power_corr_dict["Hj"]['nodes']

        # TODO
        # AFter the data re-implementation the name of the variables
        # in the commondata table will change as indicated in the metadata.
        # When this happens, this part must be updated.
        eta = cd_table['kin1'].to_numpy()
        pT = cd_table['kin2'].to_numpy()
        q2 = pT * pT

        pc_func = JET_pc(pc_jet_nodes, pT, eta)
        for pars_pc in pars_combs:
            deltas[pars_pc['label']] = pc_func(pars_pc['comb']['Hj'])

    elif process_type == 'DIJET':
        raise RuntimeError(f"No implementation for {exp_name} yet.")

    else:
        raise RuntimeError(f"{process_type} has not been implemented.")

    return deltas
