"""MSR + VSR sum rule enforcement for perturbed PDFs.

Applies the n3fit-style multiplicative normalization to perturbed PDF sets.
Sum rules are always applied in the evolution basis. Unlike an n3fit network,
the perturbed grids need not obey V = V24 = V35 before normalization, so each
valence channel uses its own integral.

FK_FLAVOURS ordering:
  0=photon, 1=Sigma, 2=g, 3=V, 4=V3, 5=V8, 6=V15, 7=V24, 8=V35, 9=T3, 10=T8, 11=T15, 12=T24, 13=T35

Momentum sum rule (MSR):
  int_0^1 dx x [g(x,Q) + Sigma(x,Q)] = 1
  Since the grid stores xf(x):  int xf(x) dx = momentum fraction.

Valence sum rules (VSR):
  int_0^1 dx V(x,Q) = int_0^1 dx V8(x,Q) = 3
  int_0^1 dx V3(x,Q) = 1
  int_0^1 dx V15(x,Q) = 3
  Since the grid stores xf(x):  int xf(x)/x dx = int f(x) dx = number.
"""

import numpy as np


def gen_integration_input(nx=2000):
    """Generate x-grid and trapezoidal weights for sum rule integrals.

    Same grid as n3fit.msr.gen_integration_input: nx/2 log-spaced points from 1e-9 to 0.1, then nx/2 linearly-spaced from 0.1 to 1.

    Parameters
    ----------
    nx : int
        Total number of grid points (default 2000).

    Returns
    -------
    xgrid : np.ndarray, shape (nx,)
    weights : np.ndarray, shape (nx,)
    """
    lognx = nx // 2
    linnx = nx - lognx
    xgrid_log = np.logspace(-9, -1, lognx + 1)
    xgrid_lin = np.linspace(0.1, 1, linnx)
    xgrid = np.concatenate([xgrid_log[:-1], xgrid_lin])
    # Trapezoidal weights
    spacing = np.zeros(nx + 1)
    for i in range(1, nx):
        spacing[i] = abs(xgrid[i - 1] - xgrid[i])
    weights = np.array([(spacing[i] + spacing[i + 1]) / 2.0
                        for i in range(nx)])
    return xgrid, weights


def compute_sumrule_normalization(gv_evol14, xgrid, weights):
    """Compute MSR + VSR normalization constants from evolution-basis PDF.

    Parameters
    ----------
    gv_evol14 : np.ndarray, shape (nrep, 14, nx)
        PDF xf(x) on the integration grid, in FK_FLAVOURS order.
    xgrid : np.ndarray, shape (nx,)
    weights : np.ndarray, shape (nx,)

    Returns
    -------
    norm : np.ndarray, shape (nrep, 14)
        Multiplicative normalization constant per replica per flavour.
        Non-constrained channels get 1.0.

    Raises
    ------
    ValueError
        If a required integral is non-finite or numerically zero. A finite
        multiplicative correction cannot restore a nonzero target then.

    Notes
    -----
    The same channels are normalized for every coalition, including the
    empty coalition. Coalition membership does not select compensators.
    Negative integrals require signed normalization factors.
    """
    nrep = gv_evol14.shape[0]
    norm = np.ones((nrep, 14))

    # Momentum integrals: int xf(x) dx  (grid stores xf(x))
    mom = np.einsum('rfx,x->rf', gv_evol14, weights)  # (nrep, 14)

    # Number integrals: int xf(x)/x dx = int f(x) dx  (for valence)
    inv_x = 1.0 / np.clip(xgrid, 1e-30, None)
    num = np.einsum('rfx,x,x->rf', gv_evol14, inv_x, weights)  # (nrep, 14)

    # MSR: gluon normalised so photon + Sigma + g momentum = 1
    sigma_mom = mom[:, 1]   # Sigma
    photon_mom = mom[:, 0]  # photon
    gluon_mom = mom[:, 2]   # g
    denominators = np.column_stack((gluon_mom, num[:, 3:9]))
    targets = np.empty_like(denominators)
    targets[:, 0] = 1.0 - sigma_mom - photon_mom
    targets[:, 1:] = [3.0, 1.0, 3.0, 3.0, 3.0, 3.0]
    if not np.all(np.isfinite(denominators)) or not np.all(np.isfinite(targets)):
        raise ValueError("Non-finite integral in sum-rule normalization.")
    zero = np.abs(denominators) < 1e-30
    if np.any(zero & (targets != 0.0)):
        raise ValueError("Cannot restore a nonzero sum-rule target from a zero integral.")
    # If both the integral and its target vanish, leave that channel alone.
    norm[:, 2:9] = np.divide(
        targets, denominators, out=np.ones_like(targets), where=~zero
    )

    return norm
