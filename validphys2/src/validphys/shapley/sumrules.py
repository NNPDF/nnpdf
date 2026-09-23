"""MSR + VSR sum rule enforcement for perturbed PDFs.

Reproduces the normalization logic of n3fit.layers.msr_normalization in pure numpy for applying to perturbated PDF sets. Sum rules are always applied in the evolution basis.

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
    norm[:, 2] = (1.0 - sigma_mom - photon_mom) / np.clip(
        np.abs(gluon_mom), 1e-30, None
    )

    # VSR: valence channels normalised to quark number targets
    v_num = num[:, 3]   # V integral
    norm[:, 3] = 3.0 / np.clip(np.abs(v_num), 1e-30, None)    # V
    norm[:, 4] = 1.0 / np.clip(np.abs(num[:, 4]), 1e-30, None) # V3
    norm[:, 5] = 3.0 / np.clip(np.abs(num[:, 5]), 1e-30, None) # V8
    norm[:, 6] = 3.0 / np.clip(np.abs(num[:, 6]), 1e-30, None) # V15
    norm[:, 7] = 3.0 / np.clip(np.abs(v_num), 1e-30, None)    # V24 (same as V)
    norm[:, 8] = 3.0 / np.clip(np.abs(v_num), 1e-30, None)    # V35 (same as V)

    return norm
