"""Gaussian perturbation of PDF grid values.

Provides tools for perturbing selected flavour channels with a Gaussian bump. Two independent parameters control the perturbation:
  - mode: 'additive', 'multiplicative', 'calibrated', or 'ablation'
  - xspace: 'linear' or 'logx'

Calibrated mode
---------------
In 'calibrated' mode the Gaussian amplitude is not a global constant but is
scaled per-flavor by the replica standard deviation evaluated at x ~ mu:

    amplitude_j = alpha * std_rep[ f_j(x_mu) ]

so the perturbation probes a shift of *alpha* times what the fit itself is
uncertain about at that x value.  The Gaussian shape is preserved; only its
peak height varies by flavor.  This removes the bias of both additive (gluon
is much larger than sea quarks) and multiplicative (gluon is much more
tightly constrained in relative terms) approaches.

Ablation mode
-------------
In 'ablation' mode the selected flavour channels are set to zero over the
entire x grid:

    xf_j^pert(x) = 0  for all x

"""

import numpy as np

# Valid perturbation parameter choices
PERTURBATION_MODES = ('additive', 'multiplicative', 'calibrated', 'ablation')
PERTURBATION_XSPACES = ('linear', 'logx')


def gaussian_profile(xgrid, mu, sigma, amplitude, xspace='linear'):
    """Build the perturbation profile G(x) for a given x-space.

    Returns an array of shape (nx,).

    Parameters
    ----------
    xgrid : array-like
        x values.
    mu : float
        Centre of the Gaussian.
    sigma : float
        Width. For xspace='logx', sigma is in decades of log10(x).
    amplitude : float
        Peak height of the Gaussian.
    xspace : str
        'linear': G(x) = A * exp(-0.5 * ((x - mu) / sigma)^2)
        'logx':   G(x) = A * exp(-0.5 * ((log10(x) - log10(mu)) / sigma)^2)
    """
    x = np.asarray(xgrid)
    if xspace == 'logx':
        log_x = np.log10(np.clip(x, 1e-30, None))
        log_mu = np.log10(max(mu, 1e-30))
        return amplitude * np.exp(-0.5 * ((log_x - log_mu) / sigma) ** 2)
    else:
        return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def apply_gaussian_perturbation(gv, local_flavor_idx, mu, sigma, amplitude,
                                xgrid, mode='additive', xspace='linear'):
    """Perturb selected flavour channels with a Gaussian bump.

    Parameters
    ----------
    gv : np.ndarray, shape (nrep, nfl, nx)
        Grid values to perturb.
    local_flavor_idx : list of int
        Indices into the flavour axis to perturb.
    mu, sigma, amplitude : float
        Gaussian parameters.
    xgrid : array-like
        x values corresponding to the last axis of gv.
    mode : str
        'additive':      f_j -> f_j + A * G(x)
        'multiplicative': f_j -> f_j * (1 + A * G(x))
        'calibrated':    f_j -> f_j + (A * sigma_rep_j) * G(x)
            where sigma_rep_j = std of replicas at x closest to mu.
            A acts as a dimensionless scaling factor (alpha).
        'ablation':      f_j -> 0 for all x.
            mu, sigma, amplitude are ignored.
    xspace : str
        'linear' or 'logx' (ignored for mode='ablation')

    Returns
    -------
    gv_pert : np.ndarray
        Perturbed copy of gv.

    Notes
    -----
    Positivity is always enforced on perturbed flavour channels:
    values that would become negative after perturbation are clipped to 0.
    """
    if len(local_flavor_idx) == 0:
        return gv
    if mode not in PERTURBATION_MODES:
        raise ValueError(
            f"Unknown perturbation mode '{mode}'. "
            f"Choose from {PERTURBATION_MODES}."
        )
    if xspace not in PERTURBATION_XSPACES:
        raise ValueError(
            f"Unknown perturbation xspace '{xspace}'. "
            f"Choose from {PERTURBATION_XSPACES}."
        )
    gv_pert = gv.copy()

    if mode == 'ablation':
        # Zero out the selected flavour channels entirely over all x.
        # mu/sigma/amplitude are not used.
        for fi in local_flavor_idx:
            gv_pert[:, fi, :] = 0.0
        return gv_pert  # positivity trivially satisfied; skip clipping
    elif mode == 'calibrated':
        # Per-flavor amplitude: alpha * replica std at x closest to mu.
        xgrid_arr = np.asarray(xgrid)
        idx_mu = int(np.argmin(np.abs(xgrid_arr - mu)))
        for fi in local_flavor_idx:
            sigma_rep = float(gv[:, fi, idx_mu].std())
            gauss = gaussian_profile(xgrid, mu, sigma, amplitude * sigma_rep, xspace)
            gv_pert[:, fi, :] += gauss[np.newaxis, :]
    elif mode == 'additive':
        gauss = gaussian_profile(xgrid, mu, sigma, amplitude, xspace)
        for fi in local_flavor_idx:
            gv_pert[:, fi, :] += gauss[np.newaxis, :]
    else:  # multiplicative
        gauss = gaussian_profile(xgrid, mu, sigma, amplitude, xspace)
        for fi in local_flavor_idx:
            gv_pert[:, fi, :] *= (1.0 + gauss[np.newaxis, :])

    # Enforce positivity only on perturbed channels and only when the perturbation would drive a previously non-negative point below zero.
    for fi in local_flavor_idx:
        became_negative = (gv[:, fi, :] >= 0.0) & (gv_pert[:, fi, :] < 0.0)
        gv_pert[:, fi, :] = np.where(became_negative, 0.0, gv_pert[:, fi, :])
    return gv_pert
