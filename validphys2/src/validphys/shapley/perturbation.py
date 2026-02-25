"""Gaussian perturbation of PDF grid values.

Provides tools for perturbing selected flavour channels with a Gaussian
bump. Two independent parameters control the perturbation:
  - mode: 'additive' or 'multiplicative'
  - xspace: 'linear' or 'logx'
"""

import numpy as np

# Valid perturbation parameter choices
PERTURBATION_MODES = ('additive', 'multiplicative')
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
        'additive': f -> f + A * G(x)
        'multiplicative': f -> f * (1 + A * G(x))
    xspace : str
        'linear' or 'logx'

    Returns
    -------
    gv_pert : np.ndarray
        Perturbed copy of gv.
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
    gauss = gaussian_profile(xgrid, mu, sigma, amplitude, xspace)

    if mode == 'additive':
        for fi in local_flavor_idx:
            gv_pert[:, fi, :] += gauss[np.newaxis, :]
    else:  # multiplicative
        for fi in local_flavor_idx:
            gv_pert[:, fi, :] *= (1.0 + gauss[np.newaxis, :])
    return gv_pert
