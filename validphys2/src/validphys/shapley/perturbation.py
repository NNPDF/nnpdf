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
                                xgrid, mode='additive', xspace='linear',
                                random_sign=False, rng=None,
                                flavor_signs=None, calibration_gv=None):
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
    random_sign : bool
        When True the Gaussian amplitude is independently flipped to +1 or -1
        for each replica/flavour pair (drawn uniformly from {-1, +1}).
        Ignored for mode='ablation'.
    rng : numpy.random.Generator or None
        Optional random-number generator for reproducibility.  When None a
        fresh, non-seeded generator is created per call.
    flavor_signs : np.ndarray or None
        Optional explicit sign matrix with shape ``(nrep, n_selected_flavours)``
        matching ``local_flavor_idx`` order. When provided, these fixed signs
        override the internal random-sign draw.
    calibration_gv : np.ndarray or None
        Optional grid-value ensemble used only to determine the calibrated
        per-flavour replica spread. Useful when the evaluated members differ
        from the replica ensemble, e.g. central-only runs.

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
    nrep = gv.shape[0]

    if mode == 'ablation':
        # Zero out the selected flavour channels entirely over all x.
        # mu/sigma/amplitude are not used.
        for fi in local_flavor_idx:
            gv_pert[:, fi, :] = 0.0
        return gv_pert  # positivity trivially satisfied; skip clipping

    # Sign matrix: one sign per replica/flavour pair for the selected flavours.
    if flavor_signs is not None:
        sign_matrix = np.asarray(flavor_signs, dtype=float)
        if sign_matrix.ndim == 1:
            sign_matrix = sign_matrix[:, np.newaxis]
        expected_shape = (nrep, len(local_flavor_idx))
        if sign_matrix.shape != expected_shape:
            raise ValueError(
                "flavor_signs must have shape "
                f"{expected_shape}, got {sign_matrix.shape}."
            )
    elif random_sign:
        _rng = rng if rng is not None else np.random.default_rng()
        sign_matrix = _rng.choice(
            np.array([-1.0, 1.0]), size=(nrep, len(local_flavor_idx))
        )
    else:
        sign_matrix = np.ones((nrep, len(local_flavor_idx)), dtype=float)

    if mode == 'calibrated':
        # Per-flavor amplitude: alpha * replica std at x closest to mu.
        xgrid_arr = np.asarray(xgrid)
        idx_mu = int(np.argmin(np.abs(xgrid_arr - mu)))
        gv_sigma = gv if calibration_gv is None else np.asarray(calibration_gv, dtype=float)
        for col, fi in enumerate(local_flavor_idx):
            sigma_rep = float(gv_sigma[:, fi, idx_mu].std())
            gauss = gaussian_profile(xgrid, mu, sigma, amplitude * sigma_rep, xspace)
            signs = sign_matrix[:, col][:, np.newaxis]
            gv_pert[:, fi, :] += signs * gauss[np.newaxis, :]
    elif mode == 'additive':
        gauss = gaussian_profile(xgrid, mu, sigma, amplitude, xspace)
        for col, fi in enumerate(local_flavor_idx):
            signs = sign_matrix[:, col][:, np.newaxis]
            gv_pert[:, fi, :] += signs * gauss[np.newaxis, :]
    else:  # multiplicative
        gauss = gaussian_profile(xgrid, mu, sigma, amplitude, xspace)
        for col, fi in enumerate(local_flavor_idx):
            signs = sign_matrix[:, col][:, np.newaxis]
            gv_pert[:, fi, :] *= (1.0 + signs * gauss[np.newaxis, :])

    # Enforce positivity only on perturbed channels and only when the perturbation would drive a previously non-negative point below zero.
    for fi in local_flavor_idx:
        became_negative = (gv[:, fi, :] >= 0.0) & (gv_pert[:, fi, :] < 0.0)
        gv_pert[:, fi, :] = np.where(became_negative, 0.0, gv_pert[:, fi, :])
    return gv_pert
