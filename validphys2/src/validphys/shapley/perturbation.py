"""Gaussian perturbation of PDF grid values.

Provides tools for perturbing selected flavour channels with a Gaussian bump. Two independent parameters control the perturbation:
  - mode: 'additive', 'multiplicative', 'calibrated', or 'ablation'
  - xspace: 'linear' or 'logx'

Calibrated mode
---------------
In 'calibrated' mode the Gaussian amplitude is not a global constant but is
scaled per-flavor from calibrated +1sigma/-1sigma envelopes evaluated at
x ~ mu. The up/down shifts are built independently and can be asymmetric:

    amplitude_j_plus  = alpha * (q84_j(x_mu) - c_j(x_mu))
    amplitude_j_minus = alpha * (c_j(x_mu) - q16_j(x_mu))

where c_j is the central calibrated reference (replica mean), and q84/q16 are
the 84th/16th percentiles over replicas. The Gaussian shape is preserved; only
its peak height varies by flavor and by sign. This avoids assuming a symmetric
down-shift equal to the negative of the up-shift.

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


def _calibrated_amplitudes(calib_vals, amplitude, calibration_stats=None):
    """Per-flavour up/down calibrated amplitudes from a member ensemble.

    ``calib_vals`` is a 1-D array of the calibration members evaluated at a
    single (flavour, x) point. Returns ``(amp_plus, amp_minus)``.

    When ``calibration_stats`` (a ``pdf.stats_class`` callable) is given, the
    central value and 68% band come from validphys' Stats classes, which are
    correct for Monte-Carlo (percentiles) and Hessian (rescaled quadrature)
    sets alike; ``calib_vals`` must then start with the central member (row 0).
    Otherwise mean + 16/84 percentile envelope is used.
    """
    calib_vals = np.asarray(calib_vals, dtype=float)
    if calibration_stats is not None:
        stats = calibration_stats(calib_vals[:, np.newaxis])
        central = float(np.asarray(stats.central_value()).ravel()[0])
        lo, hi = stats.errorbar68()
        lo = float(np.asarray(lo).ravel()[0])
        hi = float(np.asarray(hi).ravel()[0])
    else:
        central = float(np.mean(calib_vals))
        lo = float(np.percentile(calib_vals, 16.0))
        hi = float(np.percentile(calib_vals, 84.0))
    amp_plus = float(amplitude) * max(hi - central, 0.0)
    amp_minus = float(amplitude) * max(central - lo, 0.0)
    return amp_plus, amp_minus


def apply_gaussian_perturbation(gv, local_flavor_idx, mu, sigma, amplitude,
                                xgrid, mode='additive', xspace='linear',
                                random_sign=False, rng=None,
                                flavor_signs=None, calibration_gv=None,
                                calibration_stats=None):
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
        'calibrated':    f_j -> f_j + delta_j^(sign) * G(x)
            where delta_j^(+)  = A * (q84_j - c_j),
                  delta_j^(-)  = A * (c_j - q16_j),
            at x closest to mu, with c_j the calibrated replica mean.
            This builds distinct calibrated up/down templates and does not
            assume delta_j^- = -delta_j^+.
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
        per-flavour spread. The full member array (central at row 0 followed by
        every error member) is expected when ``calibration_stats`` is provided.
        Useful when the evaluated members differ from the calibration ensemble,
        e.g. central-only runs.
    calibration_stats : callable or None
        Optional ``pdf.stats_class`` callable used to turn the calibration
        member array into a central value and 68% band that is correct for both
        Monte-Carlo and Hessian sets. When provided (calibrated mode), the
        up/down amplitudes come from ``stats.central_value()`` and
        ``stats.errorbar68()`` instead of raw mean/percentiles. When None, the
        legacy mean + 16/84 percentile envelope is used (Monte-Carlo only).

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
        # Per-flavor amplitudes from calibrated +1sigma / -1sigma envelopes
        # at x closest to mu, allowing asymmetric up/down shifts.
        xgrid_arr = np.asarray(xgrid)
        idx_mu = int(np.argmin(np.abs(xgrid_arr - mu)))
        gv_sigma = gv if calibration_gv is None else np.asarray(calibration_gv, dtype=float)
        for col, fi in enumerate(local_flavor_idx):
            signs = sign_matrix[:, col][:, np.newaxis]
            calib_vals = np.asarray(gv_sigma[:, fi, idx_mu], dtype=float)
            amp_plus, amp_minus = _calibrated_amplitudes(
                calib_vals, amplitude, calibration_stats=calibration_stats
            )

            gauss_plus = gaussian_profile(xgrid, mu, sigma, amp_plus, xspace)
            gauss_minus = gaussian_profile(xgrid, mu, sigma, amp_minus, xspace)

            # Build replica-wise signed perturbations using fixed up/down
            # templates for this flavor.
            delta = np.where(signs >= 0.0, gauss_plus[np.newaxis, :], -gauss_minus[np.newaxis, :])
            gv_pert[:, fi, :] += delta
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


def apply_multi_gaussian_perturbation(gv, specs, sigma, amplitude,
                                      xgrid, mode='calibrated', xspace='logx',
                                      flavor_signs=None, calibration_gv=None,
                                      calibration_stats=None):
    """Perturb flavours with a superposition of per-player Gaussian bumps.

    Each entry in *specs* describes one (flavor, x) player in a coalition.
    Multiple specs sharing the same local_fi accumulate their bumps additively
    on that flavour channel, then positivity is enforced once at the end.

    Parameters
    ----------
    gv : np.ndarray, shape (nrep, nfl, nx)
        Grid values to perturb.
    specs : list of (local_fi, mu_k, sign_col)
        One tuple per player in the coalition:
        - local_fi  : index into the flavour axis of *gv*
        - mu_k      : Gaussian centre (x value for this player)
        - sign_col  : column index in *flavor_signs* carrying this player's signs
    sigma : float
        Gaussian width (shared across all players).
    amplitude : float
        Peak amplitude multiplier (shared across all players).
    xgrid : array-like
        x values corresponding to the last axis of *gv*.
    mode : str
        'additive', 'multiplicative', 'calibrated', or 'ablation'.
    xspace : str
        'linear' or 'logx'.
    flavor_signs : np.ndarray or None
        Shape ``(nrep, n_players)`` -- column *sign_col* carries the per-replica
        sign for the corresponding player.  None -> all signs +1.
    calibration_gv : np.ndarray or None
        Member ensemble used only in calibrated mode to compute the up/down
        envelope. Full member array (central at row 0) when calibration_stats
        is given.
    calibration_stats : callable or None
        Optional ``pdf.stats_class`` callable; see apply_gaussian_perturbation.

    Returns
    -------
    gv_pert : np.ndarray
        Perturbed copy of *gv* with positivity enforced on modified channels.
    """
    if not specs:
        return gv

    if mode not in PERTURBATION_MODES:
        raise ValueError(
            f"Unknown perturbation mode '{mode}'. Choose from {PERTURBATION_MODES}."
        )
    if xspace not in PERTURBATION_XSPACES:
        raise ValueError(
            f"Unknown perturbation xspace '{xspace}'. Choose from {PERTURBATION_XSPACES}."
        )

    gv_pert = gv.copy()
    nrep = gv.shape[0]
    xgrid_arr = np.asarray(xgrid)

    if mode == 'ablation':
        unique_fi = set(fi for fi, _, _ in specs)
        for fi in unique_fi:
            gv_pert[:, fi, :] = 0.0
        return gv_pert

    # Group specs by flavour index so we can accumulate bumps per flavour.
    from collections import defaultdict
    fi_to_specs = defaultdict(list)
    for fi, mu_k, sign_col in specs:
        fi_to_specs[fi].append((mu_k, sign_col))

    gv_sigma = gv if calibration_gv is None else np.asarray(calibration_gv, dtype=float)

    for fi, mu_sign_list in fi_to_specs.items():
        # Accumulate delta across all (mu_k, sign_col) entries for this flavour.
        delta = np.zeros((nrep, xgrid_arr.shape[0]), dtype=float)

        for mu_k, sign_col in mu_sign_list:
            signs = (
                flavor_signs[:, sign_col:sign_col + 1].astype(float)
                if flavor_signs is not None
                else np.ones((nrep, 1), dtype=float)
            )

            if mode == 'calibrated':
                idx_mu = int(np.argmin(np.abs(xgrid_arr - mu_k)))
                calib_vals = gv_sigma[:, fi, idx_mu].astype(float)
                amp_plus, amp_minus = _calibrated_amplitudes(
                    calib_vals, amplitude, calibration_stats=calibration_stats
                )
                gauss_plus = gaussian_profile(xgrid_arr, mu_k, sigma, amp_plus, xspace)
                gauss_minus = gaussian_profile(xgrid_arr, mu_k, sigma, amp_minus, xspace)
                d = np.where(
                    signs >= 0.0,
                    gauss_plus[np.newaxis, :],
                    -gauss_minus[np.newaxis, :],
                )
            elif mode == 'additive':
                gauss = gaussian_profile(xgrid_arr, mu_k, sigma, amplitude, xspace)
                d = signs * gauss[np.newaxis, :]
            else:  # multiplicative -- accumulate additive deltas; (1+G1)(1+G2) ~ 1+G1+G2
                gauss = gaussian_profile(xgrid_arr, mu_k, sigma, amplitude, xspace)
                d = signs * gauss[np.newaxis, :]

            delta += d

        if mode in ('additive', 'calibrated'):
            gv_pert[:, fi, :] += delta
        else:  # multiplicative: apply compound factor (1 + sum_k G_k)
            gv_pert[:, fi, :] *= (1.0 + delta)

    # Enforce positivity on all modified channels.
    for fi in fi_to_specs:
        became_negative = (gv[:, fi, :] >= 0.0) & (gv_pert[:, fi, :] < 0.0)
        gv_pert[:, fi, :] = np.where(became_negative, 0.0, gv_pert[:, fi, :])

    return gv_pert
