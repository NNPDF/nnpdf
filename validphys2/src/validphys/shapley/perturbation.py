"""Gaussian perturbation of PDF grid values.

Provides tools for perturbing selected flavour channels with a Gaussian bump. Two independent parameters control the perturbation:
  - mode: 'additive', 'multiplicative', 'calibrated', 'calibrated_centre',
    'calibrated_symmetric', or 'ablation'
  - xspace: 'linear' or 'logx'

Calibrated mode (default)
-------------------------
In 'calibrated' mode the Gaussian is scaled by the LOCAL 68% band of the
calibration set at every x:

    delta_j^(+)(x) = alpha * G(x) * (q84_j(x) - c_j(x))
    delta_j^(-)(x) = alpha * G(x) * (c_j(x) - q16_j(x))

where c_j is the central member and q16/q84 the 68% band (validphys Stats:
replica percentiles for Monte-Carlo sets, rescaled quadrature for Hessian
sets). Since G <= 1 the bump never leaves alpha times the local band, and it
equals that band at x = mu. In units of the local band the bump is exactly
G(x) centred at mu; in absolute terms its weight sits below mu where the band
falls steeply (high x). The band is taken on the grid nodes of
``calibration_gv`` (or ``gv``). Up and down shifts are built independently.

Calibrated-centre mode (legacy)
-------------------------------
'calibrated_centre' is the definition used until October 2026: the band is
taken once, at the centre, and copied over the whole Gaussian,

    amplitude_j_plus  = alpha * (q84_j(x_mu) - c_j(x_mu))
    amplitude_j_minus = alpha * (c_j(x_mu) - q16_j(x_mu))

Where the PDF falls steeply (high x) the tails of this bump exceed the local
band by orders of magnitude, so it overstates large-x sensitivity. Kept to
reproduce earlier results.

Calibrated-symmetric mode (legacy)
----------------------------------
'calibrated_symmetric' restores the pre-March-2026 definition: one amplitude
per flavour, at the centre, used for both signs,

    amplitude_j_plus = amplitude_j_minus = alpha * std_j(x_mu)

where std_j is the 1-sigma uncertainty of the calibration set at x ~ mu
(replica standard deviation for Monte-Carlo sets, rescaled quadrature for
Hessian sets, via validphys' Stats.std_error).

Ablation mode
-------------
In 'ablation' mode the selected flavour channels are set to zero over the
entire x grid:

    xf_j^pert(x) = 0  for all x

"""

import numpy as np

# Valid perturbation parameter choices
PERTURBATION_MODES = (
    'additive', 'multiplicative', 'calibrated', 'calibrated_centre',
    'calibrated_symmetric', 'ablation'
)
# Modes whose amplitude is the band at the centre, copied over the Gaussian.
CENTRE_CALIBRATED_MODES = ('calibrated_centre', 'calibrated_symmetric')
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


def _calibrated_amplitudes(calib_vals, amplitude, calibration_stats=None,
                           symmetric=False):
    """Per-flavour up/down calibrated amplitudes from a member ensemble.

    ``calib_vals`` is a 1-D array of the calibration members evaluated at a
    single (flavour, x) point. Returns ``(amp_plus, amp_minus)``.

    With ``symmetric=True`` both amplitudes equal ``amplitude * std``, the
    1-sigma uncertainty (``stats.std_error()``, or the replica standard
    deviation without ``calibration_stats``).

    When ``calibration_stats`` (a ``pdf.stats_class`` callable) is given, the
    central value and 68% band come from validphys' Stats classes, which are
    correct for Monte-Carlo (percentiles) and Hessian (rescaled quadrature)
    sets alike; ``calib_vals`` must then start with the central member (row 0).
    Otherwise mean + 16/84 percentile envelope is used.
    """
    calib_vals = np.asarray(calib_vals, dtype=float)
    if symmetric:
        if calibration_stats is not None:
            std = calibration_stats(calib_vals[:, np.newaxis]).std_error()
            std = float(np.asarray(std).ravel()[0])
        else:
            std = float(np.std(calib_vals))
        amp = float(amplitude) * max(std, 0.0)
        return amp, amp
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


def _local_band(gv_sigma, fi, calibration_stats=None, amplitude_cache=None):
    """Up/down 1-sigma band of flavour ``fi`` at every grid node: ``(b_plus, b_minus)``.

    ``gv_sigma`` has shape (ncalib, nfl, nx); with ``calibration_stats`` the
    central member is row 0 and the band comes from the Stats class, otherwise
    mean and 16/84 percentiles. Cached in ``amplitude_cache`` (tied to one grid).
    """
    key = ('local_band', int(fi))
    if amplitude_cache is not None and key in amplitude_cache:
        return amplitude_cache[key]
    vals = np.asarray(gv_sigma[:, fi, :], dtype=float)
    if calibration_stats is not None:
        stats = calibration_stats(vals)
        central = np.asarray(stats.central_value(), dtype=float).ravel()
        lo, hi = stats.errorbar68()
        lo = np.asarray(lo, dtype=float).ravel()
        hi = np.asarray(hi, dtype=float).ravel()
    else:
        central = vals.mean(axis=0)
        lo = np.percentile(vals, 16.0, axis=0)
        hi = np.percentile(vals, 84.0, axis=0)
    band = (np.clip(hi - central, 0.0, None), np.clip(central - lo, 0.0, None))
    if amplitude_cache is not None:
        amplitude_cache[key] = band
    return band


def _cached_calibrated_amplitudes(gv_sigma, fi, idx_mu, amplitude,
                                  calibration_stats, amplitude_cache,
                                  symmetric=False, calibration_at_mu=None, mu=None):
    """``_calibrated_amplitudes`` for ``gv_sigma[:, fi, idx_mu]``, optionally memoised.

    ``amplitude_cache`` is a dict owned by the caller and tied to ONE
    calibration array (one FK grid). The amplitudes depend only on that
    array's column, the amplitude multiplier and the stats class, so caching
    returns bit-identical values without recomputing percentiles.

    With ``calibration_at_mu`` (members evaluated at the bump centre ``mu``)
    the amplitudes come from ``calibration_at_mu[:, fi]`` and are keyed by
    ``mu`` instead of the grid node, so they no longer depend on the FK grid.
    """
    if calibration_at_mu is None:
        key = (int(fi), int(idx_mu), float(amplitude), bool(symmetric))
    else:
        key = ('mu', int(fi), float(mu), float(amplitude), bool(symmetric))
    if amplitude_cache is not None:
        cached = amplitude_cache.get(key)
        if cached is not None:
            return cached
    calib_vals = np.asarray(
        gv_sigma[:, fi, idx_mu] if calibration_at_mu is None
        else np.asarray(calibration_at_mu)[:, fi], dtype=float,
    )
    amps = _calibrated_amplitudes(
        calib_vals, amplitude, calibration_stats=calibration_stats,
        symmetric=symmetric,
    )
    if amplitude_cache is not None:
        amplitude_cache[key] = amps
    return amps


def apply_gaussian_perturbation(gv, local_flavor_idx, mu, sigma, amplitude,
                                xgrid, mode='additive', xspace='linear',
                                random_sign=False, rng=None,
                                flavor_signs=None, calibration_gv=None,
                                calibration_stats=None, amplitude_cache=None,
                                calibration_at_mu=None):
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
        'calibrated':    f_j -> f_j + sign * A * G(x) * band_j^(sign)(x)
            with band_j^(+)(x) = q84_j(x) - c_j(x) and
            band_j^(-)(x) = c_j(x) - q16_j(x), the local 68% band at every
            grid node. The bump never leaves A times the local band.
        'calibrated_centre': f_j -> f_j + delta_j^(sign) * G(x)
            where delta_j^(+)  = A * (q84_j - c_j),
                  delta_j^(-)  = A * (c_j - q16_j),
            at mu when calibration_at_mu is supplied, otherwise at the
            closest grid point (legacy; tails can exceed the local band).
        'calibrated_symmetric': as 'calibrated_centre' with A * std_j(mu)
            for both signs.
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
    amplitude_cache : dict or None
        Optional memo for calibrated amplitudes, owned by the caller and tied
        to this ``calibration_gv`` (or ``gv``). None recomputes every call.
    calibration_at_mu : np.ndarray or None, shape (ncalib, nfl)
        Centre-calibrated modes only. Calibration members evaluated at the
        actual bump centre. Overrides
        calibration_gv for the amplitude calculation, so the same bump can
        be evaluated on different FK and sum-rule grids. The flavour axis
        must match gv; member conventions are those of calibration_stats.

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

    if mode in CENTRE_CALIBRATED_MODES:
        # Per-flavor amplitudes from calibrated +1sigma / -1sigma envelopes
        # at the supplied centre (or legacy nearest grid point).
        xgrid_arr = np.asarray(xgrid)
        idx_mu = int(np.argmin(np.abs(xgrid_arr - mu)))
        gv_sigma = gv if calibration_gv is None else np.asarray(calibration_gv, dtype=float)
        for col, fi in enumerate(local_flavor_idx):
            signs = sign_matrix[:, col][:, np.newaxis]
            amp_plus, amp_minus = _cached_calibrated_amplitudes(
                gv_sigma, fi, idx_mu, amplitude, calibration_stats, amplitude_cache,
                symmetric=(mode == 'calibrated_symmetric'),
                calibration_at_mu=calibration_at_mu, mu=mu,
            )

            gauss_plus = gaussian_profile(xgrid, mu, sigma, amp_plus, xspace)
            gauss_minus = gaussian_profile(xgrid, mu, sigma, amp_minus, xspace)

            # Build replica-wise signed perturbations using fixed up/down
            # templates for this flavor.
            delta = np.where(signs >= 0.0, gauss_plus[np.newaxis, :], -gauss_minus[np.newaxis, :])
            gv_pert[:, fi, :] += delta
    elif mode == 'calibrated':
        # Gaussian shape times the local up/down band at every grid node.
        gv_sigma = gv if calibration_gv is None else np.asarray(calibration_gv, dtype=float)
        gauss = gaussian_profile(xgrid, mu, sigma, amplitude, xspace)
        for col, fi in enumerate(local_flavor_idx):
            signs = sign_matrix[:, col][:, np.newaxis]
            band_plus, band_minus = _local_band(gv_sigma, fi, calibration_stats, amplitude_cache)
            delta = np.where(signs >= 0.0, (gauss * band_plus)[np.newaxis, :],
                             -(gauss * band_minus)[np.newaxis, :])
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
                                      calibration_stats=None, amplitude_cache=None,
                                      calibration_at_mu=None):
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
        As in apply_gaussian_perturbation. In 'calibrated' mode bumps of one
        flavour superpose as band(x) * sum_k sign_k G_k(x), with the up band
        for positive and the down band for negative signs.
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
    amplitude_cache : dict or None
        Optional calibrated-amplitude memo; see apply_gaussian_perturbation.
    calibration_at_mu : mapping or None
        Centre-calibrated modes only. Maps each bump centre to a (ncalib, nfl) member array evaluated at
        that centre, with the same flavour ordering as gv. Overrides the
        nearest-grid-point calibration for every spec when supplied.

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
                band_plus, band_minus = _local_band(
                    gv_sigma, fi, calibration_stats, amplitude_cache
                )
                gauss = gaussian_profile(xgrid_arr, mu_k, sigma, amplitude, xspace)
                d = np.where(
                    signs >= 0.0,
                    (gauss * band_plus)[np.newaxis, :],
                    -(gauss * band_minus)[np.newaxis, :],
                )
            elif mode in CENTRE_CALIBRATED_MODES:
                idx_mu = int(np.argmin(np.abs(xgrid_arr - mu_k)))
                amp_plus, amp_minus = _cached_calibrated_amplitudes(
                    gv_sigma, fi, idx_mu, amplitude, calibration_stats, amplitude_cache,
                    symmetric=(mode == 'calibrated_symmetric'),
                    calibration_at_mu=(
                        None if calibration_at_mu is None else calibration_at_mu[mu_k]
                    ),
                    mu=mu_k,
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

        if mode in ('additive', 'calibrated') + CENTRE_CALIBRATED_MODES:
            gv_pert[:, fi, :] += delta
        else:  # multiplicative: apply compound factor (1 + sum_k G_k)
            gv_pert[:, fi, :] *= (1.0 + delta)

    # Enforce positivity on all modified channels.
    for fi in fi_to_specs:
        became_negative = (gv[:, fi, :] >= 0.0) & (gv_pert[:, fi, :] < 0.0)
        gv_pert[:, fi, :] = np.where(became_negative, 0.0, gv_pert[:, fi, :])

    return gv_pert
