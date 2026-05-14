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
import pandas as pd
from validphys.pdfbases import ALL_FLAVOURS

# Valid perturbation parameter choices
PERTURBATION_MODES = ('calibrated', 'ablation', 'additive', 'multiplicative', 'physical')
PERTURBATION_XSPACES = ('linear', 'logx')

#path_to_ceilings = "/home/daksh/anaconda3/envs/environment_nnpdf/lib/python3.12/site-packages/nnpdf/validphys2/src/validphys/shapley"
flavor_info = {
    'indices': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 
    'pdg_codes': [-6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6], 
    'names': ['tbar', 'bbar', 'cbar', 'sbar', 'ubar', 'dbar', 'g', 'd', 'u', 's', 'c', 'b', 't'], 
    'n_flavors': 13
    }

flavor_9_info = {
    'indices': [2, 3, 4, 5, 6, 7, 8, 9, 10],
    'pdg_codes': [-4, -3, -2, -1, 21, 1, 2, 3, 4],
    'names': ['cbar', 'sbar', 'ubar', 'dbar', 'g', 'd', 'u', 's', 'c'],
    'n_flavors': 9,
}
# Quick initialization (done for largest possible basis to avoid any error)
ceilings_plus, ceilings_minus = float('inf')*np.ones(14), float('inf')*np.ones(14)
    
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
                                flavor_signs=None, calibration_gv=None, 
                                enforce_SR=False, nf=6, flavor_pdg_codes=None):
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
        per-flavour replica spread. Useful when the evaluated members differ
        from the replica ensemble, e.g. central-only runs.
    enforce_SR: bool
        True when mode='physical'
    nf: int
        No. of quark flavors

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
        
        # shift required to access correct parton 
        # from flavor_9_info dict, depending on nf
        shift = 6-nf

        # build once outside the loop
        all_flavours_list = list(ALL_FLAVOURS)
        pdg_to_name = dict(zip(flavor_info['pdg_codes'], flavor_info['names']))

        if flavor_pdg_codes is not None:
            # flavor_pdg_codes[col] is the PDG code of local_flavor_idx[col].
            fi_to_name = {
                fi: pdg_to_name.get(flavor_pdg_codes[col], f'unknown_pdg_{flavor_pdg_codes[col]}')
                for col, fi in enumerate(local_flavor_idx)
            }
        else:
            # Fallback: assume fi is an absolute index into ALL_FLAVOURS (14-entry list).
            fi_to_name = {
                fi: pdg_to_name.get(all_flavours_list[fi], f'unknown_pdg_{all_flavours_list[fi]}')
                for fi in local_flavor_idx
            }

        ceiling_by_name_plus = dict(zip(flavor_info['names'], ceilings_plus))  
        ceiling_by_name_minus = dict(zip(flavor_info['names'], ceilings_minus)) 

        with open("/tmp/dict.log", "a") as f:
            f.write(f"\nThe flavor dict was {flavor_pdg_codes} and \nfi_to_name is {fi_to_name}. Ceilings are:\n cp: {ceiling_by_name_plus}")

        for col, fi in enumerate(local_flavor_idx):
            signs = sign_matrix[:, col][:, np.newaxis]
            calib_vals = np.asarray(gv_sigma[:, fi, idx_mu], dtype=float)
            c_ref = float(np.mean(calib_vals))
            q16 = float(np.percentile(calib_vals, 16.0))
            q84 = float(np.percentile(calib_vals, 84.0))

            amp_plus = float(amplitude) * max(q84 - c_ref, 0.0)
            amp_minus = float(amplitude) * max(c_ref - q16, 0.0)
            
            if enforce_SR:
                old_amp_plus = amp_plus
                old_amp_minus = amp_minus
                parton = fi_to_name[fi]
                with open("/tmp/ceiling_debug.log", "a") as f:
                    f.write(f"\nparton={parton}")
                    f.write(f"\nceilings_plus length={len(ceilings_plus)}")
                    f.write(f"\nceiling_by_name_plus keys={list(ceiling_by_name_plus.keys())}")
                    f.write(f"\nflavor_info names={flavor_info['names']}")
                amp_plus = min(amp_plus, ceiling_by_name_plus[parton])
                amp_minus = min(amp_minus, ceiling_by_name_minus[parton])
                with open("/tmp/sr.log", "a") as f:
                    f.write(f"\nparton is {parton}, and ceiling is {ceiling_by_name_plus[parton]}")
                if old_amp_plus != amp_plus:
                    with open("/tmp/ceilings_applied.log", "a") as f:
                        f.write(f"\n[DEBUG] fi={fi}, shift={shift}, nf={nf}, fi-shift={fi-shift}, names[fi-shift]={flavor_info['names'][fi-shift]}, parton={parton}")
                        f.write(f"\n[INFO]: For mu: {mu} and sigma: {sigma}\nOld {parton} amp_plus: {old_amp_plus} --- New {parton} amp_plus: {amp_plus}\n")
                if old_amp_minus != amp_minus:
                    with open("/tmp/ceilings_applied.log", "a") as f:
                        f.write(f"\n[DEBUG] fi={fi}, shift={shift}, nf={nf}, fi-shift={fi-shift}, names[fi-shift]={flavor_info['names'][fi-shift]}, parton={parton}")
                        f.write(f"\n[INFO]: For mu: {mu} and sigma: {sigma}\nOld {parton} amp_plus: {old_amp_plus} --- New {parton} amp_plus: {amp_plus}\n")
                
                
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
