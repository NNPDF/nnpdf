"""NNPDF Shapley value analyzer.

Wraps NNPDF chi2 evaluation into a value function for the generic
ExactShapley solver from the ``shapley_values`` package. Handles:
  - Evolution and physical-flavor perturbation bases
  - Multi-FK datasets with operations (ASY, RATIO, ADD)
  - MSR + VSR sum rule enforcement
  - Grid-value caching for performance
"""

import functools
import math
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations

import numpy as np
import matplotlib.pyplot as plt

from validphys.convolution import FK_FLAVOURS
from validphys.pdfbases import evolution as evol_basis, ALL_FLAVOURS

from shapley_values import ExactShapley, plot_shapley_bar

from .setup import (
    get_pdf_grid_values,
    get_pdf_flavor_grid_values,
    get_pdf_grid_values_all14,
    FLAVOR_PDG_NAMES,
)
from .perturbation import (
    apply_gaussian_perturbation,
    PERTURBATION_MODES,
    PERTURBATION_XSPACES,
)
from .sumrules import (
    gen_integration_input,
    compute_sumrule_normalization,
)


class NNPDFShapleyAnalyzer:
    """Compute Shapley values for PDF flavour importance.

    Supports two perturbation bases:
    - evolution : players are evolution-basis flavours (Sigma, g, V, T3, ...).
    - flavor : players are physical quarks/gluon rotated before convolution, Warning: if sumrule_enforced=true, the perturbation will be modify in the evolution basis for normalized sumrules.

    Parameters
    ----------
    pdf : validphys.core.PDF
    observables : list of NNPDFObservable
    flavor_info : dict
        Flavor metadata from setup_observables.
    n_replicas : int or None
    basis : str
        'evolution' or 'flavor'.
    enforce_sumrules : bool
        If True, re-impose momentum and valence sum rules on the
        perturbed PDF before convolution.
    """

    EVOL_LABELS = {
        "photon": "gamma (photon)",
        "\\Sigma": "Sigma (singlet)",
        "g": "g (gluon)",
        "V": "V (valence)",
        "V3": "V3", "V8": "V8", "V15": "V15", "V24": "V24", "V35": "V35",
        "T3": "T3", "T8": "T8", "T15": "T15", "T24": "T24", "T35": "T35",
    }

    EVOL_SHORT = {
        "photon": "gamma",
        "\\Sigma": "Sigma",
        "g": "g",
        "V": "V",
        "V3": "V3", "V8": "V8", "V15": "V15", "V24": "V24", "V35": "V35",
        "T3": "T3", "T8": "T8", "T15": "T15", "T24": "T24", "T35": "T35",
    }

    def __init__(self, pdf, observables, flavor_info, n_replicas=None,
                 basis='evolution', enforce_sumrules=False,
                 member_mode='replicas'):
        self.pdf = pdf
        self.observables = observables
        self.flavor_info = flavor_info
        self.n_replicas = n_replicas
        self.basis = basis
        self.enforce_sumrules = enforce_sumrules
        self.member_mode = str(member_mode).strip().lower()
        if self.member_mode not in {"replicas", "central"}:
            raise ValueError(
                f"Unknown member_mode '{member_mode}'. Use 'replicas' or 'central'."
            )

        self._gv_cache: dict = {}
        self._gv_calib_cache: dict = {}
        self._gv_cache_lock = threading.Lock()
        self._sumrule_lock = threading.Lock()

        # Sum rule integration grid (lazily initialized)
        self._sr_xgrid = None
        self._sr_weights = None
        self._sr_gv_evol = None
        self._sr_gv_flav = None
        self._sr_calib_gv_evol = None
        self._sr_calib_gv_flav = None
        self._sr_rotation = None

        if basis == 'evolution':
            self.flavor_names = flavor_info["names"]
            self.flavor_indices = flavor_info["indices"]
            self.n_flavors = flavor_info["n_flavors"]
            self.flavor_labels = [
                self.EVOL_LABELS.get(n, n) for n in self.flavor_names
            ]
            self.flavor_short = [
                self.EVOL_SHORT.get(n, n) for n in self.flavor_names
            ]
        elif basis == 'flavor':
            self.flavor_names = flavor_info["names"]
            self.flavor_indices = flavor_info["indices"]
            self.n_flavors = flavor_info["n_flavors"]
            self.flavor_labels = self.flavor_names
            self.flavor_short = self.flavor_names
        else:
            raise ValueError(
                f"Unknown basis '{basis}'. Use 'evolution' or 'flavor'."
            )

    # -- Sum rule setup and computation ------------------------------------

    def _setup_sumrule_grid(self):
        """Lazily initialise the integration grid for sum rules."""
        sr_xgrid, sr_weights = gen_integration_input(2000)
        Q0 = self.observables[0].Q0

        sr_gv_evol = get_pdf_grid_values_all14(
            self.pdf,
            type("SumruleTarget", (), {"Q0": Q0, "xgrid": sr_xgrid})(),
            n_replicas=self.n_replicas,
            member_mode=self.member_mode,
        )
        if self.member_mode == "central":
            sr_gv_calib_evol = get_pdf_grid_values_all14(
                self.pdf,
                type("SumruleTarget", (), {"Q0": Q0, "xgrid": sr_xgrid})(),
                n_replicas=self.n_replicas,
                member_mode="replicas",
            )
        else:
            sr_gv_calib_evol = sr_gv_evol

        sr_gv_flav = None
        sr_gv_calib_flav = None
        sr_rotation = None

        if self.basis == 'flavor':
            sr_gv_flav = get_pdf_flavor_grid_values(
                self.pdf,
                type("SumruleTarget", (), {"Q0": Q0, "xgrid": sr_xgrid})(),
                n_replicas=self.n_replicas,
                member_mode=self.member_mode,
            )
            if self.member_mode == "central":
                sr_gv_calib_flav = get_pdf_flavor_grid_values(
                    self.pdf,
                    type("SumruleTarget", (), {"Q0": Q0, "xgrid": sr_xgrid})(),
                    n_replicas=self.n_replicas,
                    member_mode="replicas",
                )
            else:
                sr_gv_calib_flav = sr_gv_flav

            evol_row_inds = evol_basis._to_indexes(FK_FLAVOURS)
            sr_rotation = evol_basis.from_flavour_mat[evol_row_inds, :]

        self._sr_xgrid = sr_xgrid
        self._sr_weights = sr_weights
        self._sr_gv_evol = sr_gv_evol
        self._sr_gv_flav = sr_gv_flav
        self._sr_calib_gv_evol = sr_gv_calib_evol
        self._sr_calib_gv_flav = sr_gv_calib_flav
        self._sr_rotation = sr_rotation

    def _compute_sumrule_norm(self, flavor_subset, mu, sigma, amplitude,
                              mode, xspace, random_sign_matrix=None):
        """Compute normalization constants for a perturbed coalition."""
        ready = (
            self._sr_xgrid is not None
            and self._sr_weights is not None
            and self._sr_gv_evol is not None
            and (
                self.basis != 'flavor'
                or (
                    self._sr_gv_flav is not None
                    and self._sr_rotation is not None
                )
            )
        )
        if not ready:
            with self._sumrule_lock:
                ready = (
                    self._sr_xgrid is not None
                    and self._sr_weights is not None
                    and self._sr_gv_evol is not None
                    and (
                        self.basis != 'flavor'
                        or (
                            self._sr_gv_flav is not None
                            and self._sr_rotation is not None
                        )
                    )
                )
                if not ready:
                    self._setup_sumrule_grid()

        if self.basis == 'evolution':
            gv = self._sr_gv_evol.copy()
            perturb_idx = [self.flavor_indices[p] for p in flavor_subset]
            perturb_signs = (
                random_sign_matrix[:, flavor_subset]
                if random_sign_matrix is not None else None
            )
            gv_pert = apply_gaussian_perturbation(
                gv, perturb_idx, mu, sigma, amplitude,
                self._sr_xgrid, mode=mode, xspace=xspace,
                flavor_signs=perturb_signs,
                calibration_gv=self._sr_calib_gv_evol,
            )
        else:
            gv_flav = self._sr_gv_flav.copy()
            perturb_idx = [self.flavor_indices[p] for p in flavor_subset]
            perturb_signs = (
                random_sign_matrix[:, flavor_subset]
                if random_sign_matrix is not None else None
            )
            gv_flav_pert = apply_gaussian_perturbation(
                gv_flav, perturb_idx, mu, sigma, amplitude,
                self._sr_xgrid, mode=mode, xspace=xspace,
                flavor_signs=perturb_signs,
                calibration_gv=self._sr_calib_gv_flav,
            )
            gv_pert = np.einsum(
                'ef,rfx->rex', self._sr_rotation, gv_flav_pert
            )

        return compute_sumrule_normalization(
            gv_pert, self._sr_xgrid, self._sr_weights
        )

    @staticmethod
    def _apply_norm_to_gv(gv, norm, flavor_indices):
        """Multiply grid values by per-flavour normalization constants."""
        local_norm = norm[:, flavor_indices]
        return gv * local_norm[:, :, np.newaxis]

    # -- Grid-value caching ------------------------------------------------

    def clear_cache(self):
        """Drop cached grid values."""
        self._gv_cache.clear()
        self._gv_calib_cache.clear()

    def _get_gv_for_entry(self, obs, entry_idx):
        """Cached evolution-basis grid values (FK-subset flavours)."""
        key = (obs.name, entry_idx, 'evol')
        if key not in self._gv_cache:
            with self._gv_cache_lock:
                if key not in self._gv_cache:
                    entry = obs.fk_entries[entry_idx]
                    self._gv_cache[key] = get_pdf_grid_values(
                        self.pdf, entry, n_replicas=self.n_replicas,
                        member_mode=self.member_mode,
                    )
        return self._gv_cache[key]

    def _get_calibration_gv_for_entry(self, obs, entry_idx):
        """Cached evolution-basis replica ensemble used for calibration."""
        key = (obs.name, entry_idx, 'evol_calib')
        if key not in self._gv_calib_cache:
            with self._gv_cache_lock:
                if key not in self._gv_calib_cache:
                    entry = obs.fk_entries[entry_idx]
                    self._gv_calib_cache[key] = get_pdf_grid_values(
                        self.pdf, entry, n_replicas=self.n_replicas,
                        member_mode='replicas',
                    )
        return self._gv_calib_cache[key]

    def _get_flavor_gv_for_entry(self, obs, entry_idx):
        """Cached physical flavor-basis grid values."""
        key = (obs.name, entry_idx, 'flavor')
        if key not in self._gv_cache:
            with self._gv_cache_lock:
                if key not in self._gv_cache:
                    entry = obs.fk_entries[entry_idx]
                    self._gv_cache[key] = get_pdf_flavor_grid_values(
                        self.pdf, entry, n_replicas=self.n_replicas,
                        member_mode=self.member_mode,
                    )
        return self._gv_cache[key]

    def _get_calibration_flavor_gv_for_entry(self, obs, entry_idx):
        """Cached flavor-basis replica ensemble used for calibration."""
        key = (obs.name, entry_idx, 'flavor_calib')
        if key not in self._gv_calib_cache:
            with self._gv_cache_lock:
                if key not in self._gv_calib_cache:
                    entry = obs.fk_entries[entry_idx]
                    self._gv_calib_cache[key] = get_pdf_flavor_grid_values(
                        self.pdf, entry, n_replicas=self.n_replicas,
                        member_mode='replicas',
                    )
        return self._gv_calib_cache[key]

    def _get_gv_all14_for_entry(self, obs, entry_idx):
        """Cached full 14-flavour evolution-basis grid values."""
        key = (obs.name, entry_idx, 'evol_all14')
        if key not in self._gv_cache:
            with self._gv_cache_lock:
                if key not in self._gv_cache:
                    entry = obs.fk_entries[entry_idx]
                    self._gv_cache[key] = get_pdf_grid_values_all14(
                        self.pdf, entry, n_replicas=self.n_replicas,
                        member_mode=self.member_mode,
                    )
        return self._gv_cache[key]

    def _get_calibration_gv_all14_for_entry(self, obs, entry_idx):
        """Cached full 14-flavour replica ensemble used for calibration."""
        key = (obs.name, entry_idx, 'evol_all14_calib')
        if key not in self._gv_calib_cache:
            with self._gv_cache_lock:
                if key not in self._gv_calib_cache:
                    entry = obs.fk_entries[entry_idx]
                    self._gv_calib_cache[key] = get_pdf_grid_values_all14(
                        self.pdf, entry, n_replicas=self.n_replicas,
                        member_mode='replicas',
                    )
        return self._gv_calib_cache[key]

    def _get_gv_list(self, obs):
        """Get list of baseline grid values for all FK entries."""
        result = []
        for idx, entry in enumerate(obs.fk_entries):
            if entry.hadronic:
                result.append(self._get_gv_all14_for_entry(obs, idx))
            else:
                result.append(self._get_gv_for_entry(obs, idx))
        return result

    def _local_flavor_indices_for_entry(self, entry, global_flavor_subset):
        """Map global Shapley-player indices to local FK column indices."""
        local = []
        fi_list = entry.flavor_indices.tolist()
        for player in global_flavor_subset:
            global_fi = self.flavor_indices[player]
            if global_fi in fi_list:
                local.append(fi_list.index(global_fi))
        return local

    def _local_flavor_indices_and_players_for_entry(self, entry, global_flavor_subset):
        """Map global player indices to local FK columns, preserving order."""
        local = []
        players = []
        fi_list = entry.flavor_indices.tolist()
        for player in global_flavor_subset:
            global_fi = self.flavor_indices[player]
            if global_fi in fi_list:
                local.append(fi_list.index(global_fi))
                players.append(player)
        return local, players

    def _infer_n_members(self):
        """Infer the number of PDF members currently loaded by the analyzer."""
        obs0 = self.observables[0]
        if self.basis == 'flavor':
            return int(self._get_flavor_gv_for_entry(obs0, 0).shape[0])

        entry0 = obs0.fk_entries[0]
        if entry0.hadronic:
            return int(self._get_gv_all14_for_entry(obs0, 0).shape[0])
        return int(self._get_gv_for_entry(obs0, 0).shape[0])

    @staticmethod
    def _build_sign_matrices(
        n_members,
        n_flavors,
        n_sign_samples,
        random_seed=None,
        unique_rows=False,
    ):
        """Build one or more fixed sign tables for sign-mask sampling.

        Parameters
        ----------
        unique_rows : bool
            When True and n_sign_samples=1, draw per-member sign masks
            without replacement so every member gets a distinct mask.
        """
        if int(n_sign_samples) < 1:
            raise ValueError(f"n_sign_samples must be >= 1, got {n_sign_samples}")

        rng = np.random.default_rng(random_seed)
        n_sign_samples = int(n_sign_samples)
        if n_sign_samples == 1:
            if unique_rows:
                n_unique_max = 2 ** int(n_flavors)
                if int(n_members) > n_unique_max:
                    raise ValueError(
                        "Requested unique per-member sign masks but n_members "
                        f"({n_members}) exceeds available unique masks "
                        f"(2**n_flavors = {n_unique_max})."
                    )

                chosen = rng.choice(n_unique_max, size=int(n_members), replace=False)
                bit_pos = np.arange(int(n_flavors), dtype=np.int64)
                bits = ((chosen[:, None] >> bit_pos[None, :]) & 1).astype(float)
                signs = np.where(bits > 0.5, 1.0, -1.0)
                return [signs]

            return [
                rng.choice(np.array([-1.0, 1.0]), size=(n_members, n_flavors))
            ]

        if n_sign_samples % 2 != 0:
            raise ValueError(
                "n_sign_samples must be even when random_sign is enabled "
                "and more than one sign sample is requested."
            )

        n_pairs = n_sign_samples // 2
        base = [
            rng.choice(np.array([-1.0, 1.0]), size=(n_members, n_flavors))
            for _ in range(n_pairs)
        ]
        return base + [-m for m in base]

    # -- Chi2 evaluation with perturbation ----------------------------------

    def _evaluate_chi2(self, flavor_subset, mu, sigma, amplitude,
                       mode='additive', xspace='linear',
                       per_replica=False, random_sign=False,
                       random_sign_matrix=None):
        """Chi2/Ndata across all observables for a given perturbation.

        Parameters
        ----------
        flavor_subset : list of int
            Global Shapley-player indices to perturb.
        mu, sigma, amplitude : float
        mode : str
        xspace : str
        per_replica : bool
            When False (default) return the mean over replicas as a float.
            When True return the per-replica array of shape (nrep,) so that
            Shapley values can be computed independently for every replica.
        random_sign : bool
            When True use signed perturbations. When ``random_sign_matrix`` is
            provided, those fixed replica/flavour signs are reused for all
            coalition evaluations in the run.

        Returns
        -------
        chi2 : float or np.ndarray
            Float when per_replica=False; shape (nrep,) when per_replica=True.
        """
        sr_norm = None
        if self.enforce_sumrules:
            sr_norm = self._compute_sumrule_norm(
                flavor_subset, mu, sigma, amplitude, mode, xspace,
                random_sign_matrix=random_sign_matrix,
            )

        rng = (
            np.random.default_rng()
            if random_sign and random_sign_matrix is None else None
        )

        total_chi2 = 0.0          # used when per_replica=False
        total_chi2_rep = None     # used when per_replica=True  shape (nrep,)
        total_ndata = 0

        for obs in self.observables:
            if self.basis == 'flavor':
                gv_pert_list = []
                perturb_players = list(flavor_subset)
                perturb_signs = (
                    random_sign_matrix[:, perturb_players]
                    if random_sign_matrix is not None else None
                )
                for idx, entry in enumerate(obs.fk_entries):
                    gv_flav = self._get_flavor_gv_for_entry(obs, idx)
                    gv_flav_calib = self._get_calibration_flavor_gv_for_entry(obs, idx)
                    perturb_idx = [
                        self.flavor_indices[p] for p in flavor_subset
                    ]
                    gv_pert = apply_gaussian_perturbation(
                        gv_flav, perturb_idx, mu, sigma, amplitude,
                        entry.xgrid, mode=mode, xspace=xspace,
                        random_sign=random_sign, rng=rng,
                        flavor_signs=perturb_signs,
                        calibration_gv=gv_flav_calib,
                    )
                    gv_pert_list.append(gv_pert)

                if sr_norm is not None:
                    gv_evol_list = obs.rotate_to_evolution(gv_pert_list)
                    gv_evol_list = [
                        self._apply_norm_to_gv(
                            gv, sr_norm,
                            range(14) if entry.hadronic
                            else entry.flavor_indices
                        )
                        for gv, entry in zip(gv_evol_list, obs.fk_entries)
                    ]
                    chi2_arr = obs.chi2(gv_evol_list)
                else:
                    chi2_arr = obs.chi2_from_flavor(gv_pert_list)
            else:
                gv_pert_list = []
                for idx, entry in enumerate(obs.fk_entries):
                    if entry.hadronic:
                        gv = self._get_gv_all14_for_entry(obs, idx)
                        gv_calib = self._get_calibration_gv_all14_for_entry(obs, idx)
                        perturb_players = list(flavor_subset)
                        perturb_idx = [
                            self.flavor_indices[p] for p in flavor_subset
                        ]
                    else:
                        gv = self._get_gv_for_entry(obs, idx)
                        gv_calib = self._get_calibration_gv_for_entry(obs, idx)
                        perturb_idx, perturb_players = (
                            self._local_flavor_indices_and_players_for_entry(
                                entry, flavor_subset
                            )
                        )
                    perturb_signs = (
                        random_sign_matrix[:, perturb_players]
                        if random_sign_matrix is not None else None
                    )
                    gv_pert = apply_gaussian_perturbation(
                        gv, perturb_idx, mu, sigma, amplitude,
                        entry.xgrid, mode=mode, xspace=xspace,
                        random_sign=random_sign, rng=rng,
                        flavor_signs=perturb_signs,
                        calibration_gv=gv_calib,
                    )
                    if sr_norm is not None:
                        fi = (range(14) if entry.hadronic else entry.flavor_indices)
                        gv_pert = self._apply_norm_to_gv(
                            gv_pert, sr_norm, fi
                        )
                    gv_pert_list.append(gv_pert)
                chi2_arr = obs.chi2(gv_pert_list)

            if per_replica:
                if total_chi2_rep is None:
                    total_chi2_rep = chi2_arr.copy()
                else:
                    total_chi2_rep += chi2_arr
            else:
                total_chi2 += np.mean(chi2_arr)
            total_ndata += obs.ndata

        if per_replica:
            return total_chi2_rep / total_ndata
        return total_chi2 / total_ndata

    # -- Build value function for ExactShapley ------------------------------

    def build_value_function(self, mu, sigma, amplitude,
                             mode='additive', xspace='linear',
                             _coalition_log=None,
                             random_sign=False,
                             random_sign_matrix=None):
        """Return a callable v(coalition) -> float for ExactShapley.

        The wrapper tracks progress: coalition count, elapsed time, and estimated time remaining (ETA).

        Parameters
        ----------
        mu, sigma, amplitude : float
            Gaussian perturbation parameters.
        mode : str
        xspace : str
        _coalition_log : list or None
            If provided, every (coalition_tuple, chi2) pair is appended to
            this list during evaluation.  Use this to build the full
            per-coalition chi2 record for diagnostics.
        random_sign : bool
            Forwarded to _evaluate_chi2.
        random_sign_matrix : np.ndarray or None
            Fixed sign table with shape ``(n_members, n_flavors)`` reused for
            every coalition evaluation in the run.

        Returns
        -------
        value_function : callable
            f(List[int]) -> float
        """
        import time
        import sys

        total_coalitions = 2 ** self.n_flavors
        state = {"count": 0, "t_start": None, "last_player_line": False}
        interactive_stderr = sys.stderr.isatty()
        log_every = max(1, total_coalitions // 100)

        def v(coalition):
            if state["t_start"] is None:
                state["t_start"] = time.time()

            result = self._evaluate_chi2(
                coalition, mu, sigma, amplitude,
                mode=mode, xspace=xspace,
                random_sign=random_sign,
                random_sign_matrix=random_sign_matrix,
            )

            # Record coalition -> chi2 for diagnostics if requested.
            if _coalition_log is not None:
                _coalition_log.append((tuple(sorted(coalition)), result))

            state["count"] += 1
            n = state["count"]
            elapsed = time.time() - state["t_start"]

            if n == 1:
                msg = (
                    f"    [{n}/{total_coalitions}] "
                    f"elapsed {elapsed:.0f}s ..."
                )
            else:
                rate = elapsed / n
                remaining = rate * (total_coalitions - n)
                mins, secs = divmod(int(remaining), 60)
                msg = (
                    f"    [{n}/{total_coalitions}] "
                    f"elapsed {elapsed:.0f}s | "
                    f"~{rate:.1f}s/eval | "
                    f"ETA {mins}m{secs:02d}s   "
                )

            if interactive_stderr:
                sys.stderr.write(f"\r{msg}")
                sys.stderr.flush()
            elif (
                n == 1
                or n == total_coalitions
                or n % log_every == 0
            ):
                sys.stderr.write(f"{msg}\n")
                sys.stderr.flush()

            if n == total_coalitions:
                sys.stderr.write("\n")
                sys.stderr.flush()

            return result

        return v

    # -- Coalition diagnostics ----------------------------------------------

    def _compute_diagnostics(self, coalition_log, outlier_n_sigma=3.0):
        """Compute per-coalition chi2 statistics and marginal contribution stats.

        Parameters
        ----------
        coalition_log : list of (coalition_tuple, chi2)
            Raw log produced by build_value_function when _coalition_log is set.
        outlier_n_sigma : float
            Flag a coalition as an outlier when its chi2 exceeds
            mean + outlier_n_sigma * std.  Default is 3.0.

        Returns
        -------
        diag : dict
            chi2_stats, outliers, per_player_marginals, marginal_contributions.
        """
        if not coalition_log:
            return {}

        # Build O(1) lookup dict (coalitions are already sorted tuples).
        chi2_dict = {c: v for c, v in coalition_log}
        coalitions = list(chi2_dict.keys())
        chi2_arr = np.array([chi2_dict[c] for c in coalitions])

        mean_chi2 = float(np.mean(chi2_arr))
        std_chi2 = float(np.std(chi2_arr))
        threshold = mean_chi2 + outlier_n_sigma * std_chi2

        outlier_coalitions = [
            {
                "coalition": list(c),
                "coalition_labels": [self.flavor_short[i] for i in c],
                "size": len(c),
                "chi2": float(chi2_dict[c]),
                "z_score": float((chi2_dict[c] - mean_chi2) / std_chi2)
                           if std_chi2 > 0 else 0.0,
            }
            for c in coalitions
            if chi2_dict[c] > threshold
        ]
        outlier_coalitions.sort(key=lambda x: -x["chi2"])

        # Per-player marginal contribution analysis:
        # delta_v(i, S) = v(S U {i}) - v(S) for all S not containing i.
        per_player_marginals = {}
        all_marginals = []   # flat list used for the contributions CSV

        for player_idx in range(self.n_flavors):
            label = self.flavor_short[player_idx]
            deltas = []
            for coalition in coalitions:
                if player_idx not in coalition:
                    coal_with = tuple(sorted(coalition + (player_idx,)))
                    if coal_with in chi2_dict:
                        v_without = chi2_dict[coalition]
                        v_with = chi2_dict[coal_with]
                        delta = v_with - v_without
                        deltas.append(delta)
                        all_marginals.append({
                            "player_idx": player_idx,
                            "player": label,
                            "coalition_without": list(coalition),
                            "v_without": float(v_without),
                            "v_with": float(v_with),
                            "delta_v": float(delta),
                        })

            if deltas:
                deltas_arr = np.array(deltas)
                mean_d = float(np.mean(deltas_arr))
                std_d = float(np.std(deltas_arr))
                thr_d = mean_d + outlier_n_sigma * std_d
                thr_d_low = mean_d - outlier_n_sigma * std_d
                per_player_marginals[label] = {
                    "n_marginals": len(deltas),
                    "mean": mean_d,
                    "std": std_d,
                    "min": float(np.min(deltas_arr)),
                    "max": float(np.max(deltas_arr)),
                    "p05": float(np.percentile(deltas_arr, 5)),
                    "p95": float(np.percentile(deltas_arr, 95)),
                    "n_outliers_high": int(np.sum(deltas_arr > thr_d)),
                    "n_outliers_low": int(np.sum(deltas_arr < thr_d_low)),
                }

        # Flag outlier marginals in the flat list.
        for entry in all_marginals:
            lbl = entry["player"]
            stats = per_player_marginals.get(lbl, {})
            mean_d = stats.get("mean", 0.0)
            std_d = stats.get("std", 0.0)
            thr_hi = mean_d + outlier_n_sigma * std_d
            thr_lo = mean_d - outlier_n_sigma * std_d
            dv = entry["delta_v"]
            entry["is_outlier"] = int(dv > thr_hi or dv < thr_lo)

        diag = {
            "outlier_n_sigma": float(outlier_n_sigma),
            "chi2_stats": {
                "n_coalitions": len(chi2_arr),
                "mean": mean_chi2,
                "std": std_chi2,
                "min": float(np.min(chi2_arr)),
                "max": float(np.max(chi2_arr)),
                "median": float(np.median(chi2_arr)),
                "p95": float(np.percentile(chi2_arr, 95)),
                "p99": float(np.percentile(chi2_arr, 99)),
            },
            "outlier_chi2_threshold": float(threshold),
            "n_outlier_coalitions": len(outlier_coalitions),
            "outlier_coalitions": outlier_coalitions,
            "per_player_marginals": per_player_marginals,
            "_marginal_contributions": all_marginals,  # used internally for CSV
        }
        return diag

    # -- Convenience: run full Shapley analysis -----------------------------

    def _iter_all_coalitions(self):
        """Yield all coalition tuples in deterministic order."""
        players = tuple(range(self.n_flavors))
        for size in range(self.n_flavors + 1):
            for coalition in combinations(players, size):
                yield coalition

    @staticmethod
    def _coalition_with_player(coalition, player):
        """Return sorted coalition tuple with one extra player."""
        return tuple(sorted(coalition + (player,)))

    @staticmethod
    def _coalition_to_bitmask(coalition):
        """Encode a coalition tuple as an integer bitmask."""
        mask = 0
        for player in coalition:
            mask |= (1 << int(player))
        return int(mask)

    @staticmethod
    def _compute_shapley_from_cache(value_cache, all_coalitions, n_flavors):
        """Compute exact Shapley values from a coalition->value cache."""
        factorial = [math.factorial(k) for k in range(n_flavors + 1)]
        factorial_n = factorial[n_flavors]

        first_value = next(iter(value_cache.values()))
        is_vector = np.ndim(first_value) > 0
        if is_vector:
            n_members = int(np.asarray(first_value).shape[0])
            shapley_vals = np.zeros((n_members, n_flavors), dtype=float)
        else:
            shapley_vals = np.zeros(n_flavors, dtype=float)

        for player in range(n_flavors):
            for coalition in all_coalitions:
                if player in coalition:
                    continue
                size = len(coalition)
                weight = (
                    factorial[size] * factorial[n_flavors - size - 1]
                ) / factorial_n
                coalition_with = tuple(sorted(coalition + (player,)))
                delta = value_cache[coalition_with] - value_cache[coalition]
                if is_vector:
                    shapley_vals[:, player] += weight * delta
                else:
                    shapley_vals[player] += weight * delta

        return shapley_vals

    def _print_even_odd_validation_checks(self, checks, per_replica):
        """Print numerical consistency checks for the deterministic dual-game construction."""
        def _fmt(entry):
            return (
                f"max|diff|={entry['max_abs_diff']:.3e}, "
                f"tol={entry['tol']:.1e}, pass={entry['pass']}"
            )

        print("\n--- Deterministic calibrated up/down checks ---")
        print(f"  completeness even : {_fmt(checks['completeness_even'])}")
        print(f"  completeness odd  : {_fmt(checks['completeness_odd'])}")
        print(f"  v_even(empty)=0   : {_fmt(checks['empty_even_zero'])}")
        print(f"  v_odd(empty)=0    : {_fmt(checks['empty_odd_zero'])}")
        print(
            "  odd from up/down : "
            f"chi2_plus-minus residual={checks['chi2_symmetry_residual']:.3e}, "
            f"max|phi_odd|={checks['odd_zero_if_symmetric']['phi_odd_max_abs']:.3e}, "
            f"pass={checks['odd_zero_if_symmetric']['pass']}"
        )
        print(f"  phi_even relation : {_fmt(checks['even_vs_pm'])}")
        print(f"  phi_odd relation  : {_fmt(checks['odd_vs_pm'])}")
        if per_replica:
            print("  (checks evaluated over replica vectors using max-abs residuals)")
        print()

    def _compute_exact_shap_even_odd_parallel(
        self,
        mu,
        sigma,
        amplitude,
        mode,
        xspace,
        n_jobs,
        verbose=True,
        per_replica=False,
    ):
        """Compute deterministic calibrated up/down exact Shapley values.

        For every coalition S, evaluate once:
          chi2_plus(S)  = chi2(p0 + Delta_S_plus)
          chi2_minus(S) = chi2(p0 + Delta_S_minus)

        where Delta_S_plus and Delta_S_minus are built from the calibrated
        +1sigma and -1sigma flavour templates, respectively, and are not
        assumed to be opposite perturbations.

        Then define the two games:
          v_even(S) = 0.5 * (chi2_plus + chi2_minus) - chi2_baseline
          v_odd(S)  = 0.5 * (chi2_plus - chi2_minus)
        """
        n = self.n_flavors
        all_coalitions = list(self._iter_all_coalitions())
        total_coalitions = len(all_coalitions)
        non_empty = [c for c in all_coalitions if len(c) > 0]
        coalition_bitmasks = {
            coalition: self._coalition_to_bitmask(coalition)
            for coalition in all_coalitions
        }

        if verbose:
            print(
                f"Computing deterministic calibrated even/odd Shapley for {n} players "
                f"({total_coalitions} coalitions) with {n_jobs} worker(s)..."
            )

        n_members = self._infer_n_members()
        plus_sign_matrix = np.ones((n_members, n), dtype=float)
        minus_sign_matrix = -np.ones((n_members, n), dtype=float)

        t0 = time.time()
        progress_start = t0
        interactive_stderr = sys.stderr.isatty()
        log_every = max(1, total_coalitions // 100)

        # Empty coalition is shared: plus and minus are both the baseline.
        chi2_baseline = self._evaluate_chi2(
            [], mu, sigma, amplitude,
            mode=mode, xspace=xspace,
            per_replica=per_replica,
            random_sign=False,
            random_sign_matrix=plus_sign_matrix,
        )

        chi2_plus_cache = {tuple(): chi2_baseline}
        chi2_minus_cache = {tuple(): chi2_baseline}

        def _evaluate_one(coalition):
            chi2_plus = self._evaluate_chi2(
                coalition, mu, sigma, amplitude,
                mode=mode, xspace=xspace,
                per_replica=per_replica,
                random_sign=False,
                random_sign_matrix=plus_sign_matrix,
            )
            chi2_minus = self._evaluate_chi2(
                coalition, mu, sigma, amplitude,
                mode=mode, xspace=xspace,
                per_replica=per_replica,
                random_sign=False,
                random_sign_matrix=minus_sign_matrix,
            )
            return coalition, chi2_plus, chi2_minus

        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            future_map = {
                executor.submit(_evaluate_one, coalition): coalition
                for coalition in non_empty
            }
            completed = 1  # empty coalition already done
            for future in as_completed(future_map):
                coalition, chi2_plus, chi2_minus = future.result()
                chi2_plus_cache[coalition] = chi2_plus
                chi2_minus_cache[coalition] = chi2_minus
                completed += 1

                if verbose:
                    elapsed = time.time() - progress_start
                    if completed == 1:
                        msg = (
                            f"    [{completed}/{total_coalitions}] "
                            f"elapsed {elapsed:.0f}s ..."
                        )
                    else:
                        rate = elapsed / completed
                        remaining = rate * (total_coalitions - completed)
                        mins, secs = divmod(int(remaining), 60)
                        msg = (
                            f"    [{completed}/{total_coalitions}] "
                            f"elapsed {elapsed:.0f}s | "
                            f"~{rate:.1f}s/coalition | ETA {mins}m{secs:02d}s   "
                        )
                    if interactive_stderr:
                        sys.stderr.write(f"\r{msg}")
                        sys.stderr.flush()
                    elif (
                        completed == 1
                        or completed == total_coalitions
                        or completed % log_every == 0
                    ):
                        sys.stderr.write(f"{msg}\n")
                        sys.stderr.flush()

        if verbose:
            sys.stderr.write("\n")
            sys.stderr.flush()

        v_plus_cache = {}
        v_minus_cache = {}
        v_even_cache = {}
        v_odd_cache = {}
        for coalition in all_coalitions:
            chi2_plus = chi2_plus_cache[coalition]
            chi2_minus = chi2_minus_cache[coalition]
            v_plus_cache[coalition] = chi2_plus - chi2_baseline
            v_minus_cache[coalition] = chi2_minus - chi2_baseline
            v_even_cache[coalition] = 0.5 * (chi2_plus + chi2_minus) - chi2_baseline
            v_odd_cache[coalition] = 0.5 * (chi2_plus - chi2_minus)

        phi_even = self._compute_shapley_from_cache(v_even_cache, all_coalitions, n)
        phi_odd = self._compute_shapley_from_cache(v_odd_cache, all_coalitions, n)
        phi_plus = self._compute_shapley_from_cache(v_plus_cache, all_coalitions, n)
        phi_minus = self._compute_shapley_from_cache(v_minus_cache, all_coalitions, n)

        coalition_full = tuple(range(n))
        tol = 1e-10

        def _max_abs(x):
            return float(np.max(np.abs(np.asarray(x, dtype=float))))

        lhs_even = np.sum(phi_even, axis=0) if np.ndim(phi_even) == 1 else np.sum(phi_even, axis=1)
        rhs_even = v_even_cache[coalition_full] - v_even_cache[tuple()]
        lhs_odd = np.sum(phi_odd, axis=0) if np.ndim(phi_odd) == 1 else np.sum(phi_odd, axis=1)
        rhs_odd = v_odd_cache[coalition_full] - v_odd_cache[tuple()]

        chi2_symmetry_residual = max(
            _max_abs(chi2_plus_cache[c] - chi2_minus_cache[c])
            for c in all_coalitions
        )
        odd_zero_tol = 1e-10
        odd_max_abs = _max_abs(phi_odd)
        odd_zero_pass = (chi2_symmetry_residual > odd_zero_tol) or (odd_max_abs <= odd_zero_tol)

        even_pm_diff = phi_even - 0.5 * (phi_plus + phi_minus)
        odd_pm_diff = phi_odd - 0.5 * (phi_plus - phi_minus)

        checks = {
            "completeness_even": {
                "max_abs_diff": _max_abs(lhs_even - rhs_even),
                "tol": tol,
                "pass": bool(_max_abs(lhs_even - rhs_even) <= tol),
            },
            "completeness_odd": {
                "max_abs_diff": _max_abs(lhs_odd - rhs_odd),
                "tol": tol,
                "pass": bool(_max_abs(lhs_odd - rhs_odd) <= tol),
            },
            "empty_even_zero": {
                "max_abs_diff": _max_abs(v_even_cache[tuple()]),
                "tol": tol,
                "pass": bool(_max_abs(v_even_cache[tuple()]) <= tol),
            },
            "empty_odd_zero": {
                "max_abs_diff": _max_abs(v_odd_cache[tuple()]),
                "tol": tol,
                "pass": bool(_max_abs(v_odd_cache[tuple()]) <= tol),
            },
            "chi2_symmetry_residual": chi2_symmetry_residual,
            "odd_zero_if_symmetric": {
                "phi_odd_max_abs": odd_max_abs,
                "tol": odd_zero_tol,
                "pass": bool(odd_zero_pass),
            },
            "even_vs_pm": {
                "max_abs_diff": _max_abs(even_pm_diff),
                "tol": tol,
                "pass": bool(_max_abs(even_pm_diff) <= tol),
            },
            "odd_vs_pm": {
                "max_abs_diff": _max_abs(odd_pm_diff),
                "tol": tol,
                "pass": bool(_max_abs(odd_pm_diff) <= tol),
            },
        }

        elapsed = time.time() - t0

        if per_replica:
            nrep = int(np.asarray(phi_even).shape[0])
            mean_even = np.mean(phi_even, axis=0)
            std_even_rep = (
                np.std(phi_even, axis=0, ddof=1) if nrep > 1
                else np.zeros(n, dtype=float)
            )
            err_even_rep = std_even_rep / np.sqrt(max(nrep, 1))

            # Sign dependence should not cancel between replicas:
            # S_abs_j = < |phi_odd_j| >_replicas
            phi_odd_abs = np.abs(phi_odd)
            mean_odd = np.mean(phi_odd_abs, axis=0)

            if verbose:
                print(f"Baseline (mean) : {float(np.mean(chi2_baseline)):.6f}")
                print(f"Max |phi_even|  : {float(np.max(np.abs(mean_even))):.6f}")
                print(f"Max S_abs       : {float(np.max(mean_odd)):.6f}")
                print(f"Max even std    : {float(np.max(std_even_rep)):.6f}")
                print(f"Elapsed         : {elapsed:.1f}s")
                self._print_even_odd_validation_checks(checks, per_replica=True)

            return {
                "shapley_values": mean_even,
                "shapley_values_odd": mean_odd,
                "shapley_values_odd_signed": np.mean(phi_odd, axis=0),
                "shapley_values_per_replica": phi_even,
                "shapley_values_odd_per_replica": phi_odd_abs,
                "shapley_values_odd_signed_per_replica": phi_odd,
                "shapley_std": std_even_rep,
                "shapley_err": err_even_rep,
                "shapley_odd_std": None,
                "shapley_odd_err": None,
                "mean_even": mean_even,
                "std_even_rep": std_even_rep,
                "mean_odd": mean_odd,
                "std_odd_rep": None,
                "baseline": float(np.mean(chi2_baseline)),
                "baseline_per_replica": chi2_baseline,
                "player_labels": self.flavor_labels,
                "player_short": self.flavor_short,
                "coalitions_evaluated": total_coalitions,
                "total_coalitions": total_coalitions,
                "theory_evaluations": 1 + 2 * (total_coalitions - 1),
                "elapsed_seconds": elapsed,
                "n_replicas": nrep,
                "checks": checks,
                "value_cache_even": v_even_cache,
                "value_cache_odd": v_odd_cache,
                "value_cache_plus": v_plus_cache,
                "value_cache_minus": v_minus_cache,
                "chi2_plus_cache": chi2_plus_cache,
                "chi2_minus_cache": chi2_minus_cache,
                "_mean_value_cache": {
                    coalition: float(np.mean(v_even_cache[coalition]))
                    for coalition in all_coalitions
                },
                "_mean_value_cache_even": {
                    coalition: float(np.mean(v_even_cache[coalition]))
                    for coalition in all_coalitions
                },
                "_mean_value_cache_odd": {
                    coalition: float(np.mean(v_odd_cache[coalition]))
                    for coalition in all_coalitions
                },
                "coalition_bitmasks": coalition_bitmasks,
                "_sign_matrices": [plus_sign_matrix, minus_sign_matrix],
            }

        if verbose:
            print(f"Baseline        : {float(chi2_baseline):.6f}")
            print(f"Max |phi_even|  : {float(np.max(np.abs(phi_even))):.6f}")
            print(f"Max |phi_odd|   : {float(np.max(np.abs(phi_odd))):.6f}")
            print(f"Sum phi_even    : {float(np.sum(phi_even)):.6f}")
            print(f"Sum phi_odd     : {float(np.sum(phi_odd)):.6f}")
            print(f"Elapsed         : {elapsed:.1f}s")
            self._print_even_odd_validation_checks(checks, per_replica=False)

        return {
            "shapley_values": phi_even,
            "shapley_values_odd": np.abs(phi_odd),
            "shapley_values_odd_signed": phi_odd,
            "shapley_values_per_replica": None,
            "shapley_values_odd_per_replica": None,
            "shapley_std": None,
            "shapley_err": None,
            "shapley_odd_std": None,
            "shapley_odd_err": None,
            "mean_even": phi_even,
            "std_even_rep": None,
            "mean_odd": np.abs(phi_odd),
            "std_odd_rep": None,
            "baseline": float(chi2_baseline),
            "player_labels": self.flavor_labels,
            "player_short": self.flavor_short,
            "coalitions_evaluated": total_coalitions,
            "total_coalitions": total_coalitions,
            "theory_evaluations": 1 + 2 * (total_coalitions - 1),
            "elapsed_seconds": elapsed,
            "checks": checks,
            "value_cache_even": v_even_cache,
            "value_cache_odd": v_odd_cache,
            "value_cache_plus": v_plus_cache,
            "value_cache_minus": v_minus_cache,
            "chi2_plus_cache": chi2_plus_cache,
            "chi2_minus_cache": chi2_minus_cache,
            "_mean_value_cache": {
                coalition: float(v_even_cache[coalition])
                for coalition in all_coalitions
            },
            "_mean_value_cache_even": {
                coalition: float(v_even_cache[coalition])
                for coalition in all_coalitions
            },
            "_mean_value_cache_odd": {
                coalition: float(v_odd_cache[coalition])
                for coalition in all_coalitions
            },
            "coalition_bitmasks": coalition_bitmasks,
            "_sign_matrices": [plus_sign_matrix, minus_sign_matrix],
        }

    def _compute_exact_shap_parallel(self, mu, sigma, amplitude, mode, xspace,
                                     n_jobs, verbose=True,
                                     per_replica=False, random_sign=False,
                                     sign_matrices=None):
        """Compute exact Shapley values from a parallel coalition cache.

        When per_replica=True every coalition is evaluated independently for
        each replica, yielding a (nrep, n_flavors) Shapley matrix from which
        the ensemble mean phi_j = mean_k phi_j^(k) and standard deviation
        are derived.
        """
        n = self.n_flavors
        factorial = [math.factorial(k) for k in range(n + 1)]
        factorial_n = factorial[n]
        all_coalitions = list(self._iter_all_coalitions())
        total_coalitions = len(all_coalitions)
        value_cache = {}
        sign_matrices = [None] if sign_matrices is None else list(sign_matrices)
        n_sign_samples = len(sign_matrices)

        if verbose:
            print(
                f"Computing exact Shapley values for {n} players "
                f"({total_coalitions} coalitions) with {n_jobs} worker(s)..."
            )

        t0 = time.time()
        progress_start = t0
        interactive_stderr = sys.stderr.isatty()
        log_every = max(1, total_coalitions // 100)

        def _evaluate_one(coalition):
            if n_sign_samples == 1:
                value = self._evaluate_chi2(
                    coalition, mu, sigma, amplitude, mode=mode, xspace=xspace,
                    per_replica=per_replica, random_sign=random_sign,
                    random_sign_matrix=sign_matrices[0],
                )
            else:
                value = np.stack([
                    self._evaluate_chi2(
                        coalition, mu, sigma, amplitude, mode=mode, xspace=xspace,
                        per_replica=per_replica, random_sign=random_sign,
                        random_sign_matrix=sign_matrix,
                    )
                    for sign_matrix in sign_matrices
                ], axis=0)
            return coalition, value

        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            future_map = {
                executor.submit(_evaluate_one, coalition): coalition
                for coalition in all_coalitions
            }
            completed = 0
            for future in as_completed(future_map):
                coalition, value = future.result()
                value_cache[coalition] = value
                completed += 1

                if verbose:
                    elapsed = time.time() - progress_start
                    if completed == 1:
                        msg = (
                            f"    [{completed}/{total_coalitions}] "
                            f"elapsed {elapsed:.0f}s ..."
                        )
                    else:
                        rate = elapsed / completed
                        remaining = rate * (total_coalitions - completed)
                        mins, secs = divmod(int(remaining), 60)
                        msg = (
                            f"    [{completed}/{total_coalitions}] "
                            f"elapsed {elapsed:.0f}s | "
                            f"~{rate:.1f}s/eval | ETA {mins}m{secs:02d}s   "
                        )
                    if interactive_stderr:
                        sys.stderr.write(f"\r{msg}")
                        sys.stderr.flush()
                    elif (
                        completed == 1
                        or completed == total_coalitions
                        or completed % log_every == 0
                    ):
                        sys.stderr.write(f"{msg}\n")
                        sys.stderr.flush()

        if verbose:
            sys.stderr.write("\n")
            sys.stderr.flush()

        baseline_raw = value_cache[tuple()]
        mean_value_cache = {
            coalition: float(np.mean(value))
            for coalition, value in value_cache.items()
        }

        # per-replica path 
        if per_replica:
            first_value = next(iter(value_cache.values()))
            if n_sign_samples == 1:
                nrep = len(first_value)
                baseline = float(np.mean(baseline_raw))

                # shapley_per_replica[k, j] = phi_j^(k)
                shapley_per_replica = np.zeros((nrep, n), dtype=float)

                for player in range(n):
                    if verbose:
                        print(f"  Player {player}: {self.flavor_labels[player]}",
                              end="", flush=True)

                    for coalition in all_coalitions:
                        if player in coalition:
                            continue
                        s = len(coalition)
                        weight = (
                            factorial[s] * factorial[n - s - 1]
                        ) / factorial_n
                        coalition_with_player = self._coalition_with_player(
                            coalition, player
                        )
                        v_with = value_cache[coalition_with_player]
                        v_without = value_cache[coalition]
                        shapley_per_replica[:, player] += weight * (v_with - v_without)

                    if verbose:
                        sv_mean = float(np.mean(shapley_per_replica[:, player]))
                        sv_std = float(np.std(shapley_per_replica[:, player], ddof=1))
                        print(f"  ->  SV = {sv_mean:+.6f} +/- {sv_std:.6f}")

                shapley_vals = shapley_per_replica.mean(axis=0)
                shapley_std = (
                    shapley_per_replica.std(axis=0, ddof=1) if nrep > 1
                    else np.zeros(n, dtype=float)
                )
                shapley_err = shapley_std / np.sqrt(nrep)
                shapley_per_sign_sample = None
                shapley_sign_std = None
                shapley_sign_err = None
            else:
                nrep = int(first_value.shape[1])
                baseline = float(np.mean(baseline_raw))
                shapley_sample_replica = np.zeros((n_sign_samples, nrep, n), dtype=float)

                for player in range(n):
                    if verbose:
                        print(f"  Player {player}: {self.flavor_labels[player]}",
                              end="", flush=True)

                    for coalition in all_coalitions:
                        if player in coalition:
                            continue
                        s = len(coalition)
                        weight = (
                            factorial[s] * factorial[n - s - 1]
                        ) / factorial_n
                        coalition_with_player = self._coalition_with_player(
                            coalition, player
                        )
                        v_with = value_cache[coalition_with_player]
                        v_without = value_cache[coalition]
                        shapley_sample_replica[:, :, player] += (
                            weight * (v_with - v_without)
                        )

                    if verbose:
                        sv_mean = float(np.mean(shapley_sample_replica[:, :, player]))
                        sv_std = float(
                            np.std(
                                shapley_sample_replica.mean(axis=1)[:, player],
                                ddof=1,
                            )
                        ) if n_sign_samples > 1 else 0.0
                        print(f"  ->  SV = {sv_mean:+.6f} +/- {sv_std:.6f}")

                shapley_per_replica = shapley_sample_replica.mean(axis=0)
                shapley_per_sign_sample = shapley_sample_replica.mean(axis=1)
                shapley_vals = shapley_per_replica.mean(axis=0)
                shapley_std = (
                    shapley_per_replica.std(axis=0, ddof=1) if nrep > 1
                    else np.zeros(n, dtype=float)
                )
                shapley_err = shapley_std / np.sqrt(nrep)
                shapley_sign_std = (
                    shapley_per_sign_sample.std(axis=0, ddof=1)
                    if n_sign_samples > 1 else np.zeros(n, dtype=float)
                )
                shapley_sign_err = shapley_sign_std / np.sqrt(n_sign_samples)

            elapsed = time.time() - t0

            if verbose:
                print(f"Baseline (mean) : {baseline:.6f}")
                print(f"Max |SV|        : {np.max(np.abs(shapley_vals)):.6f}")
                print(f"Mean |SV|       : {np.mean(np.abs(shapley_vals)):.6f}")
                print(f"Max SV std      : {np.max(shapley_std):.6f}")
                print(f"Elapsed         : {elapsed:.1f}s")

            return {
                "shapley_values": shapley_vals,
                "shapley_values_per_replica": shapley_per_replica,
                "shapley_values_per_sign_sample": shapley_per_sign_sample,
                "shapley_std": shapley_std,
                "shapley_err": shapley_err,
                "shapley_sign_std": shapley_sign_std,
                "shapley_sign_err": shapley_sign_err,
                "baseline": baseline,
                "baseline_per_replica": baseline_raw,
                "player_labels": self.flavor_labels,
                "player_short": self.flavor_short,
                "coalitions_evaluated": len(value_cache),
                "total_coalitions": total_coalitions,
                "elapsed_seconds": elapsed,
                "n_replicas": nrep,
                "n_sign_samples": n_sign_samples,
                "_mean_value_cache": mean_value_cache,
            }

        # scalar (mean-chi2) path 
        if n_sign_samples == 1:
            baseline = float(baseline_raw)
            shapley_vals = np.zeros(n, dtype=float)

            for player in range(n):
                if verbose:
                    print(f"  Player {player}: {self.flavor_labels[player]}", end="",
                          flush=True)

                for coalition in all_coalitions:
                    if player in coalition:
                        continue

                    s = len(coalition)
                    weight = (
                        factorial[s] * factorial[n - s - 1]
                    ) / factorial_n
                    coalition_with_player = self._coalition_with_player(
                        coalition, player
                    )
                    v_with = value_cache[coalition_with_player]
                    v_without = value_cache[coalition]
                    shapley_vals[player] += weight * (v_with - v_without)

                if verbose:
                    print(f"  ->  SV = {shapley_vals[player]:+.6f}")

            shapley_per_sign_sample = None
            shapley_std = None
            shapley_err = None
            shapley_sign_std = None
            shapley_sign_err = None
        else:
            baseline = float(np.mean(baseline_raw))
            shapley_per_sign_sample = np.zeros((n_sign_samples, n), dtype=float)

            for player in range(n):
                if verbose:
                    print(f"  Player {player}: {self.flavor_labels[player]}", end="",
                          flush=True)

                for coalition in all_coalitions:
                    if player in coalition:
                        continue

                    s = len(coalition)
                    weight = (
                        factorial[s] * factorial[n - s - 1]
                    ) / factorial_n
                    coalition_with_player = self._coalition_with_player(
                        coalition, player
                    )
                    v_with = value_cache[coalition_with_player]
                    v_without = value_cache[coalition]
                    shapley_per_sign_sample[:, player] += weight * (v_with - v_without)

                if verbose:
                    print(
                        f"  ->  SV = {float(np.mean(shapley_per_sign_sample[:, player])):+.6f}"
                    )

            shapley_vals = shapley_per_sign_sample.mean(axis=0)
            shapley_sign_std = (
                shapley_per_sign_sample.std(axis=0, ddof=1)
                if n_sign_samples > 1 else np.zeros(n, dtype=float)
            )
            shapley_sign_err = shapley_sign_std / np.sqrt(n_sign_samples)
            shapley_std = shapley_sign_std
            shapley_err = shapley_sign_err

        elapsed = time.time() - t0

        if verbose:
            print(f"Baseline        : {baseline:.6f}")
            print(f"Max |SV|        : {np.max(np.abs(shapley_vals)):.6f}")
            print(f"Mean |SV|       : {np.mean(np.abs(shapley_vals)):.6f}")
            print(f"Sum SV          : {np.sum(shapley_vals):.6f}")
            print(f"Sum |SV|        : {np.sum(np.abs(shapley_vals)):.6f}")
            print(f"Elapsed         : {elapsed:.1f}s")

        return {
            "shapley_values": shapley_vals,
            "shapley_values_per_replica": None,
            "shapley_values_per_sign_sample": shapley_per_sign_sample,
            "shapley_std": shapley_std,
            "shapley_err": shapley_err,
            "shapley_sign_std": shapley_sign_std,
            "shapley_sign_err": shapley_sign_err,
            "baseline": baseline,
            "player_labels": self.flavor_labels,
            "player_short": self.flavor_short,
            "coalitions_evaluated": len(value_cache),
            "total_coalitions": total_coalitions,
            "elapsed_seconds": elapsed,
            "n_sign_samples": n_sign_samples,
            "value_cache": {
                str(k): v for k, v in value_cache.items()
            },
            "_mean_value_cache": mean_value_cache,
        }

    def exact_shap(self, mu, sigma, amplitude, mode='additive',
                   xspace='linear', plot=True, n_jobs=1,
                   diagnostic=False, outlier_n_sigma=3.0,
                   per_replica=False, random_sign=False,
                   n_sign_samples=1, random_seed=None,
                   deterministic_sign_symmetrized=True):
        """Compute exact Shapley values for all flavour players.

        Uses the ``shapley_values.ExactShapley`` solver with the NNPDF
        chi2 value function.

        Parameters
        ----------
        mu, sigma, amplitude : float
        mode : str
        xspace : str
        plot : bool
            Whether to generate plots.
        n_jobs : int
            Number of worker threads for coalition evaluation.
            If n_jobs=1, use the serial ExactShapley solver (scalar path)
            or the single-worker parallel path (per_replica=True).
        diagnostic : bool
            When True, record chi2 for every coalition evaluated and compute
            per-coalition and per-player-marginal statistics.  The results
            dict will contain ``coalition_log`` (raw list of
            ``(coalition_tuple, chi2)`` pairs) and ``diagnostic`` (statistics
            dict).  Useful for spotting extreme coalitions that corrupt the
            Shapley values.
        outlier_n_sigma : float
            Z-score threshold used to flag outlier coalitions/marginals.
            Default is 3.0 (mean +/- 3 sigma).
        per_replica : bool
            When True compute phi_j^(k) for every replica k and return the
            full (nrep, n_flavors) matrix together with the ensemble mean,
            standard deviation, and standard error.  Forces use of the
            parallel evaluation path (n_jobs>=1).  This implements
            Eq. (shapley_replicas) from the paper.
        random_sign : bool
            When True draw a fixed sign table indexed by replica and flavour,
            and reuse it for every coalition evaluation in the run.
        n_sign_samples : int
            Number of independent sign-table games to average when
            ``random_sign=True``. Defaults to 1.
        random_seed : int or None
            Optional seed for reproducible sign-mask sampling.
        deterministic_sign_symmetrized : bool
            When True (default), use deterministic calibrated up/down dual
            games. When False, use the legacy sign-sampling path controlled
            by ``random_sign`` and ``n_sign_samples``.

        Returns
        -------
        results : dict
        """
        if mode not in PERTURBATION_MODES:
            raise ValueError(
                f"Unknown mode '{mode}'. Choose from {PERTURBATION_MODES}."
            )
        if xspace not in PERTURBATION_XSPACES:
            raise ValueError(
                f"Unknown xspace '{xspace}'. Choose from {PERTURBATION_XSPACES}."
            )
        if int(n_jobs) < 1:
            raise ValueError(f"n_jobs must be >= 1, got {n_jobs}")
        if int(n_sign_samples) < 1:
            raise ValueError(f"n_sign_samples must be >= 1, got {n_sign_samples}")

        deterministic_sign_symmetrized = bool(deterministic_sign_symmetrized)

        basis_name = "flavor" if self.basis == 'flavor' else "evolution"
        print(f"Perturbation basis : {basis_name}")
        print(f"Perturbation mode  : {mode}")
        print(f"Perturbation xspace: {xspace}")
        print(f"Sum rules          : {'ON' if self.enforce_sumrules else 'OFF'}")
        print(f"PDF member mode    : {self.member_mode}")
        print(f"Parallel workers   : {int(n_jobs)}")
        print(f"Per-replica SVs    : {'ON' if per_replica else 'OFF'}")
        if deterministic_sign_symmetrized:
            print("Sign game          : deterministic calibrated up/down")
            if random_sign:
                print(
                    "Note: deterministic calibrated up/down games are enabled; "
                    "ignoring random_sign=True."
                )
            if int(n_sign_samples) != 1:
                print(
                    "Note: deterministic calibrated up/down games do not use "
                    f"n_sign_samples={int(n_sign_samples)}; forcing to 1."
                )

            effective_random_sign = False
            effective_n_sign_samples = 1
            results = self._compute_exact_shap_even_odd_parallel(
                mu, sigma, amplitude, mode, xspace,
                n_jobs=int(n_jobs),
                per_replica=per_replica,
                verbose=True,
            )
            mean_value_cache_for_diag = results.get("_mean_value_cache_even", {})
        else:
            print("Sign game          : configurable random/fixed signs")
            effective_random_sign = bool(random_sign)
            if effective_random_sign:
                effective_n_sign_samples = int(n_sign_samples)
                sign_matrices = self._build_sign_matrices(
                    n_members=self._infer_n_members(),
                    n_flavors=self.n_flavors,
                    n_sign_samples=effective_n_sign_samples,
                    random_seed=random_seed,
                    unique_rows=(
                        self.member_mode == "replicas"
                        and bool(per_replica)
                        and int(n_sign_samples) == 1
                    ),
                )
            else:
                if int(n_sign_samples) != 1:
                    print(
                        "Note: random_sign=False so n_sign_samples is ignored; "
                        "forcing to 1."
                    )
                effective_n_sign_samples = 1
                sign_matrices = [None]

            results = self._compute_exact_shap_parallel(
                mu, sigma, amplitude, mode, xspace,
                n_jobs=int(n_jobs),
                verbose=True,
                per_replica=per_replica,
                random_sign=effective_random_sign,
                sign_matrices=sign_matrices,
            )
            results["_sign_matrices"] = sign_matrices
            results.setdefault("shapley_values_odd", None)
            results.setdefault("shapley_values_odd_per_replica", None)
            results.setdefault("shapley_values_odd_signed", None)
            results.setdefault("shapley_values_odd_signed_per_replica", None)
            mean_value_cache_for_diag = results.get("_mean_value_cache", {})

        coalition_log = None
        if diagnostic:
            coalition_log = [
                (coalition, mean_value_cache_for_diag[coalition])
                for coalition in self._iter_all_coalitions()
                if coalition in mean_value_cache_for_diag
            ]

        # Attach coalition log and diagnostics if requested.
        if diagnostic and coalition_log is not None:
            results["coalition_log"] = coalition_log
            results["diagnostic"] = self._compute_diagnostics(
                coalition_log, outlier_n_sigma=outlier_n_sigma
            )
            diag = results["diagnostic"]
            print("\n--- Coalition diagnostics ---")
            cs = diag.get("chi2_stats", {})
            print(
                f"  Coalitions : {cs.get('n_coalitions', '?')}  "
                f"chi2 mean={cs.get('mean', 0):.4f}  "
                f"std={cs.get('std', 0):.4f}  "
                f"min={cs.get('min', 0):.4f}  "
                f"max={cs.get('max', 0):.4f}"
            )
            print(
                f"  Outlier threshold ({outlier_n_sigma} sigma): "
                f"{diag.get('outlier_chi2_threshold', 0):.4f}  "
                f"({diag.get('n_outlier_coalitions', 0)} outlier coalitions)"
            )
            if diag.get("outlier_coalitions"):
                print("  Top outliers (chi2 > threshold):")
                for oc in diag["outlier_coalitions"][:5]:
                    print(
                        f"    coalition={oc['coalition_labels']}  "
                        f"chi2={oc['chi2']:.4f}  z={oc['z_score']:.2f}"
                    )
            print("  Per-player marginal |max|:")
            for lbl, pm in diag.get("per_player_marginals", {}).items():
                print(
                    f"    {lbl:>8s}: mean={pm['mean']:+.4f}  "
                    f"std={pm['std']:.4f}  ["
                    f"{pm['min']:+.4f}, {pm['max']:+.4f}]  "
                    f"outliers: {pm['n_outliers_high']+pm['n_outliers_low']}"
                )
            print()
        else:
            results["coalition_log"] = None
            results["diagnostic"] = None

        # Add NNPDF-specific metadata
        results["baseline_chi2"] = results["baseline"]
        results["mode"] = mode
        results["xspace"] = xspace
        results["enforce_sumrules"] = self.enforce_sumrules
        results["basis"] = self.basis
        results["flavor_labels"] = self.flavor_labels
        results["flavor_short"] = self.flavor_short
        results["n_jobs"] = int(n_jobs)
        results["per_replica"] = bool(per_replica)
        results["random_sign"] = bool(effective_random_sign)
        results["member_mode"] = self.member_mode
        results["n_sign_samples"] = int(effective_n_sign_samples)
        results["antithetic_sign"] = bool(
            effective_random_sign and int(effective_n_sign_samples) > 1
        )
        results["random_seed"] = random_seed
        results["deterministic_sign_symmetrized"] = bool(
            deterministic_sign_symmetrized
        )
        results["deterministic_calibrated_updown"] = bool(
            deterministic_sign_symmetrized
        )
        # Convenience: scalar uncertainty arrays (None when per_replica=False)
        # shapley_std[j]  = std_k phi_j^(k)  (replica-to-replica spread)
        # shapley_err[j]  = shapley_std[j] / sqrt(nrep)  (standard error of mean)
        if results.get("shapley_std") is None:
            results["shapley_std"] = None
        if results.get("shapley_err") is None:
            results["shapley_err"] = None
        if results.get("shapley_values_per_sign_sample") is None:
            results["shapley_values_per_sign_sample"] = None
        if results.get("shapley_sign_std") is None:
            results["shapley_sign_std"] = None
        if results.get("shapley_sign_err") is None:
            results["shapley_sign_err"] = None
        results.pop("_mean_value_cache", None)

        # Plot
        fig_pdfs = None
        fig_bar = None
        fig_bar_odd = None
        if plot:
            fig_pdfs = self.plot_pdfs(
                amplitude=amplitude, mu=mu, sigma=sigma,
                mode=mode, xspace=xspace,
            )
            if mode == 'ablation':
                bar_title = (
                    f"PDF Flavour Importance ({basis_name})  "
                    f"mode={mode}, xspace={xspace}"
                )
            elif mode == 'calibrated':
                bar_title = (
                    f"PDF Flavour Importance ({basis_name})  "
                    f"mu={mu}, sigma={sigma}, A={amplitude}\u03c3_rep, "
                    f"mode={mode}, xspace={xspace}"
                )
            else:
                bar_title = (
                    f"PDF Flavour Importance ({basis_name})  "
                    f"mu={mu}, sigma={sigma}, A={amplitude}, "
                    f"mode={mode}, xspace={xspace}"
                )
            # Single main bar plot from phi_even only.
            # Odd quantities are kept in outputs but not plotted.
            sv_even = np.asarray(results["shapley_values"], dtype=float)
            fig_bar, ax = plt.subplots(figsize=(12, 6))
            pos_color = "#518500"
            neg_color = "#FFBF00"
            bar_colors = [pos_color if v >= 0.0 else neg_color for v in sv_even]
            ax.bar(self.flavor_short, sv_even, color=bar_colors, alpha=0.85)
            ax.axhline(0.0, color="#444444", lw=0.9, alpha=0.9)
            ax.set_ylabel("Shapley Value (delta chi2/N)")
            ax.set_title(bar_title)
            ax.grid(axis="y", ls="--", alpha=0.35)
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

            _sv_std = results.get("shapley_std")
            if _sv_std is not None:
                ax.errorbar(
                    self.flavor_short,
                    sv_even,
                    yerr=np.asarray(_sv_std, dtype=float),
                    fmt="none",
                    ecolor="black",
                    elinewidth=1.2,
                    capsize=4,
                    capthick=1.2,
                    zorder=5,
                )
            ax.legend(
                handles=[
                    plt.Rectangle((0, 0), 1, 1, color=pos_color, alpha=0.85, label="SV > 0"),
                    plt.Rectangle((0, 0), 1, 1, color=neg_color, alpha=0.85, label="SV < 0"),
                ],
                loc="upper right",
                frameon=False,
            )
            fig_bar.tight_layout()

            fig_bar_odd = None
            plt.show()

        results["fig_pdfs"] = fig_pdfs
        results["fig_bar"] = fig_bar
        results["fig_bar_even"] = fig_bar
        results["fig_bar_odd"] = fig_bar_odd

        return results

    # -- Plotting -----------------------------------------------------------

    def plot_pdfs(self, amplitude=0.08, mu=0.02, sigma=0.1,
                  mode='additive', xspace='linear', x_points=200):
        """Plot reference vs perturbed PDFs for all Shapley-player flavours."""
        x_min = 1e-5
        x_max = 10 ** (-0.001)
        if mode == "ablation":
            x_plot = np.linspace(x_min, x_max, x_points)
            x_axis_scale = "linear"
        else:
            x_plot = np.logspace(-5, -0.001, x_points)
            x_axis_scale = "linear" if float(mu) > 0.1 else "log"
        Q0 = self.observables[0].Q0

        class PlotTarget:
            pass

        plot_target = PlotTarget()
        plot_target.Q0 = Q0
        plot_target.xgrid = x_plot
        plot_target.flavor_indices = np.asarray(self.flavor_indices)

        if self.basis == 'flavor':
            gv_ref_all = get_pdf_flavor_grid_values(
                self.pdf, plot_target, n_replicas=self.n_replicas,
                member_mode=self.member_mode,
            )
            gv_calib_all = get_pdf_flavor_grid_values(
                self.pdf, plot_target, n_replicas=self.n_replicas,
                member_mode='replicas',
            )
            # get_pdf_flavor_grid_values returns all 14 PDG flavours in
            # ALL_FLAVOURS order. Shapley players are a subset given by
            # self.flavor_indices (indices into that 14-flavour axis).
            sel = np.asarray(self.flavor_indices, dtype=int)
            gv_ref = gv_ref_all[:, sel, :]
            gv_calib = gv_calib_all[:, sel, :]
        else:
            gv_ref = get_pdf_grid_values(
                self.pdf, plot_target, n_replicas=self.n_replicas,
                member_mode=self.member_mode,
            )
            gv_calib = get_pdf_grid_values(
                self.pdf, plot_target, n_replicas=self.n_replicas,
                member_mode='replicas',
            )

        gv_pert = apply_gaussian_perturbation(
            gv_ref,
            local_flavor_idx=list(range(gv_ref.shape[1])),
            mu=mu,
            sigma=sigma,
            amplitude=amplitude,
            xgrid=x_plot,
            mode=mode,
            xspace=xspace,
            calibration_gv=gv_calib,
        )

        n = len(self.flavor_short)
        ncols = min(3, n)
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows))
        axes = np.atleast_1d(axes).ravel()

        for i in range(n):
            ax = axes[i]
            ref_mean = np.mean(gv_ref[:, i, :], axis=0)
            ref_std = np.std(gv_ref[:, i, :], axis=0)
            pert_mean = np.mean(gv_pert[:, i, :], axis=0)

            ax.plot(x_plot, ref_mean, "b-", lw=1.5, label="Reference")
            ax.fill_between(
                x_plot, ref_mean - ref_std, ref_mean + ref_std,
                alpha=0.25, color="blue"
            )
            ax.plot(x_plot, pert_mean, "r--", lw=1.5, label="Perturbed")
            ax.set_xscale(x_axis_scale)
            ax.set_title(self.flavor_labels[i])
            ax.set_xlabel("x")
            ax.set_ylabel("xf(x)")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3, ls="--")

        for j in range(n, len(axes)):
            axes[j].axis("off")

        basis_name = "Flavor" if self.basis == 'flavor' else "Evolution"
        fig.suptitle(
            f"{basis_name}-basis PDFs: A={amplitude}, mu={mu}, "
            f"sigma={sigma}, mode={mode}, xspace={xspace}, x-axis={x_axis_scale}",
            fontsize=14, fontweight="bold", y=1.01,
        )
        plt.tight_layout()
        plt.show()
        return fig

    def plot_predictions(self):
        """Plot theory predictions vs experimental data."""
        n_obs = len(self.observables)
        ncols = min(3, n_obs)
        nrows = int(np.ceil(n_obs / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows))
        axes = np.atleast_1d(axes).ravel()

        for i, obs in enumerate(self.observables):
            ax = axes[i]
            gv_list = self._get_gv_list(obs)
            preds = obs.convolve(gv_list)
            pred_mean = preds.mean(axis=1)
            pred_std = preds.std(axis=1)
            idx = np.arange(obs.ndata)

            ax.scatter(idx, obs.data_central, s=12, c="blue", alpha=0.7,
                       label="Data", zorder=3, edgecolors="darkblue", lw=0.3)
            ax.errorbar(idx, pred_mean, yerr=pred_std, fmt="x", color="red",
                        capsize=2, ms=4, lw=1, label="Prediction", alpha=0.8)
            ax.set_title(obs.name, fontsize=8, fontweight="bold")
            ax.set_xlabel("Data point")
            ax.set_ylabel("Observable")
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3, ls="--")

        for j in range(n_obs, len(axes)):
            axes[j].axis("off")

        fig.suptitle("Predictions vs Data",
                     fontsize=13, fontweight="bold")
        plt.tight_layout()
        plt.show()
        return fig
