"""NNPDF Shapley value analyzer.

Wraps NNPDF chi2 evaluation into a value function for the generic
ExactShapley solver from the ``shapley_values`` package. Handles:
  - Evolution and physical-flavor perturbation bases
  - Multi-FK datasets with operations (ASY, RATIO, ADD)
  - MSR + VSR sum rule enforcement
  - Grid-value caching for performance
"""

import functools
import time

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from validphys.convolution import FK_FLAVOURS
from validphys.pdfbases import evolution as evol_basis, ALL_FLAVOURS
from validphys.gridvalues import grid_values as _lhapdf_grid_values

from shapley_values import ExactShapley, plot_shapley_bar

from .setup import (
    get_pdf_grid_values,
    get_pdf_flavor_grid_values,
    get_pdf_grid_values_all14,
    FLAVOR_PDG_NAMES,
)
from .perturbation import (
    gaussian_profile,
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
    - flavor : players are physical quarks/gluon; rotated before convolution.

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
                 basis='evolution', enforce_sumrules=False):
        self.pdf = pdf
        self.observables = observables
        self.flavor_info = flavor_info
        self.n_replicas = n_replicas
        self.basis = basis
        self.enforce_sumrules = enforce_sumrules

        self._gv_cache: dict = {}

        # Sum rule integration grid (lazily initialised)
        self._sr_xgrid = None
        self._sr_weights = None
        self._sr_gv_evol = None
        self._sr_gv_flav = None
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
        self._sr_xgrid, self._sr_weights = gen_integration_input(2000)
        Q0 = self.observables[0].Q0

        self._sr_gv_evol = evol_basis.grid_values(
            self.pdf, qmat=[Q0], vmat=FK_FLAVOURS, xmat=self._sr_xgrid
        ).squeeze(-1)
        if self.n_replicas is not None:
            self._sr_gv_evol = self._sr_gv_evol[1: self.n_replicas + 1]

        if self.basis == 'flavor':
            self._sr_gv_flav = _lhapdf_grid_values(
                self.pdf, list(ALL_FLAVOURS), self._sr_xgrid, [Q0]
            ).squeeze(-1)
            if self.n_replicas is not None:
                self._sr_gv_flav = self._sr_gv_flav[1: self.n_replicas + 1]

            evol_row_inds = evol_basis._to_indexes(FK_FLAVOURS)
            self._sr_rotation = evol_basis.from_flavour_mat[evol_row_inds, :]

    def _compute_sumrule_norm(self, flavor_subset, mu, sigma, amplitude,
                              mode, xspace):
        """Compute normalization constants for a perturbed coalition."""
        if self._sr_xgrid is None:
            self._setup_sumrule_grid()

        if self.basis == 'evolution':
            gv = self._sr_gv_evol.copy()
            perturb_idx = [self.flavor_indices[p] for p in flavor_subset]
            gv_pert = apply_gaussian_perturbation(
                gv, perturb_idx, mu, sigma, amplitude,
                self._sr_xgrid, mode=mode, xspace=xspace
            )
        else:
            gv_flav = self._sr_gv_flav.copy()
            perturb_idx = [self.flavor_indices[p] for p in flavor_subset]
            gv_flav_pert = apply_gaussian_perturbation(
                gv_flav, perturb_idx, mu, sigma, amplitude,
                self._sr_xgrid, mode=mode, xspace=xspace
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

    def _get_gv_for_entry(self, obs, entry_idx):
        """Cached evolution-basis grid values (FK-subset flavours)."""
        key = (obs.name, entry_idx, 'evol')
        if key not in self._gv_cache:
            entry = obs.fk_entries[entry_idx]
            self._gv_cache[key] = get_pdf_grid_values(
                self.pdf, entry, n_replicas=self.n_replicas
            )
        return self._gv_cache[key]

    def _get_flavor_gv_for_entry(self, obs, entry_idx):
        """Cached physical flavor-basis grid values."""
        key = (obs.name, entry_idx, 'flavor')
        if key not in self._gv_cache:
            entry = obs.fk_entries[entry_idx]
            self._gv_cache[key] = get_pdf_flavor_grid_values(
                self.pdf, entry, n_replicas=self.n_replicas
            )
        return self._gv_cache[key]

    def _get_gv_all14_for_entry(self, obs, entry_idx):
        """Cached full 14-flavour evolution-basis grid values."""
        key = (obs.name, entry_idx, 'evol_all14')
        if key not in self._gv_cache:
            entry = obs.fk_entries[entry_idx]
            self._gv_cache[key] = get_pdf_grid_values_all14(
                self.pdf, entry, n_replicas=self.n_replicas
            )
        return self._gv_cache[key]

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

    # -- Chi2 evaluation with perturbation ----------------------------------

    def _evaluate_chi2(self, flavor_subset, mu, sigma, amplitude,
                       mode='additive', xspace='linear'):
        """Mean chi2/Ndata across all observables for a given perturbation.

        Parameters
        ----------
        flavor_subset : list of int
            Global Shapley-player indices to perturb.
        mu, sigma, amplitude : float
        mode : str
        xspace : str

        Returns
        -------
        chi2 : float
        """
        sr_norm = None
        if self.enforce_sumrules:
            sr_norm = self._compute_sumrule_norm(
                flavor_subset, mu, sigma, amplitude, mode, xspace
            )

        total_chi2 = 0.0
        total_ndata = 0

        for obs in self.observables:
            if self.basis == 'flavor':
                gv_pert_list = []
                for idx, entry in enumerate(obs.fk_entries):
                    gv_flav = self._get_flavor_gv_for_entry(obs, idx)
                    perturb_idx = [
                        self.flavor_indices[p] for p in flavor_subset
                    ]
                    gv_pert = apply_gaussian_perturbation(
                        gv_flav, perturb_idx, mu, sigma, amplitude,
                        entry.xgrid, mode=mode, xspace=xspace
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
                        perturb_idx = [
                            self.flavor_indices[p] for p in flavor_subset
                        ]
                    else:
                        gv = self._get_gv_for_entry(obs, idx)
                        perturb_idx = self._local_flavor_indices_for_entry(
                            entry, flavor_subset
                        )
                    gv_pert = apply_gaussian_perturbation(
                        gv, perturb_idx, mu, sigma, amplitude,
                        entry.xgrid, mode=mode, xspace=xspace
                    )
                    if sr_norm is not None:
                        fi = (range(14) if entry.hadronic
                              else entry.flavor_indices)
                        gv_pert = self._apply_norm_to_gv(
                            gv_pert, sr_norm, fi
                        )
                    gv_pert_list.append(gv_pert)
                chi2_arr = obs.chi2(gv_pert_list)

            total_chi2 += np.mean(chi2_arr)
            total_ndata += obs.ndata

        return total_chi2 / total_ndata

    # -- Build value function for ExactShapley ------------------------------

    def build_value_function(self, mu, sigma, amplitude,
                             mode='additive', xspace='linear'):
        """Return a callable v(coalition) -> float for ExactShapley.

        Parameters
        ----------
        mu, sigma, amplitude : float
            Gaussian perturbation parameters.
        mode : str
        xspace : str

        Returns
        -------
        value_function : callable
            f(List[int]) -> float
        """
        def v(coalition):
            return self._evaluate_chi2(
                coalition, mu, sigma, amplitude,
                mode=mode, xspace=xspace,
            )
        return v

    # -- Convenience: run full Shapley analysis -----------------------------

    def exact_shap(self, mu, sigma, amplitude, mode='additive',
                   xspace='linear', plot=True):
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

        basis_name = "flavor" if self.basis == 'flavor' else "evolution"
        print(f"Perturbation basis : {basis_name}")
        print(f"Perturbation mode  : {mode}")
        print(f"Perturbation xspace: {xspace}")
        print(f"Sum rules          : {'ON' if self.enforce_sumrules else 'OFF'}")

        # Build value function and run ExactShapley
        v = self.build_value_function(mu, sigma, amplitude, mode, xspace)

        solver = ExactShapley(
            n_players=self.n_flavors,
            value_function=v,
            player_labels=self.flavor_labels,
            player_short=self.flavor_short,
        )

        results = solver.compute(verbose=True)

        # Add NNPDF-specific metadata
        results["baseline_chi2"] = results["baseline"]
        results["mode"] = mode
        results["xspace"] = xspace
        results["enforce_sumrules"] = self.enforce_sumrules
        results["basis"] = self.basis
        results["flavor_labels"] = self.flavor_labels
        results["flavor_short"] = self.flavor_short

        # Plot
        fig_pdfs = None
        fig_bar = None
        if plot:
            fig_pdfs = self.plot_pdfs(
                amplitude=amplitude, mu=mu, sigma=sigma,
                mode=mode, xspace=xspace,
            )
            fig_bar = plot_shapley_bar(
                results["shapley_values"],
                self.flavor_short,
                title=(
                    f"PDF Flavour Importance ({basis_name})  "
                    f"mu={mu}, sigma={sigma}, A={amplitude}, "
                    f"mode={mode}, xspace={xspace}"
                ),
                ylabel="Shapley Value (delta chi2/N)",
            )
            plt.show()

        results["fig_pdfs"] = fig_pdfs
        results["fig_bar"] = fig_bar

        return results

    # -- Plotting -----------------------------------------------------------

    def plot_pdfs(self, amplitude=0.08, mu=0.02, sigma=0.1,
                  mode='additive', xspace='linear', x_points=200):
        """Plot reference vs perturbed PDFs for all Shapley-player flavours."""
        x_plot = np.logspace(-5, -0.001, x_points)
        Q0 = self.observables[0].Q0

        if self.basis == 'flavor':
            pdg_codes = [int(ALL_FLAVOURS[i]) for i in self.flavor_indices]
            gv_ref = _lhapdf_grid_values(
                self.pdf, pdg_codes, x_plot, [Q0]
            ).squeeze(-1)
        else:
            gv_func = functools.partial(evol_basis.grid_values, self.pdf)
            all_names = FK_FLAVOURS[self.flavor_indices]
            gv_ref = gv_func(
                qmat=[Q0], vmat=all_names, xmat=x_plot
            ).squeeze(-1)

        if self.n_replicas is not None:
            gv_ref = gv_ref[1: self.n_replicas + 1]

        gauss = gaussian_profile(x_plot, mu, sigma, amplitude, xspace)
        gv_pert = gv_ref.copy()
        if mode == 'additive':
            for i in range(gv_pert.shape[1]):
                gv_pert[:, i, :] += gauss[np.newaxis, :]
        else:
            for i in range(gv_pert.shape[1]):
                gv_pert[:, i, :] *= (1.0 + gauss[np.newaxis, :])

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
            ax.set_xscale("log")
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
            f"sigma={sigma}, mode={mode}, xspace={xspace}",
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
