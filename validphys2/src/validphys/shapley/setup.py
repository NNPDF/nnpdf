"""Dataset loading, FK convolution, and chi2 for Shapley analysis.

Loads datasets via validphys, builds dense FK tensors for fast convolution, and provides chi2 computation. Supports multi-FK datasets with operations (NULL, ASY, RATIO, ADD, COM, SMN, SMT).

Why not use validphys.convolution directly?
  1. Dense tensors: validphys keeps FK tables as sparse DataFrames and re-indexes via groupby each call. We pre-densify once for fast numpy einsum/matmul over O(N * 2^N) evaluations.
  2. Decoupled interface: validphys.convolution.predictions() is coupled to LHAPDF PDF objects. We need to inject perturbed numpy arrays.
  3. Flavor rotation: validphys has no flavor-to-evolution rotation (it always starts from evolution basis). We add it for flavor-basis perturbation.

The underlying math (convolution, chi2, multi-FK operations) is identical
to validphys.
"""

import functools
from typing import List, Optional

import numpy as np
import pandas as pd
from scipy.linalg import cholesky, solve_triangular

from validphys.loader import Loader
from validphys.fkparser import load_fktable
from validphys.convolution import FK_FLAVOURS, NFK
from validphys.pdfbases import evolution as evol_basis, ALL_FLAVOURS
from validphys.covmats import construct_covmat
from validphys.gridvalues import grid_values as _lhapdf_grid_values

# PDG code to human-readable label
FLAVOR_PDG_NAMES = {
    -6: "tbar", -5: "bbar", -4: "cbar", -3: "sbar",
    -2: "ubar", -1: "dbar", 21: "g",
    1: "d", 2: "u", 3: "s", 4: "c", 5: "b", 6: "t", 22: "photon",
}


# ---------------------------------------------------------------------------
# Single FK table with precomputed dense tensor
# ---------------------------------------------------------------------------

class FKEntry:
    """One FK table with precomputed dense tensor and convolution methods.

    NNPDFObservable holds one or more of these and combines their
    predictions via the dataset operation.

    Converts the validphys sparse DataFrame (fk.sigma) into dense numpy
    arrays at construction time. Also precomputes the flavor-to-evolution
    rotation matrix.
    """

    def __init__(self, sigma, flavor_indices, xgrid, Q0, ndata, hadronic):
        self.flavor_indices = flavor_indices
        self.xgrid = xgrid
        self.Q0 = Q0
        self.ndata = ndata
        self.hadronic = hadronic
        self._build_rotation()
        self._build_dense_tensor(sigma)

    def _build_rotation(self):
        """Build flavor-to-evolution rotation matrix for this FK table.

        Uses evol_basis.from_flavour_mat from validphys.pdfbases.
        """
        if self.hadronic:
            evol_row_inds = evol_basis._to_indexes(FK_FLAVOURS)
        else:
            evol_row_inds = evol_basis._to_indexes(
                FK_FLAVOURS[self.flavor_indices]
            )
        fl_col_mask = evol_basis._flmask_from_base_indexes(evol_row_inds)
        self.fl_col_inds = np.where(fl_col_mask)[0]
        self.transform_mat = evol_basis.from_flavour_mat[
            np.ix_(evol_row_inds, self.fl_col_inds)
        ]

    def _build_dense_tensor(self, sigma):
        """Convert sparse sigma DataFrame into dense numpy arrays.

        DIS:      _W_dis      shape (ndata, nfl,   nx)
        Hadronic: _W_had_flat  shape (ndata, npairs*nx*nx)
                  plus pair-index arrays _fl1, _fl2.
        """
        fm = sigma.columns
        nx = len(self.xgrid)

        if not self.hadronic:
            nfl = len(fm)
            W = np.zeros((self.ndata, nfl, nx), dtype=np.float64)
            data_idx = sigma.index.get_level_values(0)
            x_idx = sigma.index.get_level_values(1)
            vals = sigma.values
            unique_data = data_idx.unique()
            data_map = {d: i for i, d in enumerate(unique_data)}
            d_mapped = np.array([data_map[d] for d in data_idx])
            for fi in range(nfl):
                W[d_mapped, fi, x_idx] = vals[:, fi]
            self._W_dis = W
        else:
            npairs = len(fm)
            W = np.zeros((self.ndata, npairs, nx, nx), dtype=np.float64)
            data_idx = sigma.index.get_level_values(0)
            x1_idx = sigma.index.get_level_values(1)
            x2_idx = sigma.index.get_level_values(2)
            vals = sigma.values
            unique_data = data_idx.unique()
            data_map = {d: i for i, d in enumerate(unique_data)}
            d_mapped = np.array([data_map[d] for d in data_idx])
            for pi in range(npairs):
                W[d_mapped, pi, x1_idx, x2_idx] = vals[:, pi]
            self._W_had_flat = W.reshape(self.ndata, -1)
            all_fl1, all_fl2 = np.indices((NFK, NFK))
            self._fl1 = all_fl1.ravel()[fm]
            self._fl2 = all_fl2.ravel()[fm]

    # -- Convolution --------------------------------------------------------

    def convolve(self, gv):
        """FK convolution: DIS (einsum) or hadronic (matmul).

        Returns predictions with shape (ndata, nrep).
        """
        if self.hadronic:
            return self._convolve_hadronic(gv)
        return self._convolve_dis(gv)

    def _convolve_dis(self, gv):
        return np.einsum('dfx,rfx->dr', self._W_dis, gv)

    def _convolve_hadronic(self, gv):
        nrep = gv.shape[0]
        gv1 = gv[:, self._fl1, :]
        gv2 = gv[:, self._fl2, :]
        lumi_flat = (
            gv1[:, :, :, np.newaxis] * gv2[:, :, np.newaxis, :]
        ).reshape(nrep, -1)
        return self._W_had_flat @ lumi_flat.T

    # -- Rotation -----------------------------------------------------------

    def rotate_to_evolution(self, gv_flav):
        """Rotate physical flavor-basis to the evolution subset.

        Input:  (nrep, 14, nx) in PDG flavour order.
        Output: (nrep, nfl_used, nx) in the evolution columns this FK
                table actually needs.
        """
        gv_sub = gv_flav[:, self.fl_col_inds, :]
        return np.einsum('ef,rfx->rex', self.transform_mat, gv_sub)


# ---------------------------------------------------------------------------
# Dataset container with multi-FK support
# ---------------------------------------------------------------------------

class NNPDFObservable:
    """Container for one dataset: FK tables, data, covariance, chi2.

    Supports multi-FK datasets with operations: NULL, ASY, RATIO, ADD, COM, SMN, SMT.
    """

    def __init__(self, name, fk_entries, operation, data_central,
                 sqrt_covmat, ndata):
        self.name = name
        self.fk_entries = fk_entries
        self.operation = operation
        self.data_central = data_central
        self.sqrt_covmat = sqrt_covmat
        self.ndata = ndata
        self._precompute_inv_covmat()

    # -- Properties (delegate to first FK entry) ----------------------------

    @property
    def xgrid(self):
        return self.fk_entries[0].xgrid

    @property
    def Q0(self):
        return self.fk_entries[0].Q0

    @property
    def hadronic(self):
        return self.fk_entries[0].hadronic

    @property
    def flavor_indices(self):
        return self.fk_entries[0].flavor_indices

    @property
    def fl_col_inds(self):
        return self.fk_entries[0].fl_col_inds

    @property
    def transform_mat(self):
        return self.fk_entries[0].transform_mat

    @property
    def n_fk(self):
        return len(self.fk_entries)

    def _precompute_inv_covmat(self):
        """Precompute L_inv = inv(cholesky(covmat)) once."""
        n = self.sqrt_covmat.shape[0]
        self._L_inv = solve_triangular(
            self.sqrt_covmat, np.eye(n), lower=True, check_finite=False
        )

    # -- Combining predictions from multiple FK tables ----------------------

    def _combine_predictions(self, pred_list):
        """Apply the FK operation to combine predictions.

        Supports NULL, ASY, RATIO, ADD, COM, SMN, SMT operations
        matching validphys.convolution.OP.
        """
        if self.operation == 'NULL':
            return pred_list[0]
        elif self.operation == 'ASY':
            return (pred_list[0] - pred_list[1]) / (pred_list[0] + pred_list[1])
        elif self.operation == 'RATIO':
            return pred_list[0] / pred_list[1]
        elif self.operation == 'ADD':
            return pred_list[0] + pred_list[1]
        elif self.operation == 'SMN':
            # (a + b) / (c + d)
            return (pred_list[0] + pred_list[1]) / (pred_list[2] + pred_list[3])
        elif self.operation == 'COM':
            # (sum of first 10) / (sum of last 10)
            n = len(pred_list) // 2
            num = sum(pred_list[:n])
            den = sum(pred_list[n:])
            return num / den
        elif self.operation == 'SMT':
            # sum of all (up to 10)
            return sum(pred_list)
        else:
            raise ValueError(f"Unknown operation: {self.operation}")

    # -- FK convolution -----------------------------------------------------

    def convolve(self, gv_list):
        """FK convolution with operation combination.

        Parameters
        ----------
        gv_list : np.ndarray or list of np.ndarray
        """
        if not isinstance(gv_list, list):
            gv_list = [gv_list]
        pred_list = [
            entry.convolve(gv)
            for entry, gv in zip(self.fk_entries, gv_list)
        ]
        return self._combine_predictions(pred_list)

    # -- Rotation -----------------------------------------------------------

    def rotate_to_evolution(self, gv_flav_list):
        """Rotate flavor-basis grid values to evolution for each FK entry."""
        if not isinstance(gv_flav_list, list):
            gv_flav_list = [gv_flav_list]
        return [
            entry.rotate_to_evolution(gv)
            for entry, gv in zip(self.fk_entries, gv_flav_list)
        ]

    def chi2_from_flavor(self, gv_flav_list):
        """Chi2 from flavor-basis grid values: rotate, convolve, compare."""
        gv_evol_list = self.rotate_to_evolution(gv_flav_list)
        return self.chi2(gv_evol_list)

    # -- Chi2 ---------------------------------------------------------------

    def chi2(self, gv_list):
        """Chi2 per replica using precomputed L_inv.

        Returns shape (nrep,).
        """
        preds = self.convolve(gv_list)
        diffs = preds - self.data_central.reshape(-1, 1)
        vec = self._L_inv @ diffs
        return np.einsum('ir,ir->r', vec, vec)


# ---------------------------------------------------------------------------
# Main setup function
# ---------------------------------------------------------------------------

def _parse_dataset_entry(entry):
    """Parse a dataset entry from the runcard. Useful for handling the variant of the dataset.

    Accepts either a plain string or a dict with keys:
        dataset (str), variant (str, optional), cfac (list, optional)

    Returns (name, variant, cfac) tuple.
    """
    if isinstance(entry, str):
        return entry, None, ()
    elif isinstance(entry, dict):
        name = entry["dataset"]
        variant = entry.get("variant", None)
        cfac = tuple(entry.get("cfac", []))
        return name, variant, cfac
    else:
        raise TypeError(f"Dataset entry must be str or dict, got {type(entry)}")


def setup_observables(
    pdf_name: str,
    datasets: list,
    theory_id: int = 708,
    use_cuts: str = "internal",
    variant: Optional[str] = None,
):
    """Load datasets and build NNPDFObservable containers for SV analysis.

    Parameters
    ----------
    pdf_name : str
        LHAPDF set name, e.g. 'NNPDF40_nnlo_as_01180'.
    datasets : list
        Dataset entries. Each can be a plain string (dataset name) or a dict
        with keys 'dataset', 'variant' (optional), 'cfac' (optional).
    theory_id : int
        NNPDF theory ID (default 708).
    use_cuts : str
        Cut policy: 'internal', 'nocuts', or 'fromfit'.
    variant : str or None
        Global uncertainty variant fallback, e.g. 'legacy'.
        Per-dataset variants override this.

    Returns
    -------
    pdf : validphys.core.PDF
    observables : list of NNPDFObservable
    flavor_info : dict
        Evolution-basis flavor metadata.
    flavor_basis_info : dict
        Physical (PDG) flavor metadata.
    """
    loader = Loader()
    theoryid = loader.check_theoryID(theory_id)
    pdf = loader.check_pdf(pdf_name)

    print(f"PDF set     : {pdf_name}")
    print(f"Theory ID   : {theory_id}")
    print(f"Cut policy  : {use_cuts}")
    print(f"Datasets    : {len(datasets)}\n")

    observables = []
    all_flavor_indices = set()

    for entry in datasets:
        ds_name, ds_variant, ds_cfac = _parse_dataset_entry(entry)
        # Per-dataset variant overrides global variant
        effective_variant = ds_variant if ds_variant is not None else variant
        ds_kw = {}
        if effective_variant is not None:
            ds_kw['variant'] = effective_variant
        if ds_cfac:
            ds_kw['cfac'] = ds_cfac
        ds = loader.check_dataset(ds_name, theoryid=theoryid, cuts=use_cuts, **ds_kw)

        # Load all FK tables for this dataset
        fk_entries = []
        for spec in ds.fkspecs:
            fk = load_fktable(spec)
            if ds.cuts is not None:
                fk = fk.with_cuts(ds.cuts)

            sigma = fk.sigma
            fm = sigma.columns

            if fk.hadronic:
                for col in fm:
                    i, j = divmod(col, NFK)
                    all_flavor_indices.add(i)
                    all_flavor_indices.add(j)
            else:
                all_flavor_indices.update(fm.tolist())

            entry = FKEntry(
                sigma=sigma,
                flavor_indices=fm,
                xgrid=fk.xgrid,
                Q0=fk.Q0,
                ndata=fk.ndata,
                hadronic=fk.hadronic,
            )
            fk_entries.append(entry)

        # Experimental data and covariance
        loaded_cd = ds.commondata.load()
        if ds.cuts is not None:
            loaded_cd = loaded_cd.with_cuts(ds.cuts.load())

        data_central = loaded_cd.central_values.to_numpy()
        stat_errors = loaded_cd.stat_errors.to_numpy()
        sys_errors = loaded_cd.systematic_errors()
        covmat = construct_covmat(stat_errors, sys_errors)
        L = cholesky(covmat, lower=True)

        obs = NNPDFObservable(
            name=ds_name,
            fk_entries=fk_entries,
            operation=ds.op,
            data_central=data_central,
            sqrt_covmat=L,
            ndata=loaded_cd.ndata,
        )
        observables.append(obs)

        # Print summary
        fk0 = fk_entries[0]
        tag = "hadronic" if fk0.hadronic else "DIS"
        op_str = f"  op={ds.op}" if ds.op != 'NULL' else ""
        nfk_str = f"  #FK={len(fk_entries)}" if len(fk_entries) > 1 else ""
        if fk0.hadronic:
            print(f"  {ds_name:45s}  ndata={obs.ndata:4d}  "
                  f"nx={len(fk0.xgrid):3d}  [{tag}]{op_str}{nfk_str}")
        else:
            flavour_names = FK_FLAVOURS[fk0.flavor_indices].tolist()
            print(f"  {ds_name:45s}  ndata={obs.ndata:4d}  "
                  f"nx={len(fk0.xgrid):3d}  [{tag}]{op_str}{nfk_str}  "
                  f"fl={flavour_names}")

    # Evolution-basis flavor metadata
    sorted_fi = sorted(all_flavor_indices)
    flavor_info = {
        "indices": sorted_fi,
        "names": FK_FLAVOURS[sorted_fi].tolist(),
        "n_flavors": len(sorted_fi),
    }

    print(f"\nEvolution flavours used ({flavor_info['n_flavors']}):")
    for idx, name in zip(sorted_fi, flavor_info["names"]):
        print(f"  [{idx:2d}] {name}")

    # Physical flavor-basis metadata
    all_fl_col_inds = set()
    for obs in observables:
        for entry in obs.fk_entries:
            all_fl_col_inds.update(entry.fl_col_inds.tolist())
    sorted_fl = sorted(all_fl_col_inds)

    flavor_basis_info = {
        "indices": sorted_fl,
        "pdg_codes": [int(ALL_FLAVOURS[i]) for i in sorted_fl],
        "names": [FLAVOR_PDG_NAMES[int(ALL_FLAVOURS[i])] for i in sorted_fl],
        "n_flavors": len(sorted_fl),
    }

    print(f"\nPhysical flavours contributing ({flavor_basis_info['n_flavors']}):\n")
    for idx, name, pdg in zip(
        sorted_fl, flavor_basis_info["names"], flavor_basis_info["pdg_codes"]
    ):
        print(f"  [{idx:2d}] PDG {pdg:+3d}  {name}")

    return pdf, observables, flavor_info, flavor_basis_info


# ---------------------------------------------------------------------------
# Grid-value evaluation
# ---------------------------------------------------------------------------

def get_pdf_grid_values(pdf, target, n_replicas=None):
    """Evaluate PDF in evolution basis (FK-subset flavours only).

    Returns shape (nrep, nfl_used, nx).
    """
    gv_func = functools.partial(evol_basis.grid_values, pdf)
    flavour_names = FK_FLAVOURS[target.flavor_indices]
    gv = gv_func(
        qmat=[target.Q0], vmat=flavour_names, xmat=target.xgrid
    ).squeeze(-1)
    if n_replicas is not None:
        gv = gv[1: n_replicas + 1]
    return gv


def get_pdf_flavor_grid_values(pdf, target, n_replicas=None):
    """Evaluate PDF in physical flavor basis (all 14 PDG codes).

    Returns shape (nrep, 14, nx).
    """
    pdg_codes = list(ALL_FLAVOURS)
    gv = _lhapdf_grid_values(
        pdf, pdg_codes, target.xgrid, [target.Q0]
    ).squeeze(-1)
    if n_replicas is not None:
        gv = gv[1: n_replicas + 1]
    return gv


def get_pdf_grid_values_all14(pdf, target, n_replicas=None):
    """Evaluate PDF in all 14 evolution-basis flavours.

    Returns shape (nrep, 14, nx).
    """
    gv = evol_basis.grid_values(
        pdf, qmat=[target.Q0], vmat=FK_FLAVOURS, xmat=target.xgrid
    ).squeeze(-1)
    if n_replicas is not None:
        gv = gv[1: n_replicas + 1]
    return gv
