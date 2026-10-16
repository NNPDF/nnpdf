"""
Weight-correlation comparison across n3fit fits (Bayesian vs standard NNPDF).

For each runcard, loads ``<results>/<runcard>/nnfit/replica_<i>/weights.weights.h5``
for every replica, flattens the network weights into one vector per replica, and
computes ONE correlation matrix across replicas (shape p x p, p = #weights).

Outputs (in --out-dir):
    fig_heatmaps.png      layer-ordered rho heatmaps, one panel per runcard, shared colour scale
    fig_offdiag_hist.png  histogram of off-diagonal rho with the 1/sqrt(N) (and 1/sqrt(N_eff)) null overlaid
    fig_spectrum.png      eigenvalue spectrum + cumulative variance, Marchenko-Pastur edge marked
    fig_layer_blocks.png  mean |rho| per (layer, layer) block
    fig_diff_<a>_vs_<b>.png  rho_a - rho_b (only when the parameter layouts match)
    (--per-chain: correlations within BNN chains only; --last-layer: last layer only)
    fig_correlation_notation.png    notation, rules and the admissible weight correlations
    fig_cross_weight_over_replicas.png  rho^r_kk'  (pairs of weights, summed over replicas)
    fig_cross_weight_over_samples.png   rho^n_kk'  (pairs of weights, over the samples of a training)
    fig_cross_weight_over_trainings.png rho^m_kk'  (pairs of weights, over the training means)
    fig_cross_replica_over_weights.png  rho^k_rr'  (pairs of replicas, summed over weights)
    fig_cross_sample_over_weights.png   rho^k_nn'  (pairs of samples of a training, over weights)
    fig_cross_training_over_weights.png rho^k_mm'  (pairs of trainings, over weights)
    fig_mode_separation.png             T_k, distance of the training modes in their own widths
                          (see plt_corr_weight_correlations.tex for the notation and the rules)
    scalar_table.csv/.md  single-number summaries (also printed)
    <runcard>_W.npy       cached (N, p) weight matrix

Weights that never vary across replicas (e.g. the deterministic hidden layers of a
Bayesian fit with a single BNN model, where every pseudo-replica shares them) have zero
variance, so their correlation is undefined. They are drawn grey and excluded from all
statistics; the table reports p_var (number of varying weights) next to p.

Parameter layout is the same for dense and VBDense layers: kernel as (in, out) flattened,
then bias. For VBDense the replica weight is the frozen eval draw used by n3fit,
w = mu_w + sqrt(exp(logsig2_w)) * random_w (and likewise for the bias if bayesian_bias),
or just mu_w with --which mean. Layers are ordered dense layers first, then VBDense layers.

Examples:
    python plt_corr.py bay_1x1000 bay_10x100 bay_100x10 nnpdf40-like
    python plt_corr.py bay_100x10 nnpdf40-like --n-replicas 500 --reference nnpdf40-like
"""

import argparse
import glob
import math
import os
import re
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DEFAULT_RESULTS = os.path.join(
    os.path.expanduser("~"), "miniconda3/envs/environment_nnpdf/share/NNPDF/results"
)
LAYER_RE = re.compile(r"(?P<layer>(?:correlated_low_rank_)?(?:vb_)?dense(?:_(?P<idx>\d+))?)/vars/(?P<var>\d+)$")
LOGSIG_MIN, LOGSIG_MAX = -20.0, 11.0  # float32 clip bounds used by VBDense.s2_w/s2_b (and CorrelatedLowRankVBDense.s2_w/s2_b)
COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]


# --------------------------------------------------------------------------- loading
def _layer_vars(h5):
    """Return {layer_name: {var_index: dataset_path}} for dense / vb_dense layers (no flows)."""
    found, has_flow = {}, set()

    def visit(name, obj):
        if not hasattr(obj, "shape") or "preprocessing" in name:
            return
        m = LAYER_RE.search(name)
        if m is None:
            return
        if "/flow/" in name:
            has_flow.add(m.group("layer"))
            return
        # keep the full path prefix so two layers with the same short name can't collide
        key = (name[: m.start()], m.group("layer"), int(m.group("idx") or 0))
        found.setdefault(key, {})[int(m.group("var"))] = name

    h5.visititems(visit)
    # dense layers first (by suffix), then VB layers (plain or correlated-low-rank)
    ordered = sorted(found, key=lambda k: (k[1].endswith("vb_dense"), k[2]))
    return [(k[1], found[k], k[1] in has_flow) for k in ordered]


def _correlated_weight(get, n, which):
    """CorrelatedLowRankVBDense eval-time frozen draw: w = mu + sqrt(D)*random [+ U^T @ random_eps],
    b = bias [+ sqrt(D_b)*random_b] [+ U_b^T @ random_eps]. Variable order follows build() in
    base_layers.py exactly (bias, mu_w, logsig2_w, [u_w], [bias_logsig2], [u_b], random,
    [random_eps], [random_b]); n=4/6/6/9 correspond to the 4 (rank, bayesian_bias) combinations,
    the two n=6 cases disambiguated by the ndim of var 3 (u_w is rank-3, bias_logsig2 is rank-1)."""
    bias, mu, logsig2 = get(0), get(1), get(2)
    u_w = u_b = bias_logsig2 = random_eps = random_b = None
    if n == 4:
        random = get(3)
    elif n == 6 and get(3).ndim == 3:  # rank > 0, bayesian_bias=False
        u_w, random, random_eps = get(3), get(4), get(5)
    elif n == 6:  # rank == 0, bayesian_bias=True
        bias_logsig2, random, random_b = get(3), get(4), get(5)
    elif n == 9:  # rank > 0, bayesian_bias=True
        u_w, bias_logsig2, u_b = get(3), get(4), get(5)
        random, random_eps, random_b = get(6), get(7), get(8)
    else:
        raise ValueError(f"unexpected CorrelatedLowRankVBDense variable count {n}")

    if which == "mean":
        return mu.T, bias

    w = mu + np.sqrt(np.exp(np.clip(logsig2, LOGSIG_MIN, LOGSIG_MAX))) * random
    b = bias
    if bias_logsig2 is not None:
        s2b = np.exp(np.clip(bias_logsig2, LOGSIG_MIN, LOGSIG_MAX))
        b = bias + np.sqrt(s2b) * random_b
    if u_w is not None:
        w = w + np.einsum("r,roi->oi", random_eps, u_w)
        if u_b is not None:
            b = b + np.einsum("r,ro->o", random_eps, u_b)
    return w.T, b


def _flatten_replica(path, which, include_preproc):
    """One replica -> (vector, [(layer_label, n_params), ...], [bias_size per block])."""
    parts, blocks, biases = [], [], []
    with h5py.File(path, "r") as h5:
        for k, (name, v, flow) in enumerate(_layer_vars(h5)):
            get = lambda i: h5[v[i]][()].astype(np.float64)
            if name == "correlated_low_rank_vb_dense":
                w, b = _correlated_weight(get, len(v), which)
            elif not name.endswith("vb_dense"):
                w, b = get(0), get(1)  # kernel (in, out), bias (out,)
            else:
                n = len(v)
                if n not in (4, 6):
                    raise ValueError(f"{path}: unexpected VBDense variable count {n}")
                bias, mu, logsig2 = get(0), get(1), get(2)
                if which == "mean":
                    w, b = mu.T, bias
                else:
                    rand_w = get(4 if n == 6 else 3)
                    w = (mu + np.sqrt(np.exp(np.clip(logsig2, LOGSIG_MIN, LOGSIG_MAX))) * rand_w).T
                    b = bias
                    if n == 6:  # bayesian_bias: bias_logsig2 = var 3, random_b = var 5
                        s2b = np.exp(np.clip(get(3), LOGSIG_MIN, LOGSIG_MAX))
                        b = bias + np.sqrt(s2b) * get(5)
                if flow:
                    print(f"  warning: {name} has a normalizing flow; using the Gaussian part only")
            label = "corrVB" if name == "correlated_low_rank_vb_dense" else ("VB" if name.endswith("vb_dense") else "dense")
            parts += [w.ravel(), b.ravel()]
            blocks.append((f"L{k + 1} {label}", w.size + b.size))
            biases.append(b.size)
        if include_preproc:
            pre = [h5[n][()].ravel() for n in sorted(_preproc_names(h5), key=_natural)]
            if pre:
                parts.append(np.concatenate(pre).astype(np.float64))
                blocks.append(("preproc", parts[-1].size))
                biases.append(0)
    return np.concatenate(parts), blocks, biases


def _preproc_names(h5):
    out = []
    h5.visititems(lambda n, o: out.append(n) if hasattr(o, "shape") and "preprocessing/vars" in n else None)
    return out


def _natural(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def load_runcard(name, results, n_replicas, seed, which, include_preproc, out_dir, reload):
    """Return (W (N,p), blocks, biases) or None if the runcard has no saved weights."""
    files = glob.glob(os.path.join(results, name, "nnfit", "replica_*", "weights.weights.h5"))
    if not files:
        print(f"[{name}] no weights.weights.h5 found under {results}/{name}/nnfit/replica_*/ -- "
              "skipping (was the fit run with `save: weights.weights.h5`?)")
        return None
    files.sort(key=_natural)
    if n_replicas and n_replicas < len(files):
        rng = np.random.default_rng(seed)
        files = [files[i] for i in sorted(rng.choice(len(files), n_replicas, replace=False))]
    tag = f"{name}_{which}{'_pre' if include_preproc else ''}_n{len(files)}_s{seed}"
    cache = os.path.join(out_dir, f"{tag}_W.npy")
    if os.path.exists(cache) and not reload:
        d = np.load(cache, allow_pickle=True).item()
        if "biases" in d:  # older caches lack the bias sizes: re-read them
            print(f"[{name}] loaded cached matrix {d['W'].shape} from {cache}")
            return d["W"], d["blocks"], d["biases"]
    print(f"[{name}] reading {len(files)} replicas ...")
    rows, blocks, biases = [], None, None
    for f in files:
        v, b, bi = _flatten_replica(f, which, include_preproc)
        if blocks is None:
            blocks, biases = b, bi
        elif b != blocks:
            raise ValueError(f"{f}: layout {b} differs from first replica {blocks}")
        rows.append(v)
    W = np.stack(rows)
    np.save(cache, {"W": W, "blocks": blocks, "biases": biases}, allow_pickle=True)
    return W, blocks, biases


# ------------------------------------------------------------------- single replica
def load_one_replica(name, results, replica, which, include_preproc, last_layer):
    """One replica's weight vector -> (w, blocks) or None if that replica has no saved weights."""
    path = os.path.join(results, name, "nnfit", f"replica_{replica}", "weights.weights.h5")
    if not os.path.exists(path):
        print(f"[{name}] no replica_{replica}/weights.weights.h5 under {results}/{name}/nnfit/ -- skipping")
        return None
    w, blocks, biases = _flatten_replica(path, which, include_preproc)
    if last_layer:
        lo = w.shape[0] - blocks[-1][1]
        w, nb = w[lo:], biases[-1]
        blocks = [(f"{blocks[-1][0]} kernel", w.shape[0] - nb), (f"{blocks[-1][0]} bias", nb)]
    return w, blocks


def single_replica_matrix(w):
    """Standardised self-outer-product z_i*z_j of ONE replica's weight vector, z = (w-mean)/std.

    This is NOT a statistical correlation (no ensemble, no notion of sampling uncertainty): it is a
    deterministic, rank-1 structure map of a single trained/sampled network, answering "do large-
    magnitude weights in this one replica coincide in sign (red) or oppose (blue)". The diagonal
    (z_i^2) is not comparable to the off-diagonal and is not plotted; the true ensemble correlation
    (needs many replicas) is fig_heatmaps.png.
    """
    z = (w - w.mean()) / w.std()
    M = np.outer(z, z)
    np.fill_diagonal(M, np.nan)
    return M


def fig_single_replica(mats, blocks, replica, out, ncols=None):
    n = len(mats)
    fig, axs = _panel_grid(n, ncols, 5.2, 5.2)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#d9d9d9")
    vmax = np.nanpercentile(np.concatenate([np.abs(m[np.isfinite(m)]) for m in mats.values()]), 99)
    for ax, (name, M) in zip(axs, mats.items()):
        b = blocks[name]
        edges = np.cumsum([0] + [k for _, k in b])
        im = ax.imshow(M, cmap=cmap, vmin=-vmax, vmax=vmax, interpolation="nearest")
        for e in edges[1:-1]:
            ax.axhline(e - 0.5, c="k", lw=0.5)
            ax.axvline(e - 0.5, c="k", lw=0.5)
        mid = (edges[:-1] + edges[1:]) / 2
        ax.set_xticks(mid, [x for x, _ in b], rotation=45, ha="right", fontsize=7)
        ax.set_yticks(mid, [x for x, _ in b], fontsize=7)
        ax.set_title(f"{name}\nreplica {replica}, p={M.shape[0]}", fontsize=9)
    fig.suptitle(r"single-replica structure map $z_i z_j$, $z=(w-\bar{w})/\sigma_w$ "
                 "(NOT an across-replica correlation -- see docstring)", fontsize=8.5, y=1.02)
    fig.colorbar(im, ax=[a for a in axs if a.get_visible()], label=r"$z_i z_j$ (diagonal omitted; colour clipped at 99th pct)",
                shrink=0.8)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _panel_grid(n, ncols, width, height):
    """n panels in rows of ncols (default: one row); unused axes are hidden"""
    ncols = min(ncols or n, n)
    nrows = -(-n // ncols)
    fig, axs = plt.subplots(nrows, ncols, figsize=(width * ncols + 0.8, height * nrows), squeeze=False)
    flat = axs.ravel()
    for ax in flat[n:]:
        ax.set_visible(False)
    return fig, flat


# ------------------------------------------------------------------------ statistics
def correlation(W):
    """Full (p, p) correlation with NaN rows/cols for zero-variance weights; plus the varying mask."""
    std = W.std(axis=0)
    var = std > 1e-12 * max(np.abs(W).max(), 1e-30)
    rho = np.full((W.shape[1],) * 2, np.nan)
    idx = np.flatnonzero(var)
    rho[np.ix_(idx, idx)] = np.corrcoef(W[:, idx], rowvar=False)
    return rho, var


def chain_size(name, results):
    """Samples per BNN chain (n_bnn_samples from the fit's filter.yml); 1 for a standard fit."""
    import yaml

    path = os.path.join(results, name, "filter.yml")
    if not os.path.exists(path):
        return 1
    with open(path) as f:
        cfg = yaml.safe_load(f)
    params = cfg.get("parameters") or {}
    bayes = params.get("bayesian") or {}  # parameters::bayesian, or older flat runcards
    return int(bayes.get("n_bnn_samples", params.get("n_bnn_samples", cfg.get("n_bnn_samples", 1))))


def effective_n(W, blocks):
    """Per-weight effective sample size: number of distinct values the weight's layer block takes
    across replicas. Equals N when every replica has its own weights (standard fit, sampled VB
    layer); equals the number of BNN models for layers shared by all samples of a model."""
    edges = np.cumsum([0] + [n for _, n in blocks])
    neff = np.empty(W.shape[1])
    for lo, hi in zip(edges[:-1], edges[1:]):
        neff[lo:hi] = len(np.unique(W[:, lo:hi], axis=0))
    return neff


def spectrum_stats(rho_v, n, neff_v=None):
    """Metrics of the (p_var x p_var) correlation matrix of the varying weights.
    neff_v: per-weight effective sample size, used for the noise-floor of |rho|."""
    p = rho_v.shape[0]
    if p < 2:
        return dict(p_var=p)
    lam = np.clip(np.linalg.eigvalsh(rho_v)[::-1], 0, None)
    q = lam / lam.sum()
    q = q[q > 0]
    off = rho_v[~np.eye(p, dtype=bool)]
    # E|r| for independent Gaussian variables is sqrt(2/(pi n)); for a pair the smaller effective
    # sample size limits the estimate
    neff_v = np.full(p, float(n)) if neff_v is None else neff_v
    pair_neff = np.minimum.outer(neff_v, neff_v)
    null_mean_abs = np.sqrt(2 / (np.pi * pair_neff))[~np.eye(p, dtype=bool)].mean()
    mp_edge = (1 + np.sqrt(p / n)) ** 2
    return dict(
        p_var=p,
        mean_abs_rho=np.abs(off).mean(),
        null_mean_abs_rho=null_mean_abs,
        excess_mean_abs_rho=np.abs(off).mean() - null_mean_abs,
        n_eff_min=neff_v.min(),
        max_abs_rho=np.abs(off).max(),
        frob_dev_over_p=np.linalg.norm(rho_v - np.eye(p)) / p,
        participation_ratio_over_p=lam.sum() ** 2 / (lam**2).sum() / p,
        entropy_eff_rank_over_p=np.exp(-(q * np.log(q)).sum()) / p,
        frac_eig_above_MP=(lam > mp_edge).sum() / p,
        top_eig_share=lam[0] / lam.sum(),
    )


def block_matrix(rho, blocks):
    edges = np.cumsum([0] + [n for _, n in blocks])
    L = len(blocks)
    M = np.full((L, L), np.nan)
    for i in range(L):
        for j in range(L):
            sub = np.abs(rho[edges[i] : edges[i + 1], edges[j] : edges[j + 1]]).copy()
            if i == j:
                np.fill_diagonal(sub, np.nan)
            if np.isfinite(sub).any():
                M[i, j] = np.nanmean(sub)
    return M, edges


# --------------------------------------------------------------------------- plotting
def fig_heatmaps(data, out, ncols=None):
    n = len(data)
    fig, axs = _panel_grid(n, ncols, 5.2, 5.2)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#d9d9d9")
    for ax, (name, d) in zip(axs, data.items()):
        im = ax.imshow(d["rho"], cmap=cmap, vmin=-1, vmax=1, interpolation="nearest")
        for e in d["edges"][1:-1]:
            ax.axhline(e - 0.5, c="k", lw=0.5)
            ax.axvline(e - 0.5, c="k", lw=0.5)
        mid = (d["edges"][:-1] + d["edges"][1:]) / 2
        ax.set_xticks(mid, [b for b, _ in d["blocks"]], rotation=45, ha="right", fontsize=7)
        ax.set_yticks(mid, [b for b, _ in d["blocks"]], fontsize=7)
        ax.set_title(f"{name}\nN={d['N']}, p_var={d['stats']['p_var']}/{d['rho'].shape[0]}, min N_eff={d['neff_min']}", fontsize=9)
    fig.colorbar(im, ax=[a for a in axs if a.get_visible()], label=r"$\rho$ (grey = constant across replicas)", shrink=0.8)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_hist(data, out):
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    bins = np.linspace(-1, 1, 101)
    for c, (name, d) in zip(COLORS, data.items()):
        v = d["rho_v"]
        off = v[np.triu_indices_from(v, 1)]
        ax.hist(off, bins=bins, density=True, histtype="step", color=c, lw=1.5, label=f"{name} (N={d['N']})")
        s = 1 / np.sqrt(d["N"])
        x = np.linspace(-1, 1, 400)
        ax.plot(x, np.exp(-0.5 * (x / s) ** 2) / (s * np.sqrt(2 * np.pi)), "--", c=c, lw=0.8)
        if d["neff_min"] < d["N"]:  # weights shared by many replicas: noisier, wider null
            s = 1 / np.sqrt(d["neff_min"])
            ax.plot(x, np.exp(-0.5 * (x / s) ** 2) / (s * np.sqrt(2 * np.pi)), ":", c=c, lw=1.2)
    ax.set_yscale("log")
    ax.set_xlabel(r"off-diagonal $\rho_{ij}$")
    ax.set_ylabel("density")
    ax.set_title("dashed: null $N(0,1/N)$; dotted: null with $N_{eff}$ = distinct values of shared layers", fontsize=8)
    ax.legend(fontsize=8)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_spectrum(data, out):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for c, (name, d) in zip(COLORS, data.items()):
        lam = np.clip(np.linalg.eigvalsh(d["rho_v"])[::-1], 1e-12, None)
        p = len(lam)
        a1.loglog(np.arange(1, p + 1), lam, c=c, label=name)
        a1.axhline((1 + np.sqrt(p / d["N"])) ** 2, c=c, ls="--", lw=0.8)
        a2.plot(np.arange(1, p + 1) / p, np.cumsum(lam) / lam.sum(), c=c, label=name)
    a1.set_xlabel("eigenvalue rank")
    a1.set_ylabel(r"eigenvalue of $\rho$")
    a1.set_title("dashed: Marchenko-Pastur upper edge", fontsize=9)
    a2.set_xlabel("fraction of modes")
    a2.set_ylabel("cumulative variance fraction")
    a1.legend(fontsize=8)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_blocks(data, out, ncols=None):
    n = len(data)
    fig, axs = _panel_grid(n, ncols, 3.6, 3.6)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#d9d9d9")
    vmax = max(np.nanmax(d["blockM"]) if np.isfinite(d["blockM"]).any() else 0 for d in data.values())
    for ax, (name, d) in zip(axs, data.items()):
        M = d["blockM"]
        im = ax.imshow(M, cmap=cmap, vmin=0, vmax=vmax)
        labels = [b for b, _ in d["blocks"]]
        ax.set_xticks(range(len(labels)), labels, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(labels)), labels, fontsize=7)
        for i in range(len(labels)):
            for j in range(len(labels)):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7, color="w")
        ax.set_title(name, fontsize=9)
    fig.colorbar(im, ax=[a for a in axs if a.get_visible()], label=r"mean $|\rho|$", shrink=0.8)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_diff(a, b, da, db, out):
    diff = da["rho"] - db["rho"]
    cmap = plt.get_cmap("PuOr_r").copy()
    cmap.set_bad("#d9d9d9")
    fig, ax = plt.subplots(figsize=(6, 5.2))
    im = ax.imshow(diff, cmap=cmap, vmin=-1, vmax=1, interpolation="nearest")
    for e in da["edges"][1:-1]:
        ax.axhline(e - 0.5, c="k", lw=0.5)
        ax.axvline(e - 0.5, c="k", lw=0.5)
    mid = (da["edges"][:-1] + da["edges"][1:]) / 2
    ax.set_xticks(mid, [x for x, _ in da["blocks"]], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(mid, [x for x, _ in da["blocks"]], fontsize=7)
    ax.set_title(rf"$\rho$({a}) - $\rho$({b})   (grey: undefined in either)", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------- the full set of weight correlations
# Notation (details in plt_corr_weight_correlations.tex): rho^c_ab = sum_c (w_a^d - mu_a^e)(w_b^d - mu_b^e) /
# (sigma_a^f sigma_b^f), labels m (training), n (sample), k (weight), r = (m, n) (replica); named
# 'cross-<a> correlation over <c>'. Admissible: R1 c disjoint from a, R2 d = c, R3 no n across m,
# R4 e = c, R5 f = e (1024 -> 6 combinations), plus the mode separation T_k.

def _tol(W):
    return 1e-12 * max(np.abs(W).max(), 1e-30)


def _pearson_cols(X):
    """Correlation between the columns of X (observations in rows), NaN for constant columns"""
    std = X.std(0)
    ok = std > _tol(X)
    out = np.full((X.shape[1],) * 2, np.nan)
    idx = np.flatnonzero(ok)
    if len(idx) > 1:
        out[np.ix_(idx, idx)] = np.corrcoef(X[:, idx], rowvar=False)
    return out


def _pearson_rows(X):
    """Correlation between the rows of X over its columns"""
    return np.corrcoef(X)


def _offdiag(M):
    v = M[~np.eye(M.shape[0], dtype=bool)]
    return v[np.isfinite(v)]


def _null(L):
    """Mean |rho| of the Pearson correlation of two independent Gaussian variables from L observations,
    exactly: Gamma((L-1)/2) / (sqrt(pi) Gamma(L/2)) (-> sqrt(2 / (pi L)) for large L; 0.273 vs 0.252
    at L = 10)"""
    if L < 3:
        return np.nan
    return math.exp(math.lgamma((L - 1) / 2) - math.lgamma(L / 2)) / math.sqrt(math.pi)


def _null_avg(L, n):
    """Mean |rho| of the average of n independent such correlations (each of variance 1/(L-1)):
    sqrt(2 / (pi n (L-1)))"""
    return np.sqrt(2 / (np.pi * n * (L - 1))) if L > 1 and n > 0 else np.nan


def _trainings(W, K):
    """Training view of W (N, p): weights (N_t, K, p), training means m, within-training spreads s and the
    mask of the weights that vary within every training"""
    n_t = W.shape[0] // K
    Wt = W[: n_t * K].reshape(n_t, K, W.shape[1])
    m, s = Wt.mean(1), Wt.std(1)
    return Wt, m, s, (s > _tol(W)).all(axis=0)


def cross_weight_over_replicas(W):
    """rho^r_kk': Pearson correlation of every pair of weights over all replicas (the total correlation,
    a mixture of the within- and between-training parts)."""
    return _pearson_cols(W)


def cross_weight_over_samples(W, K, training=0, average=False):
    """rho^n_kk': correlation of every pair of weights over the K samples of one training (m fixed;
    default the first): the correlation inside one model's posterior. average=True: the correlation of
    every training, averaged over the trainings (the weights themselves are not averaged). Only weights
    that vary within the training(s) used (the Bayesian layer); NaN elsewhere."""
    Wt, _, _, within_all = _trainings(W, K)
    p = W.shape[1]
    out = np.full((p, p), np.nan)
    within = within_all if average else Wt[training].std(0) > _tol(W)
    if within.any():
        idx = np.flatnonzero(within)
        ms = range(Wt.shape[0]) if average else [training]
        out[np.ix_(idx, idx)] = np.mean([_pearson_cols(Wt[m][:, within]) for m in ms], axis=0)
    return out


def cross_weight_over_trainings(W, K, sample=None):
    """rho^m_kk': correlation of every pair of weights over the training means (n averaged first):
    how the mode locations of different weights move together across trainings. With sample=n, one
    sample per training (w_k^{m,n}) instead of the mean, i.e. with the posterior sampling noise."""
    Wt, m, _, _ = _trainings(W, K)
    return _pearson_cols(m if sample is None else Wt[:, sample])


def cross_replica_over_weights(W, K=None):
    """rho^k_rr': correlation of every pair of replicas over the weights. Returns a dict: 'raw (e=f=k)'
    compares the raw weight vectors; after a per-weight pre-standardisation, "pre-std e'=f'=r" uses the
    ensemble mean and spread (includes the training offset) and "pre-std e'=f'=n" the training mean
    and width (the posterior fluctuation only; needs K)."""
    var = W.std(0) > _tol(W)
    Wv = W[:, var]
    out = {"raw (e=f=k)": _pearson_rows(Wv), "pre-std e'=f'=r": _pearson_rows((Wv - Wv.mean(0)) / Wv.std(0))}
    if K and K > 1 and W.shape[0] % K == 0:
        Wt, m, s, within = _trainings(W, K)
        if within.any():
            resid = (Wt[:, :, within] - m[:, None, within]) / s[:, None, within]
            out["pre-std e'=f'=n"] = _pearson_rows(resid.reshape(W.shape[0], -1))
    return out


def cross_sample_over_weights(R, K, training=0, average=False):
    """rho^k_nn': correlation of every pair of samples of one training (m fixed; default the first)
    over the weights, i.e. that training's diagonal block of cross_replica_over_weights. average=True:
    the blocks of all trainings averaged entry by entry (sample labels are arbitrary, so only the level
    of the off-diagonal entries is meaningful). Returns, per variant, the block, its mean |off-diagonal|
    and its mean off-diagonal (signed)."""
    out = {}
    for v, M in R.items():
        ms = range(M.shape[0] // K) if average else [training]
        B = np.mean([M[m * K:(m + 1) * K, m * K:(m + 1) * K] for m in ms], axis=0)
        off = B[~np.eye(K, dtype=bool)]
        out[v] = (B, np.mean(np.abs(off)), np.mean(off))
    return out


def cross_training_over_weights(W, K, sample=None):
    """rho^k_mm': correlation of every pair of trainings over the weights, using the training means
    (n averaged first), or with sample=n one sample per training (w_k^{m,n}, sampling noise kept).
    Returns a dict: 'raw (e=f=k)', and "pre-std e'=m, f'=n" with u_k^m = (x_k^m - <x_k>_m) / s_k^m
    (offset from the consensus in units of the training's own width; x = mean or the sample)."""
    Wt, m, s, within = _trainings(W, K)
    x = m if sample is None else Wt[:, sample]
    out = {"raw (e=f=k)": _pearson_rows(x[:, x.std(0) > _tol(W)])}
    if within.any():
        u = (x[:, within] - x[:, within].mean(0)) / s[:, within]
        B = u @ u.T / u.shape[1]
        d = np.sqrt(np.diag(B))
        out["pre-std e'=m, f'=n"] = B / np.outer(d, d)
    return out


def mode_separation(W, K):
    """T_k = <((m_k^m - mu_k) / s_k^m)^2>_m per weight (NaN where a weight does not vary within the
    trainings) and its average over the weights per training: how far the trainings' modes sit from
    the consensus in units of their own posterior width (~1: one shared mode; >> 1: several)."""
    _, m, s, within = _trainings(W, K)
    T = np.full(W.shape[1], np.nan)
    u = (m[:, within] - m[:, within].mean(0)) / s[:, within]
    T[within] = (u**2).mean(0)
    return T, (u**2).mean(1)


def correlation_set(W, K, trainings_ok):
    """All weight correlations of one fit. W: (N, p) weights, K samples per training (replica r ->
    training r // K). Training-based quantities need trainings_ok (all replicas, in order). A standard
    fit (K = 1) is a set of trainings with one sample each: sets 2 and 5 and the mode separation are
    then undefined, set 3 equals set 1 and set 6 (raw) equals set 4 (raw)."""
    N, p = W.shape
    n_t = N // K
    out = dict(K=K, N=N, p=p, n_t=n_t, set1=cross_weight_over_replicas(W), set1_null=_null(N))
    if not trainings_ok or n_t * K != N:
        return out
    _, _, _, within = _trainings(W, K)
    out.update(within=within, p_within=int(within.sum()))
    if K > 1:
        out["set2"], out["set2_null"] = cross_weight_over_samples(W, K), _null(K)
        out["set2_avg"], out["set2_avg_null"] = cross_weight_over_samples(W, K, average=True), _null_avg(K, n_t)
    if n_t >= 3:
        out["set3"], out["set3_null"] = cross_weight_over_trainings(W, K), _null(n_t)
        # the same with one sample per training (the first): sampling noise kept (= set 3 for K = 1)
        out["set3_sample"], out["set3_sample_null"] = cross_weight_over_trainings(W, K, sample=0), _null(n_t)
    out["set4"] = cross_replica_over_weights(W, K)
    out["set4_null"] = _null(int((W.std(0) > _tol(W)).sum()))
    same = (np.arange(N)[:, None] // K) == (np.arange(N)[None, :] // K)
    off = ~np.eye(N, dtype=bool)
    if K > 1:
        out["set5"] = cross_sample_over_weights(out["set4"], K)
        out["set5_avg"] = cross_sample_over_weights(out["set4"], K, average=True)
    out["set4_within"] = {v: np.mean(np.abs(M[same & off])) if (same & off).any() else np.nan for v, M in out["set4"].items()}
    out["set4_between"] = {v: np.mean(np.abs(M[~same])) if (~same).any() else np.nan for v, M in out["set4"].items()}
    if n_t >= 2:
        out["set6"] = cross_training_over_weights(W, K)
        out["set6_sample"] = cross_training_over_weights(W, K, sample=0)
        _, m, _, _ = _trainings(W, K)
        out["set6_null"] = _null(int((m.std(0) > _tol(W)).sum()))
        if within.any():
            out["ratio_T"], out["ratio_Btt"] = mode_separation(W, K)
    return out


def correlation_set_summary(st):
    """Scalars for the table: mean |off-diagonal| of each set and its null"""
    row = {"set1_mean_abs": np.mean(np.abs(_offdiag(st["set1"]))), "set1_null": st["set1_null"]}
    if "set2" in st:
        row.update(set2_mean_abs=np.mean(np.abs(_offdiag(st["set2"]))), set2_null=st["set2_null"])
    if "set3" in st:
        row.update(set3_mean_abs=np.mean(np.abs(_offdiag(st["set3"]))), set3_null=st["set3_null"])
    if "set2_avg" in st:
        row.update(set2_avg_mean_abs=np.mean(np.abs(_offdiag(st["set2_avg"]))), set2_avg_null=st["set2_avg_null"])
    if "set3_sample" in st:
        row.update(set3_sample_mean_abs=np.mean(np.abs(_offdiag(st["set3_sample"]))))
    for key in st.get("set4", {}):
        tag = key.split()[0] if key.startswith("raw") else key.split()[1]
        row[f"set4_within_{tag}"] = st["set4_within"][key]
        if key in st.get("set5", {}):
            row[f"set5_first_training_{tag}"] = st["set5"][key][1]
            row[f"set5_avg_trainings_{tag}"] = st["set5_avg"][key][1]
        row[f"set4_between_{tag}"] = st["set4_between"][key]
    if "set4_null" in st:
        row["set4_null"] = st["set4_null"]
    for key, M in st.get("set6", {}).items():
        tag = key.split()[0] if key.startswith("raw") else "pre"
        row[f"set6_mean_{tag}"] = np.mean(_offdiag(M))
        if key in st.get("set6_sample", {}):
            row[f"set6_sample_mean_{tag}"] = np.mean(_offdiag(st["set6_sample"][key]))
    if "set6" in st:
        row["set6_independent"] = -1 / (st["n_t"] - 1)
    if "ratio_T" in st:
        row["ratio_T_median"] = np.nanmedian(st["ratio_T"])
    return row


# --------------------------------------------------------------------------- figures
CMAP = plt.get_cmap("RdBu_r").copy()
CMAP.set_bad("#d9d9d9")


def _caption(fig, notation, isolates, expect, y=-0.03):
    fig.text(0.5, y, f"{notation}\nIsolates: {isolates}\nHealthy mean-field BNN: {expect}",
             ha="center", va="top", fontsize=8, wrap=True)


def _grid(n_rows, n_cols, w=5.0, h=4.8):
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(w * n_cols + 1.0, h * n_rows), squeeze=False,
                            gridspec_kw={"wspace": 0.35, "hspace": 0.45})
    return fig, axs


def _weight_heat(ax, M, d, title):
    im = ax.imshow(M, cmap=CMAP, vmin=-1, vmax=1, interpolation="nearest")
    for e in d["edges"][1:-1]:
        ax.axhline(e - 0.5, c="k", lw=0.5)
        ax.axvline(e - 0.5, c="k", lw=0.5)
    mid = (d["edges"][:-1] + d["edges"][1:]) / 2
    ax.set_xticks(mid, [b for b, _ in d["blocks"]], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(mid, [b for b, _ in d["blocks"]], fontsize=7)
    ax.set_title(title, fontsize=8.5)
    return im


def _lines(ax, K, N):
    if 1 < K < N and N // K <= 30:
        for e in np.arange(K, N, K):
            ax.axhline(e - 0.5, c="k", lw=0.3)
            ax.axvline(e - 0.5, c="k", lw=0.3)


def _why(st, key):
    """why a correlation is missing for a fit"""
    if "within" not in st:
        return "not computed\n(replicas incomplete or out of order)"
    if st["K"] == 1 and key.split("_")[0] in ("set2", "set5", "ratio"):
        return "not defined:\none network per training (K = 1)"
    if st["n_t"] < 3:
        return "not defined:\none training (N_t = 1)" if st["n_t"] == 1 else f"not defined:\n{st['n_t']} trainings"
    return "not defined for this fit"


def _undefined(ax, name=None, reason=None):
    """a cell of a figure grid for which the correlation does not exist: left blank"""
    ax.set_axis_off()


def fig_weight_set(data, key, out, notation, isolates, expect, title):
    """sets 1-3 and their variants: one weight-weight heatmap per fit (all fits; undefined ones say why)"""
    if not any(key in d["set"] for d in data.values()):
        print(f"{os.path.basename(out)}: no fit qualifies; skipped")
        return
    items = {n: d for n, d in data.items() if key in d["set"]}
    fig, axs = _grid(1, len(items))
    for ax, (name, d) in zip(axs[0], items.items()):
        st = d["set"]
        M = st[key]
        im = _weight_heat(ax, M, d, f"{name}\nmean |off-diag| {np.mean(np.abs(_offdiag(M))):.3f} "
                                    f"(null {st[key + '_null']:.3f})")
    fig.colorbar(im, ax=list(axs[0]), shrink=0.8, label=r"$\rho_{kk'}$ (grey: undefined)")
    fig.suptitle(title, fontsize=10, y=1.03)
    _caption(fig, notation, isolates, expect)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_set4(data, out, shuffle_seed=None):
    """shuffle_seed: show the replicas in a random order (the same permutation for rows and columns; the
    statistics in the titles do not depend on it), so that the training blocks are not lined up"""
    items = {n: d for n, d in data.items() if "set4" in d["set"]}
    if not items:
        print(f"{os.path.basename(out)}: no BNN fit with trainings; skipped")
        return
    variants = max((list(d["set"]["set4"]) for d in items.values()), key=len)
    fig, axs = _grid(len(variants), len(items))
    for j, (name, d) in enumerate(items.items()):
        st = d["set"]
        for i, v in enumerate(variants):
            ax = axs[i, j]
            if v not in st["set4"]:
                _undefined(ax, f"{name}  [{v}]", "not defined:\none network per training (K = 1)")
                continue
            M = st["set4"][v]
            if shuffle_seed is None:
                im = ax.imshow(M, cmap=CMAP, vmin=-1, vmax=1, interpolation="nearest")
                _lines(ax, st["K"], st["N"])
            else:
                order = np.random.default_rng(shuffle_seed).permutation(st["N"])
                im = ax.imshow(M[np.ix_(order, order)], cmap=CMAP, vmin=-1, vmax=1, interpolation="nearest")
            w, b = st["set4_within"][v], st["set4_between"][v]
            parts = ([f"within training {w:.3f}"] if np.isfinite(w) else []) + ([f"between {b:.3f}"] if np.isfinite(b) else [])
            ax.set_title(f"{name}  [{v}]\nmean |R| " + ", ".join(parts), fontsize=8.5)
            order_txt = "" if shuffle_seed is None else " (random order)"
            ax.set_xlabel(f"replica r'{order_txt}")
            ax.set_ylabel(f"replica r{order_txt}")
    fig.colorbar(im, ax=list(axs.ravel()), shrink=0.6, label=r"$\rho_{rr'}$")
    shuffled = "" if shuffle_seed is None else f", replicas in a random order (seed {shuffle_seed})"
    fig.suptitle(r"cross-replica correlation over weights, $\rho^{k}_{rr'}$" + shuffled, fontsize=10, y=1.0)
    _caption(fig, "a,b = r;  c = k;  d = k;  e = f = k (raw).  Pre-standardised rows: z_k = (w_k - mu_k^e')/sigma_k^f' first.",
             "raw: similarity of two networks' weight vectors (large-spread weights dominate); e'=r: what two replicas "
             "share relative to the ensemble, incl. their training's offset; e'=n: the posterior fluctuations only",
             "e'=n: -1/(K-1) within a training (K residuals around their own mean sum to 0), ~0 between trainings; "
             "|.| at the null sqrt(2/(pi p)) once that is larger", y=-0.01)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


def fig_set5(data, out, key="set5"):
    """set 5: the K x K sample-sample block of the first training (key 'set5') or the blocks of all
    trainings averaged entry by entry ('set5_avg'); all fits, undefined ones say why"""
    avg = key.endswith("_avg")
    items = {n: d for n, d in data.items() if key in d["set"]}
    if not items:
        print(f"{os.path.basename(out)}: no BNN fit with trainings; skipped")
        return
    variants = max((list(d["set"][key]) for d in items.values()), key=len)
    fig, axs = _grid(len(variants), len(items), 5.0, 4.6)
    for j, (name, d) in enumerate(items.items()):
        st = d["set"]
        for i, v in enumerate(variants):
            ax = axs[i, j]
            if key not in st or v not in st[key]:
                _undefined(ax, f"{name}\n[{v}]", _why(st, key))
                continue
            block, mean_abs, mean = st[key][v]
            im = ax.imshow(block, cmap=CMAP, vmin=-1, vmax=1, interpolation="nearest")
            which = ("the only training" if st["n_t"] == 1 else f"average of the {st['n_t']} trainings") if avg else "first training"
            ax.set_title(f"{name}\n[{v}] {which} ({st['K']} samples)\n"
                         f"mean |off-diag| {mean_abs:.3f}, mean off-diag {mean:.3f}", fontsize=8)
            ax.set_xlabel("sample n'")
            ax.set_ylabel("sample n")
    fig.colorbar(im, ax=list(axs.ravel()), shrink=0.6, label=r"$\rho_{nn'}$")
    if avg:
        fig.suptitle(r"cross-sample correlation over weights, $\rho^{k}_{nn'}$, averaged over the trainings "
                     "(entry by entry)", fontsize=10, y=1.0)
        notation = ("a,b = n;  c = k;  one correlation per training, then averaged over m (the weights are not averaged);  "
                    "d = k;  e = f = k (raw), pre-standardised rows as in set 4.  Sample labels are arbitrary: only the "
                    "level of the off-diagonal entries is meaningful, not their pattern.")
    else:
        fig.suptitle(r"cross-sample correlation over weights of the first training, $\rho^{k}_{nn'}$", fontsize=10, y=1.0)
        notation = "a,b = n;  c = k;  m fixed (the first training);  d = k;  e = f = k (raw), pre-standardised rows as in set 4."
    p_w = next(iter(items.values()))["set"].get("p_within", 0)
    _caption(fig, notation,
             "whether the posterior draws of one model are independent of each other (e'=n), or how much they share "
             "the training's offset (e'=r)",
             "e'=n: -1/(K-1) (K residuals around their own mean sum to 0) plus noise of size sqrt(2/(pi p)) "
             f"(= {_null(p_w):.3f} here{', / sqrt(N_t) after averaging' if avg else ''}); e'=r: ~1 when the modes are far apart",
             y=-0.01)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


def fig_set6(data, out, key="set6"):
    """set 6 from the training means (key 'set6') or from one sample per training ('set6_sample'); all
    fits, undefined ones say why"""
    one = key.endswith("_sample")
    items = {n: d for n, d in data.items() if key in d["set"]}
    if not items:
        print(f"{os.path.basename(out)}: no fit with >= 2 trainings; skipped")
        return
    variants = max((list(d["set"][key]) for d in items.values()), key=len)
    fig, axs = _grid(len(variants), len(items), 4.8, 4.4)
    for j, (name, d) in enumerate(items.items()):
        st = d["set"]
        for i, v in enumerate(variants):
            ax = axs[i, j]
            if key not in st or v not in st[key]:
                reason = _why(st, key) if key not in st else "not defined:\none network per training (K = 1)"
                _undefined(ax, f"{name}  [{v}]", reason)
                continue
            M = st[key][v].copy()
            np.fill_diagonal(M, np.nan)
            im = ax.imshow(M, cmap=CMAP, vmin=-1, vmax=1, interpolation="nearest")
            ref = f" (independent {-1 / (st['n_t'] - 1):.3f})" if v.startswith("pre") else ""
            ax.set_title(f"{name}  [{v}]\n{st['n_t']} trainings: mean off-diag {np.nanmean(M):.3f}{ref}", fontsize=8.5)
            ax.set_xlabel("training m'")
            ax.set_ylabel("training m")
    fig.colorbar(im, ax=list(axs.ravel()), shrink=0.6, label=r"$\rho_{mm'}$ (diagonal not shown)")
    what = "one sample per training ($w_k^{m,n=0}$)" if one else "training means ($n$ averaged)"
    fig.suptitle(r"cross-training correlation over weights, $\rho^{k}_{mm'}$, " + what, fontsize=10, y=1.0)
    x = "w_k^{m,n=0} (one sample, no averaging over n)" if one else "m_k^m = <w_k^{m,n}>_n (n averaged first)"
    _caption(fig, f"a,b = m;  c = k;  x_k^m = {x};  d = k;  e = f = k (raw).  Pre-std: u_k^m = (x_k^m - <x_k>_m)/s_k^m "
             "(consensus over trainings, own width), then normalised over k.",
             "whether different trainings land in related modes (deviate from the consensus in the same direction)"
             + ("; with the posterior sampling noise kept" if one else ""),
             "pre-std: about -1/(N_t - 1) -> 0 for independent trainings (the offsets sum to 0)", y=-0.01)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_ratio(data, out):
    items = {n: d for n, d in data.items() if "ratio_T" in d["set"]}
    if not items:
        print(f"{os.path.basename(out)}: no fit with >= 2 trainings; skipped")
        return
    fig, axs = _grid(2, len(items), 5.4, 3.4)
    for j, (name, d) in enumerate(items.items()):
        st = d["set"]
        T, Btt = st["ratio_T"], st["ratio_Btt"]
        ax = axs[0, j]
        ok = np.isfinite(T)
        ax.plot(np.flatnonzero(ok), T[ok], ".", ms=3, color="#1f77b4")
        for e in d["edges"][1:-1]:
            ax.axvline(e - 0.5, c="grey", lw=0.5)
        ax.axhline(1, c="k", ls="--", lw=0.8)
        ax.set_yscale("log")
        mid = (d["edges"][:-1] + d["edges"][1:]) / 2
        ax.set_xticks(mid, [b for b, _ in d["blocks"]], rotation=30, ha="right", fontsize=7)
        ax.set_ylabel(r"$T_k$")
        ax.set_title(f"{name}: per weight, median {np.nanmedian(T):.3g}", fontsize=8.5)
        ax = axs[1, j]
        ax.plot(np.arange(len(Btt)), Btt, "o", ms=3, color="#1f77b4")
        ax.axhline(1, c="k", ls="--", lw=0.8)
        ax.set_yscale("log")
        ax.set_xlabel("training m")
        ax.set_ylabel(r"$\langle (u_k^{m})^2 \rangle_k$")
        ax.set_title(f"per training, median {np.median(Btt):.3g}", fontsize=8.5)
    fig.suptitle(r"mode separation $T_k$ (a variance ratio, not a correlation)", fontsize=10, y=1.02)
    _caption(fig, "a = b = k;  e = m (consensus over the training means);  f = n (spread within a training):  "
             r"T_k = <((m_k^m - mu_k) / s_k^m)^2>_m",
             "how far the trainings' modes sit from each other in units of each mode's own posterior width "
             "(cf. Gelman-Rubin)",
             "~1 if all trainings sample one shared mode; >> 1 if they sit in different modes (dashed: 1)")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_correlation_notation(out):
    """The generic formula, the index labels, the rules and the sensible combinations, as a figure"""
    fig = plt.figure(figsize=(12.5, 9.5))
    fig.text(0.02, 0.97, "Weight correlations: notation", fontsize=14, weight="bold", va="top")
    fig.text(0.02, 0.92, r"$\rho_{ab} \;=\; \sum_c \dfrac{(w_a^d - \mu_a^e)(w_b^d - \mu_b^e)}{\sigma_a^f\,\sigma_b^f},"
             r"\qquad \mu_a^e = \langle w_a^d \rangle_e,\qquad \sigma_a^f = \sqrt{\langle (w_a^d - \mu_a^e)^2 \rangle_f}$",
             fontsize=15, va="top")
    lines = [
        r"Index labels:  $m$ = training (one BNN model),  $n$ = sample of a training,  $k$ = weight,  $r = (m, n)$ = replica.",
        r"$a, b$: the pair that is correlated;  $c$: the summed index (the observations);  $d$: the label $w$ is resolved along;",
        r"$e$, $f$: what the mean and the spread average over (all standard deviations are population ones, $1/N$).",
        "",
        r"Rules (each a condition on $(a, c, d, e, f)$; count of surviving combinations after it):",
        r"  start: every label from $\{m, n, k, r\}$:  $4^5 = 1024$",
        r"  R1  $c \cap a = \emptyset$ (no label is both paired and summed; $r = (m,n)$):  $\to 512$",
        r"  R2  $d = c$ ($w$ resolved along exactly the summed label; a third label is held fixed ($m$) or averaged first ($n$)):  $\to 128$",
        r"  R3  never compare samples $n$ across trainings $m$ (removes $(a,c) = (m,n), (n,m)$):  $\to 96$",
        r"  R4  $e = c$ (subtract the mean of exactly the summed set; else the training offset leaks in):  $\to 24$",
        r"  R5  $f = e$ (a correlation, $\rho_{aa} = 1$):  $\to 6$.  Of the 18 removed, $(k, m, m, m, n)$ is a meaningful variance ratio.",
    ]
    for i, l in enumerate(lines):
        fig.text(0.02, 0.84 - 0.033 * i, l, fontsize=10.5, va="top")
    rows = [
        ("set", "name", "symbol", "a,b", "c", "d", "e = f", "isolates", "healthy mean-field BNN"),
        ("1", "cross-weight over replicas", r"$\rho^{r}_{kk'}$", "k", "r", "r", "r", "total (mixture of 2 and 3)", "what the PDFs see"),
        ("2", "cross-weight over samples", r"$\rho^{n}_{kk'}$", "k", "n (m fixed)", "n", "n", "one model's posterior", "~0 (null)"),
        ("3", "cross-weight over trainings", r"$\rho^{m}_{kk'}$", "k", "m (n avg.)", "m", "m", "mode locations", "structure possible"),
        ("4", "cross-replica over weights", r"$\rho^{k}_{rr'}$", "r", "k", "k", "k", "two replicas' weights", "see pre-std."),
        ("5", "cross-sample over weights", r"$\rho^{k}_{nn'}$", "n", "k (m fixed)", "k", "k", "two samples of a training", "-1/(K-1) (e'=n)"),
        ("6", "cross-training over weights", r"$\rho^{k}_{mm'}$", "m", "k (n avg.)", "k", "k", "two trainings' modes", "~ -1/(N_t-1) -> 0"),
        ("ratio", "mode separation", r"$T_k$", "k", "m", "m", "e=m, f=n", "modes in own widths", "~1 one mode, >>1 several"),
    ]
    y0 = 0.42
    xs = [0.02, 0.07, 0.25, 0.32, 0.37, 0.47, 0.52, 0.59, 0.80]
    for i, row in enumerate(rows):
        for x, cell in zip(xs, row):
            fig.text(x, y0 - 0.04 * i, cell, fontsize=9.5, va="top", weight="bold" if i == 0 else "normal")
    fig.text(0.02, 0.105, r"Naming: every quantity is a weight correlation ($w$ enters $\rho$), named cross-$a$ correlation over $c$, "
             r"written $\rho^{c}_{ab}$ (summed index as superscript).", fontsize=9.5, va="top")
    fig.text(0.02, 0.06, "Sets 4-6 (cross-replica, cross-sample, cross-training over weights): mu, sigma with the lower index "
             "of w (e = f = k) are averages over the weights of one network,\nso they are also shown after a per-weight pre-standardisation z_k = (w_k - mu_k^e')/sigma_k^f': "
             "e'=f'=r (ensemble: includes the training offset), e'=f'=n (posterior fluctuation only), "
             "e'=m, f'=n (offset in own widths).", fontsize=9, va="top")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _pair_matrices(st, key_a, key_b, variant):
    a, b = st.get(key_a), st.get(key_b)
    if a is None or b is None:
        return None, None
    if variant is not None:
        if variant not in a or variant not in b:
            return None, None
        a, b = a[variant], b[variant]
    return a, b


def fig_variant_difference(data, rows, out, title, notation):
    """difference a - b of two variants of a correlation (e.g. training means - one extraction), per fit,
    on its own colour scale. rows: [(key_a, key_b, variant or None, label)]"""
    diffs = {}
    for name, d in data.items():
        for i, (ka, kb, v, _) in enumerate(rows):
            A, B = _pair_matrices(d["set"], ka, kb, v)
            if A is not None:
                D = A - B
                np.fill_diagonal(D, np.nan)
                diffs[(i, name)] = (D, A, B)
    if not diffs:
        print(f"{os.path.basename(out)}: nothing to compare; skipped")
        return
    vmax = {}
    for i in range(len(rows)):  # own colour scale per row
        vals = [D[np.isfinite(D)] for (r, _), (D, _, _) in diffs.items() if r == i and np.nanmax(np.abs(D)) > 0]
        vals = np.abs(np.concatenate(vals)) if vals else np.array([0.0])
        vmax[i] = np.percentile(vals, 99) if vals.max() > 0 else 1.0
    names = [n for n in data if any((i, n) in diffs and np.nanmax(np.abs(diffs[(i, n)][0])) > 0 for i in range(len(rows)))]
    if not names:
        print(f"{os.path.basename(out)}: no fit with a difference; skipped")
        return
    fig, axs = _grid(len(rows), len(names), 4.8, 4.4)
    ims = {}
    for j, name in enumerate(names):
        d = data[name]
        st = d["set"]
        for i, (ka, kb, v, label) in enumerate(rows):
            ax = axs[i, j]
            if (i, name) not in diffs:
                reason = _why(st, ka) if ka not in st else "not defined:\none network per training (K = 1)"
                _undefined(ax, f"{name}\n{label}", reason)
                continue
            D, A, B = diffs[(i, name)]
            ok = np.isfinite(D)
            if not ok.any() or np.nanmax(np.abs(D)) == 0:
                _undefined(ax, f"{name}\n{label}", "identical by construction\n(K = 1: one extraction = the mean)")
                continue
            ims[i] = ax.imshow(D, cmap=CMAP, vmin=-vmax[i], vmax=vmax[i], interpolation="nearest")
            a, b = A[ok], B[ok]
            ax.set_title(f"{name}\n{label}\nmean |a-b| {np.abs(D[ok]).mean():.3f}, max {np.abs(D[ok]).max():.2f}, "
                         f"corr(a,b) {np.corrcoef(a, b)[0, 1]:.3f}", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
    for i, im in ims.items():
        fig.colorbar(im, ax=list(axs[i]), shrink=0.8, label="a - b (diagonal not shown)")
    fig.suptitle(title, fontsize=10, y=1.0)
    fig.text(0.5, -0.01, notation, ha="center", va="top", fontsize=8, wrap=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_relative_to_standard(data, key, out, title):
    """relative difference (rho_BNN - rho_std) / rho_std of a weight-weight correlation, for every BNN fit
    against the standard fit(s) (K = 1) with the same weights"""
    refs = {n: d for n, d in data.items() if d["set"].get("K") == 1 and key in d["set"]}
    bnns = {n: d for n, d in data.items() if d["set"].get("K", 1) > 1 and key in d["set"]}
    pairs = [(b, r) for r in refs for b in bnns if bnns[b]["set"][key].shape == refs[r]["set"][key].shape]
    if not pairs:
        print(f"{os.path.basename(out)}: no BNN / standard pair; skipped")
        return
    from matplotlib.colors import SymLogNorm
    fig, axs = _grid(1, len(pairs), 5.0, 4.8)
    for ax, (b, r) in zip(axs[0], pairs):
        A, B = bnns[b]["set"][key], refs[r]["set"][key]
        rel = (A - B) / B
        np.fill_diagonal(rel, np.nan)
        ok = np.isfinite(rel)
        small = np.mean(np.abs(B[ok]) < refs[r]["set"][key + "_null"])
        im = ax.imshow(rel, cmap=CMAP, norm=SymLogNorm(linthresh=1, vmin=-100, vmax=100), interpolation="nearest")
        ax.set_title(f"({b} - {r}) / {r}\nmedian |rel| {np.median(np.abs(rel[ok])):.2f};  "
                     f"|rho_std| < its null in {100 * small:.0f}% of the entries", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(im, ax=list(axs[0]), shrink=0.8, label="relative difference (symmetric log scale)")
    fig.suptitle(title, fontsize=10, y=1.03)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def correlation_set_figures(data, out_dir):
    fig_correlation_notation(os.path.join(out_dir, "fig_correlation_notation.png"))
    fig_weight_set(
        data, "set1", os.path.join(out_dir, "fig_cross_weight_over_replicas.png"),
        "a,b = k;  c = r;  d = r;  e = f = r", "the total weight correlation over the ensemble: a mixture of sets 2 "
        "(within) and 3 (between), weighted by their variances", "set 2 when there is one training; dominated by set 3 "
        "when the modes are far apart", r"cross-weight correlation over replicas, $\rho^{r}_{kk'}$")
    fig_weight_set(
        data, "set2", os.path.join(out_dir, "fig_cross_weight_over_samples.png"),
        "a,b = k;  c = n, m fixed (the first training);  d = n;  e = f = n", "the correlation inside one model's "
        "posterior (weights varying within the training only)", "~0 (null shown) in a mean-field Bayesian layer: the "
        "weights are independent in q by construction", r"cross-weight correlation over samples of the first training, $\rho^{n}_{kk'}$")
    fig_weight_set(
        data, "set2_avg", os.path.join(out_dir, "fig_cross_weight_over_samples_avg_trainings.png"),
        "a,b = k;  c = n, one correlation per training, then averaged over m (the weights are not averaged);  d = n;  "
        "e = f = n", "the correlation inside the posteriors, with the noise reduced by averaging over the trainings",
        "~0 (null sqrt(2/(pi N_t (K-1))) shown) in a mean-field Bayesian layer",
        r"cross-weight correlation over samples, averaged over the trainings, $\langle\rho^{n}_{kk'}\rangle_m$")
    fig_weight_set(
        data, "set3", os.path.join(out_dir, "fig_cross_weight_over_trainings.png"),
        "a,b = k;  c = m (n averaged first: training means);  d = m;  e = f = m", "how the mode locations of different "
        "weights move together across trainings (degenerate directions, hidden-unit symmetries)",
        "no fixed expectation; compare with the null 1/sqrt(N_t)", r"cross-weight correlation over trainings, $\rho^{m}_{kk'}$")
    fig_weight_set(
        data, "set3_sample", os.path.join(out_dir, "fig_cross_weight_over_trainings_one_sample.png"),
        "a,b = k;  c = m (one sample per training, n = 0: no averaging over n);  d = m;  e = f = m", "set 3 with "
        "the posterior sampling noise kept: mode locations plus each training's within-posterior fluctuation",
        "the set-3 value, diluted by the share of the spread that is posterior noise, var_within / (var_within + var_between)",
        r"cross-weight correlation over trainings, one sample per training, $\rho^{m}_{kk'}$ ($w_k^{m,n=0}$)")
    fig_set4(data, os.path.join(out_dir, "fig_cross_replica_over_weights.png"))
    fig_set4(data, os.path.join(out_dir, "fig_cross_replica_over_weights_shuffled.png"), shuffle_seed=0)
    fig_set5(data, os.path.join(out_dir, "fig_cross_sample_over_weights.png"))
    fig_set5(data, os.path.join(out_dir, "fig_cross_sample_over_weights_avg_trainings.png"), key="set5_avg")
    fig_set6(data, os.path.join(out_dir, "fig_cross_training_over_weights.png"))
    fig_set6(data, os.path.join(out_dir, "fig_cross_training_over_weights_one_sample.png"), key="set6_sample")
    fig_ratio(data, os.path.join(out_dir, "fig_mode_separation.png"))
    fig_relative_to_standard(
        data, "set3_sample", os.path.join(out_dir, "fig_cross_weight_over_trainings_one_sample_rel_nnpdf.png"),
        r"relative difference of $\rho^{m}_{kk'}$ (one extraction per training) to the standard fit")
    fig_variant_difference(
        data, [("set3", "set3_sample", None, "set 3: means - one extraction")],
        os.path.join(out_dir, "fig_diff_set3_means_minus_one_sample.png"),
        r"$\rho^{m}_{kk'}$: training means (3a) minus one extraction (3b)",
        "a = correlation from the training means, b = from one extraction (n = 1) per training; own colour scale")
    fig_variant_difference(
        data, [("set6", "set6_sample", "raw (e=f=k)", "set 6 raw: means - one extraction"),
               ("set6", "set6_sample", "pre-std e'=m, f'=n", "set 6 pre-std: means - one extraction")],
        os.path.join(out_dir, "fig_diff_set6_means_minus_one_sample.png"),
        r"$\rho^{k}_{mm'}$: training means (6a) minus one extraction (6b)",
        "a = correlation from the training means, b = from one extraction (n = 1) per training; own colour scale")


# ------------------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runcards", nargs="+", help="fit names under the results dir, e.g. bay_1x1000 nnpdf40-like")
    ap.add_argument("--results-dir", default=DEFAULT_RESULTS)
    ap.add_argument("--out-dir", default="plt_corr_out")
    ap.add_argument("--n-replicas", type=int, default=None, help="random subset size (default: all replicas)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--which", choices=["sampled", "mean"], default="sampled",
                    help="VBDense weights: frozen eval draw (default) or posterior mean mu_w")
    ap.add_argument("--include-preproc", action="store_true", help="also include preprocessing exponents")
    ap.add_argument("--last-layer", action="store_true",
                    help="use only the last layer (kernel + bias). This is the like-for-like comparison: "
                         "every replica has its own last-layer weights in every fit, whereas the hidden "
                         "layers of a multi-sample BNN model are shared by all samples of that model")
    ap.add_argument("--single-replica", type=int, default=None, metavar="REPLICA",
                    help="ALSO make a single-replica structure-map heatmap (fig_single_replica_<REPLICA>.png) "
                         "using replica_<REPLICA> of each runcard (each fit's own indexing: BNN pseudo-replicas "
                         "are 0-indexed, a standard fit's replicas are 1-indexed). This is not an across-replica "
                         "correlation (undefined for one replica) but a deterministic single-network structure "
                         "map -- see single_replica_matrix()'s docstring. Independent of --n-replicas/--per-chain.")
    ap.add_argument("--per-chain", action="store_true",
                    help="correlate WITHIN BNN chains only: subtract each chain's mean over its samples "
                         "(n_bnn_samples from filter.yml; replica r belongs to chain r // n_bnn_samples) and "
                         "correlate the residuals. Weights shared by all samples of a chain become undefined "
                         "(grey). A standard fit has no chains (every replica is a separate network), so it is "
                         "skipped. Needs all replicas (no --n-replicas)")
    ap.add_argument("--reference", default=None, help="runcard to subtract in the difference heatmaps")
    ap.add_argument("--reload", action="store_true", help="ignore cached weight matrices")
    ap.add_argument("--ncols", type=int, default=None,
                    help="panels per row in the per-fit figures (heatmaps, blocks, single replica); "
                         "default: all in one row. E.g. 6 fits with --ncols 3 give two rows of three")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    data = {}
    for name in args.runcards:
        loaded = load_runcard(name, args.results_dir, args.n_replicas, args.seed, args.which,
                              args.include_preproc, args.out_dir, args.reload)
        if loaded is None:
            continue
        W, blocks, biases = loaded
        if args.last_layer:
            lo = W.shape[1] - blocks[-1][1]
            W = W[:, lo:]
            nb = biases[-1]
            blocks = [(f"{blocks[-1][0]} kernel", W.shape[1] - nb), (f"{blocks[-1][0]} bias", nb)]
        # the full set of weight correlations, on the weights as loaded (before --per-chain)
        K = chain_size(name, args.results_dir)
        corr_set = correlation_set(W, K, trainings_ok=not args.n_replicas)
        n_chains = 0
        if args.per_chain:
            K = chain_size(name, args.results_dir)
            if K < 2:
                print(f"[{name}] --per-chain: no BNN chains (chain size 1); a standard fit has one network per "
                      "replica, so within-chain weight correlation is undefined -- skipping")
                continue
            if args.n_replicas or W.shape[0] % K:
                sys.exit(f"[{name}] --per-chain needs all replicas (N={W.shape[0]} must be a multiple of {K})")
            n_chains = W.shape[0] // K
            grp = np.arange(W.shape[0]) // K
            means = np.stack([W[grp == g].mean(0) for g in range(n_chains)])
            W = W - means[grp]
            name = f"{name} (within {n_chains} chain{'s' * (n_chains > 1)} of {K})"
        n_dof = W.shape[0] - n_chains  # degrees of freedom left after removing chain means
        neff = np.minimum(effective_n(W, blocks), n_dof)
        rho, var = correlation(W)
        rho_v = rho[np.ix_(var, var)]
        blockM, edges = block_matrix(rho, blocks)
        stats = spectrum_stats(rho_v, n_dof, neff[var])
        if neff[var].min() < n_dof:
            k = int(neff[var].min())
            print(f"[{name}] note: some weights take only {k} distinct values over the {W.shape[0]} replicas "
                  f"(shared by all samples of a BNN model), so their correlations have effective N={k} "
                  f"(noise floor of |rho| ~ {np.sqrt(2 / (np.pi * k)):.2f}). Statistics over all weights are "
                  "dominated by this; use --last-layer for a like-for-like comparison.")
        # last layer alone: the only block that varies in every fit (shared hidden layers are constant
        # when there is a single BNN model)
        lo, hi = edges[-2], edges[-1]
        lv = var[lo:hi]
        last = {}
        if not args.last_layer and lv.sum() > 1:
            last = spectrum_stats(rho[lo:hi, lo:hi][np.ix_(lv, lv)], n_dof, neff[lo:hi][lv])
        stats.update({f"lastlayer_{k}": v for k, v in last.items() if k != "p_var"})
        stats.update(N=n_dof, p=W.shape[1], **correlation_set_summary(corr_set))
        data[name] = dict(neff_min=int(neff[var].min()), rho=rho, rho_v=rho_v, blocks=blocks, edges=edges, blockM=blockM, N=n_dof, stats=stats, set=corr_set)

    if args.single_replica is not None:
        mats, blk = {}, {}
        for name in args.runcards:
            got = load_one_replica(name, args.results_dir, args.single_replica, args.which,
                                    args.include_preproc, args.last_layer)
            if got is not None:
                w, b = got
                mats[name], blk[name] = single_replica_matrix(w), b
        if mats:
            fig_single_replica(mats, blk, args.single_replica,
                               os.path.join(args.out_dir, f"fig_single_replica_{args.single_replica}.png"),
                               args.ncols)
        else:
            print("--single-replica: no runcard had that replica; figure skipped")

    if not data:
        sys.exit("no runcard had saved weights; nothing to do")

    fig_heatmaps(data, os.path.join(args.out_dir, "fig_heatmaps.png"), args.ncols)
    fig_hist(data, os.path.join(args.out_dir, "fig_offdiag_hist.png"))
    fig_spectrum(data, os.path.join(args.out_dir, "fig_spectrum.png"))
    fig_blocks(data, os.path.join(args.out_dir, "fig_layer_blocks.png"), args.ncols)
    correlation_set_figures(data, args.out_dir)
    ref = args.reference
    if ref is not None and ref not in data:
        print(f"reference {ref} not available; skipping difference plots")
        ref = None
    if ref is not None:
        for name, d in data.items():
            if name != ref and d["rho"].shape == data[ref]["rho"].shape:
                fig_diff(name, ref, d, data[ref], os.path.join(args.out_dir, f"fig_diff_{name}_vs_{ref}.png"))
            elif name != ref:
                print(f"[{name}] parameter count differs from reference; no difference plot")

    import pandas as pd

    table = pd.DataFrame({n: d["stats"] for n, d in data.items()}).T
    front = ["N", "p", "p_var"]
    table = table[front + [c for c in table.columns if c not in front]]
    table.to_csv(os.path.join(args.out_dir, "scalar_table.csv"))
    with open(os.path.join(args.out_dir, "scalar_table.md"), "w") as f:
        cols = ["runcard"] + list(table.columns)
        f.write("| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n")
        for n, row in table.iterrows():
            f.write("| " + " | ".join([n] + [f"{v:.4g}" for v in row]) + " |\n")
    print()
    print(table.to_string(float_format=lambda x: f"{x:.4g}"))
    print(f"\noutputs in {os.path.abspath(args.out_dir)}")


if __name__ == "__main__":
    main()
