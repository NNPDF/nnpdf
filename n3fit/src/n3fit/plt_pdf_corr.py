"""
PDF-space correlation comparison across n3fit fits (Bayesian vs standard NNPDF).
(DOOB-ed [Documented-Optimized-Organized-Beautified] by Claude Code (Sonnet 5))

Weight-space correlations are not comparable between a mean-field BNN and a standard fit: the
weights are not identifiable (hidden-unit permutation/sign symmetry), and a mean-field posterior
has zero within-chain weight correlation by construction. The PDFs, however, are the same object
for every fit: each replica stores xf(x) at Q0 in its ``<fit>.exportgrid``. The correlation of
xf_a(x_i) with xf_b(x_j) across replicas is the standard way to compare NNPDF-style ensembles,
independent of architecture, chain structure and how the ensemble was produced.

For every fit, flavour and pair of x points this computes rho(x_i, x_j) across replicas, using the
evolution-like combinations built from the flavour-basis exportgrid:
    g (gluon), Sigma = sum_q (q + qbar) [u,d,s,c], V = sum_q (q - qbar) [u,d,s,c],
    T3 = (u + ubar) - (d + dbar)

Outputs (in --out-dir):
    fig_pdf_corr_heatmaps.png   rho(x, x') per flavour (rows) and fit (columns), common +-1 scale
    fig_pdf_sigma.png           PDF uncertainty sigma[xf](x) per flavour, fits overlaid
    fig_pdf_spectrum.png        cumulative variance of the PDF-ensemble eigenmodes per flavour
    pdf_scalar_table.csv/.md    single numbers, mean +- std over random equal-size subsamples

Every fit's own LHAPDF setting ForcePositive (``nnfit/<fit>.info``) is applied: with ForcePositive 1
each flavour's xf < 0 is set to 0 before the combinations are built (what LHAPDF returns at Q0), with
2 to 1e-10. Likewise a fit whose runcard sets ``parameters::bayesian::per_replica::clip_negative``
(postfit writes its LHAPDF members with xf < 0 -> 0 for u, ubar, d, dbar, s, sbar, g) gets the same
clipping of those flavours in its exportgrids. --postfit-only uses only the replicas postfit kept
(``<fit>/postfit/replica_*``).

Single-number comparison (all computed from replica subsets of the SAME size n, so that every fit
has the same sampling noise; `n` defaults to half the reference fit's replicas so two disjoint
subsets of the reference can measure the noise floor):
    mean_abs_rho    mean |rho(x_i, x_j)| over i != j
    n_modes_PR      participation ratio (sum lam)^2 / sum lam^2 of the rho eigenvalues: the effective
                    number of independent PDF shapes in the ensemble
    sigma_ratio     geometric mean over x of sigma_fit(x) / sigma_ref(x)
    dist_to_ref     ||rho_fit - rho_ref||_F / n_x  (Frobenius distance of the correlation matrices)
    dist_null       the same distance between two disjoint subsets of the reference fit itself:
                    the value expected if the two ensembles had identical true correlations

Example:
    python plt_pdf_corr.py bay_1x1000 bay_10x100 bay_100x10 nnpdf40-like --reference nnpdf40-like
"""

import argparse
import glob
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

DEFAULT_RESULTS = os.path.join(
    os.path.expanduser("~"), "miniconda3/envs/environment_nnpdf/share/NNPDF/results"
)
COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]
FLAVOURS = ["g", "Sigma", "V", "T3"]
LOADER = getattr(yaml, "CSafeLoader", yaml.SafeLoader)


def _natural(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def _combos(grid, labels):
    """(nx, 14) flavour-basis xf -> dict of (nx,) arrays for the analysed combinations."""
    i = {l: k for k, l in enumerate(labels)}
    q = lambda n: grid[:, i[n]]
    qb = lambda n: grid[:, i[n + "BAR"]]
    quarks = ["U", "D", "S", "C"]
    return {
        "g": grid[:, i["GLUON"]],
        "Sigma": sum(q(n) + qb(n) for n in quarks),
        "V": sum(q(n) - qb(n) for n in quarks),
        "T3": (q("U") + qb("U")) - (q("D") + qb("D")),
    }


def force_positive(name, results):
    """The fit's LHAPDF ForcePositive setting (0: none, 1: xf < 0 -> 0, 2: xf < 1e-10 -> 1e-10)"""
    info = os.path.join(results, name, "nnfit", f"{name}.info")
    if not os.path.exists(info):
        return 0
    with open(info) as fh:
        m = re.search(r"^ForcePositive:\s*(\d)", fh.read(), flags=re.M)
    return int(m.group(1)) if m else 0


# exportgrid labels of the flavours postfit's clip_negative sets to 0 where negative
# (validphys.replica_checks.CLIP_NEGATIVE_PDG: sbar, ubar, dbar, g, d, u, s)
CLIP_NEGATIVE_LABELS = ("SBAR", "UBAR", "DBAR", "GLUON", "D", "U", "S")


def clip_negative(name, results):
    """True if the fit runcard (``<fit>/filter.yml``) sets parameters::bayesian::per_replica::clip_negative"""
    runcard = os.path.join(results, name, "filter.yml")
    if not os.path.exists(runcard):
        return False
    with open(runcard) as fh:
        params = (yaml.load(fh, Loader=LOADER) or {}).get("parameters") or {}
    return bool(((params.get("bayesian") or {}).get("per_replica") or {}).get("clip_negative", False))


def load_fit(name, results, out_dir, reload, postfit_only=False):
    """-> (x (nx,), {flavour: (N, nx)}), with the fit's ForcePositive applied. Cached as npz in out_dir."""
    files = glob.glob(os.path.join(results, name, "nnfit", "replica_*", f"{name}.exportgrid"))
    if not files:
        print(f"[{name}] no exportgrid files under {results}/{name}/nnfit/replica_*/ -- skipping")
        return None
    files.sort(key=_natural)
    if postfit_only:
        kept = {os.path.basename(os.path.realpath(p)) for p in glob.glob(os.path.join(results, name, "postfit", "replica_*"))}
        if not kept:
            sys.exit(f"[{name}] --postfit-only: no postfit replicas under {results}/{name}/postfit/")
        files = [f for f in files if os.path.basename(os.path.dirname(f)) in kept]
        print(f"[{name}] --postfit-only: {len(files)} postfit replicas")
    fpos = force_positive(name, results)
    if fpos:
        print(f"[{name}] ForcePositive {fpos}: xf < {0 if fpos == 1 else 1e-10:g} set to {0 if fpos == 1 else 1e-10:g} per flavour")
    clip = clip_negative(name, results)
    if clip:
        print(f"[{name}] clip_negative: xf < 0 set to 0 for {', '.join(CLIP_NEGATIVE_LABELS)}")
    cache = os.path.join(out_dir, f"{name}_pdfs{'_pf' if postfit_only else ''}{f'_fp{fpos}' if fpos else ''}{'_cn' if clip else ''}_n{len(files)}.npz")
    if os.path.exists(cache) and not reload:
        z = np.load(cache)
        print(f"[{name}] loaded cached PDFs from {cache}")
        return z["x"], {f: z[f] for f in FLAVOURS}
    print(f"[{name}] reading {len(files)} exportgrids ...")
    rows = {f: [] for f in FLAVOURS}
    x = None
    for path in files:
        with open(path) as fh:
            d = yaml.load(fh, Loader=LOADER)
        x = np.array([float(v) for v in d["xgrid"]])  # yaml parses '1e-09' as a string
        grid = np.array(d["pdfgrid"], dtype=float)
        if fpos:
            grid = np.maximum(grid, 0.0 if fpos == 1 else 1e-10)
        if clip:
            cols = [d["labels"].index(label) for label in CLIP_NEGATIVE_LABELS]
            grid[:, cols] = np.maximum(grid[:, cols], 0.0)
        for f, v in _combos(grid, d["labels"]).items():
            rows[f].append(v)
    pdfs = {f: np.stack(v) for f, v in rows.items()}
    np.savez(cache, x=x, **pdfs)
    return x, pdfs


def corr(F):
    """Correlation over replicas of (N, nx); NaN where the PDF is constant."""
    std = F.std(axis=0)
    ok = std > 1e-12 * max(np.abs(F).max(), 1e-30)
    rho = np.full((F.shape[1],) * 2, np.nan)
    idx = np.flatnonzero(ok)
    rho[np.ix_(idx, idx)] = np.corrcoef(F[:, idx], rowvar=False)
    return rho


def eig_stats(rho):
    ok = np.isfinite(np.diag(rho))
    r = rho[np.ix_(ok, ok)]
    lam = np.clip(np.linalg.eigvalsh(r)[::-1], 0, None)
    n = r.shape[0]
    return dict(
        mean_abs_rho=np.abs(r[~np.eye(n, dtype=bool)]).mean(),
        n_modes_PR=lam.sum() ** 2 / (lam**2).sum(),
    ), lam


def frob(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return np.sqrt(((a - b)[m] ** 2).sum()) / a.shape[0]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runcards", nargs="+")
    ap.add_argument("--reference", required=True, help="fit the others are compared to (e.g. the standard fit)")
    ap.add_argument("--results-dir", default=DEFAULT_RESULTS)
    ap.add_argument("--out-dir", default="plt_pdf_corr_out")
    ap.add_argument("--xmin", type=float, default=1e-3)
    ap.add_argument("--xmax", type=float, default=0.7)
    ap.add_argument("--n-sub", type=int, default=None,
                    help="replicas per subsample for the table (default: half the reference's replicas)")
    ap.add_argument("--n-draws", type=int, default=20, help="random subsamples for the table")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--drop-outliers", type=float, default=None, metavar="K",
                    help="drop replicas whose xf deviates from the ensemble median by more than K robust "
                         "(MAD) sigmas at any x-point in any flavour. Use K=50: ordinary replicas of multi-chain fits reach 12-15 sigma, a failed chain sits at ~100. These are failed replicas/chains "
                         "that a postfit-style filter would remove; one such chain can inflate sigma by 4x")
    ap.add_argument("--reload", action="store_true")
    ap.add_argument("--postfit-only", action="store_true", help="only the replicas postfit kept, not every nnfit replica")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    fits, x = {}, None
    for name in dict.fromkeys(args.runcards + [args.reference]):
        got = load_fit(name, args.results_dir, args.out_dir, args.reload, args.postfit_only)
        if got is None:
            if name == args.reference:
                sys.exit("reference fit unavailable")
            continue
        x, pdfs = got
        sel = (x >= args.xmin) & (x <= args.xmax)
        fits[name] = {f: v[:, sel] for f, v in pdfs.items()}
        # flag failed replicas (robust z-score per x-point) and, on request, drop them
        bad = np.zeros(len(fits[name]["g"]), bool)
        for v in fits[name].values():
            med = np.median(v, axis=0)
            mad = 1.4826 * np.median(np.abs(v - med), axis=0)
            bad |= ((np.abs(v - med) / np.maximum(mad, 1e-12 * np.abs(v).max())) > (args.drop_outliers or 50)).any(axis=1)
        if bad.any():
            print(f"[{name}] {bad.sum()}/{len(bad)} replicas have an >{args.drop_outliers or 50} MAD-sigma excursion "
                  f"(replica indices {np.flatnonzero(bad)[:12].tolist()}{'...' if bad.sum() > 12 else ''})"
                  + ("; dropping them" if args.drop_outliers else "; use --drop-outliers K to remove"))
            if args.drop_outliers:
                fits[name] = {f: v[~bad] for f, v in fits[name].items()}
    xs = x[(x >= args.xmin) & (x <= args.xmax)]
    ref = args.reference
    N = {n: fits[n]["g"].shape[0] for n in fits}
    n_sub = args.n_sub or min(min(N.values()), N[ref] // 2)
    if 2 * n_sub > N[ref]:
        sys.exit(f"reference has {N[ref]} replicas; need >= {2 * n_sub} for two disjoint subsets")
    print(f"replicas: {N}; x range [{xs[0]:.1e}, {xs[-1]:.2f}] ({len(xs)} points); n_sub = {n_sub}")
    rng = np.random.default_rng(args.seed)

    # ---------------------------------------------------------------- table
    stats = {}  # (fit, flavour) -> {metric: [values over draws]}
    for _ in range(args.n_draws):
        perm = rng.permutation(N[ref])
        ref_a, ref_b = perm[:n_sub], perm[n_sub : 2 * n_sub]
        for f in FLAVOURS:
            rho_a = corr(fits[ref][f][ref_a])
            rho_b = corr(fits[ref][f][ref_b])
            sig_ref = fits[ref][f][ref_a].std(axis=0)
            null = frob(rho_a, rho_b)
            for name in fits:
                sub = ref_b if name == ref else rng.choice(N[name], n_sub, replace=False)
                F = fits[name][f][sub]
                rho = corr(F)
                st, _ = eig_stats(rho)
                st["sigma_ratio"] = np.exp(np.mean(np.log(F.std(axis=0) / sig_ref)))
                st["dist_to_ref"] = frob(rho, rho_a)
                st["dist_null"] = null
                for k, v in st.items():
                    stats.setdefault((name, f), {}).setdefault(k, []).append(v)
    cols = ["mean_abs_rho", "n_modes_PR", "sigma_ratio", "dist_to_ref", "dist_null"]
    rows = []
    for (name, f), d in stats.items():
        rows.append([name, f, N[name]] + [f"{np.mean(d[c]):.3g} +- {np.std(d[c]):.2g}" for c in cols])
    header = ["fit", "flavour", "N_total"] + cols
    with open(os.path.join(args.out_dir, "pdf_scalar_table.md"), "w") as fh:
        fh.write(f"n_sub = {n_sub} replicas per subsample, {args.n_draws} draws (mean +- std)\n\n")
        fh.write("| " + " | ".join(header) + " |\n|" + "---|" * len(header) + "\n")
        for r in rows:
            fh.write("| " + " | ".join(map(str, r)) + " |\n")
    with open(os.path.join(args.out_dir, "pdf_scalar_table.csv"), "w") as fh:
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join(str(v).replace(",", ";") for v in r) + "\n")
    print(f"\nn_sub = {n_sub}, {args.n_draws} draws (mean +- std)")
    print("  ".join(f"{h:>17s}" for h in header))
    for r in rows:
        print("  ".join(f"{str(v):>17s}" for v in r))

    # ------------------------------------------------------------- figures
    # heatmaps and spectra use the same number of replicas for every fit (the smallest ensemble), so
    # their sampling noise is identical
    n_plot = min(N.values())
    full = {}
    for name in fits:
        sub = np.sort(rng.choice(N[name], n_plot, replace=False))
        full[name] = {f: fits[name][f][sub] for f in FLAVOURS}
    names = list(fits)
    fig, axs = plt.subplots(len(FLAVOURS), len(names), figsize=(3.4 * len(names) + 0.8, 3.2 * len(FLAVOURS)),
                            squeeze=False)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#d9d9d9")
    for r, f in enumerate(FLAVOURS):
        for c, name in enumerate(names):
            ax = axs[r, c]
            im = ax.pcolormesh(xs, xs, corr(full[name][f]), cmap=cmap, vmin=-1, vmax=1, shading="auto")
            ax.set_xscale("log")
            ax.set_yscale("log")
            if r == 0:
                ax.set_title(f"{name}", fontsize=9)
            if c == 0:
                ax.set_ylabel(f"{f}\n$x'$")
            if r == len(FLAVOURS) - 1:
                ax.set_xlabel("$x$")
    fig.suptitle(rf"$\rho(x,x')$ of $xf(x, Q_0)$ across replicas, N = {n_plot} replicas per fit", y=0.995, fontsize=10)
    fig.colorbar(im, ax=axs.ravel().tolist(), label=r"$\rho$", shrink=0.6)
    fig.savefig(os.path.join(args.out_dir, "fig_pdf_corr_heatmaps.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axs = plt.subplots(1, len(FLAVOURS), figsize=(4 * len(FLAVOURS), 3.6))
    for f, ax in zip(FLAVOURS, axs):
        for col, name in zip(COLORS, names):
            ax.plot(xs, fits[name][f].std(axis=0), c=col, label=f"{name} (N={N[name]})")
        ax.set_xscale("log")
        ax.set_xlabel("$x$")
        ax.set_title(f"$\\sigma[x{f}]$" if f != "g" else r"$\sigma[xg]$", fontsize=10)
    axs[0].legend(fontsize=7)
    fig.savefig(os.path.join(args.out_dir, "fig_pdf_sigma.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axs = plt.subplots(1, len(FLAVOURS), figsize=(4 * len(FLAVOURS), 3.6))
    for f, ax in zip(FLAVOURS, axs):
        for col, name in zip(COLORS, names):
            _, lam = eig_stats(corr(full[name][f]))
            ax.plot(np.arange(1, len(lam) + 1), np.cumsum(lam) / lam.sum(), c=col, label=name)
        ax.set_xscale("log")
        ax.set_xlabel("number of modes")
        ax.set_title(f, fontsize=10)
    axs[0].set_ylabel("cumulative variance fraction")
    axs[0].legend(fontsize=7)
    fig.savefig(os.path.join(args.out_dir, "fig_pdf_spectrum.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"\noutputs in {os.path.abspath(args.out_dir)}")


if __name__ == "__main__":
    main()
