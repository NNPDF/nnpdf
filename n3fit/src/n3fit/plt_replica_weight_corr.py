"""
Weight structure-map heatmap for one replica from each of two (or more) saved fits.

Not an across-replica correlation (that needs many replicas of the SAME fit -- see
plt_corr.py's fig_heatmaps.png); with a single weights.weights.h5 per fit there is no
ensemble to correlate over. What is plotted is plt_corr.py's single_replica_matrix(): the
standardised self-outer-product z_i*z_j of that replica's flattened weight vector,
z = (w - mean(w)) / std(w) -- a deterministic structure map of that one trained/sampled
network (do large-magnitude weights coincide in sign?), not a statistical correlation.
See single_replica_matrix()'s docstring in plt_corr.py for the full caveat.

Reuses plt_corr.py's _flatten_replica / single_replica_matrix / fig_single_replica so the
parameter layout (dense vs VBDense, --which sampled/mean) stays identical to plt_corr.py.

Example:
    python plt_replica_weight_corr.py --replica 42 \\
        bay_1x1000=/mnt/c/Users/daksh/nnpdf_reports/weights/bay_1x1000/weights.weights.h5 \\
        nnpdf40-like=/mnt/c/Users/daksh/nnpdf_reports/weights/nnpdf40-like/weights.weights.h5
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plt_corr import _flatten_replica, fig_single_replica, single_replica_matrix


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("weights", nargs="+", metavar="NAME=PATH",
                    help="label=path to a weights.weights.h5 file, one per fit")
    ap.add_argument("--replica", type=int, required=True, help="replica index, used only in the figure title")
    ap.add_argument("--which", choices=["sampled", "mean"], default="sampled",
                    help="VBDense weights: frozen eval draw (default) or posterior mean mu_w")
    ap.add_argument("--include-preproc", action="store_true", help="also include preprocessing exponents")
    ap.add_argument("--last-layer", action="store_true", help="use only the last layer (kernel + bias)")
    ap.add_argument("--out-dir", default="plt_replica_corr_out")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    mats, blocks = {}, {}
    for spec in args.weights:
        name, sep, path = spec.partition("=")
        if not sep:
            sys.exit(f"expected NAME=PATH, got {spec!r}")
        if not os.path.exists(path):
            sys.exit(f"[{name}] no such file: {path}")
        w, b, biases = _flatten_replica(path, args.which, args.include_preproc)
        if args.last_layer:
            lo = w.shape[0] - b[-1][1]
            w, nb = w[lo:], biases[-1]
            b = [(f"{b[-1][0]} kernel", w.shape[0] - nb), (f"{b[-1][0]} bias", nb)]
        mats[name], blocks[name] = single_replica_matrix(w), b
        print(f"[{name}] {path}: p={w.shape[0]} weights")

    out = os.path.join(args.out_dir, f"fig_single_replica_{args.replica}.png")
    fig_single_replica(mats, blocks, args.replica, out)
    print(f"wrote {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
