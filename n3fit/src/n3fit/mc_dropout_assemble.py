"""
mc_dropout_assemble.py
======================

Converts per-replica MC Dropout inference files (.npz) produced by
``mc_dropout_inference.py`` into a new NNPDF fit directory that is ready
for DGLAP evolution with ``evolven3fit``.

Two assembly modes are supported:

  --mode mean  (default)
    For every trained replica N write ONE exportgrid using the mean over all
    MC Dropout samples.  Gives 1 PDF member per trained replica.
    Uncertainty comes from inter-replica diversity (same as standard NNPDF).

  --mode samples
    For every trained replica N write ONE exportgrid per MC Dropout sample.
    Gives n_samples PDF members per trained replica (sequential numbering).
    Uncertainty comes from the MC Dropout spread within a single trained
    replica.  This is the mode you want to study MC Dropout uncertainty.
    Example: 1 trained replica x 100 samples -> 100 PDF members.

Usage
-----
    # mean mode (120 trained replicas, 1 member each)
    python -m n3fit.mc_dropout_assemble \\
        --source-fit nnpdf40-like-dropout-cluster \\
        --mc-dir     output/mc_dropout_replicas \\
        --out-fit    nnpdf40-like-dropout-mc

    # samples mode (1 trained replica, 100 members from MC Dropout)
    python -m n3fit.mc_dropout_assemble \\
        --source-fit nnpdf40-like-dropout-cluster \\
        --mc-dir     output/mc_dropout_replicas \\
        --out-fit    nnpdf40-like-dropout-mc-1rep-samples \\
        --n-replicas 1 \\
        --mode samples
"""

import argparse
import pathlib
import shutil
import sys

import numpy as np
from validphys.utils import yaml_safe
from n3fit.io.writer import evln2lha

# exportgrid label order - identical to what storefit() in n3fit/io/writer.py writes
EXPORTGRID_LABELS = [
    "TBAR", "BBAR", "CBAR", "SBAR", "UBAR", "DBAR",
    "GLUON", "D", "U", "S", "C", "B", "T", "PHT",
]

# ---------------------------------------------------------------------------
# YAML helpers - reuse yaml_safe exactly as storefit() in n3fit/io/writer.py
# ---------------------------------------------------------------------------

def _write_exportgrid(path, pdfgrid, xgrid, q20, replica_idx):
    """
    Write an exportgrid YAML file using the same yaml_safe.dump that
    n3fit's storefit() uses, ensuring byte-level compatibility with
    evolven3fit's yaml_safe.load.

    Parameters
    ----------
    path        : pathlib.Path
    pdfgrid     : np.ndarray, shape (n_x, 14), exportgrid label order, x*PDF values
    xgrid       : list of float, n_x values
    q20         : float, Q0^2 in GeV^2
    replica_idx : int
    """
    data = {
        "replica": int(replica_idx),
        "q20":      float(q20),
        "xgrid":    [float(x) for x in xgrid],
        "labels":   EXPORTGRID_LABELS,
        "pdfgrid":  pdfgrid.tolist(),
    }
    with open(path, "w", encoding="utf-8") as fh:
        yaml_safe.dump(data, fh)


def _load_ref_meta(src_exportgrid_path):
    """
    Read q20 and xgrid from an existing .exportgrid file using yaml_safe.load,
    matching how evolven3fit reads it.
    """
    with open(src_exportgrid_path, encoding="utf-8") as fh:
        data = yaml_safe.load(fh)
    return float(data["q20"]), list(data["xgrid"])


def _get_active_flavours(src_fit, q0):
    """
    Determine the number of active flavours at Q0 using the same logic as
    storefit() in n3fit/io/writer.py.  Reads the theoryID from filter.yml
    and queries the validphys theory database.

    Falls back to nf=3 if the theory cannot be loaded.
    """
    filter_path = src_fit / "filter.yml"
    if not filter_path.exists():
        filter_path = src_fit / "filter.yaml"
    if not filter_path.exists():
        return 3
    with open(filter_path, encoding="utf-8") as fh:
        fc = yaml_safe.load(fh)
    theory_block = fc.get("theory", {})
    # filter.yml may use "theoryid" (lowercase) or "theoryID"
    theory_id = theory_block.get("theoryid") or theory_block.get("theoryID")
    if theory_id is None:
        return 3
    try:
        from validphys.loader import Loader
        theory = Loader().check_theoryID(theory_id)
        desc = theory.get_description()
        nf = 3
        for quark in ["c", "b", "t"]:
            mass      = desc.get(f"m{quark}")
            threshold = desc.get(f"k{quark}Thr")
            if q0 < mass * threshold:
                break
            nf += 1
        if desc.get("IC") == 1:
            nf = max(4, nf)
        return nf
    except Exception:
        return 3


# ---------------------------------------------------------------------------
# Main assembly function
# ---------------------------------------------------------------------------

def assemble(src_fit, mc_dir, out_fit, n_replicas=None, mode="mean"):
    """
    Parameters
    ----------
    src_fit    : path to the original trained fit directory
    mc_dir     : directory containing replica_N.npz files
    out_fit    : path/name for the new fit directory to create
    n_replicas : only process the first N trained replicas (None = all)
    mode       : "mean"    - 1 exportgrid per trained replica (mean of samples)
                 "samples" - 1 exportgrid per MC Dropout sample
    """
    if mode not in ("mean", "samples"):
        sys.exit(f"ERROR: --mode must be 'mean' or 'samples', got '{mode}'")

    src_fit = pathlib.Path(src_fit)
    mc_dir  = pathlib.Path(mc_dir)
    out_fit = pathlib.Path(out_fit)

    src_name = src_fit.name
    out_name = out_fit.name

    print(f"Source fit       : {src_fit}")
    print(f"MC Dropout npz   : {mc_dir}")
    print(f"Assembly mode    : {mode}")
    if n_replicas is not None:
        print(f"Replicas limit   : {n_replicas}")
    print(f"New fit dir      : {out_fit}")
    print()

    # ------------------------------------------------------------------
    # Validate inputs
    # ------------------------------------------------------------------
    if not mc_dir.is_dir():
        sys.exit(f"ERROR: MC Dropout directory not found: {mc_dir}")
    if not (src_fit / "nnfit").is_dir():
        sys.exit(f"ERROR: Source fit does not contain nnfit/: {src_fit}")

    # Reference exportgrid for q20 and xgrid
    ref_eg_path = src_fit / "nnfit" / "replica_1" / f"{src_name}.exportgrid"
    if not ref_eg_path.is_file():
        sys.exit(f"ERROR: Reference exportgrid not found: {ref_eg_path}")

    q20, xgrid = _load_ref_meta(ref_eg_path)
    q0 = float(np.sqrt(q20))
    nf = _get_active_flavours(src_fit, q0)
    print(f"Q0^2      = {q20:.6f} GeV^2  (Q0 = {q0:.4f} GeV, nf = {nf})")
    print(f"n_x       = {len(xgrid)}")
    print()

    # ------------------------------------------------------------------
    # Create output directory skeleton
    # ------------------------------------------------------------------
    out_fit.mkdir(parents=True, exist_ok=True)
    (out_fit / "nnfit").mkdir(exist_ok=True)

    # Copy filter.yml (needed by evolven3fit to look up the theoryID)
    for fname in ("filter.yml", "filter.yaml"):
        src_filter = src_fit / fname
        if src_filter.is_file():
            shutil.copy2(src_filter, out_fit / fname)
            print(f"Copied {fname}")
            break

    # Symlink filter/ directory (needed by vp-comparefits to locate cut datasets)
    src_filter_dir = src_fit / "filter"
    dst_filter_dir = out_fit / "filter"
    if src_filter_dir.is_dir() and not dst_filter_dir.exists():
        dst_filter_dir.symlink_to(src_filter_dir.resolve())
        print("Symlinked filter/ directory")

    # ------------------------------------------------------------------
    # Process replicas
    # ------------------------------------------------------------------
    src_rep_dirs = sorted(
        src_fit.glob("nnfit/replica_*/"),
        key=lambda p: int(p.name.split("_")[1]),
    )
    if n_replicas is not None:
        src_rep_dirs = [p for p in src_rep_dirs if int(p.name.split("_")[1]) <= n_replicas]

    found   = []   # list of trained replica indices successfully processed
    missing = []
    n_samples_per_replica = None
    out_replica_counter = 0   # sequential exportgrid index (used in samples mode)

    for src_rep_dir in src_rep_dirs:
        rep_idx  = int(src_rep_dir.name.split("_")[1])
        npz_path = mc_dir / f"replica_{rep_idx}.npz"

        if not npz_path.is_file():
            missing.append(rep_idx)
            continue

        data = np.load(npz_path)

        # Record n_samples from the first npz seen
        if n_samples_per_replica is None and "samples" in data:
            n_samples_per_replica = int(data["samples"].shape[0])

        # ----------------------------------------------------------------
        # Build the list of (exportgrid_index, pdf_array) pairs to write
        # ----------------------------------------------------------------
        if mode == "mean":
            # 1 exportgrid = mean over all MC samples, shape (n_x, 14) EVOL
            evol_list = [(rep_idx, data["mean"])]
        else:  # samples
            if "samples" not in data:
                sys.exit(
                    f"ERROR: replica_{rep_idx}.npz has no 'samples' array. "
                    "Re-run inference or use --mode mean."
                )
            samples = data["samples"]  # (n_samples, n_x, 14) EVOL
            evol_list = [
                (out_replica_counter + s + 1, samples[s])
                for s in range(samples.shape[0])
            ]
            out_replica_counter += samples.shape[0]

        for out_rep_idx, evol in evol_list:
            # Convert EVOL -> LHAPDF exportgrid label order
            eg = evln2lha(evol.T, nf=nf).T  # (n_x, 14)

            out_rep_dir = out_fit / "nnfit" / f"replica_{out_rep_idx}"
            out_rep_dir.mkdir(parents=True, exist_ok=True)

            _write_exportgrid(
                out_rep_dir / f"{out_name}.exportgrid",
                eg, xgrid, q20, out_rep_idx,
            )

            # Copy .json and chi2exps.log (required by postfit)
            src_json = src_rep_dir / f"{src_name}.json"
            if src_json.is_file():
                shutil.copy2(src_json, out_rep_dir / f"{out_name}.json")
            else:
                print(f"  WARNING: missing {src_json.name} for trained replica {rep_idx}")

            src_chi2 = src_rep_dir / "chi2exps.log"
            if src_chi2.is_file():
                shutil.copy2(src_chi2, out_rep_dir / "chi2exps.log")
            else:
                print(f"  WARNING: missing chi2exps.log for trained replica {rep_idx}")

        found.append(rep_idx)

    # ------------------------------------------------------------------
    # Summary + metadata
    # ------------------------------------------------------------------
    import datetime
    n_total_members = out_replica_counter if mode == "samples" else len(found)
    metadata = {
        "mode":                    mode,
        "source_fit":              str(src_fit),
        "n_trained_replicas":      len(found),
        "n_samples_per_replica":   n_samples_per_replica if found else None,
        "n_total_members":         n_total_members,
        # keep legacy key so comparefits script still works
        "n_replicas_assembled":    n_total_members,
        "q0_gev":                  round(q0, 4),
        "n_active_flavours":       nf,
        "assembled_at":            datetime.datetime.now().isoformat(timespec="seconds"),
    }
    meta_path = out_fit / "mc_dropout_metadata.yaml"
    with open(meta_path, "w", encoding="utf-8") as _mf:
        for k, v in metadata.items():
            _mf.write(f"{k}: {v}\n")

    print()
    print(f"Mode                : {mode}")
    print(f"Trained replicas    : {len(found)}")
    print(f"MC samples/replica  : {n_samples_per_replica}")
    print(f"Total PDF members   : {n_total_members}")
    print(f"Metadata written to : {meta_path}")
    if missing:
        print(f"Missing .npz files  : {missing}")
        print("  (run mc_dropout_inference.py for these replicas first)")

    print()
    print("Next steps:")
    print(f"  evolven3fit evolve -f {out_fit}")
    print(f"  postfit {n_total_members} {out_fit}")

    return len(found), missing


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Assemble MC Dropout exportgrids from per-replica .npz files "
            "into a new NNPDF fit directory ready for evolven3fit."
        )
    )
    parser.add_argument(
        "--source-fit",
        default="nnpdf40-like-dropout-cluster",
        help="Path to the original trained fit directory (default: nnpdf40-like-dropout-cluster)",
    )
    parser.add_argument(
        "--mc-dir",
        default="output/mc_dropout_replicas",
        help="Directory containing replica_N.npz files (default: output/mc_dropout_replicas)",
    )
    parser.add_argument(
        "--out-fit",
        default="nnpdf40-like-dropout-mc",
        help="Path / name for the new fit directory to create (default: nnpdf40-like-dropout-mc)",
    )
    parser.add_argument(
        "--n-replicas",
        type=int,
        default=None,
        help="Only assemble the first N replicas (default: all)",
    )
    parser.add_argument(
        "--mode",
        choices=["mean", "samples"],
        default="mean",
        help=(
            "Assembly mode: "
            "'mean' (default) writes 1 exportgrid per trained replica using the mean of MC samples; "
            "'samples' writes 1 exportgrid per MC sample (captures MC Dropout uncertainty)."
        ),
    )
    args = parser.parse_args()
    n_ok, missing = assemble(
        args.source_fit, args.mc_dir, args.out_fit,
        n_replicas=args.n_replicas, mode=args.mode,
    )
    if not n_ok:
        sys.exit("ERROR: No replicas were assembled. Aborting.")


if __name__ == "__main__":
    main()
