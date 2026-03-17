#!/bin/bash
# submit_mc_dropout_pipeline.sh
# ==============================
# Master submission script for the full MC Dropout inference pipeline.
#
# Pipeline stages:
#   1. Inference (array job)      : run N MC Dropout passes on every trained replica
#                                   -> output/mc_dropout_replicas/replica_N.npz
#   2. Assemble + evolve + postfit: build new fit dir, DGLAP-evolve, run postfit
#                                   -> nnpdf40-like-dropout-mc/postfit/
#   3. Compare fits               : vp-comparefits MC fit vs NNPDF4.0 baseline
#                                   -> output/vp_compare_mc_vs_baseline/
#
# Usage:
#   bash submit_mc_dropout_pipeline.sh                        # full run, all 120 replicas
#   bash submit_mc_dropout_pipeline.sh --n-replicas 10        # run replicas 1-10
#   bash submit_mc_dropout_pipeline.sh --replica 5            # single replica test only
#
# Optional flags:
#   --n-replicas N  number of replicas to run (default: 120, runs 1-N)
#   --replica N     run only this one replica (test mode, skips stages 2+3)
#   --fit-dir DIR   trained fit directory (default: nnpdf40-like-dropout-cluster)
#   --n-samples N   MC Dropout passes     (default: 100)
#   --out-fit DIR   new fit name          (default: nnpdf40-like-dropout-mc)
#   --ref-fit DIR   reference fit         (default: 240701-02-rs-nnpdf40-baseline)
#   --mode MODE     assembly mode: 'mean' or 'samples' (default: mean)
#                   'mean':    1 PDF member per trained replica (inter-replica uncertainty)
#                   'samples': 1 PDF member per MC sample (MC Dropout uncertainty)

set -euo pipefail

REPLICA=""
N_REPLICAS=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --replica)    REPLICA="$2";    shift 2 ;;
    --n-replicas) N_REPLICAS="$2"; shift 2 ;;
    --fit-dir)    FIT_DIR="$2";    shift 2 ;;
    --n-samples)  N_SAMPLES="$2";  shift 2 ;;
    --out-fit)    OUT_FIT="$2";    shift 2 ;;
    --ref-fit)    REF_FIT="$2";    shift 2 ;;
    --mode)       MODE="$2";       shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

: "${FIT_DIR:=nnpdf40-like-dropout-cluster}"
: "${N_SAMPLES:=100}"
: "${OUT_FIT:=nnpdf40-like-dropout-mc}"
: "${REF_FIT:=240701-02-rs-nnpdf40-baseline}"
: "${N_REPLICAS:=120}"
: "${MODE:=mean}"
MC_DIR="output/mc_dropout_replicas"

echo "======================================================"
echo " MC Dropout pipeline submission"
echo "======================================================"
echo " Fit dir       : ${FIT_DIR}"
echo " Replicas      : ${REPLICA:-1-${N_REPLICAS}}"
echo " MC samples    : ${N_SAMPLES}"
echo " New fit       : ${OUT_FIT}"
echo " Reference     : ${REF_FIT}"
echo "======================================================"
echo

# ------------------------------------------------------------------
# Stage 1: inference array job
# ------------------------------------------------------------------
if [ -n "$REPLICA" ]; then
  ARRAY_SPEC="${REPLICA}-${REPLICA}"
else
  ARRAY_SPEC="1-${N_REPLICAS}%30"
fi

JOB1=$(sbatch --parsable \
  --array="${ARRAY_SPEC}" \
  --export=ALL,FIT_DIR="${FIT_DIR}",N_SAMPLES="${N_SAMPLES}" \
  slurm/nnpdf40_mc_dropout_inference.slurm)

echo "[stage 1] Inference array submitted : ${JOB1}  (array: ${ARRAY_SPEC})"
echo "          output -> ${MC_DIR}/replica_N.npz"

if [ -n "$REPLICA" ]; then
  echo
  echo "======================================================"
  echo " Single-replica test submitted."
  echo " Monitor:  squeue -u \$USER"
  echo " Log:      tail -f slurm/mc_dropout_inf-${JOB1}_${REPLICA}.out"
  echo
  echo " To run the full pipeline:"
  echo "   bash submit_mc_dropout_pipeline.sh"
  echo "======================================================"
  exit 0
fi

# ------------------------------------------------------------------
# Stage 2: assemble + evolve + postfit  (after all array tasks)
# ------------------------------------------------------------------
JOB2=$(sbatch --parsable \
  --dependency=afterok:"${JOB1}" \
  --export=ALL,SOURCE_FIT="${FIT_DIR}",MC_DIR="${MC_DIR}",OUT_FIT="${OUT_FIT}",N_REPLICAS="${N_REPLICAS}",MODE="${MODE}" \
  slurm/nnpdf40_mc_dropout_assemble_evolve.slurm)

echo "[stage 2] Assemble+evolve+postfit   : ${JOB2}  (after ${JOB1})"

# ------------------------------------------------------------------
# Stage 3: vp-comparefits  (after stage 2)
# ------------------------------------------------------------------
JOB3=$(sbatch --parsable \
  --dependency=afterok:"${JOB2}" \
  --export=ALL,MC_FIT="${OUT_FIT}",REF_FIT="${REF_FIT}",OUT_DIR="output/vp_compare_${OUT_FIT}_vs_baseline" \
  slurm/nnpdf40_mc_dropout_comparefits.slurm)

echo "[stage 3] vp-comparefits            : ${JOB3}  (after ${JOB2})"
echo
echo "======================================================"
echo " All jobs submitted. Monitor with:  squeue -u \$USER"
echo
echo " Expected outputs:"
echo "   ${MC_DIR}/replica_N.npz        (stage 1)"
echo "   ${OUT_FIT}/postfit/             (stage 2)"
echo "   output/vp_compare_${OUT_FIT}_vs_baseline/  (stage 3)"
echo "======================================================"
