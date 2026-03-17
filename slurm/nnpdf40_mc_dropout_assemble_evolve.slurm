#!/bin/bash
#SBATCH --job-name=mc_drop_assemble
#SBATCH --account=inf26_pml4hep_0
#SBATCH --partition=dcgp_usr_prod
#SBATCH --qos=normal
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/mc_dropout_assemble-%j.out
#SBATCH --error=slurm/mc_dropout_assemble-%j.err

# Full pipeline: assemble exportgrids -> DGLAP evolution -> postfit
#
# This job is designed to run AFTER all inference array tasks have finished.
# Chain it with:
#   sbatch --dependency=afterok:<INFERENCE_JOB_ID> this_script.slurm
#
# Overridable environment variables (set before submission via --export):
#   SOURCE_FIT   original trained fit dir  (default: nnpdf40-like-dropout-cluster)
#   MC_DIR       directory with .npz files (default: output/mc_dropout_replicas)
#   OUT_FIT      name of new fit to create (default: nnpdf40-like-dropout-mc)
#   N_REPLICAS   only assemble first N trained replicas (default: all)
#   MODE         assembly mode: 'mean' or 'samples' (default: mean)

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

cd /leonardo/home/userexternal/rbonnetg/projects/nnpdf_bayesian

: "${SOURCE_FIT:=nnpdf40-like-dropout-cluster}"
: "${MC_DIR:=output/mc_dropout_replicas}"
: "${OUT_FIT:=nnpdf40-like-dropout-mc}"
: "${N_REPLICAS:=}"  # empty = all replicas
: "${MODE:=mean}"    # mean or samples

echo "======================================================"
echo "MC Dropout pipeline: assemble + evolve + postfit"
echo "======================================================"
echo "Source fit : ${SOURCE_FIT}"
echo "MC npz dir : ${MC_DIR}"
echo "Output fit : ${OUT_FIT}"
echo "Node       : $(hostname)"
echo "Started    : $(date)"
echo

# ------------------------------------------------------------------
# Step 1: Assemble exportgrids from MC Dropout .npz files
# ------------------------------------------------------------------
echo "[1/3] Assembling exportgrids ..."
N_REPS_ARG=""
[[ -n "${N_REPLICAS}" ]] && N_REPS_ARG="--n-replicas ${N_REPLICAS}"
TF_CPP_MIN_LOG_LEVEL=2 PYTHONPATH=n3fit/src:${PYTHONPATH:-} \
  python -m n3fit.mc_dropout_assemble \
    --source-fit "${SOURCE_FIT}" \
    --mc-dir     "${MC_DIR}" \
    --out-fit    "${OUT_FIT}" \
    --mode       "${MODE}" \
    ${N_REPS_ARG}
echo

# ------------------------------------------------------------------
# Step 2: DGLAP evolution
# ------------------------------------------------------------------
echo "[2/3] Running evolven3fit ..."
# Remove stale log file if it exists so evolven3fit does not abort
rm -f "${OUT_FIT}/evolven3fit.log"

evolven3fit evolve -f "${OUT_FIT}"
echo

# ------------------------------------------------------------------
# Step 3: Postfit replica selection
# ------------------------------------------------------------------
: "${POSTFIT_REPS:=}"
# auto: read n_total_members from metadata written by assemble step
if [[ -z "${POSTFIT_REPS}" ]]; then
  META_FILE="${OUT_FIT}/mc_dropout_metadata.yaml"
  if [[ -f "${META_FILE}" ]]; then
    POSTFIT_REPS=$(grep "^n_total_members:" "${META_FILE}" | awk '{print $2}')
  fi
fi
# final fallback
: "${POSTFIT_REPS:=100}"

echo "[3/3] Running postfit (selecting ${POSTFIT_REPS} replicas) ..."
postfit "${POSTFIT_REPS}" "${OUT_FIT}"
echo

echo "======================================================"
echo "Done: $(date)"
echo "New fit ready at: ${OUT_FIT}/postfit/${OUT_FIT}/"
echo
echo "Next step:"
echo "  sbatch slurm/nnpdf40_mc_dropout_comparefits.slurm"
echo "======================================================"
