#!/bin/bash
#SBATCH --job-name=mc_dropout_inf
#SBATCH --account=inf26_pml4hep_0
#SBATCH --partition=dcgp_usr_prod
#SBATCH --qos=normal
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm/mc_dropout_inf-%A_%a.out
#SBATCH --error=slurm/mc_dropout_inf-%A_%a.err

# MC Dropout inference - Option A
# Runs N_SAMPLES stochastic forward passes on each trained replica and saves
# the per-replica (mean, std, samples) to output/mc_dropout_replicas/.
#
# Submit with:
#   sbatch slurm/nnpdf40_mc_dropout_inference.slurm
#
# Optional overrides (export before submitting or pass via --export):
#   FIT_DIR    path to fit directory  (default: nnpdf40-like-dropout-cluster)
#   N_SAMPLES  number of MC passes    (default: 200)

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

cd /leonardo/home/userexternal/rbonnetg/projects/nnpdf_bayesian

: "${FIT_DIR:=nnpdf40-like-dropout-cluster}"
: "${N_SAMPLES:=200}"

REPLICA=${SLURM_ARRAY_TASK_ID}
OUTDIR="output/mc_dropout_replicas"
mkdir -p "${OUTDIR}"

echo "=== MC Dropout inference - replica ${REPLICA} / N_SAMPLES=${N_SAMPLES} ==="
echo "Fit dir : ${FIT_DIR}"
echo "Node    : $(hostname)"
echo "Started : $(date)"

TF_CPP_MIN_LOG_LEVEL=2 PYTHONPATH=n3fit/src:${PYTHONPATH:-} \
  python -m n3fit.mc_dropout_inference \
    --fit-dir  "${FIT_DIR}" \
    --replica  "${REPLICA}" \
    --n-samples "${N_SAMPLES}" \
    --output   "${OUTDIR}/replica_${REPLICA}.npz"

echo "Done    : $(date)"
