#!/bin/bash
#SBATCH -J nnpdf40_drop
#SBATCH -p boost_usr_prod
#SBATCH --array=1-120%40
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail

cd /leonardo/home/userexternal/rbonnetg/projects/nnpdf_bayesian
mkdir -p logs

source /leonardo/home/userexternal/rbonnetg/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

n3fit n3fit/runcards/examples/nnpdf40-like-dropout-cluster-larger.yml "${SLURM_ARRAY_TASK_ID}"
