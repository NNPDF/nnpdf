#!/bin/bash
#SBATCH --job-name=bnn_resample
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=resample_%j.out
#SBATCH --error=resample_%j.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

cd /path/to/runcard

python -m n3fit.bnn_inference sample \
    --runcard Basic_runcard_bayesian.yml \
    --samples 100
