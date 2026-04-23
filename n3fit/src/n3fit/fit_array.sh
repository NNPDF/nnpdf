#!/bin/bash
#SBATCH --job-name=n3fit
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=fit_bnn_%j.out
#SBATCH --error=fit_bnn_%j.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

cd /path/to/runcard

n3fit Basic_runcard_bayesian.yml 1
