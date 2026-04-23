#!/bin/bash
#SBATCH --job-name=setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=setup_%j.out
#SBATCH --error=setup_%j.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

cd /path/to/runcard

vp-setupfit Basic_runcard_bayesian.yml
