#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=n1x1000
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=fit_bnn_%j.out
#SBATCH --error=fit_bnn_%j.err
#SBATCH --export=NONE

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

cd /home/dakshanshchawda/nnpdf/n3fit/runcards/examples

n3fit bay_1x1000.yml 1
