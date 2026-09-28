#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=s1x1000
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=setup_%j.out
#SBATCH --error=setup_%j.err
#SBATCH --export=NONE

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

cd /home/dakshanshchawda/nnpdf/n3fit/runcards/examples

vp-setupfit bay_1x1000.yml
