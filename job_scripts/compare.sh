#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=c1x1000
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=comp_%j.out
#SBATCH --error=comp_%j.err
#SBATCH --export=NONE

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

vp-comparefits "bay_1x1000" "nnpdf40-like" \
    --title "Comparison Report 1 Training 1000 Extractions vs NNPDF on L2 data" \
    --author "Dakshansh" \
    --keywords "Bayesian_report_40-like" \
    -o "/home/dakshanshchawda/jobs/BNN_JOBS/bay_1x1000/output_bay_1x1000_vs_NNPDF_on_L2_data"



