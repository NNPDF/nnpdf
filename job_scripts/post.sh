#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=p1x1000
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=post_%j.out
#SBATCH --error=post_%j.err
#SBATCH --export=NONE

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

cd /home/dakshanshchawda/nnpdf/n3fit/runcards/examples

evolven3fit evolve bay_1x1000

postfit 1000 bay_1x1000

cp -r bay_1x1000 /home/dakshanshchawda/miniconda3/envs/environment_nnpdf/share/NNPDF/results

vp-comparefits "bay_1x1000" "NNPDF40_nnlo_as_01180_1000" \
    --title "Comparison Report 1 Training 1000 Extractions vs full NNPDF" \
    --author "Dakshansh" \
    --keywords "Bayesian_report_40-like" \
    -o "/home/dakshanshchawda/jobs/BNN_JOBS/bay_1x1000/output_bay_1x1000_vs_full_NNPDF_sampling_fix"

vp-comparefits "bay_1x1000" "nnpdf40-like_l1_data" \
    --title "Comparison Report 1 Training 1000 Extractions vs NNPDF on L1 data" \
    --author "Dakshansh" \
    --keywords "Bayesian_report_40-like" \
    -o "/home/dakshanshchawda/jobs/BNN_JOBS/bay_1x1000/output_bay_1x1000_vs_NNPDF_on_L1_data_sampling_fix"
