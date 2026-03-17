#!/bin/bash
#SBATCH --job-name=vp_compare_mc_drop
#SBATCH --account=inf26_pml4hep_0
#SBATCH --partition=dcgp_usr_prod
#SBATCH --qos=normal
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=slurm/mc_dropout_compare-%j.out
#SBATCH --error=slurm/mc_dropout_compare-%j.err

# vp-comparefits: MC Dropout fit vs NNPDF4.0 baseline fit
#
# Run after nnpdf40_mc_dropout_assemble_evolve.slurm completes, e.g.:
#   sbatch --dependency=afterok:<ASSEMBLE_JOB_ID> this_script.slurm
#
# Overridable via --export  (all are optional):
#   MC_FIT      MC Dropout fit name   (default: nnpdf40-like-dropout-mc)
#   REF_FIT     reference fit name    (default: 240701-02-rs-nnpdf40-baseline)
#   OUT_DIR     output directory      (default: output/vp_compare_mc_vs_baseline)

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate environment_nnpdf

cd /leonardo/home/userexternal/rbonnetg/projects/nnpdf_bayesian

: "${MC_FIT:=nnpdf40-like-dropout-mc}"
: "${REF_FIT:=240701-02-rs-nnpdf40-baseline}"
: "${OUT_DIR:=output/vp_compare_mc_vs_baseline}"

echo "======================================================"
echo "vp-comparefits: MC Dropout vs NNPDF4.0 baseline"
echo "======================================================"
echo "MC fit     : ${MC_FIT}"
echo "Reference  : ${REF_FIT}"
echo "Output     : ${OUT_DIR}"
echo "Node       : $(hostname)"
echo "Started    : $(date)"
echo

# Symlink fit into NNPDF results path so vp-comparefits can find it
RESULTS_PATH=$(python -c "from validphys.loader import Loader; print(Loader().resultspath)")
ln -sfn "$(pwd)/${MC_FIT}" "${RESULTS_PATH}/${MC_FIT}"
echo "Symlinked ${MC_FIT} -> ${RESULTS_PATH}/${MC_FIT}"

# Symlink postfit LHAPDF grid so vp-comparefits can load the PDF
# (vp-comparefits auto-detection resolves to nnfit/, which doesn't exist for assembled fits)
LHAPDF_PATH=$(python -c "import lhapdf; print(lhapdf.paths()[0])")
ln -sfn "$(pwd)/${MC_FIT}/postfit/${MC_FIT}" "${LHAPDF_PATH}/${MC_FIT}"
echo "Symlinked LHAPDF grid -> ${LHAPDF_PATH}/${MC_FIT}"

# Build title and label from metadata if available
META="${MC_FIT}/mc_dropout_metadata.yaml"
if [[ -f "${META}" ]]; then
  META_MODE=$(grep  "^mode:"                "${META}" | awk '{print $2}')
  META_TREP=$(grep  "^n_trained_replicas:"  "${META}" | awk '{print $2}')
  META_SAMP=$(grep  "^n_samples_per_replica:" "${META}" | awk '{print $2}')
  META_MEMB=$(grep  "^n_total_members:"     "${META}" | awk '{print $2}')
  # fall back to legacy key
  [[ -z "${META_MEMB}" ]] && META_MEMB=$(grep "^n_replicas_assembled:" "${META}" | awk '{print $2}')
  if [[ "${META_MODE}" == "samples" ]]; then
    REPORT_TITLE="NNPDF40 MC Dropout vs NNPDF4.0 baseline | ${META_TREP} trained rep x ${META_SAMP} MC samples (${META_MEMB} PDF members)"
    CURRENT_LABEL="MC Dropout uncertainty (${META_TREP} rep x ${META_SAMP} samples)"
  else
    REPORT_TITLE="NNPDF40 MC Dropout vs NNPDF4.0 baseline | ${META_TREP} rep x ${META_SAMP} MC samples"
    CURRENT_LABEL="MC Dropout mean (${META_TREP} rep x ${META_SAMP} samples)"
  fi
else
  REPORT_TITLE="NNPDF40 MC Dropout vs NNPDF4.0 baseline"
  CURRENT_LABEL="MC Dropout (${MC_FIT})"
fi
echo "Title      : ${REPORT_TITLE}"
echo "Fit label  : ${CURRENT_LABEL}"

vp-comparefits "${MC_FIT}" "${REF_FIT}" \
  --title       "${REPORT_TITLE}" \
  --author      "rbonnetg" \
  --keywords    nnpdf40 dropout mc-dropout comparison \
  --current_fit_label   "${CURRENT_LABEL}" \
  --reference_fit_label "NNPDF4.0 baseline (no dropout)" \
  -o "${OUT_DIR}"

echo
echo "Done: $(date)"
echo "Report written to: ${OUT_DIR}"
