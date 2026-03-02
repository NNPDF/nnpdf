#!/bin/bash
#SBATCH --job-name=vp_shapley
#SBATCH --account=inf26_pml4hep_0
#SBATCH --partition=dcgp_usr_prod
#SBATCH --time=02:00:00
#SBATCH --output=/leonardo/home/userexternal/rbonnetg/projects/nnpdf/validphys2/src/validphys/shapley/slurm/slurm-%j.out
#SBATCH --error=/leonardo/home/userexternal/rbonnetg/projects/nnpdf/validphys2/src/validphys/shapley/slurm/slurm-%j.err

set -euo pipefail

RUNCARD_ARG="${1:-${RUNCARD:-}}"
if [[ -z "${RUNCARD_ARG}" ]]; then
  echo "ERROR: Missing runcard path." >&2
  echo "Usage: sbatch ... run_vp_shapley.slurm /abs/path/to/runcard.yaml" >&2
  exit 1
fi

if [[ "${RUNCARD_ARG}" = /* && -f "${RUNCARD_ARG}" ]]; then
  RUNCARD_PATH="${RUNCARD_ARG}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/${RUNCARD_ARG}" ]]; then
  RUNCARD_PATH="${SLURM_SUBMIT_DIR}/${RUNCARD_ARG}"
elif [[ -f "${PWD}/${RUNCARD_ARG}" ]]; then
  RUNCARD_PATH="${PWD}/${RUNCARD_ARG}"
else
  echo "ERROR: runcard not found: ${RUNCARD_ARG}" >&2
  if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    echo "Tried: ${SLURM_SUBMIT_DIR}/${RUNCARD_ARG}" >&2
  fi
  echo "Tried: ${PWD}/${RUNCARD_ARG}" >&2
  exit 1
fi

if [[ ! -f "${RUNCARD_PATH}" ]]; then
  echo "ERROR: runcard not found: ${RUNCARD_PATH}" >&2
  exit 1
fi

RUNCARD_PATH="$(cd "$(dirname "${RUNCARD_PATH}")" && pwd)/$(basename "${RUNCARD_PATH}")"
if [[ "${RUNCARD_PATH}" != *"/shapley/runcards/"* ]]; then
  echo "ERROR: Runcard path must be under .../validphys/shapley/runcards/..." >&2
  echo "Got: ${RUNCARD_PATH}" >&2
  exit 1
fi

SHAPLEY_DIR="${RUNCARD_PATH%/runcards/*}"

extract_yaml_key() {
  local key="$1"
  sed -nE "s/^[[:space:]]*${key}:[[:space:]]*([^#]+).*/\\1/p" "${RUNCARD_PATH}" | head -n 1 | xargs
}

RUN_LABEL="$(extract_yaml_key run_name)"
LOG_ROOT="$(extract_yaml_key log_root)"
OUTPUT_DIR_CFG="$(extract_yaml_key output_dir)"
OUTPUT_ROOT_CFG="$(extract_yaml_key output_root)"

if [[ -z "${RUN_LABEL}" ]]; then
  echo "ERROR: Missing required top-level key 'run_name' in ${RUNCARD_PATH}" >&2
  exit 1
fi
if [[ -z "${LOG_ROOT}" ]]; then
  echo "ERROR: Missing required top-level key 'log_root' in ${RUNCARD_PATH}" >&2
  exit 1
fi
if [[ -z "${OUTPUT_DIR_CFG}" && -z "${OUTPUT_ROOT_CFG}" ]]; then
  echo "ERROR: Missing required top-level key 'output_root' or 'output_dir' in ${RUNCARD_PATH}" >&2
  exit 1
fi

safe_token() {
  echo "$1" | sed -E 's/\./p/g; s/[^A-Za-z0-9_-]+/_/g; s/^_+//; s/_+$//'
}

RUNCARD_FILE="$(basename "${RUNCARD_PATH}")"
RUNCARD_STEM="$(safe_token "${RUNCARD_FILE%.*}")"

if [[ "${LOG_ROOT}" = /* ]]; then
  LOG_DIR="${LOG_ROOT}/$(safe_token "${RUN_LABEL}")"
else
  LOG_DIR="${SHAPLEY_DIR}/${LOG_ROOT}/$(safe_token "${RUN_LABEL}")"
fi

if [[ -n "${OUTPUT_DIR_CFG}" ]]; then
  if [[ "${OUTPUT_DIR_CFG}" = /* ]]; then
    RESULTS_BASE_DIR="${OUTPUT_DIR_CFG}"
  else
    RESULTS_BASE_DIR="${SHAPLEY_DIR}/${OUTPUT_DIR_CFG}"
  fi
else
  if [[ "${OUTPUT_ROOT_CFG}" = /* ]]; then
    RESULTS_BASE_DIR="${OUTPUT_ROOT_CFG}"
  else
    RESULTS_BASE_DIR="${SHAPLEY_DIR}/${OUTPUT_ROOT_CFG}"
  fi
fi

REPO_DIR="$(cd "${SHAPLEY_DIR}/../../../../" && pwd)"

# Ensure the runtime environment provides vp-shapley on compute nodes.
CONDA_HOME="${CONDA_HOME:-/leonardo/home/userexternal/rbonnetg/miniconda3}"
CONDA_ENV="${CONDA_ENV:-environment_nnpdf}"
if [[ ! -f "${CONDA_HOME}/etc/profile.d/conda.sh" ]]; then
  echo "ERROR: conda.sh not found at ${CONDA_HOME}/etc/profile.d/conda.sh" >&2
  exit 1
fi
source "${CONDA_HOME}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

mkdir -p "${LOG_DIR}" "${RESULTS_BASE_DIR}"

TS="$(date +%Y%m%d_%H%M%S)"
LOG_OUT="${LOG_DIR}/${TS}_${RUNCARD_STEM}_job${SLURM_JOB_ID:-manual}.out.log"
LOG_ERR="${LOG_DIR}/${TS}_${RUNCARD_STEM}_job${SLURM_JOB_ID:-manual}.err.log"
exec > >(tee -a "${LOG_OUT}")
exec 2> >(tee -a "${LOG_ERR}" >&2)

echo "Repository      : ${REPO_DIR}"
echo "Shapley dir     : ${SHAPLEY_DIR}"
echo "Runcard         : ${RUNCARD_PATH}"
echo "Run label       : ${RUN_LABEL}"
echo "Results base dir: ${RESULTS_BASE_DIR}"
echo "Log out         : ${LOG_OUT}"
echo "Log err         : ${LOG_ERR}"
echo "Conda env       : ${CONDA_ENV}"

cd "${REPO_DIR}"

CMD=(vp-shapley "${RUNCARD_PATH}")
if [[ -n "${OUTPUT_DIR:-}" ]]; then
  if [[ "${OUTPUT_DIR}" = /* ]]; then
    EFFECTIVE_OUTPUT_DIR="${OUTPUT_DIR}"
  else
    EFFECTIVE_OUTPUT_DIR="${SHAPLEY_DIR}/${OUTPUT_DIR}"
  fi
  mkdir -p "${EFFECTIVE_OUTPUT_DIR}"
  CMD+=(--output "${EFFECTIVE_OUTPUT_DIR}")
fi
if [[ -n "${N_JOBS:-}" ]]; then
  CMD+=(--n-jobs "${N_JOBS}")
fi
if [[ $# -gt 1 ]]; then
  CMD+=("${@:2}")
fi

echo "Running    : ${CMD[*]}"
if [[ "${DRY_RUN:-0}" = "1" ]]; then
  echo "DRY_RUN=1 -> stopping before vp-shapley execution."
  exit 0
fi
"${CMD[@]}"
