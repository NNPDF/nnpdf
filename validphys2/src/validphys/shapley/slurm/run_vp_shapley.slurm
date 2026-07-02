#!/bin/bash
#SBATCH --job-name=vp_shapley
#SBATCH --account=inf26_pml4hep_0
#SBATCH --partition=dcgp_usr_prod
#SBATCH --time=04:00:00
#SBATCH --output=/leonardo/home/userexternal/rbonnetg/projects/nnpdf/validphys2/src/validphys/shapley/slurm/logs/slurm-%j.out
#SBATCH --error=/leonardo/home/userexternal/rbonnetg/projects/nnpdf/validphys2/src/validphys/shapley/slurm/logs/slurm-%j.err

set -euo pipefail

# Optional modes
#  - --split (or SPLIT_EXPERIMENTS=1): submit an array job, one task per experiment,
#    plus a dependent aggregation job.
#  - --aggregate-only (or AGGREGATE_ONLY=1): only rebuild cross-experiment outputs
#    from per-experiment results already on disk.
SPLIT_EXPERIMENTS="${SPLIT_EXPERIMENTS:-0}"
SPLIT_MODE="${SPLIT_MODE:-0}"
AGGREGATE_ONLY="${AGGREGATE_ONLY:-0}"

if [[ "${1:-}" == "--split" ]]; then
  SPLIT_EXPERIMENTS=1
  shift
fi
if [[ "${1:-}" == "--aggregate-only" ]]; then
  AGGREGATE_ONLY=1
  shift
fi

RUNCARD_ARG="${1:-${RUNCARD:-}}"
if [[ -z "${RUNCARD_ARG}" ]]; then
  echo "ERROR: Missing runcard path." >&2
  echo "Usage:" >&2
  echo "  sbatch ... run_vp_shapley.slurm /abs/path/to/runcard.yaml" >&2
  echo "  ./run_vp_shapley.slurm --split /abs/path/to/runcard.yaml" >&2
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

SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"

# Ensure the runtime environment provides vp-shapley on compute nodes.
CONDA_HOME="${CONDA_HOME:-/leonardo/home/userexternal/rbonnetg/miniconda3}"
CONDA_ENV="${CONDA_ENV:-environment_nnpdf}"
if [[ ! -f "${CONDA_HOME}/etc/profile.d/conda.sh" ]]; then
  echo "ERROR: conda.sh not found at ${CONDA_HOME}/etc/profile.d/conda.sh" >&2
  exit 1
fi
source "${CONDA_HOME}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

resolve_effective_output_dir() {
  local out_dir="$1"
  if [[ -z "${out_dir}" ]]; then
    echo ""
    return 0
  fi
  if [[ "${out_dir}" = /* ]]; then
    echo "${out_dir}"
  else
    echo "${SHAPLEY_DIR}/${out_dir}"
  fi
}

get_experiment_count() {
  python - "$1" <<'PY'
import sys
import yaml

runcard = sys.argv[1]
with open(runcard) as f:
    cfg = yaml.safe_load(f)

exps = cfg.get("experiments")
if isinstance(exps, dict):
    n = len(list(exps.keys()))
elif isinstance(exps, list):
    n = len(exps)
else:
    n = 0

print(n)
PY
}

# Submission helper mode: run this script directly (not under SLURM)
# to submit an array job + dependent aggregation job.
if [[ -z "${SLURM_JOB_ID:-}" && "${SPLIT_EXPERIMENTS}" == "1" ]]; then
  # Extra args (after the runcard) are forwarded to the child jobs.
  EXTRA_ARGS=()
  if [[ $# -gt 1 ]]; then
    EXTRA_ARGS=("${@:2}")
  fi

  mkdir -p "${LOG_DIR}" "${RESULTS_BASE_DIR}"

  TS="$(date +%Y%m%d_%H%M%S)"
  RUN_TOKEN="$(safe_token "${RUN_LABEL}")"

  # Compute a single shared output directory for this split run.
  EFFECTIVE_OUTPUT_DIR=""
  if [[ -n "${OUTPUT_DIR:-}" ]]; then
    EFFECTIVE_OUTPUT_DIR="$(resolve_effective_output_dir "${OUTPUT_DIR}")"
  elif [[ -n "${OUTPUT_DIR_CFG}" ]]; then
    # output_dir is already a deterministic directory (no timestamp needed).
    EFFECTIVE_OUTPUT_DIR="${RESULTS_BASE_DIR}"
  else
    # output_root requires a timestamped subdir to match vp-shapley defaults.
    EFFECTIVE_OUTPUT_DIR="${RESULTS_BASE_DIR}/${TS}_${RUN_TOKEN}"
  fi
  mkdir -p "${EFFECTIVE_OUTPUT_DIR}"

  N_EXP="$(get_experiment_count "${RUNCARD_PATH}")"
  if [[ -z "${N_EXP}" || "${N_EXP}" -lt 1 ]]; then
    echo "ERROR: Could not determine experiment count from runcard." >&2
    exit 1
  fi

  ARRAY_MAX=$((N_EXP - 1))
  CPUS_PER_TASK="${CPUS_PER_TASK:-${N_JOBS:-1}}"
  if [[ -z "${CPUS_PER_TASK}" || "${CPUS_PER_TASK}" -lt 1 ]]; then
    CPUS_PER_TASK=1
  fi
  echo "Submitting split Shapley run" >&2
  echo "  Runcard    : ${RUNCARD_PATH}" >&2
  echo "  Experiments: ${N_EXP} (array 0-${ARRAY_MAX})" >&2
  echo "  Output dir : ${EFFECTIVE_OUTPUT_DIR}" >&2
  echo "  CPUs/task  : ${CPUS_PER_TASK} (set CPUS_PER_TASK or N_JOBS)" >&2

  array_submit_out="$(
    sbatch \
      --job-name="vp_shapley_${RUN_TOKEN}" \
      --array="0-${ARRAY_MAX}" \
      --cpus-per-task="${CPUS_PER_TASK}" \
      --export=ALL,RUNCARD="${RUNCARD_PATH}",OUTPUT_DIR="${EFFECTIVE_OUTPUT_DIR}",SPLIT_MODE=1 \
      "${SCRIPT_PATH}" "${RUNCARD_PATH}" "${EXTRA_ARGS[@]}"
  )"
  echo "${array_submit_out}" >&2
  ARRAY_JOB_ID="$(echo "${array_submit_out}" | awk '{print $NF}')"
  if [[ -z "${ARRAY_JOB_ID}" ]]; then
    echo "ERROR: Failed to parse array job id from: ${array_submit_out}" >&2
    exit 1
  fi

  agg_submit_out="$(
    sbatch \
      --job-name="vp_shapley_${RUN_TOKEN}_aggregate" \
      --dependency="afterok:${ARRAY_JOB_ID}" \
      --cpus-per-task="${CPUS_PER_TASK}" \
      --export=ALL,RUNCARD="${RUNCARD_PATH}",OUTPUT_DIR="${EFFECTIVE_OUTPUT_DIR}",AGGREGATE_ONLY=1 \
      "${SCRIPT_PATH}" "${RUNCARD_PATH}" "${EXTRA_ARGS[@]}"
  )"
  echo "${agg_submit_out}" >&2

  echo "Submitted array job : ${ARRAY_JOB_ID}" >&2
  echo "Submitted aggregate : (afterok:${ARRAY_JOB_ID})" >&2
  exit 0
fi

mkdir -p "${LOG_DIR}" "${RESULTS_BASE_DIR}"

TS="$(date +%Y%m%d_%H%M%S)"
TASK_SUFFIX=""
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  TASK_SUFFIX="_task${SLURM_ARRAY_TASK_ID}"
fi
LOG_OUT="${LOG_DIR}/${TS}_${RUNCARD_STEM}_job${SLURM_JOB_ID:-manual}${TASK_SUFFIX}.out.log"
LOG_ERR="${LOG_DIR}/${TS}_${RUNCARD_STEM}_job${SLURM_JOB_ID:-manual}${TASK_SUFFIX}.err.log"
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
  EFFECTIVE_OUTPUT_DIR="$(resolve_effective_output_dir "${OUTPUT_DIR}")"
  mkdir -p "${EFFECTIVE_OUTPUT_DIR}"
  CMD+=(--output "${EFFECTIVE_OUTPUT_DIR}")
fi
if [[ -n "${N_JOBS:-}" ]]; then
  CMD+=(--n-jobs "${N_JOBS}")
fi

if [[ "${SPLIT_MODE}" == "1" ]]; then
  if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SPLIT_MODE=1 requires SLURM_ARRAY_TASK_ID" >&2
    exit 1
  fi
  CMD+=(--experiment-index "${SLURM_ARRAY_TASK_ID}" --no-compare)
fi
if [[ "${AGGREGATE_ONLY}" == "1" ]]; then
  CMD+=(--aggregate-only)
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
