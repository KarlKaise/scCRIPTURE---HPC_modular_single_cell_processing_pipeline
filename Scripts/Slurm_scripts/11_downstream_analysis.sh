#!/bin/bash
#SBATCH --job-name=downstream
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=64G
#SBATCH --time=24:00:00
#SBATCH --partition=a100
#SBATCH --qos=gpu1day
#SBATCH --gres=gpu:1
#SBATCH --output=logs/downstream/downstream_%j.out
#SBATCH --error=logs/downstream/downstream_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ==============================================================================
# STEP 10-20: DOWNSTREAM ANALYSIS PIPELINE (GPU VERSION)
# ==============================================================================
#
# ENVIRONMENT ARCHITECTURE (2026-01-14):
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This pipeline has TWO different CUDA consumers with DIFFERENT requirements:
#
#   1. PYTHON (scVI, PyTorch, etc.):
#      - Needs: CUDA module 12.1.1 + conda libs in LD_LIBRARY_PATH
#      - Reason: Python packages were compiled against CUDA 12.x
#
#   2. JULIA (CUDA.jl for scICE/scLENS):
#      - Needs: Clean environment WITHOUT CUDA module paths
#      - Reason: Julia downloads its own CUDA 13.0 runtime artifacts
#      - If Julia sees CUDA 12.1.1 libs, it gets version mismatch warnings/errors
#
# SOLUTION:
# - Shell script loads CUDA module (for Python modules)
# - R module 06_scice_subclustering.R SANITIZES the environment for Julia:
#   - Strips CUDA module paths from LD_LIBRARY_PATH
#   - Strips conda paths from LD_LIBRARY_PATH  
#   - Unsets LD_PRELOAD (conflicts with Julia's libstdc++)
#
# This way:
# - Modules 0-5, 7-11: Use Python with CUDA 12.1.1 ✓
# - Module 6 (scICE): Julia runs with its own CUDA 13.0 ✓
#
# ==============================================================================

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

if [[ -n "${PROJECT_ROOT:-}" ]]; then
    BASE_DIR="${PROJECT_ROOT}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    BASE_DIR="${SLURM_SUBMIT_DIR}"
else
    BASE_DIR="/scicore/home/doetsch/kaiser0001/Revision_NatureComm_Sex/Vandebroucke_fibroblast_paper"
fi

SAMPLESHEET="${SAMPLESHEET:-${BASE_DIR}/samplesheet.csv}"

get_dataset_name() {
    tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f12
}

if [[ -z "${PREPROCESS_DIR:-}" ]]; then
    DATASET_NAME=$(get_dataset_name)
    if [[ -n "$DATASET_NAME" ]]; then
        PREPROCESS_DIR="${BASE_DIR}/Output_dir_${DATASET_NAME}/Single_cell_preprocessed"
        DOWNSTREAM_DIR="${BASE_DIR}/Output_dir_${DATASET_NAME}/Single_cell_clustering"
    else
        PREPROCESS_DIR="${BASE_DIR}"
        DOWNSTREAM_DIR="${BASE_DIR}"
    fi
fi
if [[ -z "${DOWNSTREAM_DIR:-}" ]]; then
    DOWNSTREAM_DIR="${BASE_DIR}/Output_dir_${DATASET_NAME}/Single_cell_clustering"
fi

PIPELINE_DIR="${BASE_DIR}/Scripts/scrnaseq_pipeline"
LOG_DIR="${BASE_DIR}/logs/downstream"

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

error_exit() {
    log_msg "ERROR: $1"
    exit 1
}

get_col_index() {
    local file=$1
    local col_name=$2
    head -1 "$file" | tr ',' '\n' | grep -n "^${col_name}$" | cut -d: -f1
}

get_samples_by_ventricle() {
    local ventricle=$1
    local col_id=$(get_col_index "$SAMPLESHEET" "sample_id")
    local col_vent=$(get_col_index "$SAMPLESHEET" "ventricle")

    if [[ "$ventricle" == "ALL" ]]; then
        tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | cut -d',' -f"$col_id"
    else
        tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | awk -F',' -v col="$col_vent" -v vent="$ventricle" -v id="$col_id" '$col == vent {print $id}'
    fi
}

get_samples_by_group() {
    local group_id=$1
    local col_id=$(get_col_index "$SAMPLESHEET" "sample_id")
    local col_gid=$(get_col_index "$SAMPLESHEET" "group_id")

    if [[ -z "$col_gid" ]]; then
        tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | cut -d',' -f"$col_id"
    else
        tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | awk -F',' -v col="$col_gid" -v gid="$group_id" -v id="$col_id" '$col == gid {print $id}'
    fi
}

get_group_label_by_id() {
    local group_id=$1
    local col_gid=$(get_col_index "$SAMPLESHEET" "group_id")
    local col_glabel=$(get_col_index "$SAMPLESHEET" "group_label")

    if [[ -z "$col_gid" || -z "$col_glabel" ]]; then
        echo ""
        return
    fi

    tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | awk -F',' -v col="$col_gid" -v gid="$group_id" -v lc="$col_glabel" '$col == gid {print $lc; exit}'
}

get_group_id_by_label() {
    local group_label=$1
    local col_gid=$(get_col_index "$SAMPLESHEET" "group_id")
    local col_glabel=$(get_col_index "$SAMPLESHEET" "group_label")

    if [[ -z "$col_gid" || -z "$col_glabel" ]]; then
        echo ""
        return
    fi

    tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | awk -F',' -v col="$col_glabel" -v label="$group_label" -v idc="$col_gid" '$col == label {print $idc; exit}'
}

get_ventricle_by_group() {
    local group_id=$1
    local col_gid=$(get_col_index "$SAMPLESHEET" "group_id")
    local col_vent=$(get_col_index "$SAMPLESHEET" "ventricle")

    if [[ -z "$col_gid" || -z "$col_vent" ]]; then
        echo ""
        return
    fi

    tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | awk -F',' -v col="$col_gid" -v gid="$group_id" -v vc="$col_vent" '$col == gid {print $vc; exit}'
}

# ==============================================================================
# Parse Arguments
# ==============================================================================

START_MODULE=0
STOP_MODULE=11
MODULES_ARG=""
STEPS_ARG=""
VENTRICLE=""
GROUP_LABEL_ARG=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --start) START_MODULE="$2"; shift 2 ;;
        --stop|--to) STOP_MODULE="$2"; shift 2 ;;
        --modules) MODULES_ARG="$2"; shift 2 ;;
        --steps) STEPS_ARG="$2"; shift 2 ;;
        --ventricle) VENTRICLE="$2"; shift 2 ;;
        --group-label) GROUP_LABEL_ARG="$2"; shift 2 ;;
        *) shift ;;
    esac
done

# ==============================================================================
# Determine Output Label and Filtering
# ==============================================================================

ENV_GROUP_ID="${GROUP_ID:-}"
ENV_GROUP_LABEL="${GROUP_LABEL:-}"
ENV_VENTRICLE_FILTER="${VENTRICLE_FILTER:-}"

OUTPUT_LABEL=""
ACTIVE_GROUP_ID=""
ACTIVE_GROUP_LABEL=""
ACTIVE_VENTRICLE=""

if [[ -n "$GROUP_LABEL_ARG" ]]; then
    OUTPUT_LABEL="$GROUP_LABEL_ARG"
    ACTIVE_GROUP_LABEL="$GROUP_LABEL_ARG"
    ACTIVE_GROUP_ID=$(get_group_id_by_label "$GROUP_LABEL_ARG")
    log_msg "Using command line --group-label: $GROUP_LABEL_ARG (GROUP_ID: ${ACTIVE_GROUP_ID:-<not found>})"
elif [[ -n "$ENV_GROUP_LABEL" ]]; then
    OUTPUT_LABEL="$ENV_GROUP_LABEL"
    ACTIVE_GROUP_ID="$ENV_GROUP_ID"
    ACTIVE_GROUP_LABEL="$ENV_GROUP_LABEL"
    log_msg "Using environment GROUP_LABEL: $ENV_GROUP_LABEL (GROUP_ID: $ENV_GROUP_ID)"
elif [[ -n "$VENTRICLE" ]]; then
    OUTPUT_LABEL="$VENTRICLE"
    ACTIVE_VENTRICLE="$VENTRICLE"
    log_msg "Using command line --ventricle: $VENTRICLE"
elif [[ -n "$ENV_VENTRICLE_FILTER" ]]; then
    OUTPUT_LABEL="$ENV_VENTRICLE_FILTER"
    ACTIVE_VENTRICLE="$ENV_VENTRICLE_FILTER"
    log_msg "Using environment VENTRICLE_FILTER: $ENV_VENTRICLE_FILTER"
else
    error_exit "No filter specified. Use --group-label <label>, --ventricle <LV|4V|ALL>, or set GROUP_LABEL/VENTRICLE_FILTER environment variable"
fi

if [[ -z "$VENTRICLE" && -n "$ACTIVE_VENTRICLE" ]]; then
    VENTRICLE="$ACTIVE_VENTRICLE"
fi

if [[ -n "$ACTIVE_GROUP_ID" && -z "$VENTRICLE" ]]; then
    VENTRICLE=$(get_ventricle_by_group "$ACTIVE_GROUP_ID")
fi

export VENTRICLE_FILTER="${VENTRICLE:-}"
export GROUP_ID="${ACTIVE_GROUP_ID:-$ENV_GROUP_ID}"
export GROUP_LABEL="${ACTIVE_GROUP_LABEL:-$ENV_GROUP_LABEL}"
export SAMPLESHEET="${SAMPLESHEET}"
export PROJECT_ROOT="${BASE_DIR}"

OUTPUT_DIR="${DOWNSTREAM_DIR}/10_Downstream_Analysis_${OUTPUT_LABEL}"

# ==============================================================================
# Translate --steps (10..20) into --modules (0..11)
# ==============================================================================

if [[ -n "${STEPS_ARG}" ]]; then
    STEPS_CLEAN="$(echo "${STEPS_ARG}" | tr -d '[:space:]')"
    IFS=',' read -r -a STEPS_ARR <<< "${STEPS_CLEAN}"

    MODS_LIST=()
    for s in "${STEPS_ARR[@]}"; do
        if [[ ! "${s}" =~ ^[0-9]+$ ]]; then
            error_exit "--steps contains non-numeric: '${s}'"
        fi
        if [[ "${s}" -eq 20 ]]; then
            MODS_LIST+=("10" "11")
        else
            if [[ "${s}" -lt 10 || "${s}" -gt 19 ]]; then
                error_exit "--steps values must be 10..20 (got ${s})"
            fi
            MODS_LIST+=("$((s - 10))")
        fi
    done

    MODULES_ARG="$(printf "%s\n" "${MODS_LIST[@]}" | awk '!seen[$0]++' | sort -n | paste -sd, -)"
fi

# ==============================================================================
# INITIALIZE
# ==============================================================================

log_msg "============================================================================"
log_msg "DOWNSTREAM ANALYSIS PIPELINE - GPU Version (A100)"
log_msg "============================================================================"
log_msg "Job ID:        ${SLURM_JOB_ID:-local}"
log_msg "Date:          $(date)"
log_msg "Host:          ${SLURMD_NODENAME:-$(hostname)}"
log_msg "Partition:     ${SLURM_JOB_PARTITION:-unknown}"
log_msg "Base Dir:      ${BASE_DIR}"
log_msg "============================================================================"
echo ""

if [[ -n "$ACTIVE_GROUP_ID" || -n "$ACTIVE_GROUP_LABEL" ]]; then
    log_msg ">>> GROUP-BASED ANALYSIS <<<"
    log_msg "    Group ID:    ${ACTIVE_GROUP_ID:-<not set>}"
    log_msg "    Group Label: ${ACTIVE_GROUP_LABEL}"
else
    log_msg ">>> VENTRICLE-BASED ANALYSIS: ${VENTRICLE} <<<"
fi
log_msg ">>> OUTPUT: ${OUTPUT_DIR} <<<"
echo ""

if [[ ! -f "$SAMPLESHEET" ]]; then
    error_exit "Samplesheet not found: $SAMPLESHEET"
fi

log_msg "Samplesheet: $SAMPLESHEET"

if [[ -n "$ACTIVE_GROUP_ID" ]]; then
    SAMPLES=($(get_samples_by_group "${ACTIVE_GROUP_ID}"))
    log_msg "Samples for group ${ACTIVE_GROUP_ID} (${ACTIVE_GROUP_LABEL}):"
elif [[ -n "$VENTRICLE" ]]; then
    SAMPLES=($(get_samples_by_ventricle "${VENTRICLE}"))
    log_msg "Samples for ${VENTRICLE} ventricle:"
else
    error_exit "No valid filter determined"
fi

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    error_exit "No samples found"
fi

for s in "${SAMPLES[@]}"; do
    log_msg "  - $s"
done
echo ""

mkdir -p "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

# ==============================================================================
# GPU DETECTION
# ==============================================================================

log_msg "============================================================================"
log_msg "GPU CONFIGURATION"
log_msg "============================================================================"

if command -v nvidia-smi &> /dev/null && nvidia-smi &> /dev/null; then
    log_msg "GPU DETECTED - Using GPU acceleration"
    export USE_GPU=true
    export CUDA_VISIBLE_DEVICES=0
    nvidia-smi --query-gpu=index,name,memory.total,memory.free,driver_version --format=csv,noheader 2>/dev/null | while read line; do
        log_msg "  $line"
    done
    echo ""
    nvidia-smi
else
    log_msg "WARNING: GPU not detected!"
    export USE_GPU=false
    export CUDA_VISIBLE_DEVICES=""
fi

echo ""

# ==============================================================================
# Environment Setup
# ==============================================================================

log_msg "Setting up environment..."

# ==============================================================================
# Load CUDA module for Python packages (scVI, PyTorch, etc.)
# ==============================================================================
# This is needed for Python-based GPU modules.
# Julia (Module 6) will IGNORE these paths - the R script sanitizes the
# environment before launching Julia subprocesses.
# ==============================================================================

log_msg "Loading CUDA module for Python packages..."
module purge 2>/dev/null || true
module load CUDA/12.1.1 2>/dev/null && log_msg "  CUDA/12.1.1 loaded" || log_msg "  CUDA module not available"
module list 2>&1 | grep -v "^$" | while read line; do log_msg "  $line"; done
echo ""

# Activate conda AFTER loading modules (to ensure conda paths come first)
source /scicore/home/doetsch/kaiser0001/miniforge3/etc/profile.d/conda.sh
conda activate /scicore/home/doetsch/kaiser0001/miniforge3/envs/R_4_5

log_msg "Conda environment: $CONDA_PREFIX"
log_msg "R version: $(R --version | head -1)"
log_msg "Python: $(python --version 2>&1)"
echo ""

# ==============================================================================
# Library path setup
# ==============================================================================
# Conda libs FIRST (for scipy/numba/llvmlite compatibility)
# CUDA module paths are already in LD_LIBRARY_PATH from module load
#
# IMPORTANT for Module 6 (scICE):
# The R script 06_scice_subclustering.R will create a SANITIZED environment
# for Julia subprocesses that:
#   - Removes CUDA module paths (Julia uses its own CUDA 13.0)
#   - Removes conda paths (breaks Julia JLL artifacts)
#   - Unsets LD_PRELOAD (conflicts with Julia's libstdc++)
# ==============================================================================

# Prepend conda lib to ensure it takes priority
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# LD_PRELOAD for R/Python compatibility (scipy, etc.)
export LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6:$CONDA_PREFIX/lib/libgcc_s.so.1"

log_msg "Library configuration:"
log_msg "  LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
log_msg "  LD_PRELOAD: $LD_PRELOAD"
log_msg ""
log_msg "NOTE: Module 6 (scICE/Julia) will use a sanitized environment:"
log_msg "  - Julia LD_LIBRARY_PATH: [no CUDA module, no conda paths]"
log_msg "  - Julia LD_PRELOAD: [unset]"
echo ""

# ==============================================================================
# Julia configuration (exported for R to use)
# ==============================================================================
# These variables tell the R script where Julia lives.
# The R script will create appropriate environment for Julia subprocesses.
# ==============================================================================

export JULIA_DEPOT_PATH=/scicore/home/doetsch/kaiser0001/.julia
export JULIA_PROJECT=/scicore/home/doetsch/kaiser0001/julia/julia_envs/scICE_env/scICE
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-6}
export JULIA_CUDA_MEMORY_POOL=none
export JULIA_PKG_OFFLINE=true

log_msg "Julia configuration (for Module 6):"
log_msg "  JULIA_DEPOT_PATH: $JULIA_DEPOT_PATH"
log_msg "  JULIA_PROJECT: $JULIA_PROJECT"
log_msg "  JULIA_NUM_THREADS: $JULIA_NUM_THREADS"
log_msg "  JULIA_PKG_OFFLINE: $JULIA_PKG_OFFLINE"
log_msg "  (Julia will use its own CUDA 13.0 runtime, ignoring CUDA module)"
echo ""

# Thread settings
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-6}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-6}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-6}
export R_PARALLELLY_AVAILABLE_CORES=${SLURM_CPUS_PER_TASK:-6}

# Python for reticulate (used by some R packages)
export RETICULATE_PYTHON=/scicore/home/doetsch/kaiser0001/miniforge3/envs/kaiser_test_py3.11/bin/python

# Pipeline directory export
export PIPELINE_DIR="${PIPELINE_DIR}"

log_msg "Thread/resource configuration:"
log_msg "  CPUs: ${SLURM_CPUS_PER_TASK:-6}"
log_msg "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
log_msg "  RETICULATE_PYTHON: $RETICULATE_PYTHON"
log_msg "  USE_GPU: $USE_GPU"
echo ""

# ==============================================================================
# Verify Pipeline
# ==============================================================================

if [[ ! -d "${PIPELINE_DIR}" ]]; then
    error_exit "Pipeline directory not found: ${PIPELINE_DIR}"
fi

if [[ ! -f "${PIPELINE_DIR}/run_pipeline.R" ]]; then
    error_exit "run_pipeline.R not found"
fi

log_msg "Pipeline directory: ${PIPELINE_DIR}"
echo ""

# ==============================================================================
# Run Pipeline
# ==============================================================================

log_msg "============================================================================"
if [[ -n "$ACTIVE_GROUP_LABEL" ]]; then
    log_msg "Running Downstream Analysis - Group: ${ACTIVE_GROUP_LABEL}"
else
    log_msg "Running Downstream Analysis - ${VENTRICLE}"
fi
log_msg "============================================================================"
echo ""

RSCRIPT_ARGS="--start ${START_MODULE} --stop ${STOP_MODULE}"

if [[ -n "${MODULES_ARG}" ]]; then
    log_msg "Running selected modules: ${MODULES_ARG}"
    RSCRIPT_ARGS="--modules ${MODULES_ARG}"
else
    log_msg "Running module range: ${START_MODULE} to ${STOP_MODULE}"
fi

echo ""
log_msg "Command: Rscript ${PIPELINE_DIR}/run_pipeline.R ${RSCRIPT_ARGS}"
echo ""

cd "${PIPELINE_DIR}"
Rscript "${PIPELINE_DIR}/run_pipeline.R" ${RSCRIPT_ARGS}

EXIT_CODE=$?

# ==============================================================================
# Summary
# ==============================================================================

echo ""
log_msg "============================================================================"
log_msg "DOWNSTREAM ANALYSIS COMPLETE"
log_msg "============================================================================"
log_msg "Exit code: ${EXIT_CODE}"
log_msg "Output directory: ${OUTPUT_DIR}"
log_msg "GPU used: ${USE_GPU}"
echo ""

if [[ "${USE_GPU}" == "true" ]]; then
    log_msg "Final GPU Memory Status:"
    nvidia-smi --query-gpu=memory.used,memory.free --format=csv,noheader 2>/dev/null || true
fi

log_msg "Job completed at: $(date)"

exit ${EXIT_CODE}
