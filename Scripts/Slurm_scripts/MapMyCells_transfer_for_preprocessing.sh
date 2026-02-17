#!/bin/bash
#SBATCH --job-name=MapMyCells_transfer
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --output=logs/mapmycells/mapmycells_transfer_%A.out
#SBATCH --error=logs/mapmycells/mapmycells_transfer_%A.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ==============================================================================
# MapMyCells Label Transfer for Preprocessing Objects
# ==============================================================================
# Transfers MapMyCells hierarchical cell-type annotations to:
#   1. Per-sample DecontX-corrected Seurat objects (R)
#   2. Integrated CHOIR Seurat object (R)
#   3. Integrated CHOIR AnnData object (Python)
# ==============================================================================

set -euo pipefail

echo "=============================================================================="
echo "  MapMyCells Label Transfer - SLURM Job"
echo "=============================================================================="
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Node:         $(hostname)"
echo "Start time:   $(date)"
echo "Working dir:  $(pwd)"
echo "=============================================================================="
echo ""

# --- Configuration ------------------------------------------------------------
BASE_DIR="/scicore/home/doetsch/kaiser0001/Single_cell_paper/Datasets/Human_Covid_LV_ChP_PMID_34153974"
SCRIPT_DIR="${BASE_DIR}/Scripts"
CONDA_BASE="/scicore/home/doetsch/kaiser0001/miniforge3"

# Ensure log directory exists
mkdir -p "${BASE_DIR}/logs/mapmycells"

# --- Activate conda -----------------------------------------------------------
echo "[SETUP] Initializing conda..."
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# ==============================================================================
# PART 1: R Script - Seurat objects (DecontX + CHOIR .rds)
# ==============================================================================
echo ""
echo "=============================================================================="
echo "  PART 1: Running R label transfer script"
echo "=============================================================================="
echo ""

conda activate R_4_5
echo "Conda env: $(conda info --envs | grep '*' | awk '{print $1}')"
echo "R version: $(R --version | head -1)"
echo ""

Rscript "${SCRIPT_DIR}/R_scripts/transfer_mapmycells_labels.R"

R_EXIT=$?
echo ""
echo "R script exit code: ${R_EXIT}"

if [[ ${R_EXIT} -ne 0 ]]; then
    echo "ERROR: R script failed with exit code ${R_EXIT}"
    echo "Continuing to Python script..."
fi

# ==============================================================================
# PART 2: Python Script - AnnData object (CHOIR .h5ad)
# ==============================================================================
echo ""
echo "=============================================================================="
echo "  PART 2: Running Python label transfer script"
echo "=============================================================================="
echo ""

# R_4_5 env should have scanpy/anndata; if not, adjust env name
echo "Conda env: $(conda info --envs | grep '*' | awk '{print $1}')"
echo "Python version: $(python --version 2>&1)"
echo ""

python "${SCRIPT_DIR}/Python_scripts/transfer_mapmycells_labels.py"

PY_EXIT=$?
echo ""
echo "Python script exit code: ${PY_EXIT}"

if [[ ${PY_EXIT} -ne 0 ]]; then
    echo "ERROR: Python script failed with exit code ${PY_EXIT}"
fi

# ==============================================================================
# Final Summary
# ==============================================================================
echo ""
echo "=============================================================================="
echo "  FINAL STATUS"
echo "=============================================================================="
echo "  R script:      $([ ${R_EXIT} -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo "  Python script: $([ ${PY_EXIT} -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo "  End time:      $(date)"
echo "  Elapsed:       ${SECONDS} seconds"
echo "=============================================================================="

# Exit with failure if either script failed
if [[ ${R_EXIT} -ne 0 || ${PY_EXIT} -ne 0 ]]; then
    exit 1
fi

exit 0
