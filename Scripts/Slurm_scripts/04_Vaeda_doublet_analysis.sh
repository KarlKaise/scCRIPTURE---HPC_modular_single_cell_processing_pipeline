#!/bin/bash
#SBATCH --job-name=VAEDA
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --output=logs/vaeda/vaeda_%a_%A.out
#SBATCH --error=logs/vaeda/vaeda_%a_%A.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ============================================================================
# STEP 4: VAEDA Doublet Detection (metadata only; NO removal here)
# ============================================================================
#
# UNIVERSAL MODULAR VERSION - Reads samples from samplesheet.csv
#
# Purpose:
#   Run VAEDA and STORE results in the AnnData object, but DO NOT remove doublets.
#   Downstream consensus removal happens only in Step 6 (>=2/3 voting).
#
# Input:
#   3_qClus_empty_droplets/<SAMPLE>/<SAMPLE>_qclus_dropletqc_filtered.h5ad
#
# Output:
#   4_Vaeda_doublet_detection/<SAMPLE>/<SAMPLE>_qClus_dropletqc_vaeda.h5ad
#   + summary CSV per sample
#
# IMPORTANT for Step 6 compatibility:
#   Step 6 currently looks for meta column "vaeda_prediction" in Seurat.
#   VAEDA outputs "vaeda_calls". This script writes BOTH:
#     - adata.obs["vaeda_calls"]      (native VAEDA)
#     - adata.obs["vaeda_prediction"] (compat alias for your Step 6 code)
#
# Usage:
#   sbatch --array=1-N 04_Vaeda_doublet_analysis.sh
#
# ============================================================================

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Base directory (from environment or default)
if [[ -n "${PROJECT_ROOT:-}" ]]; then
    BASE_DIR="${PROJECT_ROOT}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    BASE_DIR="${SLURM_SUBMIT_DIR}"
else
    BASE_DIR="/scicore/home/doetsch/kaiser0001/Revision_NatureComm_Sex/Vandebroucke_fibroblast_paper"
fi

# Samplesheet
SAMPLESHEET="${SAMPLESHEET:-${BASE_DIR}/samplesheet.csv}"

# Get dataset name from samplesheet (column 12)
get_dataset_name() {
    tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f12
}

# Set PREPROCESS_DIR from environment or derive from samplesheet
if [[ -z "${PREPROCESS_DIR:-}" ]]; then
    DATASET_NAME=$(get_dataset_name)
    if [[ -n "$DATASET_NAME" ]]; then
        PREPROCESS_DIR="${BASE_DIR}/Output_dir_${DATASET_NAME}/Single_cell_preprocessed"
    else
        PREPROCESS_DIR="${BASE_DIR}"
    fi
fi

# Directories
INPUT_ROOT="${PREPROCESS_DIR}/3_qClus_empty_droplets"
OUTPUT_ROOT="${PREPROCESS_DIR}/4_Vaeda_doublet_detection"
SCRIPT_DIR="${BASE_DIR}/Scripts/Python_scripts"
LOG_DIR="${BASE_DIR}/logs/vaeda"
README_FILE="${OUTPUT_ROOT}/README.txt"

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

get_unique_samples() {
    local col=$(get_col_index "$SAMPLESHEET" "sample_id")
    tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | cut -d',' -f"$col" | sort -u
}

get_sample_count() {
    get_unique_samples | wc -l
}

get_sample_by_index() {
    local idx=$1
    get_unique_samples | sed -n "${idx}p"
}

get_sample_field() {
    local sample=$1
    local field=$2
    local col=$(get_col_index "$SAMPLESHEET" "$field")
    if [[ -n "$col" ]]; then
        grep "^${sample}," "$SAMPLESHEET" | head -1 | cut -d',' -f"$col"
    fi
}

# ==============================================================================
# INITIALIZE
# ==============================================================================

log_msg "============================================================================"
log_msg "VAEDA Doublet Detection - Modular Pipeline"
log_msg "============================================================================"
log_msg "Job ID:        ${SLURM_JOB_ID:-local}"
log_msg "Array Task:    ${SLURM_ARRAY_TASK_ID:-1}"
log_msg "Date:          $(date)"
log_msg "Host:          $(hostname)"
log_msg "Base Dir:      ${BASE_DIR}"
log_msg "============================================================================"
echo ""

# Validate samplesheet
if [[ ! -f "$SAMPLESHEET" ]]; then
    error_exit "Samplesheet not found: $SAMPLESHEET"
fi

log_msg "Samplesheet: $SAMPLESHEET"

# Get sample
N_SAMPLES=$(get_sample_count)
TASK_IDX=${SLURM_ARRAY_TASK_ID:-1}

if [[ $TASK_IDX -gt $N_SAMPLES ]]; then
    error_exit "Array task $TASK_IDX exceeds sample count $N_SAMPLES"
fi

SAMPLE=$(get_sample_by_index $TASK_IDX)

if [[ -z "$SAMPLE" ]]; then
    error_exit "Could not determine sample for task index $TASK_IDX"
fi

log_msg "Processing: $SAMPLE (task $TASK_IDX of $N_SAMPLES)"

# Get sample metadata
SAMPLE_SEX=$(get_sample_field "$SAMPLE" "sex")
SAMPLE_BATCH=$(get_sample_field "$SAMPLE" "batch")
SAMPLE_VENTRICLE=$(get_sample_field "$SAMPLE" "ventricle")

log_msg "Sample metadata: Sex=${SAMPLE_SEX:-N/A}, Batch=${SAMPLE_BATCH:-N/A}, Ventricle=${SAMPLE_VENTRICLE:-N/A}"
echo ""

# Paths
INPUT_H5AD="${INPUT_ROOT}/${SAMPLE}/${SAMPLE}_qclus_dropletqc_filtered.h5ad"
SAMPLE_OUTPUT_DIR="${OUTPUT_ROOT}/${SAMPLE}"
OUTPUT_H5AD="${SAMPLE_OUTPUT_DIR}/${SAMPLE}_qClus_dropletqc_vaeda.h5ad"

START_TIME=$(date '+%H:%M:%S')

# Create directories
mkdir -p "${SAMPLE_OUTPUT_DIR}"
mkdir -p "${SCRIPT_DIR}"
mkdir -p "${LOG_DIR}"

# ==============================================================================
# INITIALIZE README (first task only)
# ==============================================================================

if [[ ${SLURM_ARRAY_TASK_ID:-1} -eq 1 ]]; then
    cat > "${README_FILE}" << EOF
================================================================================
STEP 4: VAEDA Doublet Detection (metadata only; no removal)
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')
Pipeline: Universal Modular scRNA-seq Pipeline

DESCRIPTION:
  Runs VAEDA (Variational Autoencoder Doublet Analysis) for doublet detection
  on cells that passed QClus + DropletQC filtering.
  This step DOES NOT remove doublets; it only writes metadata used later for
  consensus voting in Step 6.

SAMPLESHEET: ${SAMPLESHEET}
INPUT: 3_qClus_empty_droplets/<SAMPLE>/<SAMPLE>_qclus_dropletqc_filtered.h5ad

SAMPLES PROCESSED:
EOF
fi

# ==============================================================================
# VERIFY INPUT
# ==============================================================================

if [[ ! -f "${INPUT_H5AD}" ]]; then
    error_exit "Input file not found: ${INPUT_H5AD}"
fi

INPUT_SIZE=$(du -h "${INPUT_H5AD}" | cut -f1)
log_msg "Input: ${INPUT_H5AD} (${INPUT_SIZE})"
log_msg "Output: ${OUTPUT_H5AD}"
echo ""

# ==============================================================================
# ACTIVATE CONDA ENVIRONMENT
# ==============================================================================

log_msg "Activating VAEDA conda environment..."

source /scicore/home/doetsch/kaiser0001/miniforge3/etc/profile.d/conda.sh
conda activate /scicore/home/doetsch/kaiser0001/miniforge3/envs/vaeda_env

export PYTHONNOUSERSITE=1

log_msg "Python: $(which python)"
python -c "import vaeda; print('VAEDA available')" 2>/dev/null || log_msg "VAEDA: checking..."
echo ""

# ==============================================================================
# CREATE PYTHON SCRIPT
# ==============================================================================

PYTHON_SCRIPT="${SCRIPT_DIR}/run_vaeda_${SAMPLE}.py"

cat > "${PYTHON_SCRIPT}" << PYTHON_EOF
# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
VAEDA Doublet Detection for ${SAMPLE}
Universal Modular Version

NOTE:
  This script does NOT remove doublets. It only writes VAEDA metadata
  used later for consensus voting in Step 6.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import vaeda
import warnings
import time
from scipy import sparse
warnings.filterwarnings('ignore')

SAMPLE = "${SAMPLE}"
INPUT_H5AD = "${INPUT_H5AD}"
OUTPUT_H5AD = "${OUTPUT_H5AD}"
OUTPUT_DIR = "${SAMPLE_OUTPUT_DIR}"

print(f"\\n{'='*70}")
print(f"VAEDA Doublet Detection: {SAMPLE}")
print(f"{'='*70}\\n")

# =========================
# LOAD DATA
# =========================
print(f"[{SAMPLE}] Loading QClus + DropletQC filtered data...")
adata = sc.read_h5ad(INPUT_H5AD)
print(f"[{SAMPLE}] Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")
input_cells = adata.shape[0]

# DropletQC metadata check
dropletqc_cols = [c for c in adata.obs.columns if c.startswith('dropletqc_')]
if dropletqc_cols:
    print(f"[{SAMPLE}] DropletQC metadata present: {dropletqc_cols}")
    if 'dropletqc_status' in adata.obs.columns:
        print(f"[{SAMPLE}] DropletQC status distribution:")
        print(adata.obs['dropletqc_status'].value_counts())
else:
    print(f"[{SAMPLE}] No DropletQC metadata found (cells already filtered)")

# QClus status check
if 'qclus' in adata.obs.columns:
    print(f"\\n[{SAMPLE}] QClus status distribution:")
    print(adata.obs['qclus'].value_counts())
print()

# =========================
# PREPROCESSING
# =========================
print(f"[{SAMPLE}] Preprocessing for VAEDA...")

# keep a raw snapshot if not present
if adata.raw is None:
    adata.raw = adata.copy()

# preserve ORIGINAL counts for VAEDA and ensure sparse CSR
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()
if not sparse.issparse(adata.layers["counts"]):
    adata.layers["counts"] = sparse.csr_matrix(adata.layers["counts"])

# Preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', layer=None, subset=False)
print(f"[{SAMPLE}] Identified {adata.var.highly_variable.sum()} highly variable genes")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

print(f"[{SAMPLE}] Preprocessing complete")
print()

# =========================
# RUN VAEDA
# =========================
print(f"[{SAMPLE}] Running VAEDA doublet detection...")
print(f"[{SAMPLE}] This may take several minutes...")

start_time = time.time()

try:
    # VAEDA expects raw counts in adata.X
    adata.X = adata.layers["counts"]

    vaeda.vaeda(adata)

    elapsed = time.time() - start_time
    print(f"[{SAMPLE}] VAEDA completed in {elapsed/60:.1f} minutes")

    vaeda_cols = [c for c in adata.obs.columns if 'vaeda' in c.lower() or 'doublet' in c.lower()]
    print(f"[{SAMPLE}] VAEDA added columns: {vaeda_cols}")

    if 'vaeda_calls' in adata.obs.columns:
        print(f"\\n[{SAMPLE}] VAEDA calls:")
        print(adata.obs['vaeda_calls'].value_counts())
        n_doublets = int((adata.obs['vaeda_calls'] == 'doublet').sum())
        n_singlets = int((adata.obs['vaeda_calls'] == 'singlet').sum())
        pct_doublets = 100 * n_doublets / len(adata)
        print(f"[{SAMPLE}] Doublet rate: {pct_doublets:.2f}%")
    else:
        n_doublets = 0
        n_singlets = len(adata)

    # Provide Step 6 compatible alias WITHOUT removing anything
    if 'vaeda_calls' in adata.obs.columns:
        adata.obs['vaeda_prediction'] = adata.obs['vaeda_calls'].astype(str)
    if 'vaeda_scores' in adata.obs.columns:
        adata.obs['vaeda_score'] = pd.to_numeric(adata.obs['vaeda_scores'], errors='coerce')

except Exception as e:
    print(f"[{SAMPLE}] VAEDA error: {e}")
    print(f"[{SAMPLE}] Attempting alternative doublet detection with Scrublet...")

    import scrublet as scr

    counts = adata.raw.X if adata.raw is not None else adata.X
    scrub = scr.Scrublet(counts)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    adata.obs['vaeda_prediction'] = ['doublet' if d else 'singlet' for d in predicted_doublets]
    adata.obs['vaeda_score'] = doublet_scores

    n_doublets = int(sum(predicted_doublets))
    n_singlets = int(len(predicted_doublets) - n_doublets)
    pct_doublets = 100 * n_doublets / len(adata)
    print(f"[{SAMPLE}] Scrublet doublet rate: {pct_doublets:.2f}%")

# =========================
# CREATE SUMMARY
# =========================
print(f"\\n[{SAMPLE}] Creating summary...")

summary = {
    'sample': SAMPLE,
    'total_cells': int(len(adata)),
    'n_genes': int(adata.shape[1]),
    'n_singlets': int(n_singlets),
    'n_doublets': int(n_doublets),
    'doublet_rate': float(n_doublets / len(adata)) if len(adata) else 0.0,
}

if 'dropletqc_empty_droplet' in adata.obs.columns:
    summary['dropletqc_empty_in_input'] = int(adata.obs['dropletqc_empty_droplet'].sum())
if 'dropletqc_damaged_cell' in adata.obs.columns:
    summary['dropletqc_damaged_in_input'] = int(adata.obs['dropletqc_damaged_cell'].sum())

summary_df = pd.DataFrame([summary])
summary_path = f"{OUTPUT_DIR}/{SAMPLE}_vaeda_summary.csv"
summary_df.to_csv(summary_path, index=False)
print(f"[{SAMPLE}] Saved summary to: {summary_path}")

# =========================
# SAVE OUTPUT
# =========================
print(f"\\n[{SAMPLE}] Saving VAEDA output (no filtering)...")
adata.write(OUTPUT_H5AD)
print(f"[{SAMPLE}] Saved to: {OUTPUT_H5AD}")

# =========================
# FINAL SUMMARY
# =========================
print(f"\\n{'='*70}")
print(f"VAEDA Analysis Complete: {SAMPLE}")
print(f"{'='*70}")
print(f"\\nInput:  {input_cells} cells (QClus + DropletQC filtered)")
print(f"VAEDA calls: {n_singlets} singlets, {n_doublets} doublets")
print(f"Doublet rate: {100*summary['doublet_rate']:.2f}%")
print(f"\\nOutput files:")
print(f"  - {OUTPUT_H5AD}")
print(f"  - {summary_path}")
print(f"\\nDone!")
print(f"{'='*70}\\n")
PYTHON_EOF

chmod +x "${PYTHON_SCRIPT}"

# ==============================================================================
# RUN VAEDA
# ==============================================================================

log_msg "============================================================================"
log_msg "Running VAEDA Doublet Detection"
log_msg "============================================================================"
echo ""

python "${PYTHON_SCRIPT}"

PYTHON_EXIT=$?
END_TIME=$(date '+%H:%M:%S')

# ==============================================================================
# UPDATE README
# ==============================================================================

update_readme() {
    local sample=$1
    local status=$2

    (
        flock -x 200
        if [[ -f "${SAMPLE_OUTPUT_DIR}/${sample}_vaeda_summary.csv" ]]; then
            INPUT_CELLS=$(awk -F',' 'NR==2 {print $2}' "${SAMPLE_OUTPUT_DIR}/${sample}_vaeda_summary.csv")
            SINGLETS=$(awk -F',' 'NR==2 {print $4}' "${SAMPLE_OUTPUT_DIR}/${sample}_vaeda_summary.csv")
            DOUBLETS=$(awk -F',' 'NR==2 {print $5}' "${SAMPLE_OUTPUT_DIR}/${sample}_vaeda_summary.csv")
            echo "  ${sample}: ${status} | Cells: ${INPUT_CELLS}, Singlets: ${SINGLETS}, Doublets: ${DOUBLETS}" >> "${README_FILE}"
        else
            echo "  ${sample}: ${status}" >> "${README_FILE}"
        fi
    ) 200>"${README_FILE}.lock"
}

if [[ $PYTHON_EXIT -eq 0 ]]; then
    update_readme "${SAMPLE}" "SUCCESS"
    log_msg "Status: SUCCESS"
else
    update_readme "${SAMPLE}" "FAILED"
    error_exit "VAEDA failed"
fi

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

echo ""
log_msg "============================================================================"
log_msg "VAEDA Complete - ${SAMPLE}"
log_msg "============================================================================"
log_msg "Output: ${SAMPLE_OUTPUT_DIR}"
echo ""
ls -lh "${SAMPLE_OUTPUT_DIR}/"
echo ""
log_msg "============================================================================"
