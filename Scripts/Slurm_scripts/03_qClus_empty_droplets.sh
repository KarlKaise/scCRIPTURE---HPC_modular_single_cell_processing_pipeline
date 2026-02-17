#!/bin/bash
#SBATCH --job-name=qClus
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G
#SBATCH --time=23:55:00
#SBATCH --qos=1day
#SBATCH --output=logs/qclus/qclus_%a_%A.out
#SBATCH --error=logs/qclus/qclus_%a_%A.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ============================================================================
# STEP 3: QClus Contaminated Droplet Filtering with DropletQC Integration
# ============================================================================
#
# UNIVERSAL MODULAR VERSION - Reads samples from samplesheet.csv
#
# QClus filtering pipeline:
#   1. Initial QC filter (min genes, max genes, mito%)
#   2. Outlier filter (unspliced fraction, mito distribution)
#   3. Clustering filter (k-means on quality metrics)
#   4. Scrublet doublet filter
#
# DropletQC integration:
#   - Add empty droplet and damaged cell flags to metadata
#   - Filter out DropletQC-identified empty droplets
#
# Input:
#   - Cell Ranger output (filtered_feature_bc_matrix.h5, BAM)
#   - DropletQC results (from Step 2)
#
# Output:
#   - QClus filtered h5ad with DropletQC metadata
#
# Usage:
#   sbatch --array=1-N 03_qClus_empty_droplets.sh
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
CELLRANGER_DIR="${PREPROCESS_DIR}/1_CellRanger_output"
DROPLETQC_DIR="${PREPROCESS_DIR}/2_DropletQC_output"
OUTPUT_DIR="${PREPROCESS_DIR}/3_qClus_empty_droplets"
SCRIPT_DIR="${BASE_DIR}/Scripts/Python_scripts"
LOG_DIR="${BASE_DIR}/logs/qclus"
README_FILE="${OUTPUT_DIR}/README.txt"

# QClus settings
N_CORES=${SLURM_CPUS_PER_TASK:-8}
N_TILES=100

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
    grep -v '^#' "$file" | head -1 | tr ',' '\n' | grep -n "^${col_name}$" | cut -d: -f1
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
log_msg "QClus Contaminated Droplet Filtering - Modular Pipeline"
log_msg "============================================================================"
log_msg "Job ID:        ${SLURM_JOB_ID:-local}"
log_msg "Array Task:    ${SLURM_ARRAY_TASK_ID:-1}"
log_msg "Date:          $(date)"
log_msg "Host:          $(hostname)"
log_msg "CPUs:          ${N_CORES}"
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

# Output directory
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}"
START_TIME=$(date '+%H:%M:%S')

# Create directories
mkdir -p "${SAMPLE_OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${SCRIPT_DIR}"

cd "${BASE_DIR}"

# ==============================================================================
# INITIALIZE README
# ==============================================================================

if [[ ${SLURM_ARRAY_TASK_ID:-1} -eq 1 ]]; then
    cat > "${README_FILE}" << EOF
================================================================================
STEP 3: QClus Contaminated Droplet Filtering
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')
Pipeline: Universal Modular scRNA-seq Pipeline

SAMPLESHEET: ${SAMPLESHEET}
INPUT:
  - Cell Ranger: 1_CellRanger_output/<SAMPLE>/outs/
  - DropletQC:   2_DropletQC_output/<SAMPLE>/

SAMPLES PROCESSED:
EOF
fi

# ==============================================================================
# ACTIVATE CONDA
# ==============================================================================

log_msg "Activating QClus conda environment..."

source /scicore/home/doetsch/kaiser0001/miniforge3/etc/profile.d/conda.sh
conda activate /scicore/home/doetsch/kaiser0001/.conda/envs/qclus

export PYTHONNOUSERSITE=1

log_msg "Python: $(which python)"
python -c "import qclus; print('QClus version:', qclus.__version__)" 2>/dev/null || log_msg "QClus version: unknown"
echo ""

# ==============================================================================
# VERIFY INPUTS
# ==============================================================================

log_msg "============================================================================"
log_msg "Input Verification"
log_msg "============================================================================"

OUTS_DIR="${CELLRANGER_DIR}/${SAMPLE}/outs"

if [[ ! -d "${OUTS_DIR}" ]]; then
    error_exit "Cell Ranger output not found: ${OUTS_DIR}"
fi

H5_FILE="${OUTS_DIR}/filtered_feature_bc_matrix.h5"
BAM_FILE="${OUTS_DIR}/possorted_genome_bam.bam"
BAI_FILE="${OUTS_DIR}/possorted_genome_bam.bam.bai"
BC_FILE="${OUTS_DIR}/filtered_feature_bc_matrix/barcodes.tsv.gz"

log_msg "Cell Ranger files:"
for f in "${H5_FILE}" "${BAM_FILE}" "${BAI_FILE}" "${BC_FILE}"; do
    if [[ -f "$f" ]]; then
        log_msg "  OK: $(basename "$f")"
    else
        error_exit "Missing: $(basename "$f")"
    fi
done

# Check DropletQC files
DROPLETQC_SAMPLE_DIR="${DROPLETQC_DIR}/${SAMPLE}"
DROPLETQC_RESULTS="${DROPLETQC_SAMPLE_DIR}/${SAMPLE}_dropletqc_results.csv"

log_msg "DropletQC files:"
if [[ -f "${DROPLETQC_RESULTS}" ]]; then
    log_msg "  OK: DropletQC results found"
    DROPLETQC_AVAILABLE="True"
else
    log_msg "  WARNING: DropletQC results not found - proceeding without"
    DROPLETQC_AVAILABLE="False"
fi
echo ""

# ==============================================================================
# CREATE PYTHON SCRIPT
# ==============================================================================

PYTHON_SCRIPT="${SCRIPT_DIR}/run_qclus_${SAMPLE}.py"

cat > "${PYTHON_SCRIPT}" << PYTHON_EOF
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
QClus + DropletQC Filtering for ${SAMPLE}
Universal Modular Version
"""

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import qclus as qc
import time
import warnings
warnings.filterwarnings('ignore')

# Configuration
SAMPLE = "${SAMPLE}"
CELLRANGER_DIR = Path("${CELLRANGER_DIR}")
DROPLETQC_DIR = Path("${DROPLETQC_DIR}")
OUTPUT_DIR = Path("${SAMPLE_OUTPUT_DIR}")
N_CORES = ${N_CORES}
N_TILES = ${N_TILES}
DROPLETQC_AVAILABLE = ${DROPLETQC_AVAILABLE}

# Paths
OUTS = CELLRANGER_DIR / SAMPLE / "outs"
H5_PATH = OUTS / "filtered_feature_bc_matrix.h5"
BAM_PATH = OUTS / "possorted_genome_bam.bam"
BAI_PATH = OUTS / "possorted_genome_bam.bam.bai"
BC_PATH = OUTS / "filtered_feature_bc_matrix" / "barcodes.tsv.gz"

DROPLETQC_RESULTS = DROPLETQC_DIR / SAMPLE / f"{SAMPLE}_dropletqc_results.csv"
DROPLETQC_EMPTY = DROPLETQC_DIR / SAMPLE / f"{SAMPLE}_empty_droplet_barcodes.txt"
DROPLETQC_DAMAGED = DROPLETQC_DIR / SAMPLE / f"{SAMPLE}_damaged_cell_barcodes.txt"

print(f"\n{'='*70}")
print(f"QClus + DropletQC Analysis: {SAMPLE}")
print(f"{'='*70}\n")

# =========================
# LOAD DROPLETQC RESULTS
# =========================

print(f"[{SAMPLE}] Loading DropletQC results...")

dropletqc_df = None
empty_barcodes = set()
damaged_barcodes = set()

if DROPLETQC_AVAILABLE and DROPLETQC_RESULTS.exists():
    dropletqc_df = pd.read_csv(DROPLETQC_RESULTS, index_col=0)
    print(f"[{SAMPLE}] Loaded DropletQC for {len(dropletqc_df)} barcodes")

    if DROPLETQC_EMPTY.exists():
        with open(DROPLETQC_EMPTY, 'r') as f:
            empty_barcodes = set(line.strip() for line in f if line.strip())
        print(f"[{SAMPLE}] Empty droplets: {len(empty_barcodes)}")

    if DROPLETQC_DAMAGED.exists():
        with open(DROPLETQC_DAMAGED, 'r') as f:
            damaged_barcodes = set(line.strip() for line in f if line.strip())
        print(f"[{SAMPLE}] Damaged cells: {len(damaged_barcodes)}")
else:
    print(f"[{SAMPLE}] No DropletQC data - proceeding without")

print()

# =========================
# LOAD COUNTS
# =========================

print(f"[{SAMPLE}] Loading counts from: {H5_PATH}")
adata = sc.read_10x_h5(str(H5_PATH))
print(f"[{SAMPLE}] Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

adata.obs["barcode_cellranger"] = adata.obs_names.astype(str)
trimmed = adata.obs_names.astype(str).str.replace(r"-1\$", "", regex=True)
adata.obs["barcode_trim"] = trimmed
adata.obs_names = trimmed
adata.var_names_make_unique()

initial_cells = len(adata)

# =========================
# ADD DROPLETQC METADATA
# =========================

if dropletqc_df is not None:
    print(f"\n[{SAMPLE}] Adding DropletQC metadata...")

    adata.obs["dropletqc_nuclear_fraction"] = np.nan
    adata.obs["dropletqc_empty_droplet"] = False
    adata.obs["dropletqc_damaged_cell"] = False
    adata.obs["dropletqc_status"] = "unknown"

    dropletqc_barcodes_trim = dropletqc_df.index.astype(str).str.replace(r"-1\$", "", regex=True)
    bc_to_nf = dict(zip(dropletqc_barcodes_trim, dropletqc_df['nuclear_fraction']))
    bc_to_status = dict(zip(dropletqc_barcodes_trim, dropletqc_df['final_status']))

    for bc in adata.obs_names:
        if bc in bc_to_nf:
            adata.obs.loc[bc, "dropletqc_nuclear_fraction"] = bc_to_nf[bc]
            adata.obs.loc[bc, "dropletqc_status"] = bc_to_status.get(bc, "unknown")

    empty_barcodes_trim = set(bc.replace("-1", "") for bc in empty_barcodes)
    damaged_barcodes_trim = set(bc.replace("-1", "") for bc in damaged_barcodes)

    adata.obs["dropletqc_empty_droplet"] = adata.obs_names.isin(empty_barcodes_trim)
    adata.obs["dropletqc_damaged_cell"] = adata.obs_names.isin(damaged_barcodes_trim)

    print(f"[{SAMPLE}] Empty flagged: {adata.obs['dropletqc_empty_droplet'].sum()}")
    print(f"[{SAMPLE}] Damaged flagged: {adata.obs['dropletqc_damaged_cell'].sum()}")

# Save initial counts
counts_path = OUTPUT_DIR / f"{SAMPLE}_counts_filtered.h5ad"
adata.write(counts_path)
print(f"[{SAMPLE}] Saved: {counts_path}")

# =========================
# CALCULATE FRACTION UNSPLICED
# =========================

print(f"\n[{SAMPLE}] Calculating fraction_unspliced from BAM...")
print(f"[{SAMPLE}] This may take 1-2 hours for large BAMs...")

start_time = time.time()

fraction_unspliced = qc.utils.fraction_unspliced_from_bam(
    bam_path=str(BAM_PATH),
    bam_index_path=str(BAI_PATH),
    barcodes_path=str(BC_PATH),
    tiles=N_TILES,
    cores=N_CORES
)

elapsed = time.time() - start_time
print(f"[{SAMPLE}] Completed in {elapsed/60:.1f} minutes")

frac_path = OUTPUT_DIR / f"{SAMPLE}_fraction_unspliced.csv"
fraction_unspliced.to_csv(frac_path)

# =========================
# ALIGN BARCODES
# =========================

print(f"\n[{SAMPLE}] Aligning barcodes...")

counts_bcs = pd.Index(adata.obs_names.astype(str))
frac_bcs_orig = fraction_unspliced.index.astype(str)
frac_bcs = pd.Index(frac_bcs_orig.str.replace(r"-1\$", "", regex=True))

common = counts_bcs.intersection(frac_bcs).sort_values()
print(f"[{SAMPLE}] Common barcodes: {len(common)}")

if len(common) == 0:
    raise ValueError("No overlapping barcodes!")

adata_aligned = adata[common].copy()

frac_mapping = pd.Series(frac_bcs_orig.values, index=frac_bcs)
original_indices = frac_mapping.loc[common].values
frac_aligned = fraction_unspliced.loc[original_indices].copy()
frac_aligned.index = common

aligned_path = OUTPUT_DIR / f"{SAMPLE}_counts_aligned.h5ad"
adata_aligned.write(aligned_path)

# =========================
# RUN QCLUS
# =========================

print(f"\n[{SAMPLE}] Running QClus...")

adata_qclus = adata_aligned.copy()
adata_qclus.var_names = adata_qclus.var_names.str.upper()
adata_qclus.var_names_make_unique()

qclus_input_path = OUTPUT_DIR / f"{SAMPLE}_counts_qclus_upper.h5ad"
adata_qclus.write(qclus_input_path)

start_time = time.time()
adata_qc = qc.quickstart_qclus(
    str(qclus_input_path),
    frac_aligned,
    tissue="other"
)
elapsed = time.time() - start_time
print(f"[{SAMPLE}] QClus completed in {elapsed/60:.1f} minutes")

qclus_cells = len(adata_qc)

# =========================
# INTEGRATE DROPLETQC
# =========================

if dropletqc_df is not None:
    print(f"\n[{SAMPLE}] Integrating DropletQC into QClus output...")

    for col in ["dropletqc_nuclear_fraction", "dropletqc_empty_droplet",
                "dropletqc_damaged_cell", "dropletqc_status"]:
        if col in adata_aligned.obs.columns:
            common_bcs = adata_qc.obs_names.intersection(adata_aligned.obs_names)
            adata_qc.obs[col] = np.nan if "fraction" in col else False
            adata_qc.obs.loc[common_bcs, col] = adata_aligned.obs.loc[common_bcs, col]

qclus_path = OUTPUT_DIR / f"{SAMPLE}_qclus_filtered.h5ad"
adata_qc.write(qclus_path)

# =========================
# BUILD STATUS TABLE
# =========================

print(f"\n[{SAMPLE}] Building status table...")

status = pd.DataFrame({
    "barcode_trim": adata.obs["barcode_trim"].values
}, index=adata.obs["barcode_cellranger"].values)

input_bcs = pd.Index(adata_aligned.obs["barcode_trim"].astype(str).values)
kept_bcs = pd.Index(adata_qc.obs_names.astype(str))

status["in_qclus_input"] = status["barcode_trim"].isin(input_bcs)
status["in_qclus_kept"] = status["barcode_trim"].isin(kept_bcs)

if dropletqc_df is not None:
    status["dropletqc_empty"] = status["barcode_trim"].isin(empty_barcodes_trim)
    status["dropletqc_damaged"] = status["barcode_trim"].isin(damaged_barcodes_trim)
else:
    status["dropletqc_empty"] = False
    status["dropletqc_damaged"] = False

status["pass_qclus"] = status["in_qclus_kept"]
status["pass_dropletqc"] = ~status["dropletqc_empty"]
status["pass_combined"] = status["pass_qclus"] & status["pass_dropletqc"]

# Also add a simple 'status' column for downstream compatibility
status["status"] = "cell"
status.loc[status["dropletqc_empty"], "status"] = "empty_droplet"
status.loc[status["dropletqc_damaged"], "status"] = "damaged_cell"
status.loc[~status["pass_qclus"], "status"] = "filtered_qclus"

status_path = OUTPUT_DIR / f"{SAMPLE}_qclus_dropletqc_status.csv"
status.to_csv(status_path)

print(f"[{SAMPLE}] Status summary:")
print(f"  Total:           {len(status)}")
print(f"  QClus kept:      {status['in_qclus_kept'].sum()}")
print(f"  DropletQC empty: {status['dropletqc_empty'].sum()}")
print(f"  Pass combined:   {status['pass_combined'].sum()}")

# =========================
# FINAL FILTERED OUTPUT
# =========================

print(f"\n[{SAMPLE}] Creating final filtered output...")

final_keep_bcs = status[status["pass_combined"]]["barcode_trim"].values
adata_final = adata_qc[adata_qc.obs_names.isin(final_keep_bcs)].copy()
adata_final.obs["filter_status"] = "passed_qclus_dropletqc"

final_cells = len(adata_final)

final_path = OUTPUT_DIR / f"{SAMPLE}_qclus_dropletqc_filtered.h5ad"
adata_final.write(final_path)

# Save counts summary
counts_df = pd.DataFrame([{
    'sample': SAMPLE,
    'initial_cells': initial_cells,
    'qclus_cells': qclus_cells,
    'final_cells': final_cells
}])
counts_df.to_csv(OUTPUT_DIR / f"{SAMPLE}_cell_counts.csv", index=False)

# =========================
# SUMMARY
# =========================

print(f"\n{'='*70}")
print(f"Complete: {SAMPLE}")
print(f"{'='*70}")
print(f"  Cell Ranger: {initial_cells}")
print(f"  After QClus: {qclus_cells}")
print(f"  Final:       {final_cells}")
print(f"  Retention:   {100*final_cells/initial_cells:.1f}%")
print(f"\nUse: {final_path}")
print(f"{'='*70}\n")
PYTHON_EOF

chmod +x "${PYTHON_SCRIPT}"

# ==============================================================================
# RUN QCLUS
# ==============================================================================

log_msg "============================================================================"
log_msg "Running QClus Pipeline"
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
        if [[ -f "${SAMPLE_OUTPUT_DIR}/${sample}_cell_counts.csv" ]]; then
            INITIAL=$(awk -F',' 'NR==2 {print $2}' "${SAMPLE_OUTPUT_DIR}/${sample}_cell_counts.csv")
            QCLUS=$(awk -F',' 'NR==2 {print $3}' "${SAMPLE_OUTPUT_DIR}/${sample}_cell_counts.csv")
            FINAL=$(awk -F',' 'NR==2 {print $4}' "${SAMPLE_OUTPUT_DIR}/${sample}_cell_counts.csv")
            echo "  ${sample}: ${status} | Initial: ${INITIAL}, QClus: ${QCLUS}, Final: ${FINAL}" >> "${README_FILE}"
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
    error_exit "QClus failed"
fi

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

echo ""
log_msg "============================================================================"
log_msg "QClus Complete - ${SAMPLE}"
log_msg "============================================================================"
log_msg "Output: ${SAMPLE_OUTPUT_DIR}"
echo ""
ls -lh "${SAMPLE_OUTPUT_DIR}/"
echo ""
log_msg "Use for downstream: ${SAMPLE_OUTPUT_DIR}/${SAMPLE}_qclus_dropletqc_filtered.h5ad"
log_msg "============================================================================"
