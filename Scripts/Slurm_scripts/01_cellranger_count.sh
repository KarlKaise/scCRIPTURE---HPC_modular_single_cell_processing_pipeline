#!/bin/bash
#SBATCH --job-name=CR_count
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=64G
#SBATCH --time=23:50:00
#SBATCH --qos=1day
#SBATCH --output=logs/cellranger/cellranger_%a_%A.out
#SBATCH --error=logs/cellranger/cellranger_%a_%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ============================================================================
# STEP 1: Cell Ranger Count - Universal Modular Pipeline
# ============================================================================
#
# This script reads from two samplesheets:
#   1. samplesheet.csv       - Sample metadata (sample_id, sex, batch, etc.)
#   2. fastq_samplesheet.csv - FASTQ file locations with multi-lane support
#
# FASTQ SAMPLESHEET FORMAT:
#   sample_id,fastq_dir,fastq_prefix,lane,read1_file,read2_file,index_file,read1_renamed,read2_renamed
#
# Multiple rows with the same sample_id = multiple lanes (auto-combined)
#
# OUTPUT STRUCTURE:
#   Output_dir_<dataset_name>/Single_cell_preprocessed/1_CellRanger_output/
#
# Usage:
#   # Auto-detect sample count and submit
#   N_SAMPLES=$(tail -n +2 samplesheet.csv | grep -v '^#' | grep -v '^$' | wc -l)
#   sbatch --array=1-${N_SAMPLES} Scripts/Slurm_scripts/01_cellranger_count.sh
#
#   # Or run single sample
#   sbatch Scripts/Slurm_scripts/01_cellranger_count.sh  # Runs first sample
#
# ============================================================================

set -e

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Base directory (from environment or auto-detect)
if [[ -n "${PROJECT_ROOT:-}" ]]; then
    BASE_DIR="${PROJECT_ROOT}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    BASE_DIR="${SLURM_SUBMIT_DIR}"
else
    # Default for this dataset
    BASE_DIR="/scicore/home/doetsch/kaiser0001/Single_cell_paper/Datasets/Human_Covid_LV_ChP_PMID_34153974"
fi

# Samplesheet paths (fixed names in project root)
SAMPLESHEET="${BASE_DIR}/samplesheet.csv"
FASTQ_SAMPLESHEET="${BASE_DIR}/fastq_samplesheet.csv"

# Get dataset name from samplesheet (column 12)
get_dataset_name() {
    tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f12
}

# Set output directories based on environment or derive from samplesheet
if [[ -n "${PREPROCESS_DIR:-}" ]]; then
    OUTPUT_DIR="${PREPROCESS_DIR}/1_CellRanger_output"
else
    DATASET_NAME=$(get_dataset_name)
    if [[ -n "$DATASET_NAME" ]]; then
        OUTPUT_DIR="${BASE_DIR}/Output_dir_${DATASET_NAME}/Single_cell_preprocessed/1_CellRanger_output"
    else
        # Fallback to old structure if dataset_name not found
        OUTPUT_DIR="${BASE_DIR}/1_CellRanger_output"
    fi
fi

SCRATCH_BASE="/scratch/${USER}/cellranger_${SLURM_JOB_ID:-$$}"
LOG_DIR="${BASE_DIR}/logs/cellranger"
README_FILE="${OUTPUT_DIR}/README.txt"

# Default transcriptome (can be overridden per-sample in samplesheet.csv)
DEFAULT_TRANSCRIPTOME="/scicore/home/doetsch/kaiser0001/Single_cell_paper/Offline_resources/Human_ref_genome/refdata-gex-GRCh38-2024-A"

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

# Get column index by name (1-based for cut)
# Skips comment lines starting with #
get_col_index() {
    local file=$1
    local col_name=$2
    grep -v '^#' "$file" | head -1 | tr ',' '\n' | grep -n "^${col_name}$" | cut -d: -f1
}

# Get unique sample IDs from samplesheet
get_unique_samples() {
    local col=$(get_col_index "$SAMPLESHEET" "sample_id")
    tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | cut -d',' -f"$col" | sort -u
}

# Get sample count
get_sample_count() {
    get_unique_samples | wc -l
}

# Get nth sample (1-based)
get_sample_by_index() {
    local idx=$1
    get_unique_samples | sed -n "${idx}p"
}

# Get field from samplesheet for a sample
get_sample_field() {
    local sample=$1
    local field=$2
    local col=$(get_col_index "$SAMPLESHEET" "$field")
    if [[ -n "$col" ]]; then
        grep "^${sample}," "$SAMPLESHEET" | head -1 | cut -d',' -f"$col"
    fi
}

# Get all FASTQ rows for a sample from fastq_samplesheet
get_fastq_rows() {
    local sample=$1
    grep -v '^#' "$FASTQ_SAMPLESHEET" | grep "^${sample},"
}

# Count lanes for a sample
get_lane_count() {
    local sample=$1
    get_fastq_rows "$sample" | wc -l
}

# ==============================================================================
# VALIDATE INPUTS
# ==============================================================================

log_msg "============================================================================"
log_msg "Cell Ranger Count - Universal Modular Pipeline"
log_msg "============================================================================"
log_msg "Job ID:        ${SLURM_JOB_ID:-local}"
log_msg "Array Task:    ${SLURM_ARRAY_TASK_ID:-1}"
log_msg "Date:          $(date)"
log_msg "Host:          $(hostname)"
log_msg "CPUs:          ${SLURM_CPUS_PER_TASK:-8}"
log_msg "Base Dir:      ${BASE_DIR}"
log_msg "Output Dir:    ${OUTPUT_DIR}"
log_msg "============================================================================"
echo ""

# Validate samplesheets exist
if [[ ! -f "$SAMPLESHEET" ]]; then
    error_exit "Samplesheet not found: $SAMPLESHEET"
fi

if [[ ! -f "$FASTQ_SAMPLESHEET" ]]; then
    error_exit "FASTQ samplesheet not found: $FASTQ_SAMPLESHEET"
fi

log_msg "Samplesheet:      $SAMPLESHEET"
log_msg "FASTQ samplesheet: $FASTQ_SAMPLESHEET"
echo ""

# ==============================================================================
# DETERMINE SAMPLE TO PROCESS
# ==============================================================================

N_SAMPLES=$(get_sample_count)
log_msg "Total unique samples: $N_SAMPLES"

# Get task index
TASK_IDX=${SLURM_ARRAY_TASK_ID:-1}

if [[ $TASK_IDX -gt $N_SAMPLES ]]; then
    error_exit "Array task $TASK_IDX exceeds sample count $N_SAMPLES"
fi

# Get sample for this task
SAMPLE=$(get_sample_by_index $TASK_IDX)

if [[ -z "$SAMPLE" ]]; then
    error_exit "Could not determine sample for task index $TASK_IDX"
fi

log_msg "Processing: $SAMPLE (task $TASK_IDX of $N_SAMPLES)"
echo ""

# ==============================================================================
# GET SAMPLE METADATA
# ==============================================================================

SAMPLE_SEX=$(get_sample_field "$SAMPLE" "sex")
SAMPLE_BATCH=$(get_sample_field "$SAMPLE" "batch")
SAMPLE_VENTRICLE=$(get_sample_field "$SAMPLE" "ventricle")
SAMPLE_SPECIES=$(get_sample_field "$SAMPLE" "species")
TRANSCRIPTOME=$(get_sample_field "$SAMPLE" "cellranger_ref")

# Use default transcriptome if not specified
TRANSCRIPTOME="${TRANSCRIPTOME:-$DEFAULT_TRANSCRIPTOME}"

log_msg "Sample metadata:"
log_msg "  Sample ID:    $SAMPLE"
log_msg "  Sex:          ${SAMPLE_SEX:-N/A}"
log_msg "  Batch:        ${SAMPLE_BATCH:-N/A}"
log_msg "  Ventricle:    ${SAMPLE_VENTRICLE:-N/A}"
log_msg "  Species:      ${SAMPLE_SPECIES:-N/A}"
log_msg "  Transcriptome: $TRANSCRIPTOME"
echo ""

# Validate transcriptome
if [[ ! -d "$TRANSCRIPTOME" ]]; then
    error_exit "Transcriptome reference not found: $TRANSCRIPTOME"
fi

# ==============================================================================
# GET FASTQ INFORMATION
# ==============================================================================

N_LANES=$(get_lane_count "$SAMPLE")
log_msg "FASTQ lanes/runs for $SAMPLE: $N_LANES"

if [[ $N_LANES -eq 0 ]]; then
    error_exit "No FASTQ entries found for sample $SAMPLE in $FASTQ_SAMPLESHEET"
fi

# Get column indices for FASTQ samplesheet
FQ_SAMPLE_COL=$(get_col_index "$FASTQ_SAMPLESHEET" "sample_id")
FQ_DIR_COL=$(get_col_index "$FASTQ_SAMPLESHEET" "fastq_dir")
FQ_LANE_COL=$(get_col_index "$FASTQ_SAMPLESHEET" "lane")
FQ_R1_COL=$(get_col_index "$FASTQ_SAMPLESHEET" "read1_file")
FQ_R2_COL=$(get_col_index "$FASTQ_SAMPLESHEET" "read2_file")
FQ_R1_RENAMED_COL=$(get_col_index "$FASTQ_SAMPLESHEET" "read1_renamed")
FQ_R2_RENAMED_COL=$(get_col_index "$FASTQ_SAMPLESHEET" "read2_renamed")

# ==============================================================================
# CREATE SCRATCH DIRECTORY WITH RENAMED SYMLINKS
# ==============================================================================

SCRATCH_DIR="${SCRATCH_BASE}/${SAMPLE}"
mkdir -p "${SCRATCH_DIR}/fastqs"

log_msg "Creating Cell Ranger-compatible symlinks in scratch..."
log_msg "Scratch dir: ${SCRATCH_DIR}/fastqs"
echo ""

LANE_NUM=0
while IFS=',' read -r fq_sample fq_dir fq_prefix fq_lane fq_r1 fq_r2 fq_index fq_r1_renamed fq_r2_renamed rest; do
    # Skip if not our sample (shouldn't happen but safety check)
    [[ "$fq_sample" != "$SAMPLE" ]] && continue

    # Increment lane number (fixed: avoid ((LANE_NUM++)) which fails with set -e when LANE_NUM=0)
    LANE_NUM=$((LANE_NUM + 1))

    # Source files
    SRC_R1="${fq_dir}/${fq_r1}"
    SRC_R2="${fq_dir}/${fq_r2}"

    # Validate source files exist
    if [[ ! -f "$SRC_R1" ]]; then
        error_exit "R1 file not found: $SRC_R1"
    fi
    if [[ ! -f "$SRC_R2" ]]; then
        error_exit "R2 file not found: $SRC_R2"
    fi

    # Determine destination names
    if [[ -n "$fq_r1_renamed" && "$fq_r1_renamed" != "" ]]; then
        DEST_R1="$fq_r1_renamed"
    else
        # Auto-generate Cell Ranger compatible name
        DEST_R1="${SAMPLE}_S1_L00${LANE_NUM}_R1_001.fastq"
    fi

    if [[ -n "$fq_r2_renamed" && "$fq_r2_renamed" != "" ]]; then
        DEST_R2="$fq_r2_renamed"
    else
        DEST_R2="${SAMPLE}_S1_L00${LANE_NUM}_R2_001.fastq"
    fi

    # Handle .gz extension
    if [[ "$SRC_R1" == *.gz && "$DEST_R1" != *.gz ]]; then
        DEST_R1="${DEST_R1}.gz"
    fi
    if [[ "$SRC_R2" == *.gz && "$DEST_R2" != *.gz ]]; then
        DEST_R2="${DEST_R2}.gz"
    fi

    # Create symlinks
    ln -sf "$SRC_R1" "${SCRATCH_DIR}/fastqs/${DEST_R1}"
    ln -sf "$SRC_R2" "${SCRATCH_DIR}/fastqs/${DEST_R2}"

    log_msg "  Lane ${LANE_NUM} (${fq_lane}):"
    log_msg "    R1: $(basename $SRC_R1) -> $DEST_R1"
    log_msg "    R2: $(basename $SRC_R2) -> $DEST_R2"

done < <(get_fastq_rows "$SAMPLE")

echo ""
log_msg "Symlinks created:"
ls -la "${SCRATCH_DIR}/fastqs/"
echo ""

# ==============================================================================
# CREATE OUTPUT DIRECTORIES
# ==============================================================================

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Remove existing output if present (Cell Ranger won't overwrite)
if [[ -d "${OUTPUT_DIR}/${SAMPLE}" ]]; then
    log_msg "Removing existing output directory..."
    rm -rf "${OUTPUT_DIR}/${SAMPLE}"
fi

# ==============================================================================
# LOAD MODULES
# ==============================================================================

log_msg "Loading Cell Ranger module..."
module purge
module load CellRanger

log_msg "Cell Ranger version: $(cellranger --version 2>&1 | head -1)"
echo ""

# ==============================================================================
# RUN CELL RANGER COUNT
# ==============================================================================

# Calculate memory (leave some headroom)
MEM_GB=$(( ${SLURM_CPUS_PER_TASK:-8} * ${SLURM_MEM_PER_CPU:-64000} / 1024 - 4 ))
MEM_GB=$(( MEM_GB > 4 ? MEM_GB : 64 ))  # Minimum 64GB or calculated

log_msg "============================================================================"
log_msg "Running Cell Ranger Count"
log_msg "============================================================================"
log_msg "  Sample:       $SAMPLE"
log_msg "  FASTQs:       ${SCRATCH_DIR}/fastqs"
log_msg "  Transcriptome: $TRANSCRIPTOME"
log_msg "  Cores:        ${SLURM_CPUS_PER_TASK:-8}"
log_msg "  Memory:       ${MEM_GB} GB"
log_msg "  Output:       ${OUTPUT_DIR}/${SAMPLE}"
log_msg "  Include introns: yes (snRNA-seq)"
echo ""

START_TIME=$(date +%s)

cd "$OUTPUT_DIR"

cellranger count \
    --id="${SAMPLE}" \
    --transcriptome="${TRANSCRIPTOME}" \
    --fastqs="${SCRATCH_DIR}/fastqs" \
    --sample="${SAMPLE}" \
    --include-introns=true \
    --create-bam=true \
    --localcores=${SLURM_CPUS_PER_TASK:-8} \
    --localmem=${MEM_GB}

CR_EXIT=$?

END_TIME=$(date +%s)
DURATION=$(( (END_TIME - START_TIME) / 60 ))

# ==============================================================================
# CLEANUP AND VERIFY
# ==============================================================================

log_msg "Cleaning up scratch directory..."
rm -rf "$SCRATCH_DIR"

echo ""
log_msg "============================================================================"
log_msg "Cell Ranger Complete"
log_msg "============================================================================"
log_msg "Sample:   $SAMPLE"
log_msg "Duration: ${DURATION} minutes"
log_msg "Exit code: $CR_EXIT"
echo ""

if [[ $CR_EXIT -eq 0 ]]; then
    log_msg "Status: SUCCESS"

    OUTS_DIR="${OUTPUT_DIR}/${SAMPLE}/outs"

    if [[ -d "$OUTS_DIR" ]]; then
        log_msg "Output files:"
        for f in filtered_feature_bc_matrix.h5 raw_feature_bc_matrix.h5 \
                 web_summary.html metrics_summary.csv possorted_genome_bam.bam; do
            if [[ -f "${OUTS_DIR}/${f}" ]]; then
                SIZE=$(du -h "${OUTS_DIR}/${f}" | cut -f1)
                log_msg "  $f: $SIZE"
            fi
        done

        # Extract metrics
        if [[ -f "${OUTS_DIR}/metrics_summary.csv" ]]; then
            echo ""
            log_msg "Key metrics:"
            # Get estimated cells (first numeric column after header)
            CELLS=$(tail -1 "${OUTS_DIR}/metrics_summary.csv" | cut -d',' -f1 | tr -d '"')
            log_msg "  Estimated cells: $CELLS"
        fi
    fi
else
    log_msg "Status: FAILED"
    exit $CR_EXIT
fi

echo ""
log_msg "============================================================================"
