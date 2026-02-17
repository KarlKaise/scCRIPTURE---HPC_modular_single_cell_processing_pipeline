#!/bin/bash
#SBATCH --job-name=DropletQC
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --qos=1day
#SBATCH --output=logs/dropletqc/dropletqc_%a_%A.out
#SBATCH --error=logs/dropletqc/dropletqc_%a_%A.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ============================================================================
# STEP 2: DropletQC Nuclear Fraction & Empty Droplet Identification
# ============================================================================
#
# UNIVERSAL MODULAR VERSION - Reads samples from samplesheet.csv
#
# This script calculates nuclear fraction scores using DropletQC and
# identifies empty droplets and damaged cells in Cell Ranger filtered output.
#
# Nuclear fraction = intronic reads / (intronic + exonic reads)
#
# Cell Classification:
#   - Empty droplets: LOW nuclear fraction + LOW UMI (ambient RNA)
#   - Damaged cells:  HIGH nuclear fraction + LOW UMI (nuclear RNA leak)
#   - Intact cells:   MODERATE nuclear fraction + MODERATE-HIGH UMI
#
# Input:  Cell Ranger output (outs directory with BAM and filtered matrix)
# Output: Nuclear fraction scores, empty droplet flags, damaged cell flags
#
# Usage:
#   # Submit as array job (sample count from samplesheet)
#   sbatch --array=1-N 02_DropletQC_nuclear_fraction.sh
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

# Samplesheet (fixed name in project root)
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
OUTPUT_DIR="${PREPROCESS_DIR}/2_DropletQC_output"
SCRIPT_DIR="${BASE_DIR}/Scripts/R_scripts"
LOG_DIR="${BASE_DIR}/logs/dropletqc"
README_FILE="${OUTPUT_DIR}/README.txt"

# DropletQC settings
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

# Get column index by name (1-based for cut)
get_col_index() {
    local file=$1
    local col_name=$2
    head -1 "$file" | tr ',' '\n' | grep -n "^${col_name}$" | cut -d: -f1
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

# ==============================================================================
# INITIALIZE
# ==============================================================================

log_msg "============================================================================"
log_msg "DropletQC Nuclear Fraction Analysis - Modular Pipeline"
log_msg "============================================================================"
log_msg "Job ID:        ${SLURM_JOB_ID:-local}"
log_msg "Array Task:    ${SLURM_ARRAY_TASK_ID:-1}"
log_msg "Date:          $(date)"
log_msg "Host:          $(hostname)"
log_msg "CPUs:          ${N_CORES}"
log_msg "Tiles:         ${N_TILES}"
log_msg "Base Dir:      ${BASE_DIR}"
log_msg "============================================================================"
echo ""

# Validate samplesheet
if [[ ! -f "$SAMPLESHEET" ]]; then
    error_exit "Samplesheet not found: $SAMPLESHEET"
fi

log_msg "Samplesheet: $SAMPLESHEET"

# Get sample count and current sample
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

log_msg "Sample metadata:"
log_msg "  Sample ID:   $SAMPLE"
log_msg "  Sex:         ${SAMPLE_SEX:-N/A}"
log_msg "  Batch:       ${SAMPLE_BATCH:-N/A}"
log_msg "  Ventricle:   ${SAMPLE_VENTRICLE:-N/A}"
echo ""

# Output directory for this sample
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}"
START_TIME=$(date '+%H:%M:%S')

# Create directories
mkdir -p "${SAMPLE_OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${SCRIPT_DIR}"

# ==============================================================================
# INITIALIZE README (first task only)
# ==============================================================================

init_readme() {
    if [[ ${SLURM_ARRAY_TASK_ID:-1} -eq 1 ]]; then
        cat > "${README_FILE}" << EOF
================================================================================
STEP 2: DropletQC Nuclear Fraction & Empty Droplet Identification
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')
Pipeline: Universal Modular scRNA-seq Pipeline

DESCRIPTION:
  Calculates nuclear fraction scores and identifies empty droplets/damaged cells
  using the DropletQC R package.

  Nuclear fraction = intronic reads / (intronic + exonic reads)

  - Empty droplets: LOW nuclear fraction + LOW UMI (ambient RNA)
  - Damaged cells:  HIGH nuclear fraction + LOW UMI (nuclear RNA leak)
  - Intact cells:   MODERATE nuclear fraction + MODERATE-HIGH UMI

SAMPLESHEET: ${SAMPLESHEET}
INPUT: 1_CellRanger_output/<SAMPLE>/outs/

SAMPLES PROCESSED:
EOF
    fi
}

init_readme

# ==============================================================================
# ACTIVATE R ENVIRONMENT
# ==============================================================================

log_msg "Activating R environment..."

source /scicore/home/doetsch/kaiser0001/miniforge3/etc/profile.d/conda.sh
conda activate /scicore/home/doetsch/kaiser0001/miniforge3/envs/R_4_5

log_msg "R: $(which R)"
R --version | head -n 1
echo ""

# ==============================================================================
# VERIFY INPUTS
# ==============================================================================

log_msg "============================================================================"
log_msg "Input Verification"
log_msg "============================================================================"
echo ""

OUTS_DIR="${CELLRANGER_DIR}/${SAMPLE}/outs"

if [[ ! -d "${OUTS_DIR}" ]]; then
    error_exit "Cell Ranger output not found: ${OUTS_DIR}"
fi

# Check required files
H5_FILE="${OUTS_DIR}/filtered_feature_bc_matrix.h5"
BAM_FILE="${OUTS_DIR}/possorted_genome_bam.bam"
BAI_FILE="${OUTS_DIR}/possorted_genome_bam.bam.bai"

log_msg "Checking required files:"

for f in "${H5_FILE}" "${BAM_FILE}" "${BAI_FILE}"; do
    if [[ -f "$f" ]]; then
        size=$(du -h "$f" | cut -f1)
        log_msg "  OK: $(basename "$f") (${size})"
    else
        error_exit "Missing: $(basename "$f")"
    fi
done

BAM_SIZE=$(du -h "${BAM_FILE}" | cut -f1)
log_msg "BAM file size: ${BAM_SIZE}"
echo ""

# ==============================================================================
# CREATE R SCRIPT
# ==============================================================================

R_SCRIPT="${SCRIPT_DIR}/run_dropletqc_${SAMPLE}.R"

cat > "${R_SCRIPT}" << 'RSCRIPT_EOF'
#!/usr/bin/env Rscript
# ============================================================================
# DropletQC Nuclear Fraction & Empty Droplet Identification
# Universal Modular Version
# ============================================================================

suppressPackageStartupMessages({
    library(DropletQC)
    library(Matrix)
    library(DropletUtils)
    library(ggplot2)
    library(dplyr)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: run_dropletqc.R <sample> <outs_dir> <output_dir> <cores> <tiles>")
}

SAMPLE <- args[1]
OUTS_DIR <- args[2]
OUTPUT_DIR <- args[3]
N_CORES <- as.integer(args[4])
N_TILES <- as.integer(args[5])

cat("\n============================================================================\n")
cat("DropletQC Analysis:", SAMPLE, "\n")
cat("============================================================================\n\n")

cat("Configuration:\n")
cat("  Sample:     ", SAMPLE, "\n")
cat("  Outs dir:   ", OUTS_DIR, "\n")
cat("  Output dir: ", OUTPUT_DIR, "\n")
cat("  Cores:      ", N_CORES, "\n")
cat("  Tiles:      ", N_TILES, "\n\n")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Step 1: Calculate Nuclear Fraction
# ---------------------------------------------------------------------------

cat("Step 1: Calculating nuclear fraction from Cell Ranger BAM tags...\n")

start_time <- Sys.time()

nf <- nuclear_fraction_tags(
    outs = OUTS_DIR,
    tiles = N_TILES,
    cores = N_CORES,
    verbose = TRUE
)

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\nNuclear fraction completed in", round(elapsed, 1), "minutes\n")

cat("\nNuclear fraction statistics:\n")
print(summary(nf$nuclear_fraction))

nf_file <- file.path(OUTPUT_DIR, paste0(SAMPLE, "_nuclear_fraction.csv"))
write.csv(nf, nf_file, row.names = TRUE)
cat("Saved:", nf_file, "\n\n")

# ---------------------------------------------------------------------------
# Step 2: Get UMI counts
# ---------------------------------------------------------------------------

cat("Step 2: Loading UMI counts...\n")

h5_file <- file.path(OUTS_DIR, "filtered_feature_bc_matrix.h5")
sce <- read10xCounts(h5_file, col.names = TRUE)
umi_counts <- colSums(counts(sce))

cat("  Loaded", length(umi_counts), "cells\n")
cat("  UMI range:", min(umi_counts), "-", max(umi_counts), "\n\n")

# ---------------------------------------------------------------------------
# Step 3: Combine data
# ---------------------------------------------------------------------------

cat("Step 3: Combining nuclear fraction with UMI counts...\n")

nf_barcodes <- rownames(nf)
umi_barcodes <- names(umi_counts)

# Handle barcode suffix differences
common_barcodes <- intersect(nf_barcodes, umi_barcodes)

if (length(common_barcodes) == 0) {
    nf_barcodes_trim <- sub("-1$", "", nf_barcodes)
    umi_barcodes_trim <- sub("-1$", "", umi_barcodes)
    common_trim <- intersect(nf_barcodes_trim, umi_barcodes_trim)
    
    if (length(common_trim) > 0) {
        rownames(nf) <- nf_barcodes_trim
        names(umi_counts) <- umi_barcodes_trim
        common_barcodes <- common_trim
    } else {
        stop("No common barcodes found!")
    }
}

nf_umi <- data.frame(
    barcode = common_barcodes,
    nuclear_fraction = nf[common_barcodes, "nuclear_fraction"],
    umi_count = umi_counts[common_barcodes],
    row.names = common_barcodes
)

cat("  Combined:", nrow(nf_umi), "cells\n\n")

# ---------------------------------------------------------------------------
# Step 4: Identify empty droplets
# ---------------------------------------------------------------------------

cat("Step 4: Identifying empty droplets...\n")

nf_umi_ed <- identify_empty_drops(
    nf_umi = nf_umi[, c("nuclear_fraction", "umi_count")]
)

n_empty <- sum(nf_umi_ed$cell_status == "empty_droplet")
n_cells <- sum(nf_umi_ed$cell_status == "cell")
pct_empty <- round(100 * n_empty / nrow(nf_umi_ed), 2)

cat("  Empty droplets:", n_empty, "(", pct_empty, "%)\n")
cat("  Cells:", n_cells, "\n\n")

# ---------------------------------------------------------------------------
# Step 5: Identify damaged cells
# ---------------------------------------------------------------------------

cat("Step 5: Identifying damaged cells...\n")

nf_umi_ed$damaged_status <- "intact"

umi_threshold <- median(nf_umi_ed$umi_count) / 2
nf_threshold_high <- 0.8

damaged_idx <- which(
    nf_umi_ed$nuclear_fraction > nf_threshold_high &
    nf_umi_ed$umi_count < umi_threshold &
    nf_umi_ed$cell_status == "cell"
)

nf_umi_ed$damaged_status[damaged_idx] <- "damaged"

n_damaged <- sum(nf_umi_ed$damaged_status == "damaged")
pct_damaged <- round(100 * n_damaged / n_cells, 2)

cat("  Damaged cells:", n_damaged, "(", pct_damaged, "% of cells)\n\n")

# ---------------------------------------------------------------------------
# Step 6: Create final status
# ---------------------------------------------------------------------------

cat("Step 6: Creating final status...\n")

nf_umi_ed$final_status <- ifelse(
    nf_umi_ed$cell_status == "empty_droplet", "empty_droplet",
    ifelse(nf_umi_ed$damaged_status == "damaged", "damaged_cell", "intact_cell")
)

# Also create a simple status column for downstream compatibility
nf_umi_ed$status <- ifelse(nf_umi_ed$final_status == "intact_cell", "cell", nf_umi_ed$final_status)

print(table(nf_umi_ed$final_status))
cat("\n")

# ---------------------------------------------------------------------------
# Step 7: Save results
# ---------------------------------------------------------------------------

cat("Step 7: Saving results...\n")

# Full results
results_file <- file.path(OUTPUT_DIR, paste0(SAMPLE, "_dropletqc_results.csv"))
write.csv(nf_umi_ed, results_file, row.names = TRUE)

# Barcode lists
empty_barcodes <- rownames(nf_umi_ed)[nf_umi_ed$cell_status == "empty_droplet"]
writeLines(empty_barcodes, file.path(OUTPUT_DIR, paste0(SAMPLE, "_empty_droplet_barcodes.txt")))

damaged_barcodes <- rownames(nf_umi_ed)[nf_umi_ed$damaged_status == "damaged"]
writeLines(damaged_barcodes, file.path(OUTPUT_DIR, paste0(SAMPLE, "_damaged_cell_barcodes.txt")))

keep_barcodes <- rownames(nf_umi_ed)[nf_umi_ed$final_status == "intact_cell"]
writeLines(keep_barcodes, file.path(OUTPUT_DIR, paste0(SAMPLE, "_intact_cell_barcodes.txt")))

# Summary
summary_df <- data.frame(
    sample = SAMPLE,
    total_barcodes = nrow(nf_umi_ed),
    empty_droplets = n_empty,
    damaged_cells = n_damaged,
    intact_cells = length(keep_barcodes),
    pct_empty = pct_empty,
    pct_damaged = pct_damaged,
    median_nf_empty = median(nf_umi_ed$nuclear_fraction[nf_umi_ed$cell_status == "empty_droplet"], na.rm = TRUE),
    median_nf_cell = median(nf_umi_ed$nuclear_fraction[nf_umi_ed$cell_status == "cell"], na.rm = TRUE),
    median_umi = median(nf_umi_ed$umi_count)
)
write.csv(summary_df, file.path(OUTPUT_DIR, paste0(SAMPLE, "_dropletqc_summary.csv")), row.names = FALSE)

cat("  Saved all results\n\n")

# ---------------------------------------------------------------------------
# Step 8: Generate plots
# ---------------------------------------------------------------------------

cat("Step 8: Generating plots...\n")

# Scatter plot
p1 <- ggplot(nf_umi_ed, aes(x = nuclear_fraction, y = log10(umi_count + 1), color = final_status)) +
    geom_point(alpha = 0.5, size = 0.5) +
    scale_color_manual(values = c(
        "empty_droplet" = "#E41A1C",
        "damaged_cell" = "#FF7F00",
        "intact_cell" = "#4DAF4A"
    )) +
    labs(title = paste("DropletQC -", SAMPLE),
         x = "Nuclear Fraction", y = "log10(UMI + 1)", color = "Status") +
    theme_bw() + theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, paste0(SAMPLE, "_dropletqc_scatter.png")), p1, width = 8, height = 6, dpi = 150)

# Histogram
p2 <- ggplot(nf_umi_ed, aes(x = nuclear_fraction, fill = final_status)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c(
        "empty_droplet" = "#E41A1C",
        "damaged_cell" = "#FF7F00",
        "intact_cell" = "#4DAF4A"
    )) +
    labs(title = paste("Nuclear Fraction -", SAMPLE), x = "Nuclear Fraction", y = "Count") +
    theme_bw() + theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, paste0(SAMPLE, "_nuclear_fraction_histogram.png")), p2, width = 8, height = 6, dpi = 150)

# Violin plot
p3 <- ggplot(nf_umi_ed, aes(x = final_status, y = log10(umi_count + 1), fill = final_status)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = c(
        "empty_droplet" = "#E41A1C",
        "damaged_cell" = "#FF7F00",
        "intact_cell" = "#4DAF4A"
    )) +
    labs(title = paste("UMI by Status -", SAMPLE), x = "Status", y = "log10(UMI + 1)") +
    theme_bw() + theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, paste0(SAMPLE, "_umi_violin.png")), p3, width = 6, height = 6, dpi = 150)

cat("  Saved plots\n\n")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

cat("============================================================================\n")
cat("DropletQC Complete:", SAMPLE, "\n")
cat("============================================================================\n\n")

cat("Summary:\n")
cat("  Total barcodes: ", nrow(nf_umi_ed), "\n")
cat("  Empty droplets: ", n_empty, " (", pct_empty, "%)\n")
cat("  Damaged cells:  ", n_damaged, " (", pct_damaged, "%)\n")
cat("  Intact cells:   ", length(keep_barcodes), "\n\n")

cat("Done!\n")
RSCRIPT_EOF

chmod +x "${R_SCRIPT}"

# ==============================================================================
# RUN DROPLETQC
# ==============================================================================

log_msg "============================================================================"
log_msg "Running DropletQC"
log_msg "============================================================================"
echo ""

Rscript "${R_SCRIPT}" \
    "${SAMPLE}" \
    "${OUTS_DIR}" \
    "${SAMPLE_OUTPUT_DIR}" \
    "${N_CORES}" \
    "${N_TILES}"

R_EXIT=$?
END_TIME=$(date '+%H:%M:%S')

# ==============================================================================
# UPDATE README
# ==============================================================================

update_readme() {
    local sample=$1
    local status=$2
    
    (
        flock -x 200
        if [[ -f "${SAMPLE_OUTPUT_DIR}/${sample}_dropletqc_summary.csv" ]]; then
            TOTAL=$(awk -F',' 'NR==2 {print $2}' "${SAMPLE_OUTPUT_DIR}/${sample}_dropletqc_summary.csv")
            EMPTY=$(awk -F',' 'NR==2 {print $3}' "${SAMPLE_OUTPUT_DIR}/${sample}_dropletqc_summary.csv")
            DAMAGED=$(awk -F',' 'NR==2 {print $4}' "${SAMPLE_OUTPUT_DIR}/${sample}_dropletqc_summary.csv")
            INTACT=$(awk -F',' 'NR==2 {print $5}' "${SAMPLE_OUTPUT_DIR}/${sample}_dropletqc_summary.csv")
            echo "  ${sample}: ${status} | Total: ${TOTAL}, Empty: ${EMPTY}, Damaged: ${DAMAGED}, Intact: ${INTACT}" >> "${README_FILE}"
        else
            echo "  ${sample}: ${status}" >> "${README_FILE}"
        fi
    ) 200>"${README_FILE}.lock"
}

if [[ $R_EXIT -eq 0 ]]; then
    update_readme "${SAMPLE}" "SUCCESS"
    log_msg "Status: SUCCESS"
else
    update_readme "${SAMPLE}" "FAILED"
    error_exit "DropletQC failed"
fi

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

echo ""
log_msg "============================================================================"
log_msg "DropletQC Complete - ${SAMPLE}"
log_msg "============================================================================"
log_msg "Output: ${SAMPLE_OUTPUT_DIR}"
echo ""
ls -lh "${SAMPLE_OUTPUT_DIR}/"
echo ""
log_msg "============================================================================"
