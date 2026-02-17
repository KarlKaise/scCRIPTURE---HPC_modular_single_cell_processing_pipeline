#!/bin/bash
#SBATCH --job-name=DecontX
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --output=logs/decontx/decontx_%A_%a.out
#SBATCH --error=logs/decontx/decontx_%A_%a.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ==============================================================================
# STEP 8: DecontX Ambient RNA Correction
# ==============================================================================
#
# UNIVERSAL MODULAR VERSION - Reads samples from samplesheet.csv
#
# This step applies DecontX (from celda package) to remove ambient RNA
# contamination from scCDC-corrected counts.
#
# INPUT:
#   - 7_scCDC_correction/<SAMPLE>/<SAMPLE>_scCDC_corrected.rds
#   - Seurat v5 object with 'scCDC_corrected' layer in RNA assay
#
# OUTPUT:
#   - 8_DecontX_correction/<SAMPLE>/<SAMPLE>_decontX_corrected.rds
#   - Seurat v5 object with 'scCDC_corrected' layer UPDATED with DecontX results
#   - DecontX metadata added (contamination estimates)
#
# LAYER STRUCTURE (preserved, scCDC_corrected layer updated):
#   RNA[[]]
#   ├── counts           - Original raw counts (preserved)
#   ├── data             - Normalized data (preserved)
#   ├── scale.data       - Scaled data (preserved)
#   └── scCDC_corrected  - Now contains DecontX-corrected counts (UPDATED)
#
# NEW METADATA COLUMNS:
#   - decontX_contamination: Per-cell contamination estimate (0-1)
#
# Usage:
#   sbatch --array=1-N 08_DecontX_correction.sh
#
# ==============================================================================

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

# Input/Output directories
INPUT_DIR="${PREPROCESS_DIR}/7_scCDC_correction"
OUTPUT_DIR="${PREPROCESS_DIR}/8_DecontX_correction"
SCRIPT_DIR="${BASE_DIR}/Scripts/R_scripts"
LOG_DIR="${BASE_DIR}/logs/decontx"

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

get_sample_by_index() {
    local idx=$1
    local col=$(get_col_index "$SAMPLESHEET" "sample_id")
    tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | sed -n "${idx}p" | cut -d',' -f"$col"
}

# ==============================================================================
# DETERMINE SAMPLE
# ==============================================================================

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    SAMPLE=$(get_sample_by_index "$SLURM_ARRAY_TASK_ID")
else
    SAMPLE="${1:-}"
fi

if [[ -z "$SAMPLE" ]]; then
    error_exit "No sample specified. Use SLURM array or provide sample as argument."
fi

# ==============================================================================
# SETUP
# ==============================================================================

SAMPLE_INPUT_DIR="${INPUT_DIR}/${SAMPLE}"
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}"

INPUT_RDS="${SAMPLE_INPUT_DIR}/${SAMPLE}_scCDC_corrected.rds"
OUTPUT_RDS="${SAMPLE_OUTPUT_DIR}/${SAMPLE}_decontX_corrected.rds"
SUMMARY_FILE="${SAMPLE_OUTPUT_DIR}/${SAMPLE}_decontX_summary.csv"
PLOT_FILE="${SAMPLE_OUTPUT_DIR}/${SAMPLE}_decontX_contamination.png"

mkdir -p "${SAMPLE_OUTPUT_DIR}"
mkdir -p "${SCRIPT_DIR}"
mkdir -p "${LOG_DIR}"

# ==============================================================================
# PRINT INFO
# ==============================================================================

log_msg "============================================================================"
log_msg "DecontX Ambient RNA Correction - Step 8"
log_msg "============================================================================"
log_msg "Sample:        ${SAMPLE}"
log_msg "Job ID:        ${SLURM_JOB_ID:-local}"
log_msg "Array Task:    ${SLURM_ARRAY_TASK_ID:-N/A}"
log_msg "Date:          $(date)"
log_msg "Host:          $(hostname)"
log_msg "CPUs:          ${SLURM_CPUS_PER_TASK:-8}"
log_msg "Memory:        ${SLURM_MEM_PER_NODE:-64G}"
log_msg "============================================================================"
log_msg "Input:         ${INPUT_RDS}"
log_msg "Output:        ${OUTPUT_RDS}"
log_msg "============================================================================"
echo ""

# Validate input
if [[ ! -f "${INPUT_RDS}" ]]; then
    error_exit "Input file not found: ${INPUT_RDS}"
fi

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

log_msg "Setting up environment..."

source /scicore/home/doetsch/kaiser0001/miniforge3/etc/profile.d/conda.sh
conda activate /scicore/home/doetsch/kaiser0001/miniforge3/envs/R_4_5

# Force UTF-8 locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8

log_msg "R: $(which R)"
R --version | head -n 1
log_msg "LANG=${LANG}"
echo ""

# ==============================================================================
# CREATE R SCRIPT
# ==============================================================================

R_SCRIPT="${SCRIPT_DIR}/run_decontX_${SAMPLE}.R"

cat > "${R_SCRIPT}" << 'RSCRIPT_EOF'
#!/usr/bin/env Rscript
# ==============================================================================
# DecontX Ambient RNA Correction
# ==============================================================================
#
# Reads scCDC_corrected layer, applies DecontX via SCE, updates scCDC_corrected layer
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(decontX)
  library(ggplot2)
  library(Matrix)
})

# ------------------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: run_decontX.R <input_rds> <output_rds> <summary_csv> <plot_png>")
}

input_rds   <- args[1]
output_rds  <- args[2]
summary_csv <- args[3]
plot_png    <- args[4]

cat("==============================================================================\n")
cat("DecontX Ambient RNA Correction\n")
cat("==============================================================================\n\n")

cat("Input:  ", input_rds, "\n")
cat("Output: ", output_rds, "\n\n")

# ------------------------------------------------------------------------------
# Load Seurat object
# ------------------------------------------------------------------------------
cat("Loading Seurat object...\n")
seu <- readRDS(input_rds)

cat("  Cells:    ", ncol(seu), "\n")
cat("  Features: ", nrow(seu), "\n")
cat("  Assays:   ", paste(names(seu@assays), collapse = ", "), "\n")

# Check for RNA assay
if (!"RNA" %in% names(seu@assays)) {
  stop("ERROR: RNA assay not found in Seurat object")
}

# Check layers
rna_assay <- seu[["RNA"]]
available_layers <- Layers(rna_assay)
cat("  Available layers: ", paste(available_layers, collapse = ", "), "\n\n")

# Verify scCDC_corrected layer exists
if (!"scCDC_corrected" %in% available_layers) {
  stop("ERROR: 'scCDC_corrected' layer not found. Available layers: ", 
       paste(available_layers, collapse = ", "))
}

# ------------------------------------------------------------------------------
# Extract scCDC-corrected counts
# ------------------------------------------------------------------------------
cat("Extracting scCDC_corrected counts...\n")
counts_sccdc <- LayerData(seu, assay = "RNA", layer = "scCDC_corrected")

cat("  Class: ", class(counts_sccdc)[1], "\n")
cat("  Dimensions: ", nrow(counts_sccdc), " x ", ncol(counts_sccdc), "\n")
cat("  Total UMIs (scCDC): ", format(sum(counts_sccdc), big.mark = ","), "\n\n")

# Ensure sparse matrix format
if (!inherits(counts_sccdc, "dgCMatrix")) {
  cat("  Converting to dgCMatrix...\n")
  counts_sccdc <- as(counts_sccdc, "dgCMatrix")
}

# ------------------------------------------------------------------------------
# Create SingleCellExperiment for decontX
# ------------------------------------------------------------------------------
cat("Creating SingleCellExperiment object...\n")

sce <- SingleCellExperiment(list(counts = counts_sccdc))

# Add cluster information if available (improves decontX performance)
if ("seurat_clusters" %in% colnames(seu[[]])) {
  sce$clusters <- as.integer(as.factor(seu$seurat_clusters))
  cat("  Added 'seurat_clusters' (", length(unique(sce$clusters)), " clusters)\n")
} else if ("RNA_snn_res.0.5" %in% colnames(seu[[]])) {
  sce$clusters <- as.integer(as.factor(seu$RNA_snn_res.0.5))
  cat("  Added 'RNA_snn_res.0.5' (", length(unique(sce$clusters)), " clusters)\n")
}

cat("  SCE dimensions: ", nrow(sce), " x ", ncol(sce), "\n\n")

# ------------------------------------------------------------------------------
# Run DecontX
# ------------------------------------------------------------------------------
cat("Running decontX...\n")
cat("  This may take several minutes...\n\n")

set.seed(42)

# Run decontX - it will use clusters if present in colData
sce <- tryCatch({
  if ("clusters" %in% colnames(colData(sce))) {
    decontX(sce, z = sce$clusters, verbose = TRUE)
  } else {
    decontX(sce, verbose = TRUE)
  }
}, error = function(e) {
  cat("ERROR in decontX: ", conditionMessage(e), "\n")
  stop(e)
})

cat("\ndecontX completed.\n\n")

# ------------------------------------------------------------------------------
# Extract results
# ------------------------------------------------------------------------------
# Get decontaminated counts using the accessor function
decontx_counts <- decontXcounts(sce)

# Get per-cell contamination estimates
contamination <- sce$decontX_contamination

cat("decontX Results:\n")
cat("  Total UMIs (input):  ", format(sum(counts_sccdc), big.mark = ","), "\n")
cat("  Total UMIs (output): ", format(sum(decontx_counts), big.mark = ","), "\n")
cat("  UMI reduction: ", format(sum(counts_sccdc) - sum(decontx_counts), big.mark = ","), 
    " (", round((1 - sum(decontx_counts)/sum(counts_sccdc)) * 100, 2), "%)\n")
cat("  Mean contamination:   ", round(mean(contamination), 4), "\n")
cat("  Median contamination: ", round(median(contamination), 4), "\n")
cat("  Contamination range:  [", round(min(contamination), 4), ", ", 
    round(max(contamination), 4), "]\n\n")

# ------------------------------------------------------------------------------
# Update Seurat object
# ------------------------------------------------------------------------------
cat("Updating Seurat object...\n")

# Ensure decontx_counts is sparse and has correct dimnames
decontx_counts <- as(decontx_counts, "dgCMatrix")
rownames(decontx_counts) <- rownames(counts_sccdc)
colnames(decontx_counts) <- colnames(counts_sccdc)

# Update the scCDC_corrected layer with decontX results
LayerData(seu, assay = "RNA", layer = "scCDC_corrected") <- decontx_counts

# Add contamination estimates to metadata
seu$decontX_contamination <- contamination

# Add decontX clusters if they were generated
if ("decontX_clusters" %in% colnames(colData(sce))) {
  seu$decontX_clusters <- sce$decontX_clusters
}

# Verify the update
updated_layers <- Layers(seu[["RNA"]])
cat("  Updated layers: ", paste(updated_layers, collapse = ", "), "\n")

# Verify layer was updated
updated_counts <- LayerData(seu, assay = "RNA", layer = "scCDC_corrected")
cat("  Verification - scCDC_corrected layer sum: ", format(sum(updated_counts), big.mark = ","), "\n\n")

# ------------------------------------------------------------------------------
# Create summary statistics
# ------------------------------------------------------------------------------
cat("Creating summary statistics...\n")

summary_df <- data.frame(
  metric = c(
    "n_cells",
    "n_features",
    "total_umi_input",
    "total_umi_output",
    "umi_removed",
    "umi_reduction_pct",
    "mean_contamination",
    "median_contamination",
    "min_contamination",
    "max_contamination",
    "sd_contamination",
    "q25_contamination",
    "q75_contamination"
  ),
  value = c(
    ncol(seu),
    nrow(seu),
    sum(counts_sccdc),
    sum(decontx_counts),
    sum(counts_sccdc) - sum(decontx_counts),
    round((1 - sum(decontx_counts)/sum(counts_sccdc)) * 100, 4),
    round(mean(contamination), 6),
    round(median(contamination), 6),
    round(min(contamination), 6),
    round(max(contamination), 6),
    round(sd(contamination), 6),
    round(quantile(contamination, 0.25), 6),
    round(quantile(contamination, 0.75), 6)
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, summary_csv, row.names = FALSE)
cat("  Saved: ", summary_csv, "\n")

# ------------------------------------------------------------------------------
# Create contamination plot
# ------------------------------------------------------------------------------
cat("Creating contamination plot...\n")

plot_df <- data.frame(
  contamination = contamination,
  nCount = seu$nCount_RNA
)

p <- ggplot(plot_df, aes(x = contamination)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = mean(contamination), color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = median(contamination), color = "darkgreen", linetype = "dashed", linewidth = 1) +
  labs(
    title = paste0("decontX Contamination Estimates"),
    subtitle = paste0("Mean: ", round(mean(contamination), 4), 
                      " | Median: ", round(median(contamination), 4),
                      " | n = ", ncol(seu), " cells"),
    x = "Contamination Fraction",
    y = "Number of Cells"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  ) +
  annotate("text", x = mean(contamination), y = Inf, label = "Mean", 
           vjust = 2, hjust = -0.1, color = "red", size = 3) +
  annotate("text", x = median(contamination), y = Inf, label = "Median", 
           vjust = 3.5, hjust = -0.1, color = "darkgreen", size = 3)

ggsave(plot_png, p, width = 8, height = 6, dpi = 150)
cat("  Saved: ", plot_png, "\n\n")

# ------------------------------------------------------------------------------
# Save output
# ------------------------------------------------------------------------------
cat("Saving Seurat object...\n")
saveRDS(seu, output_rds)
cat("  Saved: ", output_rds, "\n\n")

# ------------------------------------------------------------------------------
# Final summary
# ------------------------------------------------------------------------------
cat("==============================================================================\n")
cat("decontX Correction Complete\n")
cat("==============================================================================\n")
cat("Input:  ", input_rds, "\n")
cat("Output: ", output_rds, "\n")
cat("Cells:  ", ncol(seu), "\n")
cat("UMI reduction: ", round((1 - sum(decontx_counts)/sum(counts_sccdc)) * 100, 2), "%\n")
cat("Mean contamination: ", round(mean(contamination), 4), "\n")
cat("==============================================================================\n")

cat("\nDone!\n")
RSCRIPT_EOF

chmod +x "${R_SCRIPT}"
log_msg "R script written to: ${R_SCRIPT}"

# ==============================================================================
# RUN R SCRIPT
# ==============================================================================

log_msg "============================================================================"
log_msg "Running DecontX for sample: ${SAMPLE}"
log_msg "============================================================================"
echo ""

Rscript "${R_SCRIPT}" \
    "${INPUT_RDS}" \
    "${OUTPUT_RDS}" \
    "${SUMMARY_FILE}" \
    "${PLOT_FILE}"

R_EXIT=$?

# ==============================================================================
# CLEANUP AND SUMMARY
# ==============================================================================

if [[ $R_EXIT -eq 0 ]]; then
    log_msg "============================================================================"
    log_msg "DecontX Correction Complete"
    log_msg "============================================================================"
    log_msg "Sample: ${SAMPLE}"
    log_msg "Output: ${OUTPUT_RDS}"
    
    if [[ -f "${SUMMARY_FILE}" ]]; then
        echo ""
        log_msg "Summary statistics:"
        cat "${SUMMARY_FILE}"
    fi
    
    echo ""
    log_msg "Output files:"
    ls -lh "${SAMPLE_OUTPUT_DIR}/"
    log_msg "============================================================================"
else
    error_exit "DecontX failed for sample ${SAMPLE}"
fi

exit $R_EXIT
