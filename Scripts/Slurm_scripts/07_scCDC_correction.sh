#!/bin/bash
#SBATCH --job-name=scCDC
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --output=logs/scCDC/scCDC_%a_%A.out
#SBATCH --error=logs/scCDC/scCDC_%a_%A.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ==============================================================================
# STEP 7: scCDC Contamination Correction Pipeline
# ==============================================================================
#
# UNIVERSAL MODULAR VERSION - Reads samples from samplesheet.csv
#
# Input:  6_Doublet_consensus/<SAMPLE>/<SAMPLE>_doublets_removed.rds
# Output: 7_scCDC_correction/<SAMPLE>/<SAMPLE>_scCDC_corrected.rds
#
# Usage:
#   sbatch --array=1-N 07_scCDC_correction.sh
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

# Directories
INPUT_ROOT="${PREPROCESS_DIR}/6_Doublet_consensus"
OUTPUT_ROOT="${PREPROCESS_DIR}/7_scCDC_correction"
SCRIPT_DIR="${BASE_DIR}/Scripts/R_scripts"
LOG_DIR="${BASE_DIR}/logs/scCDC"
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
log_msg "scCDC Contamination Correction - Modular Pipeline"
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
INPUT_RDS="${INPUT_ROOT}/${SAMPLE}/${SAMPLE}_doublets_removed.rds"
SAMPLE_OUTPUT_DIR="${OUTPUT_ROOT}/${SAMPLE}"
OUTPUT_RDS="${SAMPLE_OUTPUT_DIR}/${SAMPLE}_scCDC_corrected.rds"

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
STEP 7: scCDC Contamination Correction
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')
Pipeline: Universal Modular scRNA-seq Pipeline

DESCRIPTION:
  Runs scCDC (single-cell Contamination Detection and Correction) to:
    1. Detect Globally Contaminating Genes (GCGs)
    2. Quantify contamination ratios per cell
    3. Correct count matrix for contamination

SAMPLESHEET: ${SAMPLESHEET}
INPUT: 6_Doublet_consensus/<SAMPLE>/<SAMPLE>_doublets_removed.rds

SAMPLES PROCESSED:
EOF
fi

# ==============================================================================
# VERIFY INPUT
# ==============================================================================

if [[ ! -f "${INPUT_RDS}" ]]; then
    error_exit "Input file not found: ${INPUT_RDS}"
fi

INPUT_SIZE=$(du -h "${INPUT_RDS}" | cut -f1)
log_msg "Input: ${INPUT_RDS} (${INPUT_SIZE})"
log_msg "Output: ${OUTPUT_RDS}"
echo ""

# ==============================================================================
# ACTIVATE R ENVIRONMENT
# ==============================================================================

log_msg "Activating R environment..."

source /scicore/home/doetsch/kaiser0001/miniforge3/etc/profile.d/conda.sh
conda activate /scicore/home/doetsch/kaiser0001/miniforge3/envs/R_4_5

log_msg "R: $(which R)"
R --version | head -n 1

# Check scCDC
log_msg "Verifying scCDC package..."
Rscript -e "if (!requireNamespace('scCDC', quietly = TRUE)) { stop('scCDC not installed') } else { cat('scCDC:', as.character(packageVersion('scCDC')), '\n') }"
echo ""

# ==============================================================================
# CREATE R SCRIPT
# ==============================================================================

R_SCRIPT="${SCRIPT_DIR}/run_scCDC_correction_${SAMPLE}.R"

cat > "${R_SCRIPT}" << 'RSCRIPT_EOF'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(scCDC)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: run_scCDC_correction.R <input_rds> <output_rds>\n")
}

input_rds  <- normalizePath(args[[1]], mustWork = TRUE)
output_rds <- args[[2]]

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)

log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(paste0("[", timestamp, "] ", msg))
}

log_message("=== scCDC Contamination Correction Pipeline ===")
log_message(paste0("Input:  ", input_rds))
log_message(paste0("Output: ", output_rds))
log_message("")

log_message("Loading Seurat object...")
seurat_obj <- readRDS(input_rds)

input_cells <- ncol(seurat_obj)
input_genes <- nrow(seurat_obj)
log_message(paste0("Input: ", input_genes, " genes x ", input_cells, " cells"))
log_message(paste0("Default assay: ", DefaultAssay(seurat_obj)))

if (!"RNA" %in% Assays(seurat_obj)) {
  stop("No RNA assay found!")
}

DefaultAssay(seurat_obj) <- "RNA"

available_layers <- Layers(seurat_obj, assay = "RNA")
log_message(paste0("Available layers: ", paste(available_layers, collapse = ", ")))

if (!"counts" %in% available_layers) {
  stop("No counts layer found in RNA assay!")
}

log_message("Preparing data for scCDC...")

if (!"data" %in% available_layers) {
  log_message("Normalizing data (LogNormalize)...")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize",
                               scale.factor = 10000, verbose = FALSE)
}

if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  log_message("Running clustering for scCDC...")

  if (length(VariableFeatures(seurat_obj)) == 0) {
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                        nfeatures = 2000, verbose = FALSE)
  }

  available_layers <- Layers(seurat_obj, assay = "RNA")
  if (!"scale.data" %in% available_layers) {
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  }

  if (!"pca" %in% Reductions(seurat_obj)) {
    seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
  }

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
}

log_message(paste0("Number of clusters: ", length(unique(Idents(seurat_obj)))))

log_message("")
log_message("=== scCDC Contamination Detection ===")

n_GCGs <- 0

tryCatch({
  GCGs <- ContaminationDetection(seurat_obj)
  n_GCGs <- nrow(GCGs)
  log_message(paste0("Globally Contaminating Genes (GCGs): ", n_GCGs))

  if (n_GCGs == 0) {
    log_message("No contamination detected, copying counts as scCDC_corrected layer")
    seurat_obj[["RNA"]]$scCDC_corrected <- seurat_obj[["RNA"]]$counts

  } else {
    sample_base <- gsub("_scCDC_corrected\\.rds$", "", basename(output_rds))
    gcg_file <- file.path(dirname(output_rds), paste0(sample_base, "_GCGs.csv"))
    write.csv(GCGs, gcg_file, row.names = TRUE)
    log_message(paste0("GCGs saved to: ", basename(gcg_file)))

    log_message("")
    log_message("=== scCDC Contamination Quantification ===")
    contamination_ratio <- ContaminationQuantification(seurat_obj, rownames(GCGs))
    log_message("Contamination quantification complete")

    contam_file <- file.path(dirname(output_rds), paste0(sample_base, "_contamination_ratio.csv"))
    write.csv(contamination_ratio, contam_file, row.names = TRUE)
    log_message(paste0("Contamination ratios saved to: ", basename(contam_file)))

    log_message("")
    log_message("=== scCDC Contamination Correction ===")
    seurat_corrected <- ContaminationCorrection(seurat_obj, rownames(GCGs))
    log_message("Contamination correction complete")

    log_message("")
    log_message("=== Transferring corrected counts ===")

    if ("Corrected" %in% Assays(seurat_corrected)) {
      corrected_counts <- GetAssayData(seurat_corrected, assay = "Corrected", layer = "counts")
      seurat_obj[["RNA"]]$scCDC_corrected <- corrected_counts
      log_message("Added 'scCDC_corrected' layer to RNA assay")

      if ("contamination_ratio" %in% colnames(seurat_corrected@meta.data)) {
        seurat_obj$scCDC_contamination_ratio <- seurat_corrected$contamination_ratio
      }
    } else {
      log_message("WARNING: Corrected assay not found")
      seurat_obj[["RNA"]]$scCDC_corrected <- seurat_obj[["RNA"]]$counts
    }

    rm(seurat_corrected)
  }

}, error = function(e) {
  log_message(paste0("ERROR in scCDC: ", conditionMessage(e)))
  log_message("Saving original with counts copy as scCDC_corrected")
  seurat_obj[["RNA"]]$scCDC_corrected <- seurat_obj[["RNA"]]$counts
})

log_message("")
log_message("=== Final Output ===")

output_cells <- ncol(seurat_obj)
output_genes <- nrow(seurat_obj)
final_layers <- Layers(seurat_obj, assay = "RNA")

log_message(paste0("Final: ", output_genes, " genes x ", output_cells, " cells"))
log_message(paste0("RNA layers: ", paste(final_layers, collapse = ", ")))

if ("scCDC_corrected" %in% final_layers) {
  corrected_dims <- dim(seurat_obj[["RNA"]]$scCDC_corrected)
  log_message(paste0("scCDC_corrected: ", corrected_dims[1], " x ", corrected_dims[2]))
}

saveRDS(seurat_obj, file = output_rds)
log_message(paste0("Saved: ", output_rds))
log_message(paste0("Size: ", round(file.info(output_rds)$size / (1024^2), 2), " MB"))

# Save summary
summary_file <- file.path(dirname(output_rds), "scCDC_summary.txt")
cat(paste(output_cells, n_GCGs, sep=","), file=summary_file)

log_message("")
log_message(paste0("GCGs detected: ", n_GCGs))
log_message("Done!")
RSCRIPT_EOF

chmod +x "${R_SCRIPT}"

# ==============================================================================
# RUN scCDC
# ==============================================================================

log_msg "============================================================================"
log_msg "Running scCDC Correction"
log_msg "============================================================================"
echo ""

Rscript "${R_SCRIPT}" "${INPUT_RDS}" "${OUTPUT_RDS}"

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
        if [[ -f "${SAMPLE_OUTPUT_DIR}/scCDC_summary.txt" ]]; then
            CELLS=$(cut -d',' -f1 "${SAMPLE_OUTPUT_DIR}/scCDC_summary.txt")
            GCGS=$(cut -d',' -f2 "${SAMPLE_OUTPUT_DIR}/scCDC_summary.txt")
            echo "  ${sample}: ${status} | Cells: ${CELLS}, GCGs: ${GCGS}" >> "${README_FILE}"
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
    error_exit "scCDC failed"
fi

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

echo ""
log_msg "============================================================================"
log_msg "scCDC Complete - ${SAMPLE}"
log_msg "============================================================================"
log_msg "Output: ${SAMPLE_OUTPUT_DIR}"
echo ""
ls -lh "${SAMPLE_OUTPUT_DIR}/"
echo ""
log_msg "============================================================================"
