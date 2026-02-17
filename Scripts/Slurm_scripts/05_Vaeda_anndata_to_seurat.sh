#!/bin/bash
#SBATCH --job-name=h5ad2rds
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --qos=6hours
#SBATCH --output=logs/h5ad2rds/h5ad2rds_%a_%A.out
#SBATCH --error=logs/h5ad2rds/h5ad2rds_%a_%A.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ============================================================================
# STEP 5: Convert VAEDA h5ad to Seurat RDS for DoubletFinder/scDblFinder
# ============================================================================
#
# UNIVERSAL MODULAR VERSION - Reads samples from samplesheet.csv
#
# Input:  4_Vaeda_doublet_detection/<SAMPLE>/<SAMPLE>_qClus_dropletqc_vaeda.h5ad
# Output: 5_Seurat_conversion/<SAMPLE>/<SAMPLE>_qClus_dropletqc_vaeda.rds
#
# Usage:
#   sbatch --array=1-N 05_Vaeda_anndata_to_seurat.sh
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
INPUT_ROOT="${PREPROCESS_DIR}/4_Vaeda_doublet_detection"
OUTPUT_ROOT="${PREPROCESS_DIR}/5_Seurat_conversion"
SCRIPT_DIR="${BASE_DIR}/Scripts/R_scripts"
LOG_DIR="${BASE_DIR}/logs/h5ad2rds"
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
log_msg "H5AD to Seurat Conversion - Modular Pipeline"
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
INPUT_H5AD="${INPUT_ROOT}/${SAMPLE}/${SAMPLE}_qClus_dropletqc_vaeda.h5ad"
SAMPLE_OUTPUT_DIR="${OUTPUT_ROOT}/${SAMPLE}"
OUTPUT_RDS="${SAMPLE_OUTPUT_DIR}/${SAMPLE}_qClus_dropletqc_vaeda.rds"

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
STEP 5: AnnData to Seurat Conversion
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')
Pipeline: Universal Modular scRNA-seq Pipeline

DESCRIPTION:
  Converts VAEDA-annotated AnnData (h5ad) objects to Seurat RDS format
  for downstream R-based analysis (DoubletFinder, scDblFinder, scCDC, CHOIR).

  Uses capseuratconverter for robust conversion with metadata preservation.

SAMPLESHEET: ${SAMPLESHEET}
INPUT: 4_Vaeda_doublet_detection/<SAMPLE>/<SAMPLE>_qClus_dropletqc_vaeda.h5ad

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
echo ""

# ==============================================================================
# CREATE R SCRIPT
# ==============================================================================

R_SCRIPT="${SCRIPT_DIR}/convert_h5ad_to_seurat_${SAMPLE}.R"

cat > "${R_SCRIPT}" <<'REOF'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rhdf5)
  library(capseuratconverter)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: convert_h5ad_to_seurat.R <input_h5ad> <output_rds> [ignore_bad_format TRUE/FALSE]\n")
}
input_h5ad <- normalizePath(args[[1]], mustWork = TRUE)
output_rds <- args[[2]]
ignore_bad_format <- ifelse(length(args) >= 3, as.logical(args[[3]]), TRUE)

# Create output directory
dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper functions
# -----------------------------
nnz <- function(m) {
  if (is.null(m)) return(NA_integer_)
  if (inherits(m, "dgCMatrix")) return(length(m@x))
  if (inherits(m, "Matrix")) return(Matrix::nnzero(m))
  suppressWarnings(sum(m != 0))
}

is_empty_mat <- function(m) {
  is.null(m) || any(dim(m) == 0) || (inherits(m, "dgCMatrix") && length(m@x) == 0)
}

get_layer <- function(seu, layer_name, assay="RNA") {
  tryCatch(GetAssayData(seu, assay = assay, layer = layer_name), error = function(e) NULL)
}

check_h5ad_structure <- function(h5) {
  cat("\n---- H5AD structure check ----\n")
  cat("File: ", h5, "\n", sep="")
  ls <- h5ls(h5, recursive = FALSE)
  top <- unique(ls$name)
  cat("Top-level entries: ", paste(top, collapse=", "), "\n", sep="")
  has <- function(name) any(top == name)
  cat("Has X:      ", has("X"), "\n", sep="")
  cat("Has obs:    ", has("obs"), "\n", sep="")
  cat("Has var:    ", has("var"), "\n", sep="")
  cat("Has obsm:   ", has("obsm"), "\n", sep="")
  cat("Has uns:    ", has("uns"), "\n", sep="")
  cat("Has raw:    ", has("raw"), "\n", sep="")
  cat("Has layers: ", has("layers"), "\n", sep="")
  cat("---- end H5AD check ----\n")
}

summarize_layers <- function(seu, assay="RNA") {
  cat("Assay class: ", class(seu[[assay]])[1], "\n", sep="")
  lay_names <- tryCatch(Layers(seu[[assay]]), error=function(e) character())
  cat("Layers(): ", if (length(lay_names)) paste(lay_names, collapse=", ") else "(none)", "\n", sep="")

  for (ln in unique(c("counts","data","scale.data", lay_names))) {
    m <- get_layer(seu, ln, assay=assay)
    if (is.null(m)) next
    cat(sprintf("  layer %-12s dim=%dx%d nnz=%s sum=%s\n",
                ln, nrow(m), ncol(m), nnz(m), signif(sum(m), 6)))
  }
}

summarize_seurat <- function(seu, tag="") {
  cat("\n---- Seurat summary ", tag, " ----\n", sep="")
  cat("Cells: ", ncol(seu), "  Genes: ", nrow(seu), "\n", sep="")
  cat("Assays: ", paste(Assays(seu), collapse=", "), "\n", sep="")
  cat("meta.data cols: ", ncol(seu[[]]), "\n", sep="")
  cat("First meta cols: ", paste(head(colnames(seu[[]]), 16), collapse=", "), "\n", sep="")
  reds <- Reductions(seu)
  cat("Reductions: ", ifelse(length(reds) > 0, paste(reds, collapse=", "), "(none)"), "\n", sep="")
  summarize_layers(seu, assay = DefaultAssay(seu))
  cat("---- end Seurat summary ----\n")
}

rebuild_rna_counts_from_data <- function(seu) {
  DefaultAssay(seu) <- "RNA"
  m_counts <- get_layer(seu, "counts", assay="RNA")
  if (!is_empty_mat(m_counts)) return(list(seu=seu, changed=FALSE, msg="counts already non-empty"))

  m_data <- get_layer(seu, "data", assay="RNA")
  if (is_empty_mat(m_data)) stop("counts is empty and data is empty; cannot fix.")

  # Build a VALID assay with counts=m_data
  if (exists("CreateAssay5Object", where=asNamespace("SeuratObject"), inherits=FALSE)) {
    new_assay <- SeuratObject::CreateAssay5Object(counts = m_data)
  } else {
    new_assay <- SeuratObject::CreateAssayObject(counts = m_data)
  }

  # Preserve feature metadata if present
  old_assay <- seu[["RNA"]]
  old_meta <- tryCatch(old_assay@meta.data, error=function(e) NULL)
  if (!is.null(old_meta) && nrow(old_meta) > 0) {
    new_assay@meta.data <- old_meta[rownames(new_assay), , drop=FALSE]
  }

  seu[["RNA"]] <- new_assay
  DefaultAssay(seu) <- "RNA"

  # Validate
  m_counts2 <- get_layer(seu, "counts", assay="RNA")
  if (is_empty_mat(m_counts2)) stop("Fix attempt failed: counts still empty after rebuild.")

  list(seu=seu, changed=TRUE, msg="rebuilt RNA assay: counts <- data")
}

# -----------------------------
# Main conversion
# -----------------------------

cat("\n================================================================================\n")
cat("INPUT:  ", input_h5ad, "\n", sep="")
cat("OUTPUT: ", output_rds, "\n", sep="")
cat("================================================================================\n")

# 1) H5AD layout check
check_h5ad_structure(input_h5ad)

# 2) Convert
cat("\nConverting with capseuratconverter::h5ad_to_seurat() ...\n")
seu <- h5ad_to_seurat(input_h5ad, ignore_bad_format = ignore_bad_format)

# 3) Summarize before fix
DefaultAssay(seu) <- "RNA"
summarize_seurat(seu, tag="BEFORE counts-fix")

# 4) Validate & fix counts if needed
fix <- rebuild_rna_counts_from_data(seu)
seu2 <- fix$seu
changed <- fix$changed

# Validate object
valid_ok <- TRUE
tryCatch({ validObject(seu2) }, error=function(e) { valid_ok <<- FALSE })

summarize_seurat(seu2, tag="AFTER counts-fix")
cat("validObject(seurat): ", valid_ok, "\n", sep="")
cat("counts-fix: ", fix$msg, "\n", sep="")

# Check for DropletQC metadata
dropletqc_cols <- grep("dropletqc", colnames(seu2[[]]), value = TRUE, ignore.case = TRUE)
if (length(dropletqc_cols) > 0) {
  cat("\nDropletQC metadata preserved: ", paste(dropletqc_cols, collapse=", "), "\n", sep="")
} else {
  cat("\nNo DropletQC metadata found in converted object.\n")
}

# Check for VAEDA metadata
vaeda_cols <- grep("vaeda|doublet", colnames(seu2[[]]), value = TRUE, ignore.case = TRUE)
if (length(vaeda_cols) > 0) {
  cat("VAEDA metadata preserved: ", paste(vaeda_cols, collapse=", "), "\n", sep="")
} else {
  cat("No VAEDA metadata found in converted object.\n")
}

# 5) Save .rds
saveRDS(seu2, output_rds)
cat("\nSaved RDS: ", output_rds, "\n", sep="")
cat("RDS size: ", round(file.info(output_rds)$size / (1024^2), 2), " MB\n", sep="")

# 6) Reload sanity check
seu3 <- readRDS(output_rds)
DefaultAssay(seu3) <- "RNA"
m_counts3 <- get_layer(seu3, "counts", assay="RNA")
cat("Reload check -> Cells: ", ncol(seu3), " Genes: ", nrow(seu3),
    " | counts dim=", nrow(m_counts3), "x", ncol(m_counts3),
    " nnz=", nnz(m_counts3), "\n", sep="")

# Save cell/gene counts for README
cat(ncol(seu3), ",", nrow(seu3), "\n", file=paste0(dirname(output_rds), "/conversion_stats.txt"))

cat("\nDONE\n")
REOF

chmod +x "${R_SCRIPT}"

# ==============================================================================
# RUN CONVERSION
# ==============================================================================

log_msg "============================================================================"
log_msg "Running H5AD to Seurat Conversion"
log_msg "============================================================================"
echo ""

Rscript "${R_SCRIPT}" "${INPUT_H5AD}" "${OUTPUT_RDS}" TRUE

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
        if [[ -f "${SAMPLE_OUTPUT_DIR}/conversion_stats.txt" ]]; then
            CELLS=$(cut -d',' -f1 "${SAMPLE_OUTPUT_DIR}/conversion_stats.txt")
            GENES=$(cut -d',' -f2 "${SAMPLE_OUTPUT_DIR}/conversion_stats.txt")
            echo "  ${sample}: ${status} | Cells: ${CELLS}, Genes: ${GENES}" >> "${README_FILE}"
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
    error_exit "Conversion failed"
fi

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

echo ""
log_msg "============================================================================"
log_msg "Conversion Complete - ${SAMPLE}"
log_msg "============================================================================"
log_msg "Output: ${OUTPUT_RDS}"
echo ""
ls -lh "${SAMPLE_OUTPUT_DIR}/"
echo ""
log_msg "============================================================================"
