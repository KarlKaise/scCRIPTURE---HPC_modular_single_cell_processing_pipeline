#!/bin/bash
#SBATCH --job-name=CHOIR
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --qos=1day
#SBATCH --output=logs/CHOIR/CHOIR_integration_%j.out
#SBATCH --error=logs/CHOIR/CHOIR_integration_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ==============================================================================
# STEP 9: CHOIR Clustering with Harmony Integration Pipeline
# ==============================================================================
#
# DATASET: Vandebroucke Fibroblast Paper
#
# This script runs CHOIR clustering on:
#   1. Individual samples (4 samples)
#   2. LV only: M_22w_LV + F_7w_LV
#   3. 4V only: M_22w_4V + F_7w_4V
#   4. All combined (4 samples)
#
# Sample Metadata:
#   - M_22w_LV  (Male, 22 weeks, Lateral Ventricle)
#   - M_22w_4V  (Male, 22 weeks, 4th Ventricle)
#   - F_7w_LV   (Female, 7 weeks, Lateral Ventricle)
#   - F_7w_4V   (Female, 7 weeks, 4th Ventricle)
#
# Input:  8_DecontX_correction/<SAMPLE>/<SAMPLE>_decontX_corrected.rds
# Output: 9_CHOIR_integration/
#
# ==============================================================================

set -euo pipefail

# ==============================================================================
# Directory Configuration
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

# UPDATED: Input now comes from DecontX output (Step 8)
IN_DIR="${PREPROCESS_DIR}/8_DecontX_correction"
OUT_DIR="${PREPROCESS_DIR}/9_CHOIR_integration"
SCRIPT_DIR="${BASE_DIR}/Scripts/R_scripts"
LOG_DIR="${BASE_DIR}/logs/CHOIR"
README_FILE="${OUT_DIR}/README.txt"

N_CORES=${SLURM_CPUS_PER_TASK:-4}

# ==============================================================================
# Create directories
# ==============================================================================
mkdir -p "${OUT_DIR}" "${SCRIPT_DIR}" "${LOG_DIR}"
mkdir -p "${OUT_DIR}/Individual"
mkdir -p "${OUT_DIR}/LV_only"
mkdir -p "${OUT_DIR}/4V_only"
mkdir -p "${OUT_DIR}/All_combined"

# Read samples from samplesheet
if [[ -f "${SAMPLESHEET}" ]]; then
    SAMPLES=($(tail -n +2 "${SAMPLESHEET}" | cut -d',' -f1 | grep -v '^$'))
else
    SAMPLES=("M_22w_LV" "M_22w_4V" "F_7w_LV" "F_7w_4V")
fi

for SAMPLE in "${SAMPLES[@]}"; do
    mkdir -p "${OUT_DIR}/Individual/${SAMPLE}"
done

# ==============================================================================
# Initialize README
# ==============================================================================
cat > "${README_FILE}" << EOF
================================================================================
STEP 9: CHOIR Clustering with Harmony Integration
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')
Pipeline: Vandebroucke Fibroblast scRNA-seq Analysis

DESCRIPTION:
  Runs CHOIR clustering with Harmony batch integration on multiple configurations.
  Uses decontaminated counts from combined scCDC + DecontX correction.

INPUT:
  DecontX-corrected RDS: 8_DecontX_correction/<SAMPLE>/<SAMPLE>_decontX_corrected.rds
  Primary counts layer: decontX_corrected

SAMPLE METADATA:
  - M_22w_LV  (Male, 22 weeks, Lateral Ventricle)
  - M_22w_4V  (Male, 22 weeks, 4th Ventricle)
  - F_7w_LV   (Female, 7 weeks, Lateral Ventricle)
  - F_7w_4V   (Female, 7 weeks, 4th Ventricle)

INTEGRATION STRATEGIES:
  1. Individual samples (4)
  2. LV_only: M_22w_LV, F_7w_LV
  3. 4V_only: M_22w_4V, F_7w_4V
  4. All_combined: All 4 samples

PROCESSING LOG:
EOF

# ==============================================================================
# Environment Setup
# ==============================================================================
echo "============================================================"
echo "CHOIR Clustering with Harmony Integration Pipeline"
echo "Vandebroucke Fibroblast Dataset"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

source /scicore/home/doetsch/kaiser0001/miniforge3/etc/profile.d/conda.sh
conda activate /scicore/home/doetsch/kaiser0001/miniforge3/envs/R_4_5

export RETICULATE_PYTHON="/scicore/home/doetsch/kaiser0001/miniforge3/envs/R_4_5/bin/python"

echo "[env] R: $(which R)"
echo "[env] Python: ${RETICULATE_PYTHON}"
echo "[paths] BASE_DIR=${BASE_DIR}"
echo "[paths] IN_DIR=${IN_DIR}"
echo "[paths] OUT_DIR=${OUT_DIR}"
echo ""

echo "[input] Samples from samplesheet: ${SAMPLES[*]}"
for SAMPLE in "${SAMPLES[@]}"; do
    INPUT_FILE="${IN_DIR}/${SAMPLE}/${SAMPLE}_decontX_corrected.rds"
    if [[ -f "${INPUT_FILE}" ]]; then
        echo "  OK: ${SAMPLE}"
    else
        echo "  MISSING: ${SAMPLE}"
    fi
done
echo ""

# ==============================================================================
# Write R script
# ==============================================================================
R_SCRIPT="${SCRIPT_DIR}/run_CHOIR_integration.R"

cat > "${R_SCRIPT}" << 'RSCRIPT_EOF'
#!/usr/bin/env Rscript

# Setup Python
conda_python <- "/scicore/home/doetsch/kaiser0001/miniforge3/envs/R_4_5/bin/python"
if (file.exists(conda_python)) {
    Sys.setenv(RETICULATE_PYTHON = conda_python)
}

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(CHOIR)
    library(harmony)
    library(scCustomize)
    library(reticulate)
})

tryCatch({
    use_condaenv("/scicore/home/doetsch/kaiser0001/miniforge3/envs/R_4_5", required = TRUE)
}, error = function(e) {
    use_python(conda_python, required = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript run_CHOIR_integration.R <input_dir> <output_dir> <readme_file> <samplesheet> [n_cores]")

input_dir   <- normalizePath(args[[1]], mustWork = TRUE)
output_dir  <- normalizePath(args[[2]], mustWork = FALSE)
readme_file <- args[[3]]
samplesheet <- args[[4]]
n_cores     <- ifelse(length(args) >= 5, as.integer(args[[5]]), 4)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Read sample metadata from samplesheet
# ==============================================================================
ss <- read.csv(samplesheet, stringsAsFactors = FALSE)

sample_metadata <- list()
for (i in 1:nrow(ss)) {
    sample_metadata[[ss$sample_id[i]]] <- list(
        sex = ss$sex[i],
        age = ss$age[i],
        batch = ss$batch[i],
        ventricle = ss$ventricle[i]
    )
}

all_samples <- names(sample_metadata)

# Define integration groups based on ventricle
LV_samples <- ss$sample_id[ss$ventricle == "LV"]
V4_samples <- ss$sample_id[ss$ventricle == "4V"]

cat("All samples:", paste(all_samples, collapse = ", "), "\n")
cat("LV samples:", paste(LV_samples, collapse = ", "), "\n")
cat("4V samples:", paste(V4_samples, collapse = ", "), "\n")

# ==============================================================================
# Logging setup
# ==============================================================================
log_file <- file.path(output_dir, paste0("CHOIR_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

log_message <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    full_msg <- paste0("[", timestamp, "] ", msg)
    message(full_msg)
    cat(full_msg, "\n", file = log_file, append = TRUE)
}

readme_append <- function(msg) {
    cat(msg, "\n", file = readme_file, append = TRUE)
}

log_message("=== CHOIR Clustering Pipeline - Vandebroucke Dataset ===")
log_message(paste0("Input directory: ", input_dir))
log_message(paste0("Samples: ", paste(all_samples, collapse = ", ")))

# Check Python/anndata
python_available <- tryCatch({
    py_module_available("anndata")
}, error = function(e) FALSE)

if (python_available) {
    anndata <- import("anndata")
    log_message(paste0("anndata version: ", anndata$`__version__`))
}

# ==============================================================================
# Helper Functions
# ==============================================================================
add_barcode_suffix <- function(seurat_obj, suffix) {
    RenameCells(seurat_obj, new.names = paste0(colnames(seurat_obj), "_", suffix))
}

ensure_metadata <- function(seurat_obj, sample_name, metadata) {
    seurat_obj$Sample_ID <- sample_name
    seurat_obj$Sex <- metadata$sex
    seurat_obj$Age <- metadata$age
    seurat_obj$Batch <- metadata$batch
    seurat_obj$Ventricle <- metadata$ventricle
    seurat_obj
}

prepare_for_choir <- function(seurat_obj, use_decontx = TRUE) {
    # UPDATED: Use decontX_corrected layer as primary counts
    DefaultAssay(seurat_obj) <- "RNA"
    available_layers <- Layers(seurat_obj, assay = "RNA")
    
    log_message(paste0("  Available layers: ", paste(available_layers, collapse = ", ")))
    
    # Priority: decontX_corrected > scCDC_corrected > counts
    if (use_decontx && "decontX_corrected" %in% available_layers) {
        log_message("  Using decontX_corrected counts (scCDC + DecontX corrected)")
        seurat_obj[["RNA"]]$counts <- seurat_obj[["RNA"]]$decontX_corrected
    } else if ("scCDC_corrected" %in% available_layers) {
        log_message("  Using scCDC_corrected counts (DecontX layer not found)")
        seurat_obj[["RNA"]]$counts <- seurat_obj[["RNA"]]$scCDC_corrected
    } else {
        log_message("  Using original counts (no corrected layers found)")
    }
    
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
    seurat_obj
}

merge_seurat_v5 <- function(seurat_list, project_name) {
    log_message("Merging Seurat v5 objects...")
    if (length(seurat_list) == 1) return(seurat_list[[1]])
    merged <- merge(seurat_list[[1]], seurat_list[-1], project = project_name)
    log_message(paste0("After merge: ", ncol(merged), " cells"))
    merged <- JoinLayers(merged)
    merged
}

run_choir_single <- function(seurat_obj, sample_name, n_cores) {
    log_message(paste0("Running CHOIR on: ", sample_name))
    tryCatch({
        seurat_obj <- CHOIR(seurat_obj, n_cores = n_cores, alpha = 0.05, verbose = TRUE)
        seurat_obj <- runCHOIRumap(seurat_obj, reduction = "P0_reduction")
        log_message(paste0("CHOIR complete for ", sample_name))
        seurat_obj
    }, error = function(e) {
        log_message(paste0("ERROR in CHOIR: ", conditionMessage(e)))
        seurat_obj
    })
}

run_choir_integrated <- function(seurat_list, integration_name, n_cores, batch_var = "Sample_ID") {
    log_message(paste0("Running integrated CHOIR: ", integration_name))
    merged <- merge_seurat_v5(seurat_list, project_name = integration_name)
    
    log_message(paste0("Running Harmony on: ", batch_var))
    merged <- NormalizeData(merged, verbose = FALSE)
    merged <- FindVariableFeatures(merged, nfeatures = 2000, verbose = FALSE)
    merged <- ScaleData(merged, verbose = FALSE)
    merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
    merged <- RunHarmony(merged, group.by.vars = batch_var, verbose = FALSE)
    
    tryCatch({
        merged <- CHOIR(merged, n_cores = n_cores, alpha = 0.05, verbose = TRUE, reduction = "harmony")
        merged <- runCHOIRumap(merged, reduction = "P0_reduction")
        log_message(paste0("Integrated CHOIR complete: ", integration_name))
        merged
    }, error = function(e) {
        log_message(paste0("ERROR: ", conditionMessage(e)))
        merged
    })
}

convert_to_h5ad <- function(seurat_obj, out_dir, file_base) {
    if (!python_available) return(NULL)
    h5ad_path <- file.path(out_dir, paste0(file_base, ".h5ad"))
    tryCatch({
        adata <- as.anndata(x = seurat_obj, file_path = out_dir, file_name = paste0(file_base, ".h5ad"))
        log_message(paste0("Saved: ", file_base, ".h5ad"))
    }, error = function(e) {
        log_message(paste0("h5ad ERROR: ", conditionMessage(e)))
    })
}

# ==============================================================================
# Find input files - UPDATED for DecontX output
# ==============================================================================
log_message("Finding input files...")

sample_files <- list()
for (sample_name in all_samples) {
    # UPDATED: Look for decontX_corrected files
    input_file <- file.path(input_dir, sample_name, paste0(sample_name, "_decontX_corrected.rds"))
    if (file.exists(input_file)) {
        sample_files[[sample_name]] <- input_file
        log_message(paste0("  Found: ", sample_name))
    } else {
        # Fallback to scCDC output if DecontX not available
        fallback_file <- file.path(dirname(dirname(input_dir)), "7_scCDC_correction", 
                                    sample_name, paste0(sample_name, "_scCDC_corrected.rds"))
        if (file.exists(fallback_file)) {
            sample_files[[sample_name]] <- fallback_file
            log_message(paste0("  Found (fallback scCDC): ", sample_name))
        } else {
            log_message(paste0("  MISSING: ", sample_name))
        }
    }
}

if (length(sample_files) == 0) stop("No input files found!")

# ==============================================================================
# PART 1: Individual samples
# ==============================================================================
log_message("")
log_message(strrep("=", 70))
log_message("PART 1: Individual Sample CHOIR Clustering")
log_message(strrep("=", 70))

readme_append("")
readme_append("INDIVIDUAL SAMPLES:")

for (sample_name in names(sample_files)) {
    log_message(paste0("Processing: ", sample_name))
    
    seu <- readRDS(sample_files[[sample_name]])
    log_message(paste0("  Loaded: ", ncol(seu), " cells"))
    
    seu <- ensure_metadata(seu, sample_name, sample_metadata[[sample_name]])
    seu <- prepare_for_choir(seu, use_decontx = TRUE)
    seu <- add_barcode_suffix(seu, sample_name)
    seu <- run_choir_single(seu, sample_name, n_cores)
    
    out_dir <- file.path(output_dir, "Individual", sample_name)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    saveRDS(seu, file.path(out_dir, paste0(sample_name, "_CHOIR.rds")))
    readme_append(paste0("  ", sample_name, ": ", ncol(seu), " cells"))
    convert_to_h5ad(seu, out_dir, paste0(sample_name, "_CHOIR"))
    
    rm(seu); gc()
}

# ==============================================================================
# PART 2: LV only integration
# ==============================================================================
log_message("")
log_message(strrep("=", 70))
log_message("PART 2: LV Only Integration")
log_message(strrep("=", 70))

readme_append("")
readme_append("LV_ONLY INTEGRATION:")

if (length(LV_samples) >= 2) {
    lv_list <- list()
    for (sample_name in LV_samples) {
        if (sample_name %in% names(sample_files)) {
            seu <- readRDS(sample_files[[sample_name]])
            seu <- ensure_metadata(seu, sample_name, sample_metadata[[sample_name]])
            seu <- prepare_for_choir(seu, use_decontx = TRUE)
            seu <- add_barcode_suffix(seu, paste0(sample_name, "_LV"))
            lv_list[[sample_name]] <- seu
        }
    }
    
    if (length(lv_list) > 0) {
        integrated <- run_choir_integrated(lv_list, "LV_only", n_cores, batch_var = "Sample_ID")
        rm(lv_list); gc()
        
        out_dir <- file.path(output_dir, "LV_only")
        saveRDS(integrated, file.path(out_dir, "LV_only_CHOIR_integrated.rds"))
        readme_append(paste0("  LV_only: ", ncol(integrated), " cells"))
        convert_to_h5ad(integrated, out_dir, "LV_only_CHOIR_integrated")
        rm(integrated); gc()
    }
} else {
    log_message("Skipping LV integration - fewer than 2 samples")
}

# ==============================================================================
# PART 3: 4V only integration
# ==============================================================================
log_message("")
log_message(strrep("=", 70))
log_message("PART 3: 4V Only Integration")
log_message(strrep("=", 70))

readme_append("")
readme_append("4V_ONLY INTEGRATION:")

if (length(V4_samples) >= 2) {
    v4_list <- list()
    for (sample_name in V4_samples) {
        if (sample_name %in% names(sample_files)) {
            seu <- readRDS(sample_files[[sample_name]])
            seu <- ensure_metadata(seu, sample_name, sample_metadata[[sample_name]])
            seu <- prepare_for_choir(seu, use_decontx = TRUE)
            seu <- add_barcode_suffix(seu, paste0(sample_name, "_4V"))
            v4_list[[sample_name]] <- seu
        }
    }
    
    if (length(v4_list) > 0) {
        integrated <- run_choir_integrated(v4_list, "4V_only", n_cores, batch_var = "Sample_ID")
        rm(v4_list); gc()
        
        out_dir <- file.path(output_dir, "4V_only")
        saveRDS(integrated, file.path(out_dir, "4V_only_CHOIR_integrated.rds"))
        readme_append(paste0("  4V_only: ", ncol(integrated), " cells"))
        convert_to_h5ad(integrated, out_dir, "4V_only_CHOIR_integrated")
        rm(integrated); gc()
    }
} else {
    log_message("Skipping 4V integration - fewer than 2 samples")
}

# ==============================================================================
# PART 4: All combined
# ==============================================================================
log_message("")
log_message(strrep("=", 70))
log_message("PART 4: All Samples Combined")
log_message(strrep("=", 70))

readme_append("")
readme_append("ALL_COMBINED INTEGRATION:")

all_list <- list()
for (sample_name in names(sample_files)) {
    seu <- readRDS(sample_files[[sample_name]])
    seu <- ensure_metadata(seu, sample_name, sample_metadata[[sample_name]])
    seu <- prepare_for_choir(seu, use_decontx = TRUE)
    seu <- add_barcode_suffix(seu, paste0(sample_name, "_all"))
    all_list[[sample_name]] <- seu
}

if (length(all_list) > 0) {
    # For all combined, integrate by Sex (since all same batch)
    integrated <- run_choir_integrated(all_list, "All_combined", n_cores, batch_var = "Sex")
    rm(all_list); gc()
    
    out_dir <- file.path(output_dir, "All_combined")
    saveRDS(integrated, file.path(out_dir, "All_combined_CHOIR_integrated.rds"))
    readme_append(paste0("  All_combined: ", ncol(integrated), " cells"))
    convert_to_h5ad(integrated, out_dir, "All_combined_CHOIR_integrated")
    rm(integrated); gc()
}

# ==============================================================================
# Done
# ==============================================================================
log_message("")
log_message("=== PIPELINE COMPLETE ===")

readme_append("")
readme_append("NEXT STEP:")
readme_append("  Step 10: QC Visualization")
readme_append("================================================================================")

log_message("Done!")
RSCRIPT_EOF

chmod +x "${R_SCRIPT}"
echo "[script] R script written to: ${R_SCRIPT}"

# ==============================================================================
# Run the R script
# ==============================================================================
echo ""
echo "============================================================"
echo "Running CHOIR Clustering Pipeline"
echo "============================================================"
echo "Start time: $(date)"

Rscript "${R_SCRIPT}" "${IN_DIR}" "${OUT_DIR}" "${README_FILE}" "${SAMPLESHEET}" "${N_CORES}"

echo ""
echo "============================================================"
echo "Pipeline Complete"
echo "============================================================"
echo "End time: $(date)"

echo "[output] Output files:"
find "${OUT_DIR}" -type f \( -name "*.rds" -o -name "*.h5ad" \) -exec ls -lh {} \;

echo ""
echo "Job finished successfully"
