#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 01: LOAD INPUT DATA (MULTI-SAMPLE PIPELINE)
# ==============================================================================
#
# This module loads scCDC-corrected Seurat objects from preprocessing pipeline
# and merges them into a single object for downstream analysis.
#
# INPUT: scCDC-corrected RDS files from Step 7 of preprocessing
#        Located at: {base_dir}/7_scCDC_correction/{sample}/{sample}_scCDC_corrected.rds
#
# OUTPUT: Merged Seurat object with all selected samples
#         Saved to: {out_root}/objects/01_loaded_data.RData
#
# FEATURES:
# - Loads multiple samples based on params$samples_to_analyze
# - Uses scCDC_corrected layer as counts for analysis
# - Adds comprehensive metadata from sample_metadata
# - Handles batch variable assignment for downstream integration
# - Validates all input files before processing
# - Maps preprocessing QC column names to standard names
# - Optionally filters hemoglobin-high cells
# - Validates and fixes batch assignments from authoritative samplesheet
#
# UPDATES:
# - 2026-01-03: Added pct_counts_MT → percent.mt mapping
# - 2026-01-03: Added preprocessing metadata verification
# - 2026-01-03: Added hemoglobin filtering support
# - 2026-01-03: Improved Seurat v5 compatibility
# - 2026-01-12: Added batch validation and fix from authoritative samplesheet
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 01: LOAD INPUT DATA\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ==============================================================================
# Load environment from previous module
# ==============================================================================
pipeline_dir <- Sys.getenv("PIPELINE_DIR", unset = "")
if (pipeline_dir == "") {
  if (file.exists("config/params.R")) {
    pipeline_dir <- normalizePath(".")
  } else {
    stop("PIPELINE_DIR not set. Run Module 00 first or set environment variable.")
  }
}

source(file.path(pipeline_dir, "config", "params.R"))

# Load utility functions if available
utils_file <- file.path(pipeline_dir, "utils", "functions.R")
if (file.exists(utils_file)) {
  source(utils_file)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(patchwork)
})

# Load saved environment
out_base <- params$out_root
env_file <- file.path(out_base, "objects", "pipeline_environment.RData")

if (file.exists(env_file)) {
  load(env_file)
  cat("Loaded pipeline environment from:", env_file, "\n")
} else {
  stop("Environment file not found. Run 00_setup_environment.R first.")
}

cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ==============================================================================
# EXPECTED PREPROCESSING METADATA COLUMNS
# ==============================================================================
# These columns should exist from Steps 2-7 of preprocessing
EXPECTED_PREPROCESSING_COLS <- c(
  # From DropletQC (Step 2)
  "dropletqc_nuclear_fraction", "dropletqc_status", "dropletqc_empty_droplet", "dropletqc_damaged_cell",
  # From QClus (Step 3)
  "qclus", "filter_status",
  # From VAEDA (Step 4)
  "vaeda_prediction", "vaeda_calls", "vaeda_scores",
  # From Doublet Consensus (Step 6)
  "doublet_vaeda", "doublet_scDblFinder", "doublet_DoubletFinder",
  "doublet_votes", "doublet_consensus",
  "scDblFinder_score", "scDblFinder_class", "DoubletFinder_class",
  # QC metrics (various sources)
  "nCount_RNA", "nFeature_RNA", "pct_counts_MT"
)

# ==============================================================================
# VALIDATE INPUT FILES
# ==============================================================================
cat("--- Validating Input Files ---\n\n")

samples_to_load <- params$samples_to_analyze
input_paths <- params$input_paths

cat("Samples to load:", length(samples_to_load), "\n")
for (s in samples_to_load) {
  cat("  -", s, "\n")
}
cat("\n")

# Check all files exist
missing_files <- c()
for (sample in samples_to_load) {
  path <- input_paths[[sample]]
  if (!file.exists(path)) {
    missing_files <- c(missing_files, sample)
    cat("  [MISSING]", sample, ":", path, "\n")
  } else {
    size_mb <- round(file.info(path)$size / (1024^2), 2)
    cat("  [OK]", sample, ":", size_mb, "MB\n")
  }
}

if (length(missing_files) > 0) {
  stop("Missing input files for samples: ", paste(missing_files, collapse = ", "),
       "\nRun preprocessing pipeline first or update samples_to_exclude in params.R")
}

cat("\nAll input files validated.\n")

# ==============================================================================
# LOAD AUTHORITATIVE SAMPLESHEET FOR BATCH VALIDATION
# ==============================================================================
# The RDS files may contain outdated/incorrect batch metadata from preprocessing.
# We validate and override batch assignments from the authoritative samplesheet.
# ==============================================================================

authoritative_samplesheet <- NULL
samplesheet_path <- Sys.getenv("SAMPLESHEET", unset = "")

if (samplesheet_path != "" && file.exists(samplesheet_path)) {
  cat("\n--- Loading Authoritative Samplesheet for Batch Validation ---\n")
  cat("Path:", samplesheet_path, "\n")
  authoritative_samplesheet <- read.csv(samplesheet_path, stringsAsFactors = FALSE)
  
  # Standardize column name: sample_id -> sample_name
  if ("sample_id" %in% colnames(authoritative_samplesheet) && !"sample_name" %in% colnames(authoritative_samplesheet)) {
    authoritative_samplesheet$sample_name <- authoritative_samplesheet$sample_id
  }
  
  cat("Expected batch assignments from samplesheet:\n")
  for (s in samples_to_load) {
    if ("sample_name" %in% colnames(authoritative_samplesheet)) {
      row <- authoritative_samplesheet[authoritative_samplesheet$sample_name == s, ]
    } else if ("sample_id" %in% colnames(authoritative_samplesheet)) {
      row <- authoritative_samplesheet[authoritative_samplesheet$sample_id == s, ]
    } else {
      row <- data.frame()
    }
    if (nrow(row) > 0 && "batch" %in% colnames(row)) {
      cat("  ", s, "->", row$batch[1], "\n")
    } else {
      cat("  ", s, "-> [NOT FOUND IN SAMPLESHEET]\n")
    }
  }
  cat("\n")
} else {
  cat("\nNOTE: SAMPLESHEET env var not set or file not found.\n")
  cat("      Using params$sample_metadata for batch assignments.\n")
  cat("      Batch validation against authoritative source will be skipped.\n\n")
}

# ==============================================================================
# HELPER FUNCTION: Validate and fix batch assignments from samplesheet
# ==============================================================================
validate_and_fix_batch_assignments <- function(obj, sample_name, params_metadata, authoritative_sheet = NULL) {
  # Get batch from the loaded object
obj_batch <- unique(obj$batch)
  if (length(obj_batch) > 1) {
    warning("Multiple batch values in object for ", sample_name, ": ", paste(obj_batch, collapse = ", "))
    obj_batch <- obj_batch[1]
  }
  
  # Get expected batch from samplesheet (priority: authoritative > params)
  expected_batch <- NULL
  source_used <- NULL
  
  # Try authoritative samplesheet first
  if (!is.null(authoritative_sheet)) {
    if ("sample_name" %in% colnames(authoritative_sheet)) {
      row <- authoritative_sheet[authoritative_sheet$sample_name == sample_name, ]
    } else if ("sample_id" %in% colnames(authoritative_sheet)) {
      row <- authoritative_sheet[authoritative_sheet$sample_id == sample_name, ]
    } else {
      row <- data.frame()
    }
    if (nrow(row) > 0 && "batch" %in% colnames(row)) {
      expected_batch <- row$batch[1]
      source_used <- "authoritative_samplesheet"
    }
  }
  
  # Fallback to params$sample_metadata
  if (is.null(expected_batch) && !is.null(params_metadata)) {
    if ("sample_name" %in% colnames(params_metadata)) {
      row <- params_metadata[params_metadata$sample_name == sample_name, ]
    } else if ("sample_id" %in% colnames(params_metadata)) {
      row <- params_metadata[params_metadata$sample_id == sample_name, ]
    } else {
      row <- data.frame()
    }
    if (nrow(row) > 0 && "batch" %in% colnames(row)) {
      expected_batch <- row$batch[1]
      source_used <- "params_sample_metadata"
    }
  }
  
  # Compare and report
  if (is.null(expected_batch)) {
    cat("    [WARNING] No batch found in samplesheet for", sample_name, "\n")
    return(list(obj = obj, mismatch = FALSE, fixed = FALSE))
  }
  
  if (is.na(obj_batch) || obj_batch != expected_batch) {
    cat("    [BATCH MISMATCH] RDS file:", obj_batch, "-> Samplesheet:", expected_batch, "(", source_used, ")\n")
    cat("    [FIXING] Overriding batch to:", expected_batch, "\n")
    obj$batch <- expected_batch
    obj$sequencing_batch <- expected_batch
    return(list(obj = obj, mismatch = TRUE, fixed = TRUE, old_batch = obj_batch, new_batch = expected_batch))
  } else {
    cat("    [OK] Batch matches samplesheet:", expected_batch, "\n")
    return(list(obj = obj, mismatch = FALSE, fixed = FALSE))
  }
}

# ==============================================================================
# HELPER FUNCTION: Ensure QC columns exist with correct names
# ==============================================================================
ensure_qc_columns_multi <- function(obj, sample_name) {
  cat("    Mapping QC columns...\n")

  # --------------------------------------------------------------------------
  # Map MT percentage columns to percent.mt
  # --------------------------------------------------------------------------
  # Priority order: percent.mt > pct_counts_MT > percent_mt > pct_counts_mt
  mt_cols <- c("percent.mt", "pct_counts_MT", "percent_mt", "pct_counts_mt")
  mt_col_found <- intersect(mt_cols, colnames(obj@meta.data))

  if ("percent.mt" %in% colnames(obj@meta.data)) {
    # Already has percent.mt with correct name
    cat("      [OK] percent.mt already exists\n")
  } else if (length(mt_col_found) > 0) {
    # Map from alternative column name
    source_col <- mt_col_found[1]
    obj$percent.mt <- obj@meta.data[[source_col]]
    cat("      [MAPPED]", source_col, "→ percent.mt\n")
  } else {
    # Calculate MT percentage from counts
    cat("      [CALC] Computing percent.mt from counts...\n")
    mt_genes <- grep("^MT-|^mt-", rownames(obj), value = TRUE)
    if (length(mt_genes) > 0) {
      counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
      mt_counts <- Matrix::colSums(counts_mat[mt_genes, , drop = FALSE])
      total_counts <- Matrix::colSums(counts_mat)
      obj$percent.mt <- (mt_counts / total_counts) * 100
      cat("      Found", length(mt_genes), "MT genes\n")
    } else {
      obj$percent.mt <- 0
      cat("      WARNING: No MT genes found (pattern ^MT-|^mt-)\n")
    }
  }

  # --------------------------------------------------------------------------
  # Ensure nFeature_RNA exists
  # --------------------------------------------------------------------------
  if ("nFeature_RNA" %in% colnames(obj@meta.data)) {
    cat("      [OK] nFeature_RNA already exists\n")
  } else if ("n_genes_by_counts" %in% colnames(obj@meta.data)) {
    obj$nFeature_RNA <- obj$n_genes_by_counts
    cat("      [MAPPED] n_genes_by_counts → nFeature_RNA\n")
  } else {
    counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
    obj$nFeature_RNA <- Matrix::colSums(counts_mat > 0)
    cat("      [CALC] Computed nFeature_RNA from counts\n")
  }

  # --------------------------------------------------------------------------
  # Ensure nCount_RNA exists
  # --------------------------------------------------------------------------
  if ("nCount_RNA" %in% colnames(obj@meta.data)) {
    cat("      [OK] nCount_RNA already exists\n")
  } else if ("total_counts" %in% colnames(obj@meta.data)) {
    obj$nCount_RNA <- obj$total_counts
    cat("      [MAPPED] total_counts → nCount_RNA\n")
  } else {
    counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
    obj$nCount_RNA <- Matrix::colSums(counts_mat)
    cat("      [CALC] Computed nCount_RNA from counts\n")
  }

  return(obj)
}

# ==============================================================================
# HELPER FUNCTION: Calculate hemoglobin percentage
# ==============================================================================
calculate_percent_hb <- function(obj, hb_pattern = "^HB[AB]-") {
  # Get counts matrix
  counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")

  # Find hemoglobin genes (excluding HBEGF which is NOT a hemoglobin gene)
  all_hb_matches <- grep(hb_pattern, rownames(obj), value = TRUE, ignore.case = TRUE)

  # Explicitly exclude HBEGF (Heparin Binding EGF Like Growth Factor)
  hb_genes <- all_hb_matches[!grepl("HBEGF|Hbegf", all_hb_matches)]

  if (length(hb_genes) == 0) {
    cat("    No hemoglobin genes found with pattern:", hb_pattern, "\n")
    return(rep(0, ncol(obj)))
  }

  cat("    Found", length(hb_genes), "hemoglobin genes:", paste(hb_genes, collapse = ", "), "\n")

  # Calculate percentage
  hb_counts <- Matrix::colSums(counts_mat[hb_genes, , drop = FALSE])
  total_counts <- Matrix::colSums(counts_mat)
  percent_hb <- (hb_counts / total_counts) * 100
  percent_hb[is.na(percent_hb)] <- 0

  return(percent_hb)
}

# ==============================================================================
# HELPER FUNCTION: Load and prepare single sample
# ==============================================================================
load_and_prepare_sample <- function(sample_name, file_path, sample_metadata, params, authoritative_sheet = NULL) {
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("Loading sample:", sample_name, "\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")

  # Load the object
  obj <- readRDS(file_path)

  # Get sample metadata from config
  meta <- sample_metadata[sample_metadata$sample_name == sample_name, ]
  
  # If not found by sample_name, try sample_id
  if (nrow(meta) == 0 && "sample_id" %in% colnames(sample_metadata)) {
    meta <- sample_metadata[sample_metadata$sample_id == sample_name, ]
  }
  
  if (nrow(meta) == 0) {
    stop("No metadata found for sample: ", sample_name, " in params$sample_metadata")
  }

  cat("  Original dimensions:", nrow(obj), "genes x", ncol(obj), "cells\n")
  cat("  Assays:", paste(names(obj@assays), collapse = ", "), "\n")

  # Check for RNA assay
  if (!"RNA" %in% names(obj@assays)) {
    stop("No RNA assay found in ", sample_name)
  }

  # Get available layers
  available_layers <- Layers(obj[["RNA"]])
  cat("  RNA layers:", paste(available_layers, collapse = ", "), "\n")

  # Report reductions (will be discarded but good to document)
  if (length(names(obj@reductions)) > 0) {
    cat("  Reductions (will discard):", paste(names(obj@reductions), collapse = ", "), "\n")
  }

  # ===========================================================================
  # Preprocessing metadata verification
  # ===========================================================================
  present_cols <- intersect(EXPECTED_PREPROCESSING_COLS, colnames(obj@meta.data))
  missing_cols <- setdiff(EXPECTED_PREPROCESSING_COLS, colnames(obj@meta.data))

  cat("  Preprocessing metadata: ", length(present_cols), "/", length(EXPECTED_PREPROCESSING_COLS), " expected columns present\n", sep = "")

  # Report critical missing columns
  critical_cols <- c("doublet_consensus", "pct_counts_MT", "nCount_RNA", "nFeature_RNA")
  critical_missing <- intersect(missing_cols, critical_cols)
  if (length(critical_missing) > 0) {
    cat("    NOTE: Missing columns (will be computed/mapped):", paste(critical_missing, collapse = ", "), "\n")
  }

  # ===========================================================================
  # Extract scCDC-corrected counts
  # ===========================================================================
  if ("scCDC_corrected" %in% available_layers) {
    cat("  >> Using scCDC_corrected layer as counts <<\n")

    # Extract scCDC corrected counts
    sccdc_counts <- LayerData(obj, assay = "RNA", layer = "scCDC_corrected")

    # Create new Seurat object with scCDC counts as the counts layer
    new_obj <- CreateSeuratObject(
      counts = sccdc_counts,
      meta.data = obj@meta.data,
      project = sample_name
    )

    cat("  Created new object with scCDC_corrected as counts layer\n")
    obj <- new_obj
    rm(new_obj, sccdc_counts)
    gc(verbose = FALSE)

  } else if ("counts" %in% available_layers) {
    cat("  >> WARNING: scCDC_corrected layer not found <<\n")
    cat("     Using original counts - Step 7 may not have completed.\n")
  } else {
    stop("No counts layer found in ", sample_name)
  }

  # ===========================================================================
  # Add/update metadata from sample_sheet
  # ===========================================================================
  obj$sample_name <- sample_name
  obj$sex <- meta$sex
  obj$Sex <- meta$sex  # Capitalized version for compatibility
  obj$batch <- meta$batch
  obj$sequencing_batch <- meta$batch

  # ===========================================================================
  # VALIDATE AND FIX BATCH ASSIGNMENT FROM AUTHORITATIVE SAMPLESHEET
  # ===========================================================================
  batch_check <- validate_and_fix_batch_assignments(
    obj = obj,
    sample_name = sample_name,
    params_metadata = sample_metadata,
    authoritative_sheet = authoritative_sheet
  )
  obj <- batch_check$obj

  # Handle optional metadata columns safely
  # estrous_phase
  if ("estrous_phase" %in% colnames(meta) && !is.null(meta$estrous_phase) && !is.na(meta$estrous_phase)) {
    obj$estrous_phase <- meta$estrous_phase
  } else {
    obj$estrous_phase <- NA
  }

  # sequencing_round
  if ("sequencing_round" %in% colnames(meta) && !is.null(meta$sequencing_round) && !is.na(meta$sequencing_round)) {
    obj$sequencing_round <- meta$sequencing_round
  } else {
    obj$sequencing_round <- NA
  }

  # ventricle (specific to Vandebroucke dataset)
  if ("ventricle" %in% colnames(meta) && !is.null(meta$ventricle) && !is.na(meta$ventricle)) {
    obj$ventricle <- meta$ventricle
  }

  # age
  if ("age" %in% colnames(meta) && !is.null(meta$age) && !is.na(meta$age)) {
    obj$age <- meta$age
  }

  # condition
  if ("condition" %in% colnames(meta) && !is.null(meta$condition) && !is.na(meta$condition)) {
    obj$condition <- meta$condition
  }

  # Set batch variable based on params$batch_variable
  if (!is.null(params$batch_variable) && params$batch_variable != "none") {
    if (params$batch_variable %in% colnames(obj@meta.data)) {
      obj$integration_batch <- obj@meta.data[[params$batch_variable]]
    } else {
      cat("  WARNING: batch_variable '", params$batch_variable,
          "' not found in metadata. Using 'batch' column.\n", sep = "")
      obj$integration_batch <- obj$batch
    }
  }

  # ===========================================================================
  # Ensure QC columns exist with correct names
  # ===========================================================================
  obj <- ensure_qc_columns_multi(obj, sample_name)

  # ===========================================================================
  # Calculate hemoglobin percentage if filtering is enabled
  # ===========================================================================
  if (isTRUE(params$filter_hemoglobin)) {
    hb_pattern <- if (!is.null(params$hemoglobin_pattern)) params$hemoglobin_pattern else "^HB[AB]-"
    cat("    Calculating percent.hb (pattern:", hb_pattern, ")...\n")
    obj$percent.hb <- calculate_percent_hb(obj, hb_pattern)
    cat("      Median:", round(median(obj$percent.hb), 4), "%, Max:", round(max(obj$percent.hb), 2), "%\n")
  }

  # ===========================================================================
  # Summary
  # ===========================================================================
  cat("  Final dimensions:", nrow(obj), "genes x", ncol(obj), "cells\n")
  cat("  Sex:", unique(obj$sex), "\n")
  cat("  Batch:", unique(obj$batch), "\n")
  cat("  Median %MT:", round(median(obj$percent.mt, na.rm = TRUE), 2), "\n")

  return(obj)
}

# ==============================================================================
# LOAD ALL SAMPLES
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("LOADING ALL SAMPLES\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

seurat_list <- list()
batch_mismatch_summary <- data.frame(
  sample = character(),
  rds_batch = character(),
  samplesheet_batch = character(),
  fixed = logical(),
  stringsAsFactors = FALSE
)

for (sample in samples_to_load) {
  path <- input_paths[[sample]]
  seurat_list[[sample]] <- load_and_prepare_sample(
    sample_name = sample,
    file_path = path,
    sample_metadata = params$sample_metadata,
    params = params,
    authoritative_sheet = authoritative_samplesheet
  )
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ALL SAMPLES LOADED\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# ==============================================================================
# CREATE INPUT SUMMARY
# ==============================================================================
cat("\n--- Creating Input Summary ---\n")

input_summary <- data.frame(
  sample = character(),
  sex = character(),
  batch = character(),
  estrous_phase = character(),
  n_cells = integer(),
  n_genes = integer(),
  median_nFeature = numeric(),
  median_nCount = numeric(),
  median_pct_mt = numeric(),
  median_pct_hb = numeric(),
  input_file = character(),
  stringsAsFactors = FALSE
)

for (sample in names(seurat_list)) {
  obj <- seurat_list[[sample]]

  # Get median hemoglobin if available
  median_hb <- if ("percent.hb" %in% colnames(obj@meta.data)) {
    median(obj$percent.hb, na.rm = TRUE)
  } else {
    NA
  }

  # Get estrous phase (may be NA)
  estrous <- if ("estrous_phase" %in% colnames(obj@meta.data)) {
    unique(obj$estrous_phase)[1]
  } else {
    NA
  }

  input_summary <- rbind(input_summary, data.frame(
    sample = sample,
    sex = unique(obj$sex)[1],
    batch = unique(obj$batch)[1],
    estrous_phase = as.character(estrous),
    n_cells = ncol(obj),
    n_genes = nrow(obj),
    median_nFeature = median(obj$nFeature_RNA),
    median_nCount = median(obj$nCount_RNA),
    median_pct_mt = median(obj$percent.mt),
    median_pct_hb = median_hb,
    input_file = basename(input_paths[[sample]]),
    stringsAsFactors = FALSE
  ))
}

cat("\n>>> INPUT DATA SUMMARY <<<\n")
print(input_summary)

# Save summary
write.csv(input_summary, file.path(output_dirs$tables, "01_input_summary.csv"), row.names = FALSE)
cat("\nSaved:", file.path(output_dirs$tables, "01_input_summary.csv"), "\n")

# ==============================================================================
# MERGE ALL SAMPLES
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MERGING ALL SAMPLES\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

if (length(seurat_list) == 1) {
  merged_obj <- seurat_list[[1]]
  cat("Only one sample - no merge needed.\n")
} else {
  cat("Merging", length(seurat_list), "samples...\n")

  # Merge all objects
  merged_obj <- merge(
    x = seurat_list[[1]],
    y = seurat_list[2:length(seurat_list)],
    add.cell.ids = names(seurat_list),
    project = "CP_scRNAseq"
  )

  cat("Merge complete.\n")
}

# Print merged object summary
cat("\n--- Merged Object Summary ---\n")
cat("Total cells:", ncol(merged_obj), "\n")
cat("Total genes:", nrow(merged_obj), "\n")
cat("Assays:", paste(names(merged_obj@assays), collapse = ", "), "\n")
cat("RNA layers:", paste(Layers(merged_obj[["RNA"]]), collapse = ", "), "\n")

cat("\nCells per sample:\n")
print(table(merged_obj$sample_name))

cat("\nCells per sex:\n")
print(table(merged_obj$sex))

cat("\nCells per batch:\n")
print(table(merged_obj$batch))

if (params$batch_variable == "batch") {
  cat("\nCells per integration_batch (same as batch):\n")
  print(table(merged_obj$integration_batch))
}

# ==============================================================================
# BATCH ASSIGNMENT VALIDATION SUMMARY
# ==============================================================================
cat("\n--- Batch Assignment Validation Summary ---\n")
if (!is.null(authoritative_samplesheet)) {
  cat("Authoritative samplesheet used:", samplesheet_path, "\n")
  
  # Build summary of what was checked/fixed
  validation_summary <- data.frame(
    sample = character(),
    final_batch = character(),
    stringsAsFactors = FALSE
  )
  
  for (sample in names(seurat_list)) {
    obj <- seurat_list[[sample]]
    validation_summary <- rbind(validation_summary, data.frame(
      sample = sample,
      final_batch = unique(obj$batch)[1],
      stringsAsFactors = FALSE
    ))
  }
  
  cat("\nFinal batch assignments after validation:\n")
  print(validation_summary)
  
  # Cross-check with samplesheet
  cat("\nCross-check with samplesheet:\n")
  all_match <- TRUE
  for (i in 1:nrow(validation_summary)) {
    s <- validation_summary$sample[i]
    final_b <- validation_summary$final_batch[i]
    
    if ("sample_name" %in% colnames(authoritative_samplesheet)) {
      expected_row <- authoritative_samplesheet[authoritative_samplesheet$sample_name == s, ]
    } else {
      expected_row <- authoritative_samplesheet[authoritative_samplesheet$sample_id == s, ]
    }
    
    if (nrow(expected_row) > 0) {
      expected_b <- expected_row$batch[1]
      if (final_b == expected_b) {
        cat("  [OK]", s, ":", final_b, "\n")
      } else {
        cat("  [ERROR]", s, ": Final=", final_b, ", Expected=", expected_b, "\n")
        all_match <- FALSE
      }
    }
  }
  
  if (all_match) {
    cat("\n>>> ALL BATCH ASSIGNMENTS VALIDATED SUCCESSFULLY <<<\n")
  } else {
    cat("\n>>> WARNING: SOME BATCH ASSIGNMENTS DO NOT MATCH SAMPLESHEET <<<\n")
  }
} else {
  cat("No authoritative samplesheet available for validation.\n")
  cat("Batch assignments taken from params$sample_metadata.\n")
}

# ==============================================================================
# HEMOGLOBIN SUMMARY (if filtering enabled)
# ==============================================================================
if (isTRUE(params$filter_hemoglobin) && "percent.hb" %in% colnames(merged_obj@meta.data)) {
  cat("\n--- Hemoglobin Content Summary ---\n")
  cat("  Pattern used:", params$hemoglobin_pattern, "\n")
  cat("  Median %Hb:", round(median(merged_obj$percent.hb, na.rm = TRUE), 4), "\n")
  cat("  Mean %Hb:", round(mean(merged_obj$percent.hb, na.rm = TRUE), 4), "\n")
  cat("  Max %Hb:", round(max(merged_obj$percent.hb, na.rm = TRUE), 2), "\n")

  max_hb <- if (!is.null(params$max_percent_hb)) params$max_percent_hb else 2
  cells_above_threshold <- sum(merged_obj$percent.hb > max_hb, na.rm = TRUE)
  cat("  Cells >", max_hb, "% Hb:", cells_above_threshold,
      "(", round(cells_above_threshold/ncol(merged_obj)*100, 2), "%)\n", sep = "")
  cat("  Note: Filtering will be applied in Module 02 (QC Validation)\n")
}

# ==============================================================================
# CREATE QC VISUALIZATIONS (PRE-FILTERING)
# ==============================================================================
cat("\n--- Creating Pre-filtering QC Plots ---\n")

# Violin plots
p_violin <- VlnPlot(
  merged_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_name",
  pt.size = 0,
  ncol = 3
) +
  plot_annotation(title = "QC Metrics by Sample (Pre-filtering)")

ggsave(file.path(output_dirs$qc, "01_qc_violin_by_sample.png"),
       p_violin, width = 14, height = 5, dpi = 300)

# Violin plots by sex
p_violin_sex <- VlnPlot(
  merged_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sex",
  pt.size = 0,
  ncol = 3
) +
  plot_annotation(title = "QC Metrics by Sex (Pre-filtering)")

ggsave(file.path(output_dirs$qc, "01_qc_violin_by_sex.png"),
       p_violin_sex, width = 10, height = 5, dpi = 300)

# Violin plots by batch
p_violin_batch <- VlnPlot(
  merged_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "batch",
  pt.size = 0,
  ncol = 3
) +
  plot_annotation(title = "QC Metrics by Batch (Pre-filtering)")

ggsave(file.path(output_dirs$qc, "01_qc_violin_by_batch.png"),
       p_violin_batch, width = 10, height = 5, dpi = 300)

# Hemoglobin violin plot (if enabled)
if (isTRUE(params$filter_hemoglobin) && "percent.hb" %in% colnames(merged_obj@meta.data)) {
  p_hb_violin <- VlnPlot(
    merged_obj,
    features = "percent.hb",
    group.by = "sample_name",
    pt.size = 0
  ) +
    geom_hline(yintercept = params$max_percent_hb, linetype = "dashed", color = "red") +
    labs(title = "Hemoglobin Content by Sample",
         subtitle = paste("Red line: max_percent_hb =", params$max_percent_hb))

  ggsave(file.path(output_dirs$qc, "01_hemoglobin_violin_by_sample.png"),
         p_hb_violin, width = 10, height = 5, dpi = 300)
}

# Scatter plot: nFeature vs nCount colored by sample
p_scatter <- ggplot(merged_obj@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = sample_name)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "nFeature vs nCount by Sample",
       x = "nCount_RNA (log10)",
       y = "nFeature_RNA (log10)",
       color = "Sample")

ggsave(file.path(output_dirs$qc, "01_scatter_nFeature_vs_nCount.png"),
       p_scatter, width = 10, height = 8, dpi = 300)

# Scatter plot: nFeature vs percent.mt
p_scatter_mt <- ggplot(merged_obj@meta.data, aes(x = nFeature_RNA, y = percent.mt, color = sample_name)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_hline(yintercept = params$max_percent_mt, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Percent MT vs nFeature by Sample",
       subtitle = paste("Red line: max_percent_mt =", params$max_percent_mt),
       x = "nFeature_RNA",
       y = "Percent Mitochondrial",
       color = "Sample")

ggsave(file.path(output_dirs$qc, "01_scatter_pctMT_vs_nFeature.png"),
       p_scatter_mt, width = 10, height = 8, dpi = 300)

cat("QC plots saved to:", output_dirs$qc, "\n")

# ==============================================================================
# SAVE OBJECTS FOR NEXT MODULE
# ==============================================================================
cat("\n--- Saving Objects ---\n")

data_file <- file.path(output_dirs$objects, "01_loaded_data.RData")
save(merged_obj, seurat_list, input_summary, file = data_file)
cat("Data saved to:", data_file, "\n")

# Also save individual objects list for potential per-sample analysis
individual_file <- file.path(output_dirs$objects, "01_individual_samples.RData")
save(seurat_list, file = individual_file)
cat("Individual samples saved to:", individual_file, "\n")

# ==============================================================================
# WRITE README
# ==============================================================================
readme_content <- paste0(
  "================================================================================\n",
  "MODULE 01: INPUT DATA LOADING\n",
  "================================================================================\n\n",
  "Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
  "Input data from scCDC preprocessing (Step 7):\n",
  "- scCDC contamination-corrected Seurat objects\n",
  "- Using scCDC_corrected layer as counts for downstream analysis\n\n",
  "--------------------------------------------------------------------------------\n",
  "BATCH VALIDATION\n",
  "--------------------------------------------------------------------------------\n"
)

if (!is.null(authoritative_samplesheet)) {
  readme_content <- paste0(readme_content,
    "Authoritative samplesheet: ", samplesheet_path, "\n",
    "Batch assignments validated and corrected from samplesheet.\n\n"
  )
} else {
  readme_content <- paste0(readme_content,
    "No authoritative samplesheet provided (SAMPLESHEET env var not set).\n",
    "Batch assignments taken from params$sample_metadata.\n\n"
  )
}

readme_content <- paste0(readme_content,
  "--------------------------------------------------------------------------------\n",
  "SAMPLES LOADED\n",
  "--------------------------------------------------------------------------------\n"
)

for (i in 1:nrow(input_summary)) {
  readme_content <- paste0(readme_content,
    input_summary$sample[i], ":\n",
    "  Sex: ", input_summary$sex[i], "\n",
    "  Batch: ", input_summary$batch[i], "\n",
    "  Estrous: ", input_summary$estrous_phase[i], "\n",
    "  Cells: ", input_summary$n_cells[i], "\n",
    "  Genes: ", input_summary$n_genes[i], "\n",
    "  Median %MT: ", round(input_summary$median_pct_mt[i], 2), "\n"
  )
  if (!is.na(input_summary$median_pct_hb[i])) {
    readme_content <- paste0(readme_content,
      "  Median %Hb: ", round(input_summary$median_pct_hb[i], 4), "\n"
    )
  }
  readme_content <- paste0(readme_content, "\n")
}

readme_content <- paste0(readme_content,
  "--------------------------------------------------------------------------------\n",
  "MERGED OBJECT SUMMARY\n",
  "--------------------------------------------------------------------------------\n",
  "Total cells: ", ncol(merged_obj), "\n",
  "Total genes: ", nrow(merged_obj), "\n",
  "Batch variable: ", params$batch_variable, "\n\n",
  "Cells per sex:\n",
  "  Female: ", sum(merged_obj$sex == "Female"), "\n",
  "  Male: ", sum(merged_obj$sex == "Male"), "\n\n",
  "Cells per batch:\n"
)

for (b in unique(merged_obj$batch)) {
  readme_content <- paste0(readme_content, "  ", b, ": ", sum(merged_obj$batch == b), "\n")
}

if (isTRUE(params$filter_hemoglobin)) {
  readme_content <- paste0(readme_content,
    "\n--------------------------------------------------------------------------------\n",
    "HEMOGLOBIN FILTERING\n",
    "--------------------------------------------------------------------------------\n",
    "Hemoglobin filtering: ENABLED\n",
    "Pattern: ", params$hemoglobin_pattern, "\n",
    "Max threshold: ", params$max_percent_hb, "%\n",
    "Cells above threshold: ", sum(merged_obj$percent.hb > params$max_percent_hb, na.rm = TRUE), "\n",
    "Note: Actual filtering applied in Module 02\n"
  )
}

readme_content <- paste0(readme_content,
  "\n--------------------------------------------------------------------------------\n",
  "OUTPUT FILES\n",
  "--------------------------------------------------------------------------------\n",
  "01_loaded_data.RData: Merged Seurat object + individual sample list\n",
  "01_individual_samples.RData: Individual sample objects (for per-sample analysis)\n",
  "01_input_summary.csv: Summary statistics for all loaded samples\n\n",
  "QC Plots:\n",
  "  01_qc_violin_by_sample.png: QC metrics by sample\n",
  "  01_qc_violin_by_sex.png: QC metrics by sex\n",
  "  01_qc_violin_by_batch.png: QC metrics by batch\n",
  "  01_scatter_nFeature_vs_nCount.png: Feature vs count scatter\n",
  "  01_scatter_pctMT_vs_nFeature.png: MT% vs features scatter\n"
)

if (isTRUE(params$filter_hemoglobin)) {
  readme_content <- paste0(readme_content,
    "  01_hemoglobin_violin_by_sample.png: Hemoglobin content by sample\n"
  )
}

readme_content <- paste0(readme_content,
  "\n--------------------------------------------------------------------------------\n",
  "NEXT STEPS\n",
  "--------------------------------------------------------------------------------\n",
  "Module 02: QC validation and filtering\n",
  "  - Apply QC thresholds (min_features, max_features, max_percent_mt)\n"
)

if (isTRUE(params$filter_hemoglobin)) {
  readme_content <- paste0(readme_content,
    "  - Apply hemoglobin filtering (max_percent_hb = ", params$max_percent_hb, ")\n"
  )
}

readme_content <- paste0(readme_content,
  "  - Filter genes (min_cells_per_gene)\n",
  "  - Generate post-filtering QC reports\n"
)

writeLines(readme_content, file.path(output_dirs$qc, "README_01_data_loading.txt"))
cat("\nREADME written to:", file.path(output_dirs$qc, "README_01_data_loading.txt"), "\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 01 COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Summary:\n")
cat("  Samples loaded:", length(seurat_list), "\n")
cat("  Total cells:", ncol(merged_obj), "\n")
cat("  Total genes:", nrow(merged_obj), "\n")
cat("  Female cells:", sum(merged_obj$sex == "Female"), "\n")
cat("  Male cells:", sum(merged_obj$sex == "Male"), "\n")
cat("  Batches:", paste(unique(merged_obj$batch), collapse = ", "), "\n")
cat("  Batch variable for integration:", params$batch_variable, "\n")

if (isTRUE(params$filter_hemoglobin)) {
  cat("  Hemoglobin filtering: ENABLED (threshold:", params$max_percent_hb, "%)\n")
  cat("  Cells above Hb threshold:", sum(merged_obj$percent.hb > params$max_percent_hb, na.rm = TRUE), "\n")
}

cat("\nInput files used scCDC_corrected layer (ambient RNA corrected counts)\n")

cat("\n>>> MODULE 01 COMPLETE <<<\n")