#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 08: DIFFERENTIAL EXPRESSION ANALYSIS (MULTI-SAMPLE PIPELINE)
# ==============================================================================
#
# This module performs differential expression analysis between sexes using:
# 1. MAST - Single-cell hurdle model with proper component merging
# 2. DREAM - Mixed effects model (variancePartition) - RECOMMENDED for nested data
# 3. edgeR - Pseudobulk with pseudo-replicates (when true replicates unavailable)
# 4. DESeq2 - Pseudobulk with pseudo-replicates (when true replicates unavailable)
# 5. Permutation Test - Non-parametric validation
#
# Features:
# - Comprehensive MAST implementation with zlm hurdle model
# - DREAM implementation for mixed effects modeling (Hafner et al. 2025 recommendation)
# - Permutation test for non-parametric validation
# - Negative control with permuted labels to assess FPR
# - Log2FC diagnostic explaining natural log vs log2 scale
# - Pseudobulk with sample-based pseudo-replicates for variance estimation
# - Direction annotation (dynamically derived from factor levels)
# - Per-cluster analysis option
# - Control check printing after CSV saves
# - Volcano plots for each method
# - Concordance analysis: genes significant in same direction across methods
#
# INPUT: Clustered Seurat object from Module 07
# OUTPUT: DE results tables, volcano plots, and concordance tables
#
# UPDATES:
# - 2026-01-03: Improved JoinLayers handling for Seurat v5
# - 2026-01-03: Added layer verification before data extraction
# - 2026-01-04: Use "logFC" component for MAST (preferred over "C")
# - 2026-01-04: Changed gene detection threshold to 0.3% (0.003)
# - 2026-01-04: Use nFeature_RNA for cngeneson calculation
# - 2026-01-04: Added dynamic direction assignment from factor levels
# - 2026-01-04: Added comprehensive input diagnostic function
# - 2026-01-06: Added DREAM implementation (Hafner et al. 2025 recommendation)
# - 2026-01-06: Added permutation test for non-parametric validation
# - 2026-01-06: Added negative control with permuted labels
# - 2026-01-06: Added concordance analysis across methods
# - 2026-01-06: Added true pseudobulk option (one per sample)
# - 2026-01-06: Enhanced warnings for limited biological replicates
#
# REFERENCE:
# Hafner L, et al. (2025) Single-cell differential expression analysis between
# conditions within nested settings. Briefings in Bioinformatics, 26(4), bbaf397.
# https://doi.org/10.1093/bib/bbaf397
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 08: DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("*** CAUTION: Results should be validated with proper biological replicates ***\n\n")

# ==============================================================================
# Load environment and data
# ==============================================================================
pipeline_dir <- Sys.getenv("PIPELINE_DIR", unset = "")
if (pipeline_dir == "") {
  if (file.exists("config/params.R")) {
    pipeline_dir <- normalizePath(".")
  } else {
    stop("PIPELINE_DIR not set. Run Module 00 first.")
  }
}

source(file.path(pipeline_dir, "config", "params.R"))

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

# Check for ggrepel
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_ggrepel) library(ggrepel)

out_base <- params$out_root
load(file.path(out_base, "objects", "pipeline_environment.RData"))

# ==============================================================================
# PARALLEL PROCESSING SETUP FOR MAST
# ==============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SETTING UP PARALLEL PROCESSING FOR MAST\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Load BiocParallel for MAST parallelization
if (!requireNamespace("BiocParallel", quietly = TRUE)) {
  cat("WARNING: BiocParallel not available. MAST will run without parallelization.\n")
  cat("Install with: BiocManager::install('BiocParallel')\n")
  mast_parallel_enabled <- FALSE
  mast_parallel_param <- NULL
  n_cores_mast <- 1
} else {
  library(BiocParallel)
  mast_parallel_enabled <- TRUE
  
  # Determine number of cores from SLURM or system
  n_cores_available <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  if (n_cores_available == 1) {
    # Fallback to detecting cores if SLURM variable not set
    n_cores_available <- max(1, parallel::detectCores() - 1)
  }
  
  # Cap cores to avoid excessive overhead (diminishing returns above ~24 cores for MAST)
  n_cores_mast <- max(1, min(n_cores_available, 24))
  
  cat(sprintf("  SLURM_CPUS_PER_TASK: %s\n", Sys.getenv("SLURM_CPUS_PER_TASK", "<not set>")))
  cat(sprintf("  Detected available cores: %d\n", n_cores_available))
  cat(sprintf("  Cores allocated for MAST: %d\n", n_cores_mast))
  
  # Register the parallel backend
  # Using MulticoreParam for Unix/Linux systems (fork-based parallelism)
  if (.Platform$OS.type == "unix") {
    mast_parallel_param <- MulticoreParam(
      workers = n_cores_mast,
      progressbar = TRUE,
      tasks = 0L,  # Auto-determine task chunking
      stop.on.error = FALSE  # Continue even if some genes fail
    )
    cat(sprintf("  Backend: MulticoreParam (fork-based, %d workers)\n", n_cores_mast))
  } else {
    # Windows fallback (uses PSOCK clusters)
    mast_parallel_param <- SnowParam(
      workers = n_cores_mast,
      progressbar = TRUE,
      stop.on.error = FALSE
    )
    cat(sprintf("  Backend: SnowParam (PSOCK-based, %d workers)\n", n_cores_mast))
  }
  
  # Register as default
  register(mast_parallel_param)
  
  cat(sprintf("  Parallel backend registered: %s\n", class(bpparam())[1]))
}

cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ==============================================================================
# MAST MODEL SELECTION THRESHOLD
# ==============================================================================
# With >40,000 cells, random effects GLMM becomes prohibitively slow (20-30+ hours)
# In such cases, we use fixed effects for MAST and rely on pseudobulk methods
# (DESeq2, DREAM, edgeR) for proper random effects modeling

MAST_MAX_CELLS_FOR_RANDOM_EFFECTS <- 40000

cat("MAST Model Selection Threshold:\n")
cat(sprintf("  Max cells for random effects: %d\n", MAST_MAX_CELLS_FOR_RANDOM_EFFECTS))
cat("  If cell count exceeds this, fixed effects will be used for MAST\n")
cat("  (Pseudobulk methods still provide proper random effects modeling)\n\n")

# ==============================================================================
# Helper Functions (inline for module independence)
# ==============================================================================

# Print first 5 rows of saved CSV for control check
print_csv_head5_local <- function(filepath, label = "") {
  if (file.exists(filepath)) {
    cat("\n--- Control Check:", label, "---\n")
    df <- read.csv(filepath, nrow = 5)
    print(df)
    cat("... (showing first 5 rows)\n")
  }
}

# ==============================================================================
# DYNAMIC DIRECTION ASSIGNMENT FUNCTIONS
# ==============================================================================

#' Setup comparison factor with explicit levels and document the direction logic
#' @param meta Metadata data frame
#' @param var_name Name of the comparison variable
#' @param ref_level Reference level (denominator) - NULL for alphabetical
#' @param comp_level Comparison level (numerator) - NULL for alphabetical
#' @return List with updated meta, levels, and direction labels
setup_comparison_factor <- function(meta, var_name, ref_level = NULL, comp_level = NULL) {

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SETTING UP COMPARISON FACTOR\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  values <- meta[[var_name]]
  unique_vals <- unique(as.character(values[!is.na(values)]))

  cat("Comparison variable:", var_name, "\n")
  cat("Unique values found:", paste(unique_vals, collapse = ", "), "\n")

  if (length(unique_vals) != 2) {
    stop(paste("Comparison variable must have exactly 2 levels. Found:",
               length(unique_vals), "-", paste(unique_vals, collapse = ", ")))
  }

  # Determine levels
  if (is.null(ref_level) || is.null(comp_level)) {
    # Use alphabetical order
    levels_ordered <- sort(unique_vals)
    ref_level <- levels_ordered[1]
    comp_level <- levels_ordered[2]
    cat("\nUsing ALPHABETICAL order for factor levels\n")
  } else {
    # Validate provided levels
    if (!ref_level %in% unique_vals) {
      stop(paste("Reference level '", ref_level, "' not found in data. Available:",
                 paste(unique_vals, collapse = ", ")))
    }
    if (!comp_level %in% unique_vals) {
      stop(paste("Comparison level '", comp_level, "' not found in data. Available:",
                 paste(unique_vals, collapse = ", ")))
    }
    cat("\nUsing USER-SPECIFIED factor levels\n")
  }

  # Create direction labels
  up_label <- paste0(comp_level, "_Higher")
  down_label <- paste0(ref_level, "_Higher")

  cat("\n", paste(rep("-", 50), collapse = ""), "\n")
  cat("DIRECTION ASSIGNMENT LOGIC\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("  Reference level (coded as 0):", ref_level, "\n")
  cat("  Comparison level (coded as 1):", comp_level, "\n")
  cat("\n  INTERPRETATION:\n")
  cat("    Positive logFC → Higher expression in", comp_level, "→", up_label, "\n")
  cat("    Negative logFC → Higher expression in", ref_level, "→", down_label, "\n")
  cat("    Zero logFC     → No difference → None\n")
  cat(paste(rep("-", 50), collapse = ""), "\n\n")

  # Create factor with explicit levels
  meta[[var_name]] <- factor(meta[[var_name]], levels = c(ref_level, comp_level))

  # Verify factor creation
  cat("Factor levels set successfully:\n")
  cat("  levels():", paste(levels(meta[[var_name]]), collapse = " → "), "\n")
  cat("  First level (reference):", levels(meta[[var_name]])[1], "\n")
  cat("  Second level (comparison):", levels(meta[[var_name]])[2], "\n")

  # Show distribution
  cat("\nSample sizes per group:\n")
  print(table(meta[[var_name]], useNA = "always"))

  cat("\n", paste(rep("=", 70), collapse = ""), "\n\n")

  return(list(
    meta = meta,
    reference_level = ref_level,
    comparison_level = comp_level,
    up_label = up_label,
    down_label = down_label,
    contrast_name = paste0(var_name, comp_level)  # e.g., "sexMale"
  ))
}

#' Add direction column dynamically based on established factor levels
#' @param df Data frame with DE results
#' @param lfc_col Name of the log fold change column
#' @param up_label Label for positive logFC (comparison level higher)
#' @param down_label Label for negative logFC (reference level higher)
#' @param none_label Label for zero logFC
#' @param threshold Threshold for considering as "None" (default 0)
#' @return Data frame with direction column added
add_direction_column_dynamic <- function(df, lfc_col, up_label, down_label,
                                          none_label = "None", threshold = 0) {
  if (!lfc_col %in% colnames(df)) {
    warning(paste("Column", lfc_col, "not found. Skipping direction annotation."))
    return(df)
  }

  df$direction <- ifelse(df[[lfc_col]] > threshold, up_label,
                         ifelse(df[[lfc_col]] < -threshold, down_label, none_label))

  # Add numeric direction for filtering convenience
  df$direction_numeric <- sign(df[[lfc_col]])

  return(df)
}

# ==============================================================================
# INPUT DIAGNOSTIC FUNCTION
# ==============================================================================

#' Comprehensive diagnostic for DE input data
#' @param seurat_obj Seurat object to diagnose
#' @param comparison_var Variable used for comparison (e.g., "sex")
#' @param sample_var Sample variable for pseudobulk (e.g., "sample_name")
#' @param covariates Covariates to check
diagnose_de_input <- function(seurat_obj,
                               comparison_var = "sex",
                               sample_var = "sample_name",
                               covariates = c("percent.mt", "nFeature_RNA")) {

  cat("\n")
  cat(paste(rep("#", 80), collapse = ""), "\n")
  cat("#", paste(rep(" ", 25), collapse = ""), "DE INPUT DIAGNOSTIC",
      paste(rep(" ", 25), collapse = ""), "#\n")
  cat(paste(rep("#", 80), collapse = ""), "\n\n")

  meta <- seurat_obj@meta.data

  # ===========================================================================
  # 1. BASIC OBJECT INFO
  # ===========================================================================
  cat("=== 1. BASIC OBJECT INFO ===\n\n")
  cat("Dimensions:", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n")
  cat("Default assay:", DefaultAssay(seurat_obj), "\n")

  # Check layers
  assay <- seurat_obj[[DefaultAssay(seurat_obj)]]
  layers <- Layers(assay)
  cat("Available layers:", paste(layers, collapse = ", "), "\n\n")

  # ===========================================================================
  # 2. COMPARISON VARIABLE ANALYSIS
  # ===========================================================================
  cat("=== 2. COMPARISON VARIABLE: '", comparison_var, "' ===\n\n", sep = "")

  if (!comparison_var %in% colnames(meta)) {
    cat("ERROR: Variable '", comparison_var, "' not found in metadata!\n", sep = "")
    cat("Available columns:", paste(head(colnames(meta), 20), collapse = ", "), "...\n")
    return(invisible(NULL))
  }

  comp_values <- meta[[comparison_var]]

  cat("Class:", class(comp_values), "\n")
  cat("Unique values:", paste(unique(comp_values), collapse = ", "), "\n")
  cat("N per group:\n")
  print(table(comp_values, useNA = "always"))

  # Factor level check
  if (is.factor(comp_values)) {
    cat("\nFactor levels (first = reference):", paste(levels(comp_values), collapse = " -> "), "\n")
  } else {
    cat("\nNOTE: Not yet a factor. Will be converted with explicit levels.\n")
    cat("Alphabetical order would be:", paste(sort(unique(as.character(comp_values))), collapse = " -> "), "\n")
  }

  # ===========================================================================
  # 3. SAMPLE STRUCTURE (for pseudobulk)
  # ===========================================================================
  cat("\n=== 3. SAMPLE STRUCTURE ===\n\n")

  n_bio_reps <- NA

  if (!is.null(sample_var) && sample_var %in% colnames(meta)) {
    cat("Sample variable: '", sample_var, "'\n", sep = "")
    cat("\nSamples per group:\n")
    sample_table <- table(meta[[comparison_var]], meta[[sample_var]])
    print(sample_table)

    # Count biological replicates
    cat("\nBiological replicates per group:\n")
    n_reps_per_group <- c()
    for (grp in unique(meta[[comparison_var]])) {
      grp_samples <- unique(meta[[sample_var]][meta[[comparison_var]] == grp])
      n_reps_per_group[as.character(grp)] <- length(grp_samples)
      cat("  ", grp, ":", length(grp_samples), "sample(s) -",
          paste(grp_samples, collapse = ", "), "\n")
    }

    n_bio_reps <- min(n_reps_per_group)

    if (n_bio_reps < 2) {
      cat("\n")
      cat(paste(rep("!", 70), collapse = ""), "\n")
      cat("!!!", paste(rep(" ", 20), collapse = ""), "CRITICAL WARNING", paste(rep(" ", 20), collapse = ""), "!!!\n")
      cat(paste(rep("!", 70), collapse = ""), "\n")
      cat("\nMinimum biological replicates per group:", n_bio_reps, "\n")
      cat("\nBased on Hafner et al. (2025) Briefings in Bioinformatics:\n")
      cat("  - Pseudoreplication bias is a MAJOR concern with <2 replicates\n")
      cat("  - 'Cells from the same sample are not independent observations'\n")
      cat("  - 'Ignoring this dependence violates the independence assumption'\n")
      cat("  - P-values will be ARTIFICIALLY LOW (inflated Type I error)\n")
      cat("  - Results should be treated as HYPOTHESIS-GENERATING only\n")
      cat("\nPseudo-replicates will be used but DO NOT provide true variance estimates.\n")
      cat("The variance is from TECHNICAL/CELL variation, not biological variation.\n")
      cat(paste(rep("!", 70), collapse = ""), "\n")
    } else if (n_bio_reps < 3) {
      cat("\n")
      cat(paste(rep("*", 60), collapse = ""), "\n")
      cat("*** WARNING: Limited biological replicates (n=", n_bio_reps, ") ***\n", sep = "")
      cat(paste(rep("*", 60), collapse = ""), "\n")
      cat("Statistical power is limited. Consider results exploratory.\n")
    }

  } else {
    cat("Sample variable '", sample_var, "' not found.\n", sep = "")
    cat("Pseudobulk analysis will create pseudo-replicates from all cells.\n")
    n_bio_reps <- NA
  }

  # ===========================================================================
  # 4. COVARIATE CHECK
  # ===========================================================================
  cat("\n=== 4. COVARIATES FOR MODEL ===\n\n")

  for (cov in covariates) {
    if (cov %in% colnames(meta)) {
      cat(cov, ":\n")
      cat("  Class:", class(meta[[cov]]), "\n")
      cat("  Range:", round(min(meta[[cov]], na.rm = TRUE), 3), "-",
          round(max(meta[[cov]], na.rm = TRUE), 3), "\n")
      cat("  NAs:", sum(is.na(meta[[cov]])), "\n")

      # Check for group differences in covariate
      if (is.numeric(meta[[cov]])) {
        means_by_group <- tapply(meta[[cov]], meta[[comparison_var]], mean, na.rm = TRUE)
        cat("  Mean by group:", paste(names(means_by_group), "=",
                                       round(means_by_group, 2), collapse = ", "), "\n")
      }
    } else {
      cat(cov, ": NOT FOUND (will be created if needed)\n")
    }
    cat("\n")
  }

  # ===========================================================================
  # 5. EXPRESSION DATA CHECK
  # ===========================================================================
  cat("=== 5. EXPRESSION DATA ===\n\n")

  # Get counts and data layers
  tryCatch({
    counts <- GetAssayData(seurat_obj, layer = "counts")
    cat("Counts layer:\n")
    cat("  Dimensions:", nrow(counts), "x", ncol(counts), "\n")
    cat("  Sparsity:", round(100 * (1 - Matrix::nnzero(counts) / length(counts)), 1), "%\n")
    cat("  Total counts range:", min(Matrix::colSums(counts)), "-",
        max(Matrix::colSums(counts)), "\n")
  }, error = function(e) {
    cat("Counts layer: ERROR -", e$message, "\n")
  })

  tryCatch({
    data_mat <- GetAssayData(seurat_obj, layer = "data")
    cat("\nData layer (normalized):\n")
    cat("  Dimensions:", nrow(data_mat), "x", ncol(data_mat), "\n")
    cat("  Value range:", round(min(data_mat), 3), "-", round(max(data_mat), 3), "\n")

    # Check if log-transformed
    if (max(data_mat) < 20) {
      cat("  Appears to be: Log-transformed (likely natural log from LogNormalize)\n")
    } else {
      cat("  Appears to be: Raw or linear scale\n")
    }
  }, error = function(e) {
    cat("Data layer: ERROR -", e$message, "\n")
  })

  # ===========================================================================
  # 6. GENE STATISTICS
  # ===========================================================================
  cat("\n=== 6. GENE STATISTICS ===\n\n")

  tryCatch({
    data_mat <- GetAssayData(seurat_obj, layer = "data")

    detection_rate <- rowMeans(data_mat > 0)

    cat("Gene detection rate distribution:\n")
    cat("  Min:", round(min(detection_rate), 4), "\n")
    cat("  1st Qu:", round(quantile(detection_rate, 0.25), 4), "\n")
    cat("  Median:", round(median(detection_rate), 4), "\n")
    cat("  Mean:", round(mean(detection_rate), 4), "\n")
    cat("  3rd Qu:", round(quantile(detection_rate, 0.75), 4), "\n")
    cat("  Max:", round(max(detection_rate), 4), "\n")

    cat("\nGenes by detection threshold:\n")
    thresholds <- c(0, 0.001, 0.003, 0.005, 0.01, 0.05, 0.10)
    for (thresh in thresholds) {
      n_pass <- sum(detection_rate > thresh)
      marker <- if (thresh == 0.003) " <-- CURRENT THRESHOLD" else ""
      cat("  >", sprintf("%5.1f%%", thresh*100), ":", n_pass, "genes (",
          round(100*n_pass/length(detection_rate), 1), "%)", marker, "\n")
    }

  }, error = function(e) {
    cat("Gene statistics: ERROR -", e$message, "\n")
  })

  # ===========================================================================
  # 7. CONFOUNDING CHECK
  # ===========================================================================
  cat("\n=== 7. POTENTIAL CONFOUNDERS ===\n\n")

  # Check for batch/sample confounding with comparison variable
  potential_confounders <- c("batch", "sample_name", "orig.ident", "age", "timepoint")

  found_confounders <- FALSE
  for (conf in potential_confounders) {
    if (conf %in% colnames(meta) && conf != comparison_var) {
      conf_table <- table(meta[[comparison_var]], meta[[conf]])

      # Check if perfectly confounded
      if (nrow(conf_table) == ncol(conf_table) &&
          all(rowSums(conf_table > 0) == 1) &&
          all(colSums(conf_table > 0) == 1)) {
        cat("*** CONFOUNDED: '", conf, "' is perfectly confounded with '",
            comparison_var, "' ***\n", sep = "")
        print(conf_table)
        cat("\n")
        found_confounders <- TRUE
      } else if (any(rowSums(conf_table > 0) == 1) || any(colSums(conf_table > 0) == 1)) {
        cat("PARTIALLY CONFOUNDED: '", conf, "' shows imbalance with '",
            comparison_var, "'\n", sep = "")
        print(conf_table)
        cat("\n")
        found_confounders <- TRUE
      }
    }
  }

  if (!found_confounders) {
    cat("No obvious confounding detected among checked variables.\n")
  }

  # ===========================================================================
  # 8. RECOMMENDATIONS (Updated based on Hafner et al. 2025)
  # ===========================================================================
  cat("\n=== 8. RECOMMENDATIONS (Hafner et al. 2025) ===\n\n")

  recommendations <- c()

  # Replicate recommendation
  if (exists("n_bio_reps") && !is.na(n_bio_reps) && n_bio_reps < 2) {
    recommendations <- c(recommendations,
                         "CRITICAL: Only 1 sample per condition - results are EXPLORATORY only",
                         "Pseudoreplication bias cannot be avoided (Hafner et al. 2025)",
                         "Use multiple methods and focus on CONCORDANT results",
                         "Effect sizes may be more meaningful than p-values",
                         "Consider this analysis hypothesis-generating, not confirmatory")
  } else if (exists("n_bio_reps") && !is.na(n_bio_reps) && n_bio_reps < 3) {
    recommendations <- c(recommendations,
                         "Limited replicates - use DESeq2 pseudobulk as primary (Hafner recommendation)",
                         "DREAM recommended for nested/batch scenarios",
                         "Cross-validate with multiple methods")
  }

  # Method recommendations from paper
  recommendations <- c(recommendations,
                       "Per Hafner et al.: DESeq2 pseudobulk best for single dataset (AUPRC=0.93)",
                       "Per Hafner et al.: DREAM recommended for nested/atlas data",
                       "Per Hafner et al.: Permutation test robust but slow",
                       "Per Hafner et al.: MAST good for single-cell level analysis")

  # Layer recommendation
  if (length(grep("^counts\\.", layers)) > 1 || length(grep("^data\\.", layers)) > 1) {
    recommendations <- c(recommendations,
                         "Join layers before DE: seurat_obj[[\"RNA\"]] <- JoinLayers(seurat_obj[[\"RNA\"]])")
  }

  if (length(recommendations) > 0) {
    for (i in seq_along(recommendations)) {
      cat(i, ". ", recommendations[i], "\n", sep = "")
    }
  } else {
    cat("No critical issues detected.\n")
  }

  cat("\n")
  cat(paste(rep("#", 80), collapse = ""), "\n")
  cat("#", paste(rep(" ", 23), collapse = ""), "END OF DIAGNOSTIC",
      paste(rep(" ", 24), collapse = ""), "#\n")
  cat(paste(rep("#", 80), collapse = ""), "\n\n")

  # Return summary invisibly
  invisible(list(
    n_cells = ncol(seurat_obj),
    n_genes = nrow(seurat_obj),
    comparison_var = comparison_var,
    groups = unique(meta[[comparison_var]]),
    n_per_group = table(meta[[comparison_var]]),
    n_bio_reps = if(exists("n_bio_reps")) n_bio_reps else NA,
    recommendations = recommendations
  ))
}

# MAST log2FC diagnostic function
mast_log2fc_diagnostic <- function() {
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("MAST LOG2FC SCALE DIAGNOSTIC\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("
MAST uses log-transformed expression data and reports coefficients (logFC)
on the NATURAL LOG scale (ln), not log2.
To convert MAST logFC to log2FC for comparison with other methods:
  log2FC = logFC / log(2) = logFC / 0.693
The 'log2FC_if_ln' column in the output applies this conversion.
Interpretation:
- logFC > 0 means higher expression in the comparison level
- logFC < 0 means higher expression in the reference level
- The 'direction' column summarizes this dynamically based on factor levels
Note: This module uses the 'logFC' component from MAST summary when available.
This is a combined estimate from both discrete and continuous model parts.
")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
}

#' Extract MAST logFC - prefers "logFC" component, falls back to "C"
#' @param summaryDt MAST summary datatable
#' @param contrast_name Name of the contrast (e.g., "sexMale")
#' @return Data frame with gene, logFC, confidence intervals, z_score, and source
extract_mast_logfc <- function(summaryDt, contrast_name) {

  components <- unique(summaryDt$component)
  cat("\n--- MAST LogFC Extraction ---\n")
  cat("Available components:", paste(components, collapse = ", "), "\n")

  # Prefer "logFC" component for biological interpretation
  if ("logFC" %in% components) {
    logfc_data <- summaryDt[summaryDt$component == "logFC" &
                              summaryDt$contrast == contrast_name, ]
    source_comp <- "logFC"
    cat("Using 'logFC' component (RECOMMENDED: combined discrete + continuous estimate)\n")
    cat("  This captures overall fold change including expression frequency differences.\n")
  } else if ("C" %in% components) {
    logfc_data <- summaryDt[summaryDt$component == "C" &
                              summaryDt$contrast == contrast_name, ]
    source_comp <- "C"
    cat("WARNING: 'logFC' component not found, using 'C' (continuous part only)\n")
    cat("  This represents fold change among expressing cells only.\n")
    cat("  Differences in expression frequency are NOT captured.\n")
  } else {
    stop("Neither 'logFC' nor 'C' component found in MAST summary. Available: ",
         paste(components, collapse = ", "))
  }

  cat("Rows extracted:", nrow(logfc_data), "\n")

  # Build result data frame with available columns
  result <- data.frame(gene = logfc_data$primerid, stringsAsFactors = FALSE)

  if ("coef" %in% colnames(logfc_data)) {
    result$logFC <- logfc_data$coef
  } else {
    result$logFC <- NA
    warning("'coef' column not found in logFC data")
  }

  if ("ci.hi" %in% colnames(logfc_data)) {
    result$ci_high <- logfc_data$ci.hi
  }
  if ("ci.lo" %in% colnames(logfc_data)) {
    result$ci_low <- logfc_data$ci.lo
  }
  if ("z" %in% colnames(logfc_data)) {
    result$z_score <- logfc_data$z
  }

  result$logFC_source <- source_comp

  cat("LogFC extraction complete. Source:", source_comp, "\n\n")

  return(result)
}

# ==============================================================================
# PSEUDOBULK CREATION FUNCTIONS
# ==============================================================================

#' Create TRUE pseudobulk - one aggregate per sample (RECOMMENDED when possible)
#' Based on Hafner et al. 2025 recommendations
#' @param seurat_obj Seurat object
#' @param group_var Grouping variable for comparison
#' @param sample_var Sample identifier variable
#' @return List with counts matrix and sample info
create_true_pseudobulk <- function(seurat_obj, group_var = "sex",
                                    sample_var = "sample_name") {

  cat("\n--- Creating TRUE Pseudobulk (one per sample) ---\n")
  cat("This is the RECOMMENDED approach per Hafner et al. 2025\n\n")

  counts <- GetAssayData(seurat_obj, layer = "counts")
  meta <- seurat_obj@meta.data

  samples <- unique(meta[[sample_var]])

  pseudobulk_list <- list()
  sample_info <- data.frame()

  for (samp in samples) {
    samp_cells <- which(meta[[sample_var]] == samp)

    if (length(samp_cells) >= 10) {
      pseudobulk_list[[samp]] <- Matrix::rowSums(counts[, samp_cells, drop = FALSE])

      sample_info <- rbind(sample_info, data.frame(
        sample_id = samp,
        group = as.character(meta[[group_var]][samp_cells[1]]),
        n_cells = length(samp_cells),
        stringsAsFactors = FALSE
      ))

      cat("  Sample:", samp, "- Group:", meta[[group_var]][samp_cells[1]],
          "- Cells:", length(samp_cells), "\n")
    } else {
      warning(paste("Sample", samp, "has too few cells (<10), excluding"))
    }
  }

  pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
  colnames(pseudobulk_matrix) <- names(pseudobulk_list)

  cat("\nTrue pseudobulk created:", ncol(pseudobulk_matrix), "samples\n")
  cat("Samples per group:\n")
  print(table(sample_info$group))

  return(list(counts = pseudobulk_matrix, sample_info = sample_info))
}

#' Create pseudobulk with pseudo-replicates (FALLBACK when <2 samples per group)
#' WARNING: This approach has known limitations (Hafner et al. 2025)
#' @param seurat_obj Seurat object
#' @param group_var Grouping variable for comparison
#' @param sample_var Sample identifier variable
#' @param n_pseudo_reps Number of pseudo-replicates per sample
#' @param seed Random seed for reproducibility
#' @return List with counts matrix and sample info
create_pseudobulk_with_replicates_local <- function(seurat_obj, group_var = "sex",
                                                     sample_var = "sample_name",
                                                     n_pseudo_reps = 3, seed = 42) {

  cat("\n")
  cat(paste(rep("!", 70), collapse = ""), "\n")
  cat("!!! WARNING: Creating PSEUDO-REPLICATES !!!\n")
  cat(paste(rep("!", 70), collapse = ""), "\n")
  cat("
Per Hafner et al. (2025):
'Cells from the same sample are NOT independent observations'
'This can lead to underestimated standard errors and inflated type 1 error'
Pseudo-replicates provide TECHNICAL variance estimates only.
Results should be treated as EXPLORATORY.
")
  cat(paste(rep("!", 70), collapse = ""), "\n\n")

  set.seed(seed)

  counts <- GetAssayData(seurat_obj, layer = "counts")
  meta <- seurat_obj@meta.data

  groups <- unique(meta[[group_var]])

  pseudobulk_list <- list()
  sample_info <- data.frame()

  for (grp in groups) {
    grp_cells <- which(meta[[group_var]] == grp)

    if (sample_var %in% colnames(meta)) {
      samples_in_grp <- unique(meta[[sample_var]][grp_cells])

      for (samp in samples_in_grp) {
        samp_cells <- which(meta[[group_var]] == grp & meta[[sample_var]] == samp)

        if (length(samp_cells) >= 10) {
          cell_splits <- split(samp_cells, cut(seq_along(samp_cells), n_pseudo_reps, labels = FALSE))

          for (rep_idx in 1:length(cell_splits)) {
            rep_cells <- cell_splits[[rep_idx]]
            if (length(rep_cells) >= 3) {
              rep_name <- paste0(grp, "_", samp, "_rep", rep_idx)
              pseudobulk_list[[rep_name]] <- Matrix::rowSums(counts[, rep_cells, drop = FALSE])
              sample_info <- rbind(sample_info, data.frame(
                sample_id = rep_name,
                group = grp,
                original_sample = samp,
                replicate = rep_idx,
                n_cells = length(rep_cells)
              ))
            }
          }
        }
      }
    } else {
      shuffled_cells <- sample(grp_cells)
      cell_splits <- split(shuffled_cells, cut(seq_along(shuffled_cells), n_pseudo_reps, labels = FALSE))

      for (rep_idx in 1:length(cell_splits)) {
        rep_cells <- cell_splits[[rep_idx]]
        if (length(rep_cells) >= 10) {
          rep_name <- paste0(grp, "_rep", rep_idx)
          pseudobulk_list[[rep_name]] <- Matrix::rowSums(counts[, rep_cells, drop = FALSE])
          sample_info <- rbind(sample_info, data.frame(
            sample_id = rep_name,
            group = grp,
            original_sample = NA,
            replicate = rep_idx,
            n_cells = length(rep_cells)
          ))
        }
      }
    }
  }

  pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
  colnames(pseudobulk_matrix) <- names(pseudobulk_list)

  return(list(counts = pseudobulk_matrix, sample_info = sample_info))
}

#' Smart pseudobulk creation - chooses method based on available replicates
#' @param seurat_obj Seurat object
#' @param group_var Grouping variable
#' @param sample_var Sample variable
#' @param n_pseudo_reps Number of pseudo-replicates if needed
#' @param seed Random seed
#' @return List with counts, sample_info, and method used
create_pseudobulk_smart <- function(seurat_obj, group_var = "sex",
                                     sample_var = "sample_name",
                                     n_pseudo_reps = 3, seed = 42) {

  meta <- seurat_obj@meta.data

  # Count biological replicates per group
  if (sample_var %in% colnames(meta)) {
    n_reps_per_group <- c()
    for (grp in unique(meta[[group_var]])) {
      grp_samples <- unique(meta[[sample_var]][meta[[group_var]] == grp])
      n_reps_per_group[as.character(grp)] <- length(grp_samples)
    }
    min_reps <- min(n_reps_per_group)
  } else {
    min_reps <- 0
  }

  cat("\n--- Smart Pseudobulk Creation ---\n")
  cat("Minimum biological replicates per group:", min_reps, "\n")

  if (min_reps >= 2) {
    cat("Using TRUE pseudobulk (one per sample) - RECOMMENDED\n")
    pb <- create_true_pseudobulk(seurat_obj, group_var, sample_var)
    pb$method <- "true_pseudobulk"
    pb$has_true_replicates <- TRUE
  } else {
    cat("Insufficient replicates - using PSEUDO-REPLICATES (exploratory)\n")
    pb <- create_pseudobulk_with_replicates_local(seurat_obj, group_var, sample_var,
                                                   n_pseudo_reps, seed)
    pb$method <- "pseudo_replicates"
    pb$has_true_replicates <- FALSE
  }

  return(pb)
}

# ==============================================================================
# PERMUTATION TEST FUNCTION (Hafner et al. 2025)
# ==============================================================================

#' Permutation test for differential expression
#' Non-parametric method that doesn't assume specific distribution
#' @param pb_counts Pseudobulk count matrix
#' @param sample_info Sample metadata with 'sample_id' and 'group' columns
#' @param n_permutations Number of permutations
#' @param seed Random seed
#' @return Data frame with DE results
run_permutation_test <- function(pb_counts, sample_info,
                                  n_permutations = 10000, seed = 42) {

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("PERMUTATION TEST FOR DIFFERENTIAL EXPRESSION\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("
Per Hafner et al. (2025):
- Non-parametric method, doesn't assume specific distribution
- 'Permutation test performs best' in atlas scenarios with batch effects
- Computationally expensive but robust
- Minimum p-value limited to 1/n_permutations
")
  cat("Number of permutations:", n_permutations, "\n")
  cat("Minimum achievable p-value:", 1/n_permutations, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  set.seed(seed)

  groups <- sample_info$group
  group_levels <- unique(groups)

  if (length(group_levels) != 2) {
    stop("Permutation test requires exactly 2 groups")
  }

  # Get sample indices for each group
  group1_samples <- sample_info$sample_id[groups == group_levels[1]]
  group2_samples <- sample_info$sample_id[groups == group_levels[2]]

  # Normalize by library size (CPM) and log transform
  lib_sizes <- colSums(pb_counts)
  pb_cpm <- t(t(pb_counts) / lib_sizes * 1e6)
  pb_log <- log2(pb_cpm + 1)

  # Calculate observed difference in means (group2 - group1)
  observed_diff <- rowMeans(pb_log[, group2_samples, drop = FALSE]) -
                   rowMeans(pb_log[, group1_samples, drop = FALSE])

  cat("Calculating permutation distribution...\n")

  # Permutation
  n_samples <- ncol(pb_log)
  n_genes <- nrow(pb_log)

  # Initialize counter for extreme values
  extreme_count <- rep(0, n_genes)

  # Progress tracking
  progress_interval <- max(1, n_permutations %/% 10)

  for (i in 1:n_permutations) {
    if (i %% progress_interval == 0) {
      cat("  Progress:", round(100 * i / n_permutations), "%\n")
    }

    # Permute group labels
    perm_idx <- sample(1:n_samples)
    perm_groups <- groups[perm_idx]

    g1 <- colnames(pb_log)[perm_groups == group_levels[1]]
    g2 <- colnames(pb_log)[perm_groups == group_levels[2]]

    perm_diff <- rowMeans(pb_log[, g2, drop = FALSE]) -
                 rowMeans(pb_log[, g1, drop = FALSE])

    # Count permutations with more extreme difference
    extreme_count <- extreme_count + (abs(perm_diff) >= abs(observed_diff))
  }

  # Calculate two-sided p-values
  p_values <- (extreme_count + 1) / (n_permutations + 1)

  # Adjust for multiple testing
  adj_p_values <- p.adjust(p_values, method = "BH")

  results <- data.frame(
    gene = rownames(pb_counts),
    log2FC = observed_diff,
    pvalue = p_values,
    padj = adj_p_values,
    stringsAsFactors = FALSE
  )

  results <- results[order(results$padj), ]

  cat("\nPermutation test complete.\n")

  return(results)
}

# ==============================================================================
# NEGATIVE CONTROL FUNCTION (Hafner et al. 2025)
# ==============================================================================

#' Run negative control with permuted labels
#' Assesses false positive rate when no true DE exists
#' @param seurat_obj Seurat object
#' @param comparison_var Variable to permute
#' @param sample_var Sample variable
#' @param method DE method to use ("MAST", "DESeq2", "edgeR")
#' @param n_iterations Number of iterations
#' @param seed Random seed
#' @return List with FPR estimates at different thresholds
run_negative_control <- function(seurat_obj,
                                  comparison_var = "sex",
                                  sample_var = "sample_name",
                                  method = "DESeq2",
                                  n_iterations = 3,
                                  seed = 42) {

  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("NEGATIVE CONTROL: ASSESSING FALSE POSITIVE RATE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("
Per Hafner et al. (2025):
This analysis permutes condition labels to create a scenario with
NO TRUE differential expression. Any 'significant' genes are false positives.
This helps calibrate method performance and detect inflated Type I error.
Method:", method, "\n")
  cat("Iterations:", n_iterations, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  set.seed(seed)

  thresholds <- c(0.001, 0.01, 0.05, 0.1)
  fpr_results <- matrix(0, nrow = n_iterations, ncol = length(thresholds))
  colnames(fpr_results) <- paste0("FPR_", thresholds)

  for (iter in 1:n_iterations) {
    cat("Iteration", iter, "of", n_iterations, "...\n")

    # Create copy with permuted labels
    obj_perm <- seurat_obj
    original_labels <- obj_perm@meta.data[[comparison_var]]

    # Permute at sample level if possible, otherwise at cell level
    if (sample_var %in% colnames(obj_perm@meta.data)) {
      samples <- unique(obj_perm@meta.data[[sample_var]])
      sample_labels <- sapply(samples, function(s) {
        obj_perm@meta.data[[comparison_var]][obj_perm@meta.data[[sample_var]] == s][1]
      })
      perm_sample_labels <- sample(sample_labels)
      names(perm_sample_labels) <- samples

      obj_perm@meta.data[[comparison_var]] <- perm_sample_labels[obj_perm@meta.data[[sample_var]]]
    } else {
      obj_perm@meta.data[[comparison_var]] <- sample(original_labels)
    }

    # Run simplified DE based on method
    tryCatch({
      if (method == "DESeq2") {
        # Quick pseudobulk DESeq2
        pb <- create_pseudobulk_smart(obj_perm, comparison_var,
                               sample_var, n_pseudo_reps = 3,
                               seed = seed + iter)
        pb_counts_int <- round(pb$counts)
        mode(pb_counts_int) <- "integer"

        keep <- rowSums(pb_counts_int >= 10) >= 2
        pb_counts_filtered <- pb_counts_int[keep, ]

        pb$sample_info$group <- factor(pb$sample_info$group)

        suppressMessages({
          dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = pb_counts_filtered,
            colData = pb$sample_info,
            design = ~ group
          )
          dds <- DESeq2::DESeq(dds, quiet = TRUE)
          res <- DESeq2::results(dds)
        })

        pvals <- res$padj
        pvals <- pvals[!is.na(pvals)]

      } else if (method == "edgeR") {
        # Quick pseudobulk edgeR
        pb <- create_pseudobulk_with_replicates_local(obj_perm, comparison_var,
                                                       sample_var, n_pseudo_reps = 3,
                                                       seed = seed + iter)
        dge <- edgeR::DGEList(counts = pb$counts, group = pb$sample_info$group)
        keep <- edgeR::filterByExpr(dge, min.count = 10)
        dge <- dge[keep, , keep.lib.sizes = FALSE]
        dge <- edgeR::calcNormFactors(dge)

        design <- model.matrix(~ 0 + group, data = pb$sample_info)
        colnames(design) <- gsub("group", "", colnames(design))

        v <- limma::voom(dge, design, plot = FALSE)
        fit <- limma::lmFit(v, design)

        groups <- unique(pb$sample_info$group)
        contrast_formula <- paste(groups[2], "-", groups[1])
        contrast_matrix <- limma::makeContrasts(contrasts = contrast_formula, levels = design)

        fit2 <- limma::contrasts.fit(fit, contrast_matrix)
        fit2 <- limma::eBayes(fit2)

        res <- limma::topTable(fit2, number = Inf)
        pvals <- res$adj.P.Val
        pvals <- pvals[!is.na(pvals)]

      } else {
        warning("Method not implemented for negative control:", method)
        next
      }

      # Calculate FPR at each threshold
      for (t_idx in seq_along(thresholds)) {
        fpr_results[iter, t_idx] <- mean(pvals < thresholds[t_idx])
      }

    }, error = function(e) {
      cat("  Error in iteration", iter, ":", conditionMessage(e), "\n")
    })
  }

  # Summarize results
  cat("\n--- NEGATIVE CONTROL RESULTS ---\n")
  cat("\nExpected FPR (nominal) vs Observed FPR:\n\n")

  summary_df <- data.frame(
    threshold = thresholds,
    expected_fpr = thresholds,
    observed_fpr_mean = colMeans(fpr_results, na.rm = TRUE),
    observed_fpr_sd = apply(fpr_results, 2, sd, na.rm = TRUE)
  )

  summary_df$inflation_ratio <- summary_df$observed_fpr_mean / summary_df$expected_fpr
  summary_df$calibration <- ifelse(summary_df$inflation_ratio > 2, "INFLATED",
                                    ifelse(summary_df$inflation_ratio < 0.5, "Conservative", "OK"))

  print(summary_df)

  if (any(summary_df$inflation_ratio > 2)) {
    cat("\n")
    cat(paste(rep("!", 60), collapse = ""), "\n")
    cat("WARNING: Type I error appears INFLATED (>2x expected)\n")
    cat("This is expected with pseudo-replicates (Hafner et al. 2025)\n")
    cat("P-values may be artificially low. Interpret with caution.\n")
    cat(paste(rep("!", 60), collapse = ""), "\n")
  }

  return(list(
    fpr_matrix = fpr_results,
    summary = summary_df,
    method = method,
    n_iterations = n_iterations
  ))
}

# ==============================================================================
# CONCORDANCE ANALYSIS FUNCTIONS
# ==============================================================================

#' Calculate concordance between DE methods
#' Identifies genes significant in same direction across multiple methods
#' @param all_results List of DE results from different methods
#' @param fdr_threshold FDR threshold for significance
#' @param lfc_threshold Log fold change threshold (optional)
#' @param up_label Label for upregulated genes
#' @param down_label Label for downregulated genes
#' @return List with concordance tables
calculate_concordance <- function(all_results, fdr_threshold = 0.05,
                                   lfc_threshold = 0, up_label, down_label) {

  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("CONCORDANCE ANALYSIS ACROSS METHODS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("FDR threshold:", fdr_threshold, "\n")
  cat("LFC threshold:", lfc_threshold, "\n")
  cat("Methods included:", paste(names(all_results), collapse = ", "), "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  # Standardize column names and extract significant genes with direction
  method_sig_genes <- list()

  for (method_name in names(all_results)) {
    df <- all_results[[method_name]]

    # Identify relevant columns
    gene_col <- "gene"
    if ("gene" %in% colnames(df)) {
      gene_col <- "gene"
    } else if ("primerid" %in% colnames(df)) {
      gene_col <- "primerid"
    }

    # Identify FDR column
    fdr_col <- NULL
    if ("FDR" %in% colnames(df)) fdr_col <- "FDR"
    else if ("padj" %in% colnames(df)) fdr_col <- "padj"
    else if ("adj.P.Val" %in% colnames(df)) fdr_col <- "adj.P.Val"

    # Identify LFC column
    lfc_col <- NULL
    if ("logFC" %in% colnames(df)) lfc_col <- "logFC"
    else if ("log2FC" %in% colnames(df)) lfc_col <- "log2FC"
    else if ("log2FoldChange" %in% colnames(df)) lfc_col <- "log2FoldChange"

    if (is.null(fdr_col) || is.null(lfc_col)) {
      cat("Skipping", method_name, "- missing required columns\n")
      next
    }

    # Get significant genes with direction
    sig_idx <- which(df[[fdr_col]] < fdr_threshold & abs(df[[lfc_col]]) > lfc_threshold)

    if (length(sig_idx) > 0) {
      sig_df <- data.frame(
        gene = df[[gene_col]][sig_idx],
        direction = ifelse(df[[lfc_col]][sig_idx] > 0, "up", "down"),
        lfc = df[[lfc_col]][sig_idx],
        fdr = df[[fdr_col]][sig_idx],
        stringsAsFactors = FALSE
      )
      method_sig_genes[[method_name]] <- sig_df
      cat(method_name, ":", nrow(sig_df), "significant genes\n")
    } else {
      method_sig_genes[[method_name]] <- data.frame(
        gene = character(0),
        direction = character(0),
        lfc = numeric(0),
        fdr = numeric(0)
      )
      cat(method_name, ": 0 significant genes\n")
    }
  }

  if (length(method_sig_genes) < 2) {
    cat("\nInsufficient methods with results for concordance analysis.\n")
    return(NULL)
  }

  # Get all genes across all methods
  all_genes <- unique(unlist(lapply(method_sig_genes, function(x) x$gene)))
  cat("\nTotal unique significant genes across all methods:", length(all_genes), "\n")

  # Create gene-method direction matrix
  direction_matrix <- matrix("NS", nrow = length(all_genes), ncol = length(method_sig_genes))
  rownames(direction_matrix) <- all_genes
  colnames(direction_matrix) <- names(method_sig_genes)

  for (method_name in names(method_sig_genes)) {
    df <- method_sig_genes[[method_name]]
    for (i in seq_len(nrow(df))) {
      direction_matrix[df$gene[i], method_name] <- df$direction[i]
    }
  }

  # Calculate concordance for each gene
  concordance_df <- data.frame(
    gene = all_genes,
    stringsAsFactors = FALSE
  )

  for (method_name in names(method_sig_genes)) {
    concordance_df[[paste0(method_name, "_direction")]] <- direction_matrix[, method_name]
  }

  # Count methods where gene is significant
  concordance_df$n_methods_sig <- apply(direction_matrix, 1, function(x) sum(x != "NS"))

  # Count methods where gene is UP
  concordance_df$n_methods_up <- apply(direction_matrix, 1, function(x) sum(x == "up"))

  # Count methods where gene is DOWN
  concordance_df$n_methods_down <- apply(direction_matrix, 1, function(x) sum(x == "down"))

  # Determine concordant direction
  concordance_df$concordant_direction <- apply(concordance_df, 1, function(row) {
    n_up <- as.numeric(row["n_methods_up"])
    n_down <- as.numeric(row["n_methods_down"])
    n_sig <- as.numeric(row["n_methods_sig"])

    if (n_sig == 0) return("NS")
    if (n_up == n_sig) return("up")
    if (n_down == n_sig) return("down")
    return("discordant")
  })

  # Sort by number of methods
  concordance_df <- concordance_df[order(-concordance_df$n_methods_sig,
                                          concordance_df$concordant_direction), ]

  # Generate concordance subsets
  n_methods <- length(method_sig_genes)
  method_names <- names(method_sig_genes)

  concordance_subsets <- list()

  # Genes significant in ALL methods (same direction)
  all_concordant_up <- concordance_df$gene[concordance_df$n_methods_sig == n_methods &
                                             concordance_df$concordant_direction == "up"]
  all_concordant_down <- concordance_df$gene[concordance_df$n_methods_sig == n_methods &
                                               concordance_df$concordant_direction == "down"]

  concordance_subsets[["all_methods"]] <- list(
    up = all_concordant_up,
    down = all_concordant_down,
    n_methods = n_methods
  )

  cat("\n--- Genes significant in ALL", n_methods, "methods (same direction) ---\n")
  cat("  Upregulated (", up_label, "):", length(all_concordant_up), "\n")
  cat("  Downregulated (", down_label, "):", length(all_concordant_down), "\n")

  # Generate pairwise and higher-order concordances
  if (n_methods >= 2) {
    for (k in 2:(n_methods - 1)) {
      # For each combination of k methods
      method_combos <- combn(method_names, k, simplify = FALSE)

      for (combo in method_combos) {
        combo_name <- paste(combo, collapse = "_AND_")

        # Find genes significant in these specific methods (at least)
        combo_cols <- paste0(combo, "_direction")

        genes_in_combo <- concordance_df$gene[
          apply(concordance_df[, combo_cols, drop = FALSE], 1, function(row) {
            all(row != "NS")
          })
        ]

        # Check direction concordance
        if (length(genes_in_combo) > 0) {
          combo_directions <- concordance_df[concordance_df$gene %in% genes_in_combo, combo_cols, drop = FALSE]

          concordant_up <- genes_in_combo[apply(combo_directions, 1, function(row) all(row == "up"))]
          concordant_down <- genes_in_combo[apply(combo_directions, 1, function(row) all(row == "down"))]

          concordance_subsets[[combo_name]] <- list(
            up = concordant_up,
            down = concordant_down,
            n_methods = k,
            methods = combo
          )
        }
      }
    }
  }

  # Summary by number of methods
  cat("\n--- Summary by number of concordant methods ---\n")
  for (k in n_methods:2) {
    genes_k <- concordance_df$gene[concordance_df$n_methods_sig >= k &
                                     concordance_df$concordant_direction %in% c("up", "down")]
    up_k <- sum(concordance_df$concordant_direction[concordance_df$gene %in% genes_k] == "up")
    down_k <- sum(concordance_df$concordant_direction[concordance_df$gene %in% genes_k] == "down")
    cat("  >=", k, "methods:", length(genes_k), "genes (",
        up_k, "up,", down_k, "down)\n")
  }

  return(list(
    full_table = concordance_df,
    direction_matrix = direction_matrix,
    method_sig_genes = method_sig_genes,
    subsets = concordance_subsets,
    n_methods = n_methods,
    method_names = method_names
  ))
}

#' Save concordance results to CSV files
#' @param concordance_result Output from calculate_concordance()
#' @param output_dir Directory to save files
#' @param up_label Label for upregulated genes
#' @param down_label Label for downregulated genes
save_concordance_results <- function(concordance_result, output_dir,
                                      up_label, down_label) {

  if (is.null(concordance_result)) {
    cat("No concordance results to save.\n")
    return(invisible(NULL))
  }

  cat("\n--- Saving concordance results ---\n")

  # Create concordance subdirectory
  conc_dir <- file.path(output_dir, "concordance")
  dir.create(conc_dir, showWarnings = FALSE, recursive = TRUE)

  # Save full concordance table
  full_path <- file.path(conc_dir, "concordance_full_table.csv")
  write.csv(concordance_result$full_table, full_path, row.names = FALSE)
  cat("Saved:", full_path, "\n")

  # Save genes by concordance level
  n_methods <- concordance_result$n_methods

  for (k in n_methods:2) {
    # Genes concordant in >= k methods (UPREGULATED)
    genes_up <- concordance_result$full_table$gene[
      concordance_result$full_table$n_methods_sig >= k &
        concordance_result$full_table$concordant_direction == "up"
    ]

    if (length(genes_up) > 0) {
      up_df <- concordance_result$full_table[concordance_result$full_table$gene %in% genes_up, ]
      up_path <- file.path(conc_dir, paste0("genes_", up_label, "_in_", k, "plus_methods.csv"))
      write.csv(up_df, up_path, row.names = FALSE)
      cat("Saved:", up_path, "(", nrow(up_df), "genes)\n")
    }

    # Genes concordant in >= k methods (DOWNREGULATED)
    genes_down <- concordance_result$full_table$gene[
      concordance_result$full_table$n_methods_sig >= k &
        concordance_result$full_table$concordant_direction == "down"
    ]

    if (length(genes_down) > 0) {
      down_df <- concordance_result$full_table[concordance_result$full_table$gene %in% genes_down, ]
      down_path <- file.path(conc_dir, paste0("genes_", down_label, "_in_", k, "plus_methods.csv"))
      write.csv(down_df, down_path, row.names = FALSE)
      cat("Saved:", down_path, "(", nrow(down_df), "genes)\n")
    }
  }

  # Save specific method combinations
  for (subset_name in names(concordance_result$subsets)) {
    subset_data <- concordance_result$subsets[[subset_name]]

    # Up genes
    if (length(subset_data$up) > 0) {
      up_df <- concordance_result$full_table[concordance_result$full_table$gene %in% subset_data$up, ]
      filename <- paste0("genes_", up_label, "_", gsub("_AND_", "_", subset_name), ".csv")
      filepath <- file.path(conc_dir, filename)
      write.csv(up_df, filepath, row.names = FALSE)
      cat("Saved:", filepath, "(", nrow(up_df), "genes)\n")
    }

    # Down genes
    if (length(subset_data$down) > 0) {
      down_df <- concordance_result$full_table[concordance_result$full_table$gene %in% subset_data$down, ]
      filename <- paste0("genes_", down_label, "_", gsub("_AND_", "_", subset_name), ".csv")
      filepath <- file.path(conc_dir, filename)
      write.csv(down_df, filepath, row.names = FALSE)
      cat("Saved:", filepath, "(", nrow(down_df), "genes)\n")
    }
  }

  # Save summary statistics
  summary_df <- data.frame(
    n_methods_minimum = n_methods:2,
    n_genes_concordant_up = sapply(n_methods:2, function(k) {
      sum(concordance_result$full_table$n_methods_sig >= k &
            concordance_result$full_table$concordant_direction == "up")
    }),
    n_genes_concordant_down = sapply(n_methods:2, function(k) {
      sum(concordance_result$full_table$n_methods_sig >= k &
            concordance_result$full_table$concordant_direction == "down")
    })
  )
  summary_df$n_genes_total <- summary_df$n_genes_concordant_up + summary_df$n_genes_concordant_down

  summary_path <- file.path(conc_dir, "concordance_summary.csv")
  write.csv(summary_df, summary_path, row.names = FALSE)
  cat("Saved:", summary_path, "\n")

  cat("\nConcordance results saved to:", conc_dir, "\n")

  return(invisible(conc_dir))
}

# ==============================================================================
# Load clustered data
# ==============================================================================
cat("--- Loading clustered data ---\n")

leiden_file <- file.path(out_base, "objects", "07_leiden_data.RData")
scice_rds <- file.path(out_base, "objects", "scice_subclustered_object.rds")
final_rds <- file.path(out_base, "objects", "07_final_object.rds")
choir_rds <- file.path(out_base, "objects", "choir_clustered_object.rds")

if (file.exists(scice_rds)) {
  clustered_obj <- readRDS(scice_rds)
  clustering_method_used <- "scICE"
  cat("Loaded scICE subclustered object.\n")
} else if (file.exists(leiden_file)) {
  load(leiden_file)
  cat("Loaded Leiden clustered data.\n")
  cat("Clustering method:", clustering_method_used, "\n")
} else if (file.exists(final_rds)) {
  clustered_obj <- readRDS(final_rds)
  clustering_method_used <- "unknown"
  cat("Loaded final clustered object.\n")
} else if (file.exists(choir_rds)) {
  clustered_obj <- readRDS(choir_rds)
  clustering_method_used <- "CHOIR"
  cat("Loaded CHOIR clustered object.\n")
} else {
  stop("No clustered data found. Run clustering modules first.")
}

if (is.null(clustered_obj)) {
  stop("No clustered object available.")
}

# ==============================================================================
# Setup (Seurat v5 compatibility)
# ==============================================================================
DefaultAssay(clustered_obj) <- "RNA"

# Join layers for DE analysis (critical for Seurat v5)
tryCatch({
  rna_layers <- Layers(clustered_obj[["RNA"]])
  counts_layers <- grep("^counts", rna_layers, value = TRUE)
  data_layers <- grep("^data", rna_layers, value = TRUE)

  if (length(counts_layers) > 1 || length(data_layers) > 1) {
    cat("Found split layers - joining for DE analysis...\n")
    clustered_obj[["RNA"]] <- JoinLayers(clustered_obj[["RNA"]])
    cat("Successfully joined RNA layers\n")
  } else {
    cat("RNA layers already unified\n")
  }
}, error = function(e) {
  cat("Note: Layer preparation message -", e$message, "\n")
})

# Verify counts and data layers are accessible
tryCatch({
  counts_check <- GetAssayData(clustered_obj, layer = "counts")
  data_check <- GetAssayData(clustered_obj, layer = "data")
  cat("Verified: counts (", nrow(counts_check), " genes) and data layers accessible\n", sep = "")
}, error = function(e) {
  cat("WARNING: Could not verify layer access -", e$message, "\n")
})

# ==============================================================================
# Run Input Diagnostic
# ==============================================================================
cat("\n>>> RUNNING INPUT DIAGNOSTIC <<<\n")
diagnostic_result <- diagnose_de_input(
  clustered_obj,
  comparison_var = "sex",
  sample_var = "sample_name",
  covariates = c("percent.mt", "nFeature_RNA")
)

cat("\n>>> OBJECT FOR DE ANALYSIS <<<\n")
print_object_structure(clustered_obj, "DE Input")

cat("\n--- Sex distribution ---\n")
print(table(clustered_obj$sex))

# Store number of biological replicates for method decisions
N_BIO_REPS <- diagnostic_result$n_bio_reps
if (is.na(N_BIO_REPS)) N_BIO_REPS <- 0

cat("\nBiological replicates detected:", N_BIO_REPS, "\n")

# ==============================================================================
# Setup Comparison Factor with Dynamic Direction Labels
# ==============================================================================
factor_setup <- setup_comparison_factor(
  meta = clustered_obj@meta.data,
  var_name = "sex",
  ref_level = "Female",
  comp_level = "Male"
)

# Update object metadata with properly factored comparison variable
clustered_obj@meta.data <- factor_setup$meta

# Store direction labels for use throughout the script
REF_LEVEL <- factor_setup$reference_level
COMP_LEVEL <- factor_setup$comparison_level
UP_LABEL <- factor_setup$up_label
DOWN_LABEL <- factor_setup$down_label
CONTRAST_NAME <- factor_setup$contrast_name

cat("Direction labels for this analysis:\n")
cat("  Positive logFC label:", UP_LABEL, "\n")
cat("  Negative logFC label:", DOWN_LABEL, "\n")
cat("  MAST contrast name:", CONTRAST_NAME, "\n\n")

# ==============================================================================
# Setup Output Directories
# ==============================================================================
de_base <- output_dirs$de_results
de_mast <- subdirs$de_mast
de_edger <- subdirs$de_edger
de_deseq2 <- subdirs$de_deseq2
de_plots <- subdirs$plots_de

# Validate and create fallback directories if needed
if (is.null(de_base) || length(de_base) == 0 || de_base == "") {
  de_base <- file.path(out_base, "08_Differential_Expression")
  cat("WARNING: de_base not found in environment, using:", de_base, "\n")
}
if (is.null(de_mast) || length(de_mast) == 0 || de_mast == "") {
  de_mast <- file.path(de_base, "MAST")
}
if (is.null(de_edger) || length(de_edger) == 0 || de_edger == "") {
  de_edger <- file.path(de_base, "edgeR")
}
if (is.null(de_deseq2) || length(de_deseq2) == 0 || de_deseq2 == "") {
  de_deseq2 <- file.path(de_base, "DESeq2")
}
if (is.null(de_plots) || length(de_plots) == 0 || de_plots == "") {
  de_plots <- file.path(out_base, "plots", "de")
}

# Additional directories for new methods
de_dream <- file.path(de_base, "DREAM")
de_permutation <- file.path(de_base, "Permutation")
de_negative_control <- file.path(de_base, "NegativeControl")

# Ensure directories exist
dir.create(de_base, showWarnings = FALSE, recursive = TRUE)
dir.create(de_mast, showWarnings = FALSE, recursive = TRUE)
dir.create(de_edger, showWarnings = FALSE, recursive = TRUE)
dir.create(de_deseq2, showWarnings = FALSE, recursive = TRUE)
dir.create(de_dream, showWarnings = FALSE, recursive = TRUE)
dir.create(de_permutation, showWarnings = FALSE, recursive = TRUE)
dir.create(de_negative_control, showWarnings = FALSE, recursive = TRUE)
dir.create(de_plots, showWarnings = FALSE, recursive = TRUE)

cat("DE output directories:\n")
cat("  Base:", de_base, "\n")
cat("  MAST:", de_mast, "\n")
cat("  DREAM:", de_dream, "\n")
cat("  edgeR:", de_edger, "\n")
cat("  DESeq2:", de_deseq2, "\n")
cat("  Permutation:", de_permutation, "\n")
cat("  NegativeControl:", de_negative_control, "\n")
cat("  Plots:", de_plots, "\n\n")

all_de_results <- list()
de_summary <- data.frame()

# ==============================================================================
# SECTION 1: MAST ANALYSIS
# ==============================================================================
if (isTRUE(params$run_mast)) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SECTION 1: MAST DIFFERENTIAL EXPRESSION\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  mast_log2fc_diagnostic()

  if (!requireNamespace("MAST", quietly = TRUE)) {
    cat("WARNING: MAST package not available. Skipping MAST analysis.\n")
  } else {
    suppressPackageStartupMessages(library(MAST))

    tryCatch({
      cat("Preparing data for MAST...\n")

      expr_mat <- GetAssayData(clustered_obj, layer = "data")

      # Gene filtering with 0.3% detection threshold (changed from 1%)
      gene_detection_rate <- rowMeans(expr_mat > 0)
      genes_to_use <- names(gene_detection_rate[gene_detection_rate > 0.003])
      cat("Gene filtering: detection rate > 0.3% (0.003)\n")
      cat("Genes passing filter:", length(genes_to_use), "out of", nrow(expr_mat), "\n")

      expr_mat_filtered <- expr_mat[genes_to_use, ]

      fData <- data.frame(primerid = rownames(expr_mat_filtered),
                          row.names = rownames(expr_mat_filtered))

      cData <- clustered_obj@meta.data[colnames(expr_mat_filtered), , drop = FALSE]
      cData$wellKey <- rownames(cData)

      # Use nFeature_RNA for cngeneson (RECOMMENDED: stable, standard practice)
      cat("Using nFeature_RNA for cngeneson calculation (recommended)\n")
      cData$cngeneson <- scale(cData$nFeature_RNA)

      sca <- FromMatrix(
        exprsArray = as.matrix(expr_mat_filtered),
        cData = cData,
        fData = fData
      )

      cat("Created SingleCellAssay with", nrow(sca), "genes and", ncol(sca), "cells\n")

      # Check if we can use random effects
      # Criteria: (1) >1 sample per group, (2) <40K cells (otherwise too slow)
      use_random_effects <- FALSE
      random_effects_reason <- ""
      n_cells_total <- ncol(sca)
      
      cat("\n")
      cat(paste(rep("=", 70), collapse = ""), "\n")
      cat("MAST MODEL SELECTION\n")
      cat(paste(rep("=", 70), collapse = ""), "\n")
      cat(sprintf("  Total cells: %d\n", n_cells_total))
      cat(sprintf("  Cell threshold for random effects: %d\n", MAST_MAX_CELLS_FOR_RANDOM_EFFECTS))
      cat(sprintf("  Biological replicates per group: %d\n", N_BIO_REPS))
      
      if ("sample_name" %in% colnames(cData) && N_BIO_REPS >= 2) {
        n_samples <- length(unique(cData$sample_name))
        if (n_samples >= 4) {
          # Check cell count threshold
          if (n_cells_total <= MAST_MAX_CELLS_FOR_RANDOM_EFFECTS) {
            use_random_effects <- TRUE
            random_effects_reason <- "Sufficient replicates and cell count within threshold"
            cat("\n  Decision: RANDOM EFFECTS (recommended)\n")
            cat("  Reason:", random_effects_reason, "\n")
          } else {
            use_random_effects <- FALSE
            random_effects_reason <- sprintf("Cell count (%d) exceeds threshold (%d) - would take 20-30+ hours", 
                                              n_cells_total, MAST_MAX_CELLS_FOR_RANDOM_EFFECTS)
            cat("\n  Decision: FIXED EFFECTS (for computational efficiency)\n")
            cat("  Reason:", random_effects_reason, "\n")
            cat("  Note: Pseudobulk methods (DESeq2, DREAM) will still use proper random effects\n")
          }
        } else {
          random_effects_reason <- "Insufficient samples (<4 total)"
          cat("\n  Decision: FIXED EFFECTS\n")
          cat("  Reason:", random_effects_reason, "\n")
        }
      } else {
        random_effects_reason <- "Insufficient biological replicates (<2 per group)"
        cat("\n  Decision: FIXED EFFECTS\n")
        cat("  Reason:", random_effects_reason, "\n")
      }
      
      cat(paste(rep("=", 70), collapse = ""), "\n\n")

      # ==============================================================================
      # MAST MODEL FITTING WITH PARALLELIZATION
      # ==============================================================================
      
      cat(paste(rep("=", 70), collapse = ""), "\n")
      cat("RUNNING MAST WITH PARALLEL PROCESSING\n")
      cat(paste(rep("=", 70), collapse = ""), "\n")
      cat(sprintf("  Parallel enabled: %s\n", mast_parallel_enabled))
      cat(sprintf("  Parallel workers: %d\n", n_cores_mast))
      cat(sprintf("  Genes to test: %d\n", nrow(sca)))
      cat(sprintf("  Cells: %d\n", ncol(sca)))
      
      # Estimate runtime
      if (use_random_effects) {
        estimated_time_per_gene_seconds <- 3  # Conservative estimate with random effects
        cat("  Model type: Random effects (GLMM via lme4)\n")
      } else {
        estimated_time_per_gene_seconds <- 0.1  # Fixed effects are much faster
        cat("  Model type: Fixed effects (GLM)\n")
      }
      
      estimated_total_minutes <- (nrow(sca) * estimated_time_per_gene_seconds) / (n_cores_mast * 60)
      cat(sprintf("  Estimated runtime: %.1f - %.1f hours\n", 
                  estimated_total_minutes / 60 * 0.5,  # Optimistic
                  estimated_total_minutes / 60 * 1.5))  # Conservative
      cat(paste(rep("=", 70), collapse = ""), "\n\n")
      
      mast_start_time <- Sys.time()

      if (use_random_effects) {
        if ("percent.mt" %in% colnames(cData)) {
          cat("Model: ~ sex + cngeneson + percent.mt + (1|sample_name)\n")
          cat("Parallelization enabled:", mast_parallel_enabled, "\n\n")
          zlm_fit <- zlm(~ sex + cngeneson + percent.mt + (1|sample_name),
                          sca, method = "glmer", ebayes = FALSE,
                          parallel = mast_parallel_enabled)
        } else {
          cat("Model: ~ sex + cngeneson + (1|sample_name)\n")
          cat("Parallelization enabled:", mast_parallel_enabled, "\n\n")
          zlm_fit <- zlm(~ sex + cngeneson + (1|sample_name),
                          sca, method = "glmer", ebayes = FALSE,
                          parallel = mast_parallel_enabled)
        }
      } else {
        if (N_BIO_REPS < 2) {
          cat("\n")
          cat(paste(rep("!", 60), collapse = ""), "\n")
          cat("WARNING: Cannot use random effects with <2 samples per group\n")
          cat("Using fixed effects only - RESULTS ARE EXPLORATORY\n")
          cat("Pseudoreplication bias is NOT controlled (Hafner et al. 2025)\n")
          cat(paste(rep("!", 60), collapse = ""), "\n\n")
        } else if (n_cells_total > MAST_MAX_CELLS_FOR_RANDOM_EFFECTS) {
          cat("\n")
          cat(paste(rep("*", 60), collapse = ""), "\n")
          cat("NOTE: Using fixed effects due to large cell count\n")
          cat("Random effects would take 20-30+ hours with", n_cells_total, "cells\n")
          cat("DESeq2 and DREAM will provide proper random effects analysis\n")
          cat(paste(rep("*", 60), collapse = ""), "\n\n")
        }

        if ("percent.mt" %in% colnames(cData)) {
          cat("Model: ~ sex + cngeneson + percent.mt\n")
          cat("Parallelization enabled:", mast_parallel_enabled, "\n\n")
          zlm_fit <- zlm(~ sex + cngeneson + percent.mt, sca,
                          parallel = mast_parallel_enabled)
        } else {
          cat("Model: ~ sex + cngeneson\n")
          cat("Parallelization enabled:", mast_parallel_enabled, "\n\n")
          zlm_fit <- zlm(~ sex + cngeneson, sca,
                          parallel = mast_parallel_enabled)
        }
      }

      # Report timing
      mast_end_time <- Sys.time()
      mast_elapsed <- difftime(mast_end_time, mast_start_time, units = "mins")
      
      cat("\n")
      cat(paste(rep("=", 70), collapse = ""), "\n")
      cat("MAST MODEL FITTING COMPLETE\n")
      cat(paste(rep("=", 70), collapse = ""), "\n")
      cat(sprintf("  Elapsed time: %.2f minutes (%.2f hours)\n", 
                  as.numeric(mast_elapsed), 
                  as.numeric(mast_elapsed) / 60))
      cat(sprintf("  Genes processed: %d\n", nrow(sca)))
      cat(sprintf("  Average time per gene: %.3f seconds\n", 
                  as.numeric(mast_elapsed) * 60 / nrow(sca)))
      if (mast_parallel_enabled && n_cores_mast > 1) {
        cat(sprintf("  Parallel speedup factor: ~%.1fx (estimated)\n", n_cores_mast * 0.7))
      }
      cat(paste(rep("=", 70), collapse = ""), "\n\n")

      cat("\nRunning likelihood ratio test for contrast:", CONTRAST_NAME, "\n")

      summary_output <- summary(zlm_fit, doLRT = CONTRAST_NAME)
      summaryDt <- summary_output$datatable

      cat("\n--- DEBUG: MAST Summary Structure ---\n")
      cat("Columns in summaryDt:\n")
      print(colnames(summaryDt))
      cat("\nUnique components:\n")
      print(unique(summaryDt$component))
      cat("\nUnique contrasts:\n")
      print(unique(summaryDt$contrast))
      cat("\n")

      # Extract hurdle p-values (component "H")
      hurdle_pvals <- summaryDt[summaryDt$component == "H" & summaryDt$contrast == CONTRAST_NAME, ]
      hurdle_pvals <- hurdle_pvals[, c("primerid", "Pr(>Chisq)")]
      colnames(hurdle_pvals) <- c("gene", "hurdle_pval")

      # Extract logFC using the new function (prefers "logFC" component)
      logfc_data <- extract_mast_logfc(summaryDt, CONTRAST_NAME)

      # Merge hurdle p-values with logFC data
      mast_results <- merge(hurdle_pvals, logfc_data, by = "gene", all = TRUE)
      mast_results$FDR <- p.adjust(mast_results$hurdle_pval, method = "BH")
      mast_results$log2FC_if_ln <- mast_results$logFC / log(2)

      # Add dynamic direction column based on established factor levels
      mast_results <- add_direction_column_dynamic(mast_results, "logFC",
                                                    UP_LABEL, DOWN_LABEL)

      mast_results <- mast_results[order(mast_results$FDR), ]
      mast_results <- mast_results[!is.na(mast_results$gene), ]

      cat("\n>>> MAST RESULTS SUMMARY <<<\n")
      cat("Total genes tested:", nrow(mast_results), "\n")
      cat("Significant (FDR < 0.05):", sum(mast_results$FDR < 0.05, na.rm = TRUE), "\n")
      cat(UP_LABEL, "(FDR < 0.05, logFC > 0):", sum(mast_results$FDR < 0.05 & mast_results$logFC > 0, na.rm = TRUE), "\n")
      cat(DOWN_LABEL, "(FDR < 0.05, logFC < 0):", sum(mast_results$FDR < 0.05 & mast_results$logFC < 0, na.rm = TRUE), "\n")
      cat("\nDirection labels used:\n")
      cat("  Positive logFC ->", UP_LABEL, "\n")
      cat("  Negative logFC ->", DOWN_LABEL, "\n")
      cat("  LogFC source component:", unique(mast_results$logFC_source), "\n")
      if (!use_random_effects && N_BIO_REPS < 2) {
        cat("\n  *** WARNING: No random effects - p-values may be inflated ***\n")
      }

      full_results_path <- file.path(de_mast, "MAST_full_results.csv")
      write.csv(mast_results, full_results_path, row.names = FALSE)
      cat("\nFull results saved to:", full_results_path, "\n")
      print_csv_head5_local(full_results_path, "MAST Full Results")

      sig_results <- mast_results[mast_results$FDR < 0.05 & !is.na(mast_results$FDR), ]
      sig_results_path <- file.path(de_mast, "MAST_significant_FDR05.csv")
      write.csv(sig_results, sig_results_path, row.names = FALSE)
      print_csv_head5_local(sig_results_path, "MAST Significant (FDR<0.05)")

      if (!is.null(params$genes_of_interest) && length(params$genes_of_interest) > 0) {
        goi_results <- mast_results[tolower(mast_results$gene) %in% tolower(params$genes_of_interest), ]
        if (nrow(goi_results) > 0) {
          goi_path <- file.path(de_mast, "MAST_genes_of_interest.csv")
          write.csv(goi_results, goi_path, row.names = FALSE)
          print_csv_head5_local(goi_path, "MAST Genes of Interest")
        }
      }

      all_de_results[["MAST"]] <- mast_results

      cat("\nCreating MAST volcano plot...\n")

      volcano_data <- mast_results[!is.na(mast_results$logFC) & !is.na(mast_results$FDR), ]
      volcano_data$neg_log10_fdr <- -log10(volcano_data$FDR + 1e-300)
      volcano_data$significant <- volcano_data$FDR < 0.05 & abs(volcano_data$logFC) > 0.25
      volcano_data$neg_log10_fdr <- pmin(volcano_data$neg_log10_fdr, 50)

      p_volcano <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_fdr)) +
        geom_point(aes(color = significant), alpha = 0.5, size = 1) +
        scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "blue") +
        theme_minimal() +
        labs(title = "MAST Differential Expression",
             subtitle = paste0("Significant genes (FDR<0.05, |logFC|>0.25): ", sum(volcano_data$significant),
                              "\nPositive logFC = ", UP_LABEL, " | Negative logFC = ", DOWN_LABEL),
             x = "Log Fold Change (natural log scale)",
             y = "-log10(FDR)") +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5, face = "bold"))

      if (has_ggrepel) {
        top_genes <- volcano_data %>% filter(significant) %>% arrange(FDR) %>% head(10)
        if (nrow(top_genes) > 0) {
          p_volcano <- p_volcano + geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 15)
        }
      }

      save_plot_multi(p_volcano, "volcano_MAST", output_dir = de_plots, width = 10, height = 8)

      de_summary <- rbind(de_summary, data.frame(
        method = "MAST", genes_tested = nrow(mast_results),
        significant_005 = sum(mast_results$FDR < 0.05, na.rm = TRUE),
        higher_in_comparison = sum(mast_results$FDR < 0.05 & mast_results$logFC > 0, na.rm = TRUE),
        higher_in_reference = sum(mast_results$FDR < 0.05 & mast_results$logFC < 0, na.rm = TRUE),
        comparison_level = COMP_LEVEL,
        reference_level = REF_LEVEL,
        random_effects = use_random_effects
      ))

      cat("\n[SUCCESS] MAST analysis completed\n")

    }, error = function(e) {
      cat("\n[FAILED] MAST analysis failed:\n")
      cat("Error:", conditionMessage(e), "\n")
    })
  }
} else {
  cat("\nMAST analysis skipped (params$run_mast = FALSE)\n")
}

# ==============================================================================
# Clean up MAST parallel backend
# ==============================================================================
if (exists("mast_parallel_enabled") && mast_parallel_enabled && !is.null(mast_parallel_param)) {
  tryCatch({
    cat("Cleaning up MAST parallel workers...\n")
    bpstop(mast_parallel_param)
    cat("Parallel cleanup complete.\n\n")
  }, error = function(e) {
    cat("Note: Parallel cleanup message -", conditionMessage(e), "\n")
  })
}

# ==============================================================================
# SECTION 2: DREAM ANALYSIS (Hafner et al. 2025 RECOMMENDED)
# ==============================================================================
if (isTRUE(params$run_dream) || !exists("params$run_dream")) {
  # Default to TRUE if not specified - this is the recommended method
  if (!exists("params") || is.null(params$run_dream)) {
    run_dream_analysis <- TRUE
  } else {
    run_dream_analysis <- params$run_dream
  }

  if (run_dream_analysis) {
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("SECTION 2: DREAM DIFFERENTIAL EXPRESSION (RECOMMENDED)\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("
Per Hafner et al. (2025) Briefings in Bioinformatics:
- DREAM is RECOMMENDED for nested/hierarchical data
- 'DREAM provides a good balance between accuracy and runtime'
- Uses linear mixed effects models via limma/voom framework
- Properly handles variance estimation with limited replicates
")
    cat(paste(rep("=", 70), collapse = ""), "\n\n")

    # Check for required packages
    has_variance_partition <- requireNamespace("variancePartition", quietly = TRUE)
    has_edger <- requireNamespace("edgeR", quietly = TRUE)
    has_limma <- requireNamespace("limma", quietly = TRUE)

    if (!has_variance_partition) {
      cat("variancePartition not installed. Attempting installation...\n")
      tryCatch({
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install("variancePartition", ask = FALSE, update = FALSE)
        has_variance_partition <- TRUE
      }, error = function(e) {
        cat("Could not install variancePartition:", conditionMessage(e), "\n")
      })
    }

    if (!has_variance_partition || !has_edger || !has_limma) {
      cat("WARNING: Required packages not available. Skipping DREAM analysis.\n")
      cat("Required: variancePartition, edgeR, limma\n")
    } else {
      suppressPackageStartupMessages({
        library(variancePartition)
        library(edgeR)
        library(limma)
      })

      tryCatch({
        cat("Creating pseudobulk for DREAM...\n")

        # Use smart pseudobulk creation
        pb_data <- create_pseudobulk_smart(clustered_obj, group_var = "sex",
                                            sample_var = "sample_name", n_pseudo_reps = 3)

        pb_counts <- pb_data$counts
        pb_info <- pb_data$sample_info

        cat("\nPseudubulk method used:", pb_data$method, "\n")
        cat("Samples per group:\n")
        print(table(pb_info$group))

        # Create DGEList
        dge <- DGEList(counts = pb_counts)

        # Gene filtering
        keep <- filterByExpr(dge, group = pb_info$group, min.count = 5)
        dge <- dge[keep, , keep.lib.sizes = FALSE]
        cat("Genes after filtering:", nrow(dge), "\n")

        # Normalize
        dge <- calcNormFactors(dge)

        # Setup factors
        pb_info$group <- factor(pb_info$group, levels = c(REF_LEVEL, COMP_LEVEL))

        # DREAM model
        # With limited replicates, we use a simpler formula
        if (pb_data$has_true_replicates && "original_sample" %in% colnames(pb_info)) {
          # Can include sample as random effect if we have replicates
          cat("Using model with variance weighting (DREAM)\n")
        }

        # Simple formula for comparison
        formula <- ~ group

        # Calculate voom weights with DREAM
        cat("Calculating voom weights...\n")
        vobjDream <- voomWithDreamWeights(
          dge,
          formula,
          pb_info,
          plot = FALSE
        )

        # Fit the model
        cat("Fitting DREAM model...\n")
        fit <- dream(vobjDream, formula, pb_info)
        fit <- eBayes(fit)

        # Extract results
        coef_name <- paste0("group", COMP_LEVEL)

        # Check available coefficients
        cat("Available coefficients:", paste(colnames(coef(fit)), collapse = ", "), "\n")

        dream_results <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
        dream_results$gene <- rownames(dream_results)
        rownames(dream_results) <- NULL

        # Add direction
        dream_results <- add_direction_column_dynamic(
          dream_results, "logFC", UP_LABEL, DOWN_LABEL
        )

        # Reorder columns
        dream_results <- dream_results[, c("gene", "logFC", "AveExpr", "t", "P.Value",
                                            "adj.P.Val", "B", "direction", "direction_numeric")]

        cat("\n>>> DREAM RESULTS SUMMARY <<<\n")
        cat("Total genes tested:", nrow(dream_results), "\n")
        cat("Significant (adj.P.Val < 0.05):", sum(dream_results$adj.P.Val < 0.05, na.rm = TRUE), "\n")
        cat(UP_LABEL, "(adj.P.Val < 0.05, logFC > 0):",
            sum(dream_results$adj.P.Val < 0.05 & dream_results$logFC > 0, na.rm = TRUE), "\n")
        cat(DOWN_LABEL, "(adj.P.Val < 0.05, logFC < 0):",
            sum(dream_results$adj.P.Val < 0.05 & dream_results$logFC < 0, na.rm = TRUE), "\n")

        # Save results
        full_path <- file.path(de_dream, "DREAM_full_results.csv")
        write.csv(dream_results, full_path, row.names = FALSE)
        cat("\nFull results saved to:", full_path, "\n")
        print_csv_head5_local(full_path, "DREAM Full Results")

        sig_dream <- dream_results[dream_results$adj.P.Val < 0.05 & !is.na(dream_results$adj.P.Val), ]
        sig_path <- file.path(de_dream, "DREAM_significant_FDR05.csv")
        write.csv(sig_dream, sig_path, row.names = FALSE)

        all_de_results[["DREAM"]] <- dream_results

        # Volcano plot
        cat("\nCreating DREAM volcano plot...\n")

        volcano_dream <- dream_results[!is.na(dream_results$logFC) & !is.na(dream_results$adj.P.Val), ]
        volcano_dream$neg_log10_p <- pmin(-log10(volcano_dream$adj.P.Val + 1e-300), 50)
        volcano_dream$significant <- volcano_dream$adj.P.Val < 0.05 & abs(volcano_dream$logFC) > 0.5

        p_volcano_dream <- ggplot(volcano_dream, aes(x = logFC, y = neg_log10_p)) +
          geom_point(aes(color = significant), alpha = 0.5, size = 1) +
          scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "steelblue")) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
          geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
          theme_minimal() +
          labs(title = "DREAM Differential Expression (Recommended)",
               subtitle = paste0("Significant genes (adj.P.Val<0.05, |logFC|>0.5): ", sum(volcano_dream$significant),
                                "\nPositive logFC = ", UP_LABEL, " | Negative logFC = ", DOWN_LABEL),
               x = "Log2 Fold Change",
               y = "-log10(adjusted P-value)") +
          theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

        if (has_ggrepel) {
          top_genes <- volcano_dream %>% filter(significant) %>% arrange(adj.P.Val) %>% head(10)
          if (nrow(top_genes) > 0) {
            p_volcano_dream <- p_volcano_dream +
              geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 15)
          }
        }

        save_plot_multi(p_volcano_dream, "volcano_DREAM", output_dir = de_plots, width = 10, height = 8)

        de_summary <- rbind(de_summary, data.frame(
          method = "DREAM", genes_tested = nrow(dream_results),
          significant_005 = sum(dream_results$adj.P.Val < 0.05, na.rm = TRUE),
          higher_in_comparison = sum(dream_results$adj.P.Val < 0.05 & dream_results$logFC > 0, na.rm = TRUE),
          higher_in_reference = sum(dream_results$adj.P.Val < 0.05 & dream_results$logFC < 0, na.rm = TRUE),
          comparison_level = COMP_LEVEL,
          reference_level = REF_LEVEL,
          random_effects = FALSE
        ))

        cat("\n[SUCCESS] DREAM analysis completed\n")

      }, error = function(e) {
        cat("\n[FAILED] DREAM analysis failed:", conditionMessage(e), "\n")
        cat("Stack trace:\n")
        traceback()
      })
    }
  }
} else {
  cat("\nDREAM analysis skipped (params$run_dream = FALSE)\n")
}

# ==============================================================================
# SECTION 3: edgeR PSEUDOBULK ANALYSIS
# ==============================================================================
if (isTRUE(params$run_pseudobulk_edger)) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SECTION 3: edgeR PSEUDOBULK DIFFERENTIAL EXPRESSION\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  if (!requireNamespace("edgeR", quietly = TRUE) || !requireNamespace("limma", quietly = TRUE)) {
    cat("WARNING: edgeR or limma not available. Skipping.\n")
  } else {
    suppressPackageStartupMessages({ library(edgeR); library(limma) })

    tryCatch({
      cat("Creating pseudobulk for edgeR...\n")

      # Use smart pseudobulk creation
      pb_data <- create_pseudobulk_smart(clustered_obj, group_var = "sex",
                                          sample_var = "sample_name", n_pseudo_reps = 3)

      pb_counts <- pb_data$counts
      pb_info <- pb_data$sample_info

      cat("Pseudobulk method used:", pb_data$method, "\n")
      cat("Samples per group:\n")
      print(table(pb_info$group))

      dge <- DGEList(counts = pb_counts, group = pb_info$group)
      keep <- filterByExpr(dge, min.count = 10, min.total.count = 15)
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      cat("Genes after filtering:", nrow(dge), "\n")

      dge <- calcNormFactors(dge)
      design <- model.matrix(~ 0 + group, data = pb_info)
      colnames(design) <- gsub("group", "", colnames(design))

      v <- voom(dge, design, plot = FALSE)
      fit <- lmFit(v, design)

      # Create contrast dynamically
      contrast_formula <- paste(COMP_LEVEL, "-", REF_LEVEL)
      contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

      fit2 <- contrasts.fit(fit, contrast_matrix)
      fit2 <- eBayes(fit2)

      edger_results <- topTable(fit2, number = Inf, sort.by = "P")
      edger_results$gene <- rownames(edger_results)
      rownames(edger_results) <- NULL

      # Add dynamic direction column
      edger_results <- add_direction_column_dynamic(edger_results, "logFC",
                                                     UP_LABEL, DOWN_LABEL)

      edger_results <- edger_results[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "direction", "direction_numeric")]

      cat("\n>>> edgeR RESULTS SUMMARY <<<\n")
      cat("Total genes tested:", nrow(edger_results), "\n")
      cat("Significant (adj.P.Val < 0.05):", sum(edger_results$adj.P.Val < 0.05, na.rm = TRUE), "\n")
      cat(UP_LABEL, "(adj.P.Val < 0.05, logFC > 0):", sum(edger_results$adj.P.Val < 0.05 & edger_results$logFC > 0, na.rm = TRUE), "\n")
      cat(DOWN_LABEL, "(adj.P.Val < 0.05, logFC < 0):", sum(edger_results$adj.P.Val < 0.05 & edger_results$logFC < 0, na.rm = TRUE), "\n")
      if (!pb_data$has_true_replicates) {
        cat("\n  *** WARNING: Using pseudo-replicates - p-values may be inflated ***\n")
      }

      full_path <- file.path(de_edger, "edgeR_full_results.csv")
      write.csv(edger_results, full_path, row.names = FALSE)
      print_csv_head5_local(full_path, "edgeR Full Results")

      sig_edger <- edger_results[edger_results$adj.P.Val < 0.05 & !is.na(edger_results$adj.P.Val), ]
      sig_path <- file.path(de_edger, "edgeR_significant_FDR05.csv")
      write.csv(sig_edger, sig_path, row.names = FALSE)

      all_de_results[["edgeR"]] <- edger_results

      volcano_edger <- edger_results[!is.na(edger_results$logFC) & !is.na(edger_results$adj.P.Val), ]
      volcano_edger$neg_log10_p <- pmin(-log10(volcano_edger$adj.P.Val + 1e-300), 50)
      volcano_edger$significant <- volcano_edger$adj.P.Val < 0.05 & abs(volcano_edger$logFC) > 0.5

      p_volcano_edger <- ggplot(volcano_edger, aes(x = logFC, y = neg_log10_p)) +
        geom_point(aes(color = significant), alpha = 0.5, size = 1) +
        scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "darkgreen")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
        theme_minimal() +
        labs(title = "edgeR Pseudobulk DE",
             subtitle = paste0("Significant genes (adj.P.Val<0.05, |logFC|>0.5): ", sum(volcano_edger$significant),
                              "\nPositive logFC = ", UP_LABEL, " | Negative logFC = ", DOWN_LABEL),
             x = "Log2 Fold Change",
             y = "-log10(adjusted P-value)") +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

      if (has_ggrepel) {
        top_genes <- volcano_edger %>% filter(significant) %>% arrange(adj.P.Val) %>% head(10)
        if (nrow(top_genes) > 0) {
          p_volcano_edger <- p_volcano_edger +
            geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 15)
        }
      }

      save_plot_multi(p_volcano_edger, "volcano_edgeR", output_dir = de_plots, width = 10, height = 8)

      de_summary <- rbind(de_summary, data.frame(
        method = "edgeR", genes_tested = nrow(edger_results),
        significant_005 = sum(edger_results$adj.P.Val < 0.05, na.rm = TRUE),
        higher_in_comparison = sum(edger_results$adj.P.Val < 0.05 & edger_results$logFC > 0, na.rm = TRUE),
        higher_in_reference = sum(edger_results$adj.P.Val < 0.05 & edger_results$logFC < 0, na.rm = TRUE),
        comparison_level = COMP_LEVEL,
        reference_level = REF_LEVEL,
        random_effects = FALSE
      ))

      cat("\n[SUCCESS] edgeR analysis completed\n")

    }, error = function(e) {
      cat("\n[FAILED] edgeR analysis failed:", conditionMessage(e), "\n")
    })
  }
} else {
  cat("\nedgeR analysis skipped (params$run_pseudobulk_edger = FALSE)\n")
}

# ==============================================================================
# SECTION 4: DESeq2 PSEUDOBULK ANALYSIS
# ==============================================================================
if (isTRUE(params$run_pseudobulk_deseq2)) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SECTION 4: DESeq2 PSEUDOBULK DIFFERENTIAL EXPRESSION\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("
Per Hafner et al. (2025):
- DESeq2 achieved BEST performance on single dataset scenario (AUPRC=0.93)
- Recommended for non-nested data with biological replicates
- Uses negative binomial model with empirical Bayes shrinkage
")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    cat("WARNING: DESeq2 not available. Skipping.\n")
  } else {
    suppressPackageStartupMessages(library(DESeq2))

    tryCatch({
      # Use smart pseudobulk creation
      pb_data <- create_pseudobulk_smart(clustered_obj, group_var = "sex",
                                          sample_var = "sample_name", n_pseudo_reps = 3)

      pb_counts <- pb_data$counts
      pb_info <- pb_data$sample_info

      cat("Pseudobulk method used:", pb_data$method, "\n")

      pb_counts_int <- round(pb_counts)
      mode(pb_counts_int) <- "integer"

      keep <- rowSums(pb_counts_int >= 10) >= 2
      pb_counts_filtered <- pb_counts_int[keep, ]
      cat("Genes after filtering:", nrow(pb_counts_filtered), "\n")

      pb_info$group <- factor(pb_info$group, levels = c(REF_LEVEL, COMP_LEVEL))

      dds <- DESeqDataSetFromMatrix(countData = pb_counts_filtered, colData = pb_info, design = ~ group)
      dds <- DESeq(dds)

      deseq2_res <- results(dds, contrast = c("group", COMP_LEVEL, REF_LEVEL))
      deseq2_results <- as.data.frame(deseq2_res)
      deseq2_results$gene <- rownames(deseq2_results)
      rownames(deseq2_results) <- NULL

      # Add dynamic direction column
      deseq2_results <- add_direction_column_dynamic(deseq2_results, "log2FoldChange",
                                                      UP_LABEL, DOWN_LABEL)

      deseq2_results <- deseq2_results[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "direction", "direction_numeric")]
      deseq2_results <- deseq2_results[order(deseq2_results$padj), ]

      cat("\n>>> DESeq2 RESULTS SUMMARY <<<\n")
      cat("Total genes tested:", nrow(deseq2_results), "\n")
      cat("Significant (padj < 0.05):", sum(deseq2_results$padj < 0.05, na.rm = TRUE), "\n")
      cat(UP_LABEL, "(padj < 0.05, log2FC > 0):", sum(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange > 0, na.rm = TRUE), "\n")
      cat(DOWN_LABEL, "(padj < 0.05, log2FC < 0):", sum(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange < 0, na.rm = TRUE), "\n")
      if (!pb_data$has_true_replicates) {
        cat("\n  *** WARNING: Using pseudo-replicates - p-values may be inflated ***\n")
      }

      full_path <- file.path(de_deseq2, "DESeq2_full_results.csv")
      write.csv(deseq2_results, full_path, row.names = FALSE)
      print_csv_head5_local(full_path, "DESeq2 Full Results")

      sig_deseq2 <- deseq2_results[deseq2_results$padj < 0.05 & !is.na(deseq2_results$padj), ]
      sig_path <- file.path(de_deseq2, "DESeq2_significant_FDR05.csv")
      write.csv(sig_deseq2, sig_path, row.names = FALSE)

      all_de_results[["DESeq2"]] <- deseq2_results

      volcano_deseq2 <- deseq2_results[!is.na(deseq2_results$log2FoldChange) & !is.na(deseq2_results$padj), ]
      volcano_deseq2$neg_log10_p <- pmin(-log10(volcano_deseq2$padj + 1e-300), 50)
      volcano_deseq2$significant <- volcano_deseq2$padj < 0.05 & abs(volcano_deseq2$log2FoldChange) > 0.5

      p_volcano_deseq2 <- ggplot(volcano_deseq2, aes(x = log2FoldChange, y = neg_log10_p)) +
        geom_point(aes(color = significant), alpha = 0.5, size = 1) +
        scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "purple")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
        theme_minimal() +
        labs(title = "DESeq2 Pseudobulk DE",
             subtitle = paste0("Significant genes (padj<0.05, |log2FC|>0.5): ", sum(volcano_deseq2$significant),
                              "\nPositive log2FC = ", UP_LABEL, " | Negative log2FC = ", DOWN_LABEL),
             x = "Log2 Fold Change",
             y = "-log10(adjusted P-value)") +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

      if (has_ggrepel) {
        top_genes <- volcano_deseq2 %>% filter(significant) %>% arrange(padj) %>% head(10)
        if (nrow(top_genes) > 0) {
          p_volcano_deseq2 <- p_volcano_deseq2 +
            geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 15)
        }
      }

      save_plot_multi(p_volcano_deseq2, "volcano_DESeq2", output_dir = de_plots, width = 10, height = 8)

      de_summary <- rbind(de_summary, data.frame(
        method = "DESeq2", genes_tested = nrow(deseq2_results),
        significant_005 = sum(deseq2_results$padj < 0.05, na.rm = TRUE),
        higher_in_comparison = sum(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange > 0, na.rm = TRUE),
        higher_in_reference = sum(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange < 0, na.rm = TRUE),
        comparison_level = COMP_LEVEL,
        reference_level = REF_LEVEL,
        random_effects = FALSE
      ))

      cat("\n[SUCCESS] DESeq2 analysis completed\n")

    }, error = function(e) {
      cat("\n[FAILED] DESeq2 analysis failed:", conditionMessage(e), "\n")
    })
  }
} else {
  cat("\nDESeq2 analysis skipped (params$run_pseudobulk_deseq2 = FALSE)\n")
}

# ==============================================================================
# SECTION 5: PERMUTATION TEST (Hafner et al. 2025)
# ==============================================================================
run_permutation <- TRUE
if (exists("params") && !is.null(params$run_permutation)) {
  run_permutation <- params$run_permutation
}

if (run_permutation) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SECTION 5: PERMUTATION TEST DIFFERENTIAL EXPRESSION\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  tryCatch({
    # Create pseudobulk for permutation test
    pb_data <- create_pseudobulk_smart(clustered_obj, group_var = "sex",
                                        sample_var = "sample_name", n_pseudo_reps = 3)

    # Run permutation test
    perm_results <- run_permutation_test(
      pb_counts = pb_data$counts,
      sample_info = pb_data$sample_info,
      n_permutations = 10000,
      seed = 42
    )

    # Add direction column
    perm_results <- add_direction_column_dynamic(perm_results, "log2FC",
                                                  UP_LABEL, DOWN_LABEL)

    cat("\n>>> PERMUTATION TEST RESULTS SUMMARY <<<\n")
    cat("Total genes tested:", nrow(perm_results), "\n")
    cat("Significant (padj < 0.05):", sum(perm_results$padj < 0.05, na.rm = TRUE), "\n")
    cat(UP_LABEL, "(padj < 0.05, log2FC > 0):",
        sum(perm_results$padj < 0.05 & perm_results$log2FC > 0, na.rm = TRUE), "\n")
    cat(DOWN_LABEL, "(padj < 0.05, log2FC < 0):",
        sum(perm_results$padj < 0.05 & perm_results$log2FC < 0, na.rm = TRUE), "\n")

    # Save results
    full_path <- file.path(de_permutation, "Permutation_full_results.csv")
    write.csv(perm_results, full_path, row.names = FALSE)
    print_csv_head5_local(full_path, "Permutation Test Full Results")

    sig_perm <- perm_results[perm_results$padj < 0.05 & !is.na(perm_results$padj), ]
    sig_path <- file.path(de_permutation, "Permutation_significant_FDR05.csv")
    write.csv(sig_perm, sig_path, row.names = FALSE)

    all_de_results[["Permutation"]] <- perm_results

    # Volcano plot
    volcano_perm <- perm_results[!is.na(perm_results$log2FC) & !is.na(perm_results$padj), ]
    volcano_perm$neg_log10_p <- pmin(-log10(volcano_perm$padj + 1e-300), 50)
    volcano_perm$significant <- volcano_perm$padj < 0.05 & abs(volcano_perm$log2FC) > 0.5

    p_volcano_perm <- ggplot(volcano_perm, aes(x = log2FC, y = neg_log10_p)) +
      geom_point(aes(color = significant), alpha = 0.5, size = 1) +
      scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "orange")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
      theme_minimal() +
      labs(title = "Permutation Test DE",
           subtitle = paste0("Significant genes (padj<0.05, |log2FC|>0.5): ", sum(volcano_perm$significant),
                            "\nPositive log2FC = ", UP_LABEL, " | Negative log2FC = ", DOWN_LABEL),
           x = "Log2 Fold Change",
           y = "-log10(adjusted P-value)") +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

    if (has_ggrepel) {
      top_genes <- volcano_perm %>% filter(significant) %>% arrange(padj) %>% head(10)
      if (nrow(top_genes) > 0) {
        p_volcano_perm <- p_volcano_perm +
          geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 15)
      }
    }

    save_plot_multi(p_volcano_perm, "volcano_Permutation", output_dir = de_plots, width = 10, height = 8)

    de_summary <- rbind(de_summary, data.frame(
      method = "Permutation", genes_tested = nrow(perm_results),
      significant_005 = sum(perm_results$padj < 0.05, na.rm = TRUE),
      higher_in_comparison = sum(perm_results$padj < 0.05 & perm_results$log2FC > 0, na.rm = TRUE),
      higher_in_reference = sum(perm_results$padj < 0.05 & perm_results$log2FC < 0, na.rm = TRUE),
      comparison_level = COMP_LEVEL,
      reference_level = REF_LEVEL,
      random_effects = NA
    ))

    cat("\n[SUCCESS] Permutation test completed\n")

  }, error = function(e) {
    cat("\n[FAILED] Permutation test failed:", conditionMessage(e), "\n")
  })
} else {
  cat("\nPermutation test skipped\n")
}

# ==============================================================================
# SECTION 6: NEGATIVE CONTROL (Hafner et al. 2025)
# ==============================================================================
run_negative_ctrl <- TRUE
if (exists("params") && !is.null(params$run_negative_control)) {
  run_negative_ctrl <- params$run_negative_control
}

if (run_negative_ctrl) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SECTION 6: NEGATIVE CONTROL ANALYSIS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  tryCatch({
    # Run negative control with DESeq2 (fastest method)
    neg_ctrl_results <- run_negative_control(
      seurat_obj = clustered_obj,
      comparison_var = "sex",
      sample_var = "sample_name",
      method = "DESeq2",
      n_iterations = 3,
      seed = 42
    )

    # Save results
    neg_ctrl_path <- file.path(de_negative_control, "negative_control_summary.csv")
    write.csv(neg_ctrl_results$summary, neg_ctrl_path, row.names = FALSE)
    cat("\nNegative control results saved to:", neg_ctrl_path, "\n")

    cat("\n[SUCCESS] Negative control analysis completed\n")

  }, error = function(e) {
    cat("\n[FAILED] Negative control analysis failed:", conditionMessage(e), "\n")
  })
} else {
  cat("\nNegative control analysis skipped\n")
}

# ==============================================================================
# SECTION 7: CONCORDANCE ANALYSIS
# ==============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SECTION 7: CONCORDANCE ANALYSIS ACROSS METHODS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

if (length(all_de_results) >= 2) {
  tryCatch({
    # Calculate concordance
    concordance_result <- calculate_concordance(
      all_results = all_de_results,
      fdr_threshold = 0.05,
      lfc_threshold = 0,
      up_label = UP_LABEL,
      down_label = DOWN_LABEL
    )

    # Save concordance results
    if (!is.null(concordance_result)) {
      save_concordance_results(concordance_result, de_base, UP_LABEL, DOWN_LABEL)

      # Create upset plot or venn diagram if available
      if (requireNamespace("UpSetR", quietly = TRUE)) {
        cat("\nCreating UpSet plot of method overlap...\n")

        # Prepare list for UpSet
        sig_gene_list <- lapply(concordance_result$method_sig_genes, function(x) x$gene)

        if (all(sapply(sig_gene_list, length) > 0)) {
          tryCatch({
            library(UpSetR)

            pdf(file.path(de_plots, "upset_method_overlap.pdf"), width = 10, height = 6)
            print(upset(fromList(sig_gene_list),
                        order.by = "freq",
                        nsets = length(sig_gene_list),
                        sets.bar.color = "steelblue"))
            dev.off()

            png(file.path(de_plots, "upset_method_overlap.png"), width = 10, height = 6, units = "in", res = 150)
            print(upset(fromList(sig_gene_list),
                        order.by = "freq",
                        nsets = length(sig_gene_list),
                        sets.bar.color = "steelblue"))
            dev.off()

            cat("UpSet plot saved\n")
          }, error = function(e) {
            cat("Could not create UpSet plot:", conditionMessage(e), "\n")
          })
        }
      }
    }

    cat("\n[SUCCESS] Concordance analysis completed\n")

  }, error = function(e) {
    cat("\n[FAILED] Concordance analysis failed:", conditionMessage(e), "\n")
  })
} else {
  cat("Insufficient methods for concordance analysis (need at least 2).\n")
  concordance_result <- NULL
}

# ==============================================================================
# SECTION 8: COMBINED RESULTS AND SUMMARY
# ==============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SECTION 8: COMBINED RESULTS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

if (length(all_de_results) > 0) {
  cat(">>> DIFFERENTIAL EXPRESSION SUMMARY <<<\n")
  cat("\nComparison setup:\n")
  cat("  Comparison variable: sex\n")
  cat("  Reference level:", REF_LEVEL, "\n")
  cat("  Comparison level:", COMP_LEVEL, "\n")
  cat("  Positive logFC means:", UP_LABEL, "\n")
  cat("  Negative logFC means:", DOWN_LABEL, "\n")
  cat("  Biological replicates:", N_BIO_REPS, "\n\n")

  print(de_summary)
  write.csv(de_summary, file.path(de_base, "DE_methods_summary.csv"), row.names = FALSE)

  # Method agreement summary
  if (length(all_de_results) >= 2) {
    sig_genes_list <- list()
    if ("MAST" %in% names(all_de_results)) sig_genes_list[["MAST"]] <- all_de_results[["MAST"]]$gene[all_de_results[["MAST"]]$FDR < 0.05]
    if ("DREAM" %in% names(all_de_results)) sig_genes_list[["DREAM"]] <- all_de_results[["DREAM"]]$gene[all_de_results[["DREAM"]]$adj.P.Val < 0.05]
    if ("edgeR" %in% names(all_de_results)) sig_genes_list[["edgeR"]] <- all_de_results[["edgeR"]]$gene[all_de_results[["edgeR"]]$adj.P.Val < 0.05]
    if ("DESeq2" %in% names(all_de_results)) sig_genes_list[["DESeq2"]] <- all_de_results[["DESeq2"]]$gene[all_de_results[["DESeq2"]]$padj < 0.05]
    if ("Permutation" %in% names(all_de_results)) sig_genes_list[["Permutation"]] <- all_de_results[["Permutation"]]$gene[all_de_results[["Permutation"]]$padj < 0.05]

    sig_genes_list <- lapply(sig_genes_list, function(x) x[!is.na(x)])

    if (length(sig_genes_list) >= 2) {
      common_genes <- Reduce(intersect, sig_genes_list)
      cat("\nGenes significant in ALL methods:", length(common_genes), "\n")
      if (length(common_genes) > 0) {
        write.csv(data.frame(gene = common_genes),
                  file.path(de_base, "genes_significant_all_methods.csv"), row.names = FALSE)
      }

      # Pairwise overlaps
      cat("\nPairwise method overlap (number of shared significant genes):\n")
      method_names <- names(sig_genes_list)
      overlap_matrix <- matrix(0, nrow = length(method_names), ncol = length(method_names))
      rownames(overlap_matrix) <- method_names
      colnames(overlap_matrix) <- method_names

      for (i in seq_along(method_names)) {
        for (j in seq_along(method_names)) {
          overlap_matrix[i, j] <- length(intersect(sig_genes_list[[i]], sig_genes_list[[j]]))
        }
      }
      print(overlap_matrix)

      write.csv(overlap_matrix, file.path(de_base, "method_overlap_matrix.csv"))
    }
  }

  # Add warning about interpretation if limited replicates
  if (N_BIO_REPS < 2) {
    cat("\n")
    cat(paste(rep("*", 70), collapse = ""), "\n")
    cat("INTERPRETATION GUIDANCE (Limited Replicates)\n")
    cat(paste(rep("*", 70), collapse = ""), "\n")
    cat("
With only 1 sample per condition, all statistical results should be
interpreted with EXTREME CAUTION. Based on Hafner et al. (2025):
1. P-values are likely INFLATED (too many false positives expected)
2. Focus on CONCORDANT genes (significant in multiple methods, same direction)
3. Consider EFFECT SIZES (logFC) more than p-values
4. These results are HYPOTHESIS-GENERATING, not confirmatory
5. Validation with additional biological replicates is ESSENTIAL
Recommended: Use the concordance analysis results, focusing on genes
that are significant in the same direction across 3+ methods.
")
    cat(paste(rep("*", 70), collapse = ""), "\n")
  }

} else {
  cat("No DE results generated.\n")
}

# Save all results
de_file <- file.path(output_dirs$objects, "08_de_data.RData")
save(all_de_results, de_summary, concordance_result,
     REF_LEVEL, COMP_LEVEL, UP_LABEL, DOWN_LABEL, N_BIO_REPS,
     file = de_file)
cat("\nDE data saved to:", de_file, "\n")

# Create detailed README
readme_content <- paste0(
  "Differential Expression Analysis\n",
  "================================\n\n",
  "Comparison: ", COMP_LEVEL, " vs ", REF_LEVEL, " (reference)\n",
  "Biological replicates per group: ", N_BIO_REPS, "\n\n",
  "Direction Interpretation:\n",
  "  - Positive logFC: ", UP_LABEL, "\n",
  "  - Negative logFC: ", DOWN_LABEL, "\n\n",
  "Methods: MAST, DREAM, edgeR, DESeq2, Permutation Test\n",
  "Threshold: FDR < 0.05\n\n",
  "MAST Settings:\n",
  "  - Gene filter: detection rate > 0.3% (0.003)\n",
  "  - LogFC source: 'logFC' component (combined estimate)\n",
  "  - cngeneson: scaled nFeature_RNA\n",
  "  - log2FC_if_ln column: logFC / log(2) for comparison\n",
  "  - Parallelization: ", ifelse(exists("mast_parallel_enabled") && mast_parallel_enabled, 
                                   paste0("Enabled (", n_cores_mast, " cores)"), "Disabled"), "\n",
  "  - Model selection: ", ifelse(exists("use_random_effects") && use_random_effects, 
                                   "Random effects (GLMM)", "Fixed effects (GLM)"), "\n",
  "  - Cell count threshold for random effects: ", MAST_MAX_CELLS_FOR_RANDOM_EFFECTS, "\n\n",
  "DREAM (variancePartition) Settings:\n",
  "  - Recommended method per Hafner et al. (2025)\n",
  "  - Uses voom weights with limma framework\n\n",
  "Permutation Test Settings:\n",
  "  - 10,000 permutations\n",
  "  - Non-parametric, robust to distributional assumptions\n\n",
  "Concordance Analysis:\n",
  "  - See 'concordance/' subdirectory for genes significant across methods\n",
  "  - Focus on concordant genes when biological replicates are limited\n\n",
  "Note: MAST uses natural log scale. Use log2FC_if_ln for comparison with other methods.\n\n",
  "Reference:\n",
  "Hafner L, et al. (2025) Single-cell differential expression analysis between\n",
  "conditions within nested settings. Briefings in Bioinformatics, 26(4), bbaf397.\n",
  "https://doi.org/10.1093/bib/bbaf397\n"
)

write_readme(de_base, "Differential Expression Analysis", readme_content,
             list("DE_methods_summary.csv" = "Summary of all methods",
                  "method_overlap_matrix.csv" = "Pairwise method overlap",
                  "genes_significant_all_methods.csv" = "Genes significant in all methods",
                  "concordance/" = "Concordance analysis results"))

cat("\n>>> MODULE 08 COMPLETE <<<\n")
cat("Methods run:", paste(names(all_de_results), collapse = ", "), "\n")
cat("\nDirection labels used throughout:\n")
cat("  ", UP_LABEL, "(positive logFC)\n")
cat("  ", DOWN_LABEL, "(negative logFC)\n")

if (N_BIO_REPS < 2) {
  cat("\n")
  cat(paste(rep("!", 60), collapse = ""), "\n")
  cat("REMINDER: Results are EXPLORATORY due to limited replicates\n")
  cat("Focus on CONCORDANT genes across multiple methods\n")
  cat(paste(rep("!", 60), collapse = ""), "\n")
}

cat("\nKey output files:\n")
cat("  - DE_methods_summary.csv: Overview of all method results\n")
cat("  - concordance/: Genes significant in same direction across methods\n")
cat("  - NegativeControl/: FPR assessment with permuted labels\n")
cat("\n")
