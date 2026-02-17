#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 10: FINAL SUMMARY (MULTI-SAMPLE PIPELINE)
# ==============================================================================
#
# This module creates the final summary including:
# - Analysis summary table
# - Per-sample statistics
# - Method selection summary
# - scICE subclustering summary
# - Output directory structure
# - Session information
# - Final README
#
# INPUT: All previous module outputs
# OUTPUT: Summary tables and README
#
# UPDATES:
# - 2026-01-13: Fixed verification loop to properly convert Python KeysView to R
# - 2026-01-13: CRITICAL FIX - H5AD conversion now passes all data at AnnData
#               creation time to avoid reticulate/pandas assignment issues.
#               Added verification step that reloads and confirms structure.
# - 2026-01-09: Added Loupe file (.cloupe) export using loupeR package
# - 2026-01-07: Made metadata detection dynamic (no hardcoded column assumptions)
# - 2026-01-07: Fixed H5AD conversion using scipy sparse matrix approach
# - 2026-01-07: Added final_output_object.rds and .h5ad export with structure printing
# - 2026-01-06: Fixed switch() NULL return for unknown DE methods (CRITICAL FIX)
# - 2026-01-03: Added missing dplyr library (CRITICAL FIX)
# - 2026-01-03: Added Matrix library for consistency
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 10: FINAL SUMMARY\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

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
  library(dplyr)
  library(Matrix)
  library(scCustomize)
  library(reticulate)
})

# ==============================================================================
# DYNAMIC OUTPUT LOCATION
# ==============================================================================
# The output directory is determined dynamically from params$out_root
# This is set in params.R based on environment variables:
#   - DOWNSTREAM_DIR (base clustering output directory)
#   - VENTRICLE_FILTER (optional ventricle-specific subdirectory)
# The final path is constructed as: {DOWNSTREAM_DIR}/10_Downstream_Analysis_{VENTRICLE}
# ==============================================================================
out_base <- params$out_root
cat(">>> DYNAMIC OUTPUT LOCATION: ", out_base, " <<<\n\n")

load(file.path(out_base, "objects", "pipeline_environment.RData"))

# Load final clustered object
leiden_file <- file.path(out_base, "objects", "07_leiden_data.RData")
if (file.exists(leiden_file)) {
  load(leiden_file)
}

# Load normalization data
norm_file <- file.path(out_base, "objects", "03_normalization_data.RData")
if (file.exists(norm_file)) {
  load(norm_file)
}

# Load integration data
integration_file <- file.path(out_base, "objects", "04_integration_data.RData")
if (file.exists(integration_file)) {
  load(integration_file)
}

# Load DE data
de_file <- file.path(out_base, "objects", "08_de_data.RData")
if (file.exists(de_file)) {
  load(de_file)
}

# Load scICE data
scice_file <- file.path(out_base, "objects", "06_scice_data.RData")
scice_loaded <- FALSE
idclust_loaded <- FALSE
if (file.exists(scice_file)) {
  load(scice_file)
  scice_loaded <- TRUE
  cat("Loaded scICE subclustering data\n")
  if (exists("idclust_success") && isTRUE(idclust_success)) {
    idclust_loaded <- TRUE
    cat("Loaded IDclust subclustering data\n")
  }
}

# Load CHOIR data
choir_file <- file.path(out_base, "objects", "05_choir_data.RData")
if (file.exists(choir_file)) {
  load(choir_file)
  cat("Loaded CHOIR clustering data\n")
}

# ==============================================================================
# HELPER FUNCTION: Safe metadata access
# ==============================================================================
safe_get_metadata_count <- function(obj, column, value) {
  if (is.null(obj)) return("N/A")
  if (!column %in% colnames(obj@meta.data)) return("N/A")
  sum(obj@meta.data[[column]] == value, na.rm = TRUE)
}

safe_get_metadata_levels <- function(obj, column) {
  if (is.null(obj)) return(NULL)
  if (!column %in% colnames(obj@meta.data)) return(NULL)
  unique(obj@meta.data[[column]])
}

# ==============================================================================
# Detect available metadata columns
# ==============================================================================
cat("--- Detecting available metadata ---\n\n")

# Define potential grouping columns (will use whichever exist)
potential_group_cols <- c("sample_name", "sample", "orig.ident", "Sample")
potential_sex_cols <- c("sex", "Sex", "gender", "Gender")
potential_batch_cols <- c("batch", "Batch", "sequencing_batch", "library")
potential_condition_cols <- c("condition", "Condition", "treatment", "Treatment", "group", "Group")

# Find which columns actually exist in the data
find_existing_col <- function(obj, candidates) {
  if (is.null(obj)) return(NULL)
  meta_cols <- colnames(obj@meta.data)
  for (col in candidates) {
    if (col %in% meta_cols) return(col)
  }
  return(NULL)
}

# Detect columns from clustered object
sample_col <- NULL
sex_col <- NULL
batch_col <- NULL
condition_col <- NULL

if (exists("clustered_obj") && !is.null(clustered_obj)) {
  sample_col <- find_existing_col(clustered_obj, potential_group_cols)
  sex_col <- find_existing_col(clustered_obj, potential_sex_cols)
  batch_col <- find_existing_col(clustered_obj, potential_batch_cols)
  condition_col <- find_existing_col(clustered_obj, potential_condition_cols)

  cat("Detected metadata columns:\n")
  cat("  Sample column:    ", if (!is.null(sample_col)) sample_col else "not found", "\n")
  cat("  Sex column:       ", if (!is.null(sex_col)) sex_col else "not found", "\n")
  cat("  Batch column:     ", if (!is.null(batch_col)) batch_col else "not found", "\n")
  cat("  Condition column: ", if (!is.null(condition_col)) condition_col else "not found", "\n")
  cat("\n")
}

# ==============================================================================
# Create Summary Table
# ==============================================================================
cat("--- Creating Analysis Summary ---\n\n")

# Get sample information
samples_analyzed <- if (exists("params") && !is.null(params$samples_to_analyze)) {
  paste(params$samples_to_analyze, collapse = ", ")
} else {
  "N/A"
}

# Build base summary data
summary_metrics <- c(
  "Analysis name",
  "Samples analyzed",
  "Total cells",
  "Number of clusters",
  "Clustering method",
  "scICE subclustering",
  "IDclust subclustering",
  "DE analysis scope",
  "Normalization method",
  "Integration method",
  "Batch variable",
  "Batch integration enabled",
  "Analysis date",
  "R version"
)

summary_values <- c(
  "Choroid Plexus Multi-Sample Analysis",
  samples_analyzed,
  if (exists("clustered_obj") && !is.null(clustered_obj)) ncol(clustered_obj) else "N/A",
  if (exists("clustered_obj") && !is.null(clustered_obj)) length(unique(clustered_obj$seurat_clusters)) else "N/A",
  if (exists("clustering_method_used")) clustering_method_used else "N/A",
  if (scice_loaded && exists("scice_success") && scice_success) "Enabled" else "Not run",
  if (idclust_loaded) "Enabled" else "Not run",
  if (!is.null(params$de_comparison_scope)) params$de_comparison_scope else "global",
  if (exists("selected_normalization_method")) selected_normalization_method else "SCTransform",
  if (exists("best_method")) best_method else "N/A",
  if (!is.null(params$batch_variable)) params$batch_variable else "batch",
  as.character(isTRUE(params$run_batch_integration)),
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(R.version$major, R.version$minor, sep = ".")
)

# Dynamically add sex-based cell counts if sex column exists
if (!is.null(sex_col) && exists("clustered_obj") && !is.null(clustered_obj)) {
  sex_levels <- unique(clustered_obj@meta.data[[sex_col]])
  sex_levels <- sex_levels[!is.na(sex_levels)]

  for (sex_val in sex_levels) {
    count <- sum(clustered_obj@meta.data[[sex_col]] == sex_val, na.rm = TRUE)
    summary_metrics <- c(summary_metrics, paste0(sex_val, " cells"))
    summary_values <- c(summary_values, as.character(count))
  }
}

# Dynamically add condition-based counts if condition column exists
if (!is.null(condition_col) && exists("clustered_obj") && !is.null(clustered_obj)) {
  condition_levels <- unique(clustered_obj@meta.data[[condition_col]])
  condition_levels <- condition_levels[!is.na(condition_levels)]

  for (cond_val in condition_levels) {
    count <- sum(clustered_obj@meta.data[[condition_col]] == cond_val, na.rm = TRUE)
    summary_metrics <- c(summary_metrics, paste0(cond_val, " cells"))
    summary_values <- c(summary_values, as.character(count))
  }
}

summary_data <- data.frame(
  Metric = summary_metrics,
  Value = summary_values
)

cat(">>> FINAL ANALYSIS SUMMARY <<<\n\n")
print(summary_data, row.names = FALSE)

write.csv(summary_data, file.path(output_dirs$tables, "analysis_summary.csv"), row.names = FALSE)
cat("\nSaved:", file.path(output_dirs$tables, "analysis_summary.csv"), "\n")

# ==============================================================================
# Per-Sample Summary
# ==============================================================================
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("PER-SAMPLE SUMMARY\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

if (exists("clustered_obj") && !is.null(clustered_obj)) {
  meta_cols <- colnames(clustered_obj@meta.data)

  # Build grouping columns dynamically
  group_cols <- c()
  if (!is.null(sample_col) && sample_col %in% meta_cols) group_cols <- c(group_cols, sample_col)
  if (!is.null(sex_col) && sex_col %in% meta_cols) group_cols <- c(group_cols, sex_col)
  if (!is.null(batch_col) && batch_col %in% meta_cols) group_cols <- c(group_cols, batch_col)

  # If no grouping columns found, use a simple summary
 if (length(group_cols) == 0) {
    cat("No grouping columns detected. Generating simple summary.\n\n")

    simple_summary <- data.frame(
      Metric = c("Total cells", "Median nFeature_RNA", "Median nCount_RNA", "Median percent.mt", "Number of clusters"),
      Value = c(
        ncol(clustered_obj),
        if ("nFeature_RNA" %in% meta_cols) median(clustered_obj$nFeature_RNA, na.rm = TRUE) else NA,
        if ("nCount_RNA" %in% meta_cols) median(clustered_obj$nCount_RNA, na.rm = TRUE) else NA,
        if ("percent.mt" %in% meta_cols) median(clustered_obj$percent.mt, na.rm = TRUE) else NA,
        length(unique(clustered_obj$seurat_clusters))
      )
    )
    print(simple_summary, row.names = FALSE)
    write.csv(simple_summary, file.path(output_dirs$tables, "per_sample_summary.csv"), row.names = FALSE)

  } else {
    # Create dynamic grouping summary
    cat("Grouping by:", paste(group_cols, collapse = ", "), "\n\n")

    # Build the summarise call dynamically
    sample_summary <- clustered_obj@meta.data %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(
        n_cells = n(),
        median_nFeature = if ("nFeature_RNA" %in% meta_cols) median(nFeature_RNA, na.rm = TRUE) else NA,
        median_nCount = if ("nCount_RNA" %in% meta_cols) median(nCount_RNA, na.rm = TRUE) else NA,
        median_pct_mt = if ("percent.mt" %in% meta_cols) median(percent.mt, na.rm = TRUE) else NA,
        n_clusters = length(unique(seurat_clusters)),
        .groups = "drop"
      )

    print(as.data.frame(sample_summary))
    write.csv(sample_summary, file.path(output_dirs$tables, "per_sample_summary.csv"), row.names = FALSE)
  }

  cat("\nSaved:", file.path(output_dirs$tables, "per_sample_summary.csv"), "\n")
}

# ==============================================================================
# Method Selection Summary
# ==============================================================================
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("METHOD SELECTION SUMMARY\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

cat("Normalization:         ", if (exists("selected_normalization_method")) selected_normalization_method else "N/A", "\n")
cat("Integration:           ", if (exists("best_method")) best_method else "N/A", "\n")
cat("Integration reduction: ", if (exists("best_reduction")) best_reduction else "N/A", "\n")
cat("Clustering:            ", if (exists("clustering_method_used")) clustering_method_used else "N/A", "\n")
cat("Batch variable:        ", if (!is.null(params$batch_variable)) params$batch_variable else "batch", "\n")
cat("Batch integration:     ", as.character(isTRUE(params$run_batch_integration)), "\n")

# ==============================================================================
# Integration Methods Summary
# ==============================================================================
if (exists("benchmark_df") && !is.null(benchmark_df)) {
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("INTEGRATION METHODS BENCHMARKING\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")

  # Select columns that exist
  bench_cols <- intersect(c("method", "reduction", "batch_variance", "batch_asw", "lisi"),
                          colnames(benchmark_df))
  print(benchmark_df[, bench_cols, drop = FALSE], row.names = FALSE)
  cat("\nBest method selected:", if (exists("best_method")) best_method else "N/A", "\n")
}

# ==============================================================================
# scICE Subclustering Summary
# ==============================================================================
if (scice_loaded && exists("scice_results") && length(scice_results) > 0) {
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("scICE SUBCLUSTERING SUMMARY\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")

  scice_summary <- data.frame(
    Cluster = character(),
    Input_Cells = numeric(),
    Output_Cells = numeric(),
    Method = character(),
    Success = character(),
    stringsAsFactors = FALSE
  )

  for (clust_name in names(scice_results)) {
    result <- scice_results[[clust_name]]
    scice_summary <- rbind(scice_summary, data.frame(
      Cluster = clust_name,
      Input_Cells = if (!is.null(result$input_cells)) result$input_cells else NA,
      Output_Cells = if (!is.null(result$n_cells)) result$n_cells else NA,
      Method = if (!is.null(result$method)) result$method else "N/A",
      Success = if (!is.null(result$success) && result$success) "Yes" else "No",
      stringsAsFactors = FALSE
    ))
  }

  print(scice_summary, row.names = FALSE)

  n_successful <- sum(sapply(scice_results, function(x) isTRUE(x$success)))
  n_total <- length(scice_results)
  cat("\nTotal clusters processed:", n_total, "\n")
  cat("Successful:", n_successful, "\n")

  write.csv(scice_summary, file.path(output_dirs$tables, "scice_subclustering_summary.csv"), row.names = FALSE)
}
# ==============================================================================
# IDclust Subclustering Summary
# ==============================================================================
if (idclust_loaded && exists("idclust_results") && length(idclust_results) > 0) {
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("IDclust SUBCLUSTERING SUMMARY\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")

  idclust_summary <- data.frame(
    Cluster = character(),
    Input_Cells = numeric(),
    N_Subclusters = numeric(),
    Elapsed_Min = numeric(),
    Success = character(),
    stringsAsFactors = FALSE
  )

  for (clust_name in names(idclust_results)) {
    result <- idclust_results[[clust_name]]
    idclust_summary <- rbind(idclust_summary, data.frame(
      Cluster = clust_name,
      Input_Cells = if (!is.null(result$n_cells)) result$n_cells else NA,
      N_Subclusters = if (!is.null(result$n_subclusters)) result$n_subclusters else NA,
      Elapsed_Min = if (!is.null(result$elapsed_min)) round(result$elapsed_min, 1) else NA,
      Success = if (isTRUE(result$success)) "Yes" else "No",
      stringsAsFactors = FALSE
    ))
  }

  print(idclust_summary, row.names = FALSE)

  n_successful <- sum(sapply(idclust_results, function(x) isTRUE(x$success)))
  n_total <- length(idclust_results)
  cat("\nTotal clusters processed:", n_total, "\n")
  cat("Successful:", n_successful, "\n")

  write.csv(idclust_summary, file.path(output_dirs$tables, "idclust_subclustering_summary.csv"), row.names = FALSE)
}

# ==============================================================================
# DE Analysis Summary
# ==============================================================================
if (exists("all_de_results") && length(all_de_results) > 0) {
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("DIFFERENTIAL EXPRESSION SUMMARY\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")

  de_summary_table <- data.frame(
    Method = character(),
    Genes_Tested = numeric(),
    Significant = numeric(),
    stringsAsFactors = FALSE
  )

  # Define mapping of method names to their significance column names
  # Include common variations and additional methods
  sig_col_mapping <- list(
    "MAST" = "FDR",
    "mast" = "FDR",
    "edgeR" = "adj.P.Val",
    "edger" = "adj.P.Val",
    "DESeq2" = "padj",
    "deseq2" = "padj",
    "Deseq2" = "padj",
    "Wilcoxon" = "p_val_adj",
    "wilcoxon" = "p_val_adj",
    "wilcox" = "p_val_adj",
    "FindMarkers" = "p_val_adj",
    "findmarkers" = "p_val_adj",
    "t-test" = "p_val_adj",
    "bimod" = "p_val_adj",
    "LR" = "p_val_adj",
    "negbinom" = "p_val_adj",
    "poisson" = "p_val_adj",
    "MAST_seurat" = "p_val_adj"
  )

  cat("DE methods found:", paste(names(all_de_results), collapse = ", "), "\n\n")

  for (method in names(all_de_results)) {
    res <- all_de_results[[method]]

    # Skip if result is NULL or not a data frame
    if (is.null(res) || !is.data.frame(res)) {
      cat("  Warning: Skipping '", method, "' - result is NULL or not a data frame\n", sep = "")
      next
    }

    # Try to get sig_col from mapping first
    sig_col <- sig_col_mapping[[method]]

    # If not in mapping, try to auto-detect common significance column names
    if (is.null(sig_col)) {
      possible_sig_cols <- c("FDR", "adj.P.Val", "padj", "p_val_adj", "q_value",
                             "qvalue", "p.adjust", "BH", "fdr", "adjusted_pvalue")
      detected_cols <- intersect(possible_sig_cols, colnames(res))

      if (length(detected_cols) > 0) {
        sig_col <- detected_cols[1]  # Use first match
        cat("  Note: Auto-detected significance column '", sig_col, "' for method '", method, "'\n", sep = "")
      } else {
        cat("  Warning: Unknown DE method '", method, "' with columns: ",
            paste(head(colnames(res), 10), collapse = ", "),
            if (ncol(res) > 10) "..." else "", "\n", sep = "")
        cat("           Could not find a recognized significance column. Skipping.\n")
        next
      }
    }

    # Now check if the sig_col exists in the results
    if (sig_col %in% colnames(res)) {
      n_sig <- sum(res[[sig_col]] < 0.05, na.rm = TRUE)
      de_summary_table <- rbind(de_summary_table, data.frame(
        Method = method,
        Genes_Tested = nrow(res),
        Significant = n_sig,
        stringsAsFactors = FALSE
      ))
    } else {
      cat("  Warning: Expected column '", sig_col, "' not found in '", method, "' results\n", sep = "")
      cat("           Available columns: ", paste(colnames(res), collapse = ", "), "\n", sep = "")
    }
  }

  if (nrow(de_summary_table) > 0) {
    cat("\n")
    print(de_summary_table, row.names = FALSE)
    write.csv(de_summary_table, file.path(output_dirs$tables, "de_summary.csv"), row.names = FALSE)
    cat("\nSaved:", file.path(output_dirs$tables, "de_summary.csv"), "\n")
  } else {
    cat("No DE results could be summarized.\n")
  }
}

# ==============================================================================
# Save Final Clustered Object
# ==============================================================================
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("SAVING FINAL OBJECTS\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

final_obj <- NULL

# Determine which object to use as final
if (exists("choir_obj") && !is.null(choir_obj) &&
    any(c("scice_subcluster", "idclust_subcluster") %in% colnames(choir_obj@meta.data))) {
  final_obj <- choir_obj
  subcl_cols <- intersect(c("scice_subcluster", "idclust_subcluster"), colnames(choir_obj@meta.data))
  cat("Using subclustered object as final object (columns:", paste(subcl_cols, collapse = ", "), ")\n")
} else if (exists("clustered_obj") && !is.null(clustered_obj)) {
  final_obj <- clustered_obj
  cat("Using clustered object as final object\n")
}

if (!is.null(final_obj)) {
  # --------------------------------------------------------------------------
  # Save as final_clustered_object.rds (existing behavior)
  # --------------------------------------------------------------------------
  final_rds_path <- file.path(output_dirs$objects, "final_clustered_object.rds")
  saveRDS(final_obj, final_rds_path)
  cat("Saved final clustered object:", final_rds_path, "\n")
  cat("Size:", round(file.size(final_rds_path) / 1e6, 2), "MB\n")

  cat("\n>>> FINAL OBJECT STRUCTURE <<<\n")
  print_object_structure(final_obj, "Final Clustered Object")

  # --------------------------------------------------------------------------
  # Save as final_output_object.rds (NEW - to out_base root directory)
  # --------------------------------------------------------------------------
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("SAVING FINAL OUTPUT OBJECT (RDS + H5AD)\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")

  cat(">>> Output directory (out_base):", out_base, "<<<\n\n")

  # Save RDS to root output directory
  final_output_rds_path <- file.path(out_base, "final_output_object.rds")
  saveRDS(final_obj, final_output_rds_path)
  cat("Saved final_output_object.rds:", final_output_rds_path, "\n")
  cat("Size:", round(file.size(final_output_rds_path) / 1e6, 2), "MB\n")

  # --------------------------------------------------------------------------
  # Print Seurat object str()
  # --------------------------------------------------------------------------
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("SEURAT OBJECT STRUCTURE (str)\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")

  # Capture str() output to file and print
  str_output_file <- file.path(out_base, "final_output_object_seurat_str.txt")
  sink(str_output_file)
  cat("Seurat Object Structure for: final_output_object.rds\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  str(final_obj)
  sink()

  cat("Seurat str() output saved to:", str_output_file, "\n\n")

  # Also print to console (truncated for readability)
  cat(">>> Seurat str() output (console preview): <<<\n")
  str(final_obj, max.level = 2)

  # --------------------------------------------------------------------------
  # Convert to H5AD using FIXED approach - pass all data at creation time
  # --------------------------------------------------------------------------
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("CONVERTING TO H5AD FORMAT (FIXED: pass all data at creation)\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")

  h5ad_conversion_success <- FALSE
  final_output_h5ad_path <- file.path(out_base, "final_output_object.h5ad")

  tryCatch({
    # Check if Python is available
    if (!reticulate::py_available(initialize = TRUE)) {
      stop("Python is not available via reticulate")
    }

    # Check if anndata is installed
    anndata_available <- tryCatch({
      reticulate::py_module_available("anndata")
    }, error = function(e) FALSE)

    if (!anndata_available) {
      stop("Python anndata module is not available. Install with: pip install anndata")
    }

    cat("Python environment detected.\n")
    cat("Converting Seurat object to AnnData format...\n\n")

    # Import Python modules with convert = FALSE for better control
    ad <- reticulate::import("anndata", convert = FALSE)
    np <- reticulate::import("numpy", convert = FALSE)
    pd <- reticulate::import("pandas", convert = FALSE)
    sp <- reticulate::import("scipy.sparse", convert = FALSE)
    builtins <- reticulate::import_builtins()

    # Get anndata version
    ilm <- tryCatch({
      reticulate::import("importlib.metadata", convert = FALSE)
    }, error = function(e) NULL)

    anndata_ver <- if (!is.null(ilm)) {
      tryCatch(reticulate::py_to_r(ilm$version("anndata")), error = function(e) "unknown")
    } else {
      "unknown"
    }
    cat("  anndata version:", anndata_ver, "\n")

    # -------------------------------------------------------------------------
    # Helper function: Convert dgCMatrix (genes x cells) to scipy CSR (cells x genes)
    # -------------------------------------------------------------------------
    dgCMatrix_to_scipy_csr <- function(mat, sp, np) {
      stopifnot(inherits(mat, "dgCMatrix"))
      data_arr <- np$array(mat@x, dtype = "float64")
      indices_arr <- np$array(mat@i, dtype = "int32")
      indptr_arr <- np$array(mat@p, dtype = "int32")
      shape <- reticulate::tuple(as.integer(nrow(mat)), as.integer(ncol(mat)))
      csc <- sp$csc_matrix(reticulate::tuple(data_arr, indices_arr, indptr_arr), shape = shape)
      csc$transpose()$tocsr()  # Transpose to cells x genes
    }

    # -------------------------------------------------------------------------
    # Helper function: Get layer data safely
    # -------------------------------------------------------------------------
    get_layer_safe <- function(obj, assay, layer_name, fallback_slot) {
      out <- tryCatch(
        LayerData(obj, assay = assay, layer = layer_name),
        error = function(e) NULL
      )
      if (is.null(out)) {
        out <- tryCatch(
          GetAssayData(obj, assay = assay, slot = fallback_slot),
          error = function(e) NULL
        )
      }
      out
    }

    # -------------------------------------------------------------------------
    # Helper function: Safely get shape from scipy sparse matrix
    # -------------------------------------------------------------------------
    get_scipy_shape <- function(scipy_mat) {
      shape <- reticulate::py_to_r(scipy_mat$shape)
      if (is.list(shape)) {
        return(c(as.integer(shape[[1]]), as.integer(shape[[2]])))
      } else {
        return(as.integer(shape))
      }
    }

    # -------------------------------------------------------------------------
    # Helper function: Convert Python dict_keys/KeysView to R character vector
    # -------------------------------------------------------------------------
    py_keys_to_r <- function(py_keys_obj, builtins) {
      tryCatch({
        as.character(builtins$list(py_keys_obj))
      }, error = function(e) {
        tryCatch({
          as.character(reticulate::iterate(py_keys_obj))
        }, error = function(e2) {
          character(0)
        })
      })
    }

    # -------------------------------------------------------------------------
    # Extract expression data from Seurat object
    # -------------------------------------------------------------------------
    default_assay <- DefaultAssay(final_obj)
    cat("  Default assay:", default_assay, "\n")

    # Get normalized data for X
    data_mat <- get_layer_safe(final_obj, default_assay, "data", "data")
    if (is.null(data_mat)) {
      stop("Could not extract 'data' layer from Seurat object")
    }
    cat("  Extracting X (normalized data)...\n")
    X_scipy <- dgCMatrix_to_scipy_csr(data_mat, sp, np)
    X_shape <- get_scipy_shape(X_scipy)
    cat("    X shape:", X_shape[1], "x", X_shape[2], "\n")

    # Get counts
    counts_mat <- get_layer_safe(final_obj, default_assay, "counts", "counts")
    counts_scipy <- NULL
    if (!is.null(counts_mat)) {
      cat("  Extracting counts layer...\n")
      counts_scipy <- dgCMatrix_to_scipy_csr(counts_mat, sp, np)
      counts_shape <- get_scipy_shape(counts_scipy)
      cat("    counts shape:", counts_shape[1], "x", counts_shape[2], "\n")
    }

    # Get scale.data if available
    scale_mat <- get_layer_safe(final_obj, default_assay, "scale.data", "scale.data")
    scale_genes <- NULL
    if (!is.null(scale_mat) && nrow(scale_mat) > 0) {
      scale_genes <- rownames(scale_mat)
      cat("  Found scale.data for", length(scale_genes), "genes\n")
    }

    # =========================================================================
    # PREPARE ALL DATA BEFORE CREATING ANNDATA
    # This is the CRITICAL FIX - we must pass everything at creation time
    # =========================================================================
    cat("\n  Preparing all data components for AnnData...\n")

    # -------------------------------------------------------------------------
    # Prepare obs (cell metadata) as pandas DataFrame
    # -------------------------------------------------------------------------
    cat("  Preparing obs (cell metadata)...\n")
    obs_df <- final_obj@meta.data
    obs_list <- list()
    obs_cols_added <- 0
    obs_cols_failed <- 0

    for (col in colnames(obs_df)) {
      tryCatch({
        values <- obs_df[[col]]

        # Convert factors to character
        if (is.factor(values)) {
          values <- as.character(values)
        }

        # Convert logical to character
        if (is.logical(values)) {
          values <- as.character(values)
        }

        # Handle NA values
        if (any(is.na(values))) {
          if (is.numeric(values)) {
            values[is.na(values)] <- NaN
          } else {
            values[is.na(values)] <- "NA"
          }
        }

        obs_list[[col]] <- values
        obs_cols_added <- obs_cols_added + 1
      }, error = function(e) {
        obs_cols_failed <- obs_cols_failed + 1
        cat("    Warning: Failed to process obs column '", col, "': ", conditionMessage(e), "\n", sep = "")
      })
    }

    # Create pandas DataFrame for obs
    cell_names <- colnames(final_obj)
    obs_pandas <- pd$DataFrame(obs_list)
    obs_pandas$index <- np$array(cell_names)

    cat("    Prepared", obs_cols_added, "of", ncol(obs_df), "obs columns\n")
    if (obs_cols_failed > 0) {
      cat("    (", obs_cols_failed, " columns failed)\n", sep = "")
    }

    # -------------------------------------------------------------------------
    # Prepare var (gene metadata) as pandas DataFrame
    # -------------------------------------------------------------------------
    cat("  Preparing var (gene metadata)...\n")
    gene_names <- rownames(final_obj)
    var_list <- list()

    # Get any existing var metadata from Seurat
    var_df_seurat <- final_obj[[default_assay]]@meta.data
    if (nrow(var_df_seurat) > 0 && ncol(var_df_seurat) > 0) {
      for (col in colnames(var_df_seurat)) {
        tryCatch({
          values <- var_df_seurat[[col]]
          if (is.factor(values)) values <- as.character(values)
          if (is.logical(values)) values <- as.character(values)
          if (any(is.na(values))) {
            if (is.numeric(values)) {
              values[is.na(values)] <- NaN
            } else {
              values[is.na(values)] <- "NA"
            }
          }
          var_list[[col]] <- values
        }, error = function(e) {
          cat("    Warning: Failed to process var column '", col, "'\n", sep = "")
        })
      }
    }

    # Add highly variable genes indicator
    var_features <- VariableFeatures(final_obj)
    if (length(var_features) > 0) {
      var_list[["highly_variable"]] <- gene_names %in% var_features
      cat("    Added highly_variable indicator (", sum(gene_names %in% var_features), " HVGs)\n", sep = "")
    }

    # Create pandas DataFrame for var
    if (length(var_list) > 0) {
      var_pandas <- pd$DataFrame(var_list)
    } else {
      # Create empty DataFrame if no var metadata
      var_pandas <- pd$DataFrame()
    }
    var_pandas$index <- np$array(gene_names)

    cat("    Prepared", length(var_list), "var columns\n")

    # -------------------------------------------------------------------------
    # Prepare obsm (embeddings) as a Python dict
    # -------------------------------------------------------------------------
    cat("  Preparing obsm (embeddings)...\n")
    obsm_dict <- reticulate::dict()
    obsm_count <- 0

    for (red_name in names(final_obj@reductions)) {
      tryCatch({
        red <- final_obj@reductions[[red_name]]
        embeddings <- Embeddings(red)

        # Use scanpy-style naming: X_pca, X_umap, etc.
        obsm_key <- paste0("X_", tolower(gsub("[^a-zA-Z0-9]", "_", red_name)))

        obsm_dict[[obsm_key]] <- np$array(embeddings, dtype = "float64")
        cat("    obsm['", obsm_key, "']: ", nrow(embeddings), " x ", ncol(embeddings), "\n", sep = "")
        obsm_count <- obsm_count + 1
      }, error = function(e) {
        cat("    Warning: Could not prepare reduction '", red_name, "': ", conditionMessage(e), "\n", sep = "")
      })
    }

    cat("    Prepared", obsm_count, "embeddings\n")

    # -------------------------------------------------------------------------
    # Prepare varm (gene loadings) as a Python dict
    # -------------------------------------------------------------------------
    cat("  Preparing varm (gene loadings)...\n")
    varm_dict <- reticulate::dict()
    varm_count <- 0

    for (red_name in names(final_obj@reductions)) {
      tryCatch({
        red <- final_obj@reductions[[red_name]]
        loadings <- Loadings(red)

        if (nrow(loadings) > 0 && ncol(loadings) > 0) {
          # Create full loading matrix with zeros for genes without loadings
          full_loadings <- matrix(0, nrow = length(gene_names), ncol = ncol(loadings))
          rownames(full_loadings) <- gene_names

          loading_genes <- rownames(loadings)
          gene_idx <- match(loading_genes, gene_names)
          valid_idx <- !is.na(gene_idx)
          full_loadings[gene_idx[valid_idx], ] <- as.matrix(loadings[valid_idx, ])

          varm_key <- paste0("loadings_", tolower(gsub("[^a-zA-Z0-9]", "_", red_name)))
          varm_dict[[varm_key]] <- np$array(full_loadings, dtype = "float64")
          cat("    varm['", varm_key, "']: ", nrow(full_loadings), " x ", ncol(full_loadings), "\n", sep = "")
          varm_count <- varm_count + 1
        }
      }, error = function(e) {
        # Silently skip reductions without loadings
      })
    }

    cat("    Prepared", varm_count, "loading matrices\n")

    # -------------------------------------------------------------------------
    # Prepare layers as a Python dict
    # -------------------------------------------------------------------------
    cat("  Preparing layers...\n")
    layers_dict <- reticulate::dict()

    if (!is.null(counts_scipy)) {
      layers_dict[["counts"]] <- counts_scipy
      cat("    layers['counts']: scCDC-corrected raw counts\n")
    }

    # Add scale.data as dense layer if it exists and is not too large
    if (!is.null(scale_mat) && nrow(scale_mat) > 0) {
      tryCatch({
        # Create full matrix with NaN for non-scaled genes
        n_cells <- length(cell_names)
        n_genes <- length(gene_names)
        full_scale <- matrix(NaN, nrow = n_cells, ncol = n_genes)

        # Fill in scaled values for genes that have them
        gene_idx <- match(scale_genes, gene_names)
        valid_idx <- !is.na(gene_idx)
        full_scale[, gene_idx[valid_idx]] <- t(as.matrix(scale_mat[valid_idx, ]))

        layers_dict[["scale_data"]] <- np$array(full_scale, dtype = "float32")
        cat("    layers['scale_data']: scaled expression for", length(scale_genes), "genes\n")
      }, error = function(e) {
        cat("    Warning: Could not add scale_data layer:", conditionMessage(e), "\n")
      })
    }

    # -------------------------------------------------------------------------
    # Prepare uns (unstructured) metadata
    # -------------------------------------------------------------------------
    cat("  Preparing uns (unstructured metadata)...\n")
    uns_dict <- reticulate::dict()

    # Variable features
    if (length(var_features) > 0) {
      uns_dict[["var_features"]] <- var_features
      uns_dict[["n_var_features"]] <- as.integer(length(var_features))
    }

    # Scale data gene list
    if (!is.null(scale_genes)) {
      uns_dict[["scale_data_genes"]] <- scale_genes
    }

    # Pipeline metadata
    uns_dict[["pipeline"]] <- "Choroid Plexus scRNA-seq Multi-Sample Analysis"
    uns_dict[["normalization_method"]] <- if (exists("selected_normalization_method")) selected_normalization_method else "unknown"
    uns_dict[["integration_method"]] <- if (exists("best_method")) best_method else "unknown"
    uns_dict[["clustering_method"]] <- if (exists("clustering_method_used")) clustering_method_used else "unknown"
    uns_dict[["counts_source"]] <- "scCDC_corrected (contamination-corrected counts)"
    uns_dict[["conversion_date"]] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    uns_dict[["seurat_version"]] <- as.character(packageVersion("Seurat"))
    uns_dict[["R_version"]] <- paste(R.version$major, R.version$minor, sep = ".")

    cat("    Prepared pipeline metadata\n")

    # =========================================================================
    # CREATE ANNDATA WITH ALL DATA AT ONCE - CRITICAL FIX
    # =========================================================================
    cat("\n  >>> Creating AnnData object with ALL components at once <<<\n")

    adata <- ad$AnnData(
      X = X_scipy,
      obs = obs_pandas,
      var = var_pandas,
      obsm = obsm_dict,
      varm = varm_dict,
      layers = layers_dict,
      uns = uns_dict
    )

    cat("  AnnData object created successfully.\n")

    # -------------------------------------------------------------------------
    # Set raw to counts if available
    # -------------------------------------------------------------------------
    if (!is.null(counts_scipy)) {
      tryCatch({
        # Create a minimal AnnData for raw
        raw_var_pandas <- pd$DataFrame()
        raw_var_pandas$index <- np$array(gene_names)
        adata$raw <- ad$AnnData(X = counts_scipy, var = raw_var_pandas)
        cat("  Set .raw to counts matrix\n")
      }, error = function(e) {
        cat("  Note: Could not set .raw:", conditionMessage(e), "\n")
      })
    }

    # -------------------------------------------------------------------------
    # Write H5AD file
    # -------------------------------------------------------------------------
    cat("\n  Writing H5AD file...\n")
    cat("    Output:", final_output_h5ad_path, "\n")
    adata$write_h5ad(final_output_h5ad_path)

    h5ad_conversion_success <- TRUE

    if (file.exists(final_output_h5ad_path)) {
      cat("    Size:", round(file.size(final_output_h5ad_path) / 1e6, 2), "MB\n")
    }

    # =========================================================================
    # VERIFICATION: Reload and confirm structure
    # =========================================================================
    cat("\n", paste(rep("-", 60), collapse = ""), "\n")
    cat("VERIFICATION: Reloading H5AD to confirm data integrity\n")
    cat(paste(rep("-", 60), collapse = ""), "\n\n")

    # Reload the file
    adata_verify <- ad$read_h5ad(final_output_h5ad_path)

    # Get actual structure using proper Python methods
    n_obs_verify <- reticulate::py_to_r(adata_verify$n_obs)
    n_vars_verify <- reticulate::py_to_r(adata_verify$n_vars)

    # Get obs columns
    obs_cols_verify <- tryCatch({
      cols <- reticulate::py_to_r(adata_verify$obs$columns$tolist())
      if (is.null(cols)) character(0) else as.character(cols)
    }, error = function(e) character(0))

    # Get var columns
    var_cols_verify <- tryCatch({
      cols <- reticulate::py_to_r(adata_verify$var$columns$tolist())
      if (is.null(cols)) character(0) else as.character(cols)
    }, error = function(e) character(0))

    # Get obsm keys - FIXED: use builtins$list to convert KeysView
    obsm_keys_verify <- tryCatch({
      py_keys_to_r(adata_verify$obsm$keys(), builtins)
    }, error = function(e) character(0))

    # Get varm keys - FIXED
    varm_keys_verify <- tryCatch({
      py_keys_to_r(adata_verify$varm$keys(), builtins)
    }, error = function(e) character(0))

    # Get layers keys - FIXED
    layers_keys_verify <- tryCatch({
      py_keys_to_r(adata_verify$layers$keys(), builtins)
    }, error = function(e) character(0))

    # Get uns keys - FIXED
    uns_keys_verify <- tryCatch({
      py_keys_to_r(adata_verify$uns$keys(), builtins)
    }, error = function(e) character(0))

    # Check if raw exists
    has_raw_verify <- tryCatch({
      !is.null(adata_verify$raw)
    }, error = function(e) FALSE)

    # Print verification results
    cat("  >>> VERIFIED H5AD STRUCTURE <<<\n\n")
    cat("  Shape: ", n_obs_verify, " cells x ", n_vars_verify, " genes\n\n", sep = "")

    cat("  obs (cell metadata):\n")
    cat("    Columns: ", length(obs_cols_verify), "\n", sep = "")
    if (length(obs_cols_verify) > 0) {
      cat("    Names: ", paste(head(obs_cols_verify, 20), collapse = ", "),
          if (length(obs_cols_verify) > 20) "..." else "", "\n", sep = "")
    }

    cat("\n  var (gene metadata):\n")
    cat("    Columns: ", length(var_cols_verify), "\n", sep = "")
    if (length(var_cols_verify) > 0) {
      cat("    Names: ", paste(head(var_cols_verify, 20), collapse = ", "),
          if (length(var_cols_verify) > 20) "..." else "", "\n", sep = "")
    }

    cat("\n  obsm (embeddings):\n")
    cat("    Keys: ", length(obsm_keys_verify), "\n", sep = "")
    if (length(obsm_keys_verify) > 0) {
      cat("    Names: ", paste(obsm_keys_verify, collapse = ", "), "\n", sep = "")
      # Print dimensions for each embedding
      for (key in obsm_keys_verify) {
        tryCatch({
          emb_shape <- reticulate::py_to_r(adata_verify$obsm[[key]]$shape)
          cat("      ", key, ": ", emb_shape[[1]], " x ", emb_shape[[2]], "\n", sep = "")
        }, error = function(e) {})
      }
    }

    cat("\n  varm (gene loadings):\n")
    cat("    Keys: ", length(varm_keys_verify), "\n", sep = "")
    if (length(varm_keys_verify) > 0) {
      cat("    Names: ", paste(varm_keys_verify, collapse = ", "), "\n", sep = "")
    }

    cat("\n  layers:\n")
    cat("    Keys: ", length(layers_keys_verify), "\n", sep = "")
    if (length(layers_keys_verify) > 0) {
      cat("    Names: ", paste(layers_keys_verify, collapse = ", "), "\n", sep = "")
    }

    cat("\n  uns (unstructured):\n")
    cat("    Keys: ", length(uns_keys_verify), "\n", sep = "")
    if (length(uns_keys_verify) > 0) {
      cat("    Names: ", paste(uns_keys_verify, collapse = ", "), "\n", sep = "")
    }

    cat("\n  raw: ", if (has_raw_verify) "present" else "not present", "\n", sep = "")

    # -------------------------------------------------------------------------
    # Validation checks
    # -------------------------------------------------------------------------
    cat("\n  >>> VALIDATION CHECKS <<<\n\n")

    validation_passed <- TRUE

    # Check obs columns
    if (length(obs_cols_verify) == 0) {
      cat("  [FAIL] obs has no columns!\n")
      validation_passed <- FALSE
    } else if (length(obs_cols_verify) < obs_cols_added) {
      cat("  [WARN] obs has fewer columns than expected (", length(obs_cols_verify), " vs ", obs_cols_added, ")\n", sep = "")
    } else {
      cat("  [PASS] obs has ", length(obs_cols_verify), " columns\n", sep = "")
    }

    # Check obsm embeddings
    if (length(obsm_keys_verify) == 0) {
      cat("  [FAIL] obsm has no embeddings!\n")
      validation_passed <- FALSE
    } else if (length(obsm_keys_verify) < obsm_count) {
      cat("  [WARN] obsm has fewer embeddings than expected (", length(obsm_keys_verify), " vs ", obsm_count, ")\n", sep = "")
    } else {
      cat("  [PASS] obsm has ", length(obsm_keys_verify), " embeddings\n", sep = "")
    }

    # Check for UMAP specifically
    if ("X_umap" %in% obsm_keys_verify) {
      cat("  [PASS] X_umap embedding present\n")
    } else {
      cat("  [WARN] X_umap embedding not found\n")
    }

    # Check dimensions match
    if (n_obs_verify == ncol(final_obj) && n_vars_verify == nrow(final_obj)) {
      cat("  [PASS] Dimensions match Seurat object (", n_obs_verify, " x ", n_vars_verify, ")\n", sep = "")
    } else {
      cat("  [FAIL] Dimension mismatch! Seurat: ", ncol(final_obj), " x ", nrow(final_obj),
          ", H5AD: ", n_obs_verify, " x ", n_vars_verify, "\n", sep = "")
      validation_passed <- FALSE
    }

    if (validation_passed) {
      cat("\n  >>> ALL VALIDATION CHECKS PASSED <<<\n")
    } else {
      cat("\n  >>> SOME VALIDATION CHECKS FAILED - please review <<<\n")
    }

    # -------------------------------------------------------------------------
    # Save H5AD structure to file
    # -------------------------------------------------------------------------
    h5ad_str_output_file <- file.path(out_base, "final_output_object_h5ad_structure.txt")

    sink(h5ad_str_output_file)
    cat("AnnData Object Structure for: final_output_object.h5ad\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")

    cat(">>> Basic Info <<<\n")
    cat("Shape: n_obs x n_vars =", n_obs_verify, "x", n_vars_verify, "\n\n")

    cat(">>> obs (cell metadata) <<<\n")
    cat("Number of columns:", length(obs_cols_verify), "\n")
    cat("Column names:", paste(obs_cols_verify, collapse = ", "), "\n\n")

    cat(">>> var (gene metadata) <<<\n")
    cat("Number of columns:", length(var_cols_verify), "\n")
    cat("Column names:", paste(var_cols_verify, collapse = ", "), "\n\n")

    cat(">>> layers <<<\n")
    cat("Layers:", if (length(layers_keys_verify) > 0) paste(layers_keys_verify, collapse = ", ") else "none", "\n\n")

    cat(">>> obsm (embeddings) <<<\n")
    if (length(obsm_keys_verify) > 0) {
      for (key in obsm_keys_verify) {
        tryCatch({
          emb_shape <- reticulate::py_to_r(adata_verify$obsm[[key]]$shape)
          cat("  ", key, ":", emb_shape[[1]], "x", emb_shape[[2]], "\n", sep = "")
        }, error = function(e) {
          cat("  ", key, ": (shape unavailable)\n", sep = "")
        })
      }
    } else {
      cat("  none\n")
    }
    cat("\n")

    cat(">>> varm (gene loadings) <<<\n")
    if (length(varm_keys_verify) > 0) {
      for (key in varm_keys_verify) {
        tryCatch({
          load_shape <- reticulate::py_to_r(adata_verify$varm[[key]]$shape)
          cat("  ", key, ":", load_shape[[1]], "x", load_shape[[2]], "\n", sep = "")
        }, error = function(e) {
          cat("  ", key, ": (shape unavailable)\n", sep = "")
        })
      }
    } else {
      cat("  none\n")
    }
    cat("\n")

    cat(">>> uns (unstructured) <<<\n")
    cat("Keys:", if (length(uns_keys_verify) > 0) paste(uns_keys_verify, collapse = ", ") else "none", "\n\n")

    cat(">>> raw <<<\n")
    cat("Present:", if (has_raw_verify) "yes" else "no", "\n\n")

    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("END OF H5AD STRUCTURE\n")
    sink()

    cat("\n  H5AD structure saved to:", h5ad_str_output_file, "\n")

  }, error = function(e) {
    cat("\n[WARNING] H5AD conversion failed:\n")
    cat("  Error:", conditionMessage(e), "\n")
    cat("  This may be due to Python/reticulate configuration issues.\n")
    cat("  The RDS file has been saved successfully.\n")
    cat("  To manually convert, run the standalone conversion script.\n\n")
  })

  # ==========================================================================
  # Convert to Loupe file (.cloupe) using loupeR package
  # ==========================================================================
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("CONVERTING TO LOUPE FORMAT (.cloupe)\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")
  flush.console()

  loupe_conversion_success <- FALSE
  final_output_loupe_path <- file.path(out_base, "final_output_object.cloupe")

  # ==========================================================================
  # SKIP FLAG - Set to TRUE to skip Loupe conversion entirely
  # ==========================================================================
  SKIP_LOUPE_CONVERSION <- TRUE  # Change to FALSE to attempt Loupe conversion

  if (SKIP_LOUPE_CONVERSION) {
    cat("[SKIPPED] Loupe conversion disabled via SKIP_LOUPE_CONVERSION flag\n")
    cat("  To enable, set SKIP_LOUPE_CONVERSION <- FALSE in 10_final_summary.R\n")
    cat("  Note: loupeR can cause crashes on some systems\n\n")
    loupe_conversion_success <- FALSE

  } else {

    # Wrap everything in tryCatch with explicit flushing
    loupe_conversion_success <- tryCatch({

      # -------------------------------------------------------------------------
      # Check if loupeR package is available
      # -------------------------------------------------------------------------
      if (!requireNamespace("loupeR", quietly = TRUE)) {
        cat("[SKIPPED] loupeR package is not installed\n")
        cat("  Install with: remotes::install_github('10xGenomics/loupeR')\n\n")
        return(FALSE)
      }

      # Check loupeR setup (this can fail if eula not accepted)
      setup_ok <- tryCatch({
        loupeR::setup()
        TRUE
      }, error = function(e) {
        cat("[WARNING] loupeR setup check failed:", conditionMessage(e), "\n")
        cat("  You may need to run loupeR::setup() interactively first\n")
        FALSE
      })

      if (!setup_ok) {
        cat("[SKIPPED] loupeR setup not complete\n\n")
        return(FALSE)
      }

      library(loupeR)
      cat("loupeR package loaded successfully.\n")
      cat("  loupeR version:", as.character(packageVersion("loupeR")), "\n\n")
      flush.console()

      # -------------------------------------------------------------------------
      # Prepare count matrix
      # -------------------------------------------------------------------------
      cat("  Preparing count matrix...\n")
      flush.console()

      default_assay <- DefaultAssay(final_obj)

      # Get counts layer (loupeR requires raw counts)
      counts_mat <- tryCatch(
        LayerData(final_obj, assay = default_assay, layer = "counts"),
        error = function(e) {
          tryCatch(
            GetAssayData(final_obj, assay = default_assay, slot = "counts"),
            error = function(e2) NULL
          )
        }
      )

      if (is.null(counts_mat)) {
        cat("  [ERROR] Could not extract counts matrix from Seurat object\n")
        return(FALSE)
      }

      # Ensure it's a dgCMatrix (sparse matrix)
      if (!inherits(counts_mat, "dgCMatrix")) {
        cat("    Converting to sparse matrix format...\n")
        flush.console()
        counts_mat <- as(counts_mat, "dgCMatrix")
      }

      cat("    Count matrix dimensions:", nrow(counts_mat), "genes x", ncol(counts_mat), "cells\n")
      flush.console()

      # -------------------------------------------------------------------------
      # Prepare cluster information
      # -------------------------------------------------------------------------
      cat("  Preparing cluster information...\n")
      flush.console()

      clusters <- NULL
      cluster_col_used <- NULL

      # Priority order for cluster columns
      cluster_priority <- c("seurat_clusters", "CHOIR_clusters_0.05", "scice_subcluster",
                            "leiden_clusters", "RNA_snn_res.0.5")

      for (clust_col in cluster_priority) {
        if (clust_col %in% colnames(final_obj@meta.data)) {
          clust_values <- final_obj@meta.data[[clust_col]]
          if (!all(is.na(clust_values))) {
            clusters <- as.factor(clust_values)
            cluster_col_used <- clust_col
            break
          }
        }
      }

      if (!is.null(clusters)) {
        cat("    Using cluster column:", cluster_col_used, "\n")
        cat("    Number of clusters:", length(levels(clusters)), "\n")
      } else {
        cat("    Warning: No cluster information found\n")
      }
      flush.console()

      # -------------------------------------------------------------------------
      # Prepare projections (embeddings) - SIMPLIFIED
      # -------------------------------------------------------------------------
      cat("  Preparing projections (embeddings)...\n")
      flush.console()

      projections <- list()

      # Only try to get UMAP - most important for Loupe
      if ("umap" %in% names(final_obj@reductions)) {
        tryCatch({
          umap_emb <- Embeddings(final_obj, "umap")[, 1:2, drop = FALSE]
          projections[["UMAP"]] <- umap_emb
          cat("    Added UMAP projection:", nrow(umap_emb), "x 2\n")
        }, error = function(e) {
          cat("    Warning: Could not extract UMAP\n")
        })
      }

      if (length(projections) == 0) {
        cat("    Warning: No projections available - Loupe may not work well\n")
      }
      flush.console()

      # -------------------------------------------------------------------------
      # Skip metadata preparation - can cause issues
      # -------------------------------------------------------------------------
      cat("  Skipping complex metadata (stability)\n")
      flush.console()

      # -------------------------------------------------------------------------
      # Create Loupe file using loupeR - MINIMAL CALL
      # -------------------------------------------------------------------------
      cat("\n  Creating Loupe file (this may take a few minutes)...\n")
      cat("    Output:", final_output_loupe_path, "\n")
      flush.console()

      # Ensure barcodes match
      barcodes <- colnames(final_obj)

      # Name the clusters vector with cell barcodes
      if (!is.null(clusters)) {
        names(clusters) <- barcodes
      }

      # MINIMAL create_loupe call - fewer arguments = fewer crash points
      cat("    Calling loupeR::create_loupe()...\n")
      flush.console()

      loupe_result <- loupeR::create_loupe(
        count_mat = counts_mat,
        clusters = clusters,
        projections = if (length(projections) > 0) projections else NULL,
        output_dir = dirname(final_output_loupe_path),
        output_name = tools::file_path_sans_ext(basename(final_output_loupe_path)),
        force = TRUE
      )

      cat("    loupeR::create_loupe() returned\n")
      flush.console()

      # Check if file was created
      if (file.exists(final_output_loupe_path)) {
        cat("\n  Loupe file created successfully!\n")
        cat("    Size:", round(file.size(final_output_loupe_path) / 1e6, 2), "MB\n")
        flush.console()

        # Save structure summary
        loupe_str_output_file <- file.path(out_base, "final_output_object_loupe_structure.txt")

        tryCatch({
          sink(loupe_str_output_file)
          cat("Loupe File Structure for: final_output_object.cloupe\n")
          cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
          cat(paste(rep("=", 80), collapse = ""), "\n\n")
          cat("Cells:", ncol(counts_mat), "\n")
          cat("Genes:", nrow(counts_mat), "\n")
          if (!is.null(clusters)) {
            cat("Clusters:", length(levels(clusters)), "\n")
            cat("Cluster column:", cluster_col_used, "\n")
          }
          cat("Projections:", paste(names(projections), collapse = ", "), "\n")
          cat(paste(rep("=", 80), collapse = ""), "\n")
          sink()
          cat("  Loupe structure saved to:", loupe_str_output_file, "\n")
        }, error = function(e) {
          try(sink(), silent = TRUE)
          cat("  Warning: Could not save structure file\n")
        })

        return(TRUE)

      } else {
        # Check if file was created with different name
        loupe_files <- list.files(dirname(final_output_loupe_path),
                                   pattern = "\\.cloupe$", full.names = TRUE)
        if (length(loupe_files) > 0) {
          cat("    Found Loupe file(s):", paste(basename(loupe_files), collapse = ", "), "\n")
          if (length(loupe_files) == 1 && loupe_files[1] != final_output_loupe_path) {
            file.rename(loupe_files[1], final_output_loupe_path)
            cat("    Renamed to:", basename(final_output_loupe_path), "\n")
            cat("    Size:", round(file.size(final_output_loupe_path) / 1e6, 2), "MB\n")
            return(TRUE)
          }
        }
        cat("  Warning: Loupe file was not created\n")
        return(FALSE)
      }

    }, error = function(e) {
      # Make sure sink is closed if error occurred during sink
      try(sink(), silent = TRUE)

      cat("\n[WARNING] Loupe file conversion failed:\n")
      cat("  Error:", conditionMessage(e), "\n")
      cat("  This may be due to:\n")
      cat("    - loupeR external binary crash\n")
      cat("    - Missing loupeR dependencies\n")
      cat("    - Incompatible data format\n")
      cat("  The RDS and H5AD files have been saved successfully.\n")
      cat("  Set SKIP_LOUPE_CONVERSION <- TRUE to bypass this step.\n\n")
      flush.console()
      return(FALSE)
    })

  }  # end of if (!SKIP_LOUPE_CONVERSION)

  # --------------------------------------------------------------------------
  # Summary of saved files
  # --------------------------------------------------------------------------
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("FINAL OUTPUT FILES SUMMARY\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")

  cat("Output directory:", out_base, "\n\n")
  cat("Files saved:\n")
  cat("  1. final_clustered_object.rds (in objects/)\n")
  cat("     Path:", final_rds_path, "\n")
  cat("     Size:", round(file.size(final_rds_path) / 1e6, 2), "MB\n\n")

  cat("  2. final_output_object.rds (in root output directory)\n")
  cat("     Path:", final_output_rds_path, "\n")
  cat("     Size:", round(file.size(final_output_rds_path) / 1e6, 2), "MB\n\n")

  if (h5ad_conversion_success && file.exists(final_output_h5ad_path)) {
    cat("  3. final_output_object.h5ad (in root output directory)\n")
    cat("     Path:", final_output_h5ad_path, "\n")
    cat("     Size:", round(file.size(final_output_h5ad_path) / 1e6, 2), "MB\n\n")
  } else {
    cat("  3. final_output_object.h5ad - NOT CREATED (conversion failed)\n\n")
  }

  if (loupe_conversion_success && file.exists(final_output_loupe_path)) {
    cat("  4. final_output_object.cloupe (in root output directory)\n")
    cat("     Path:", final_output_loupe_path, "\n")
    cat("     Size:", round(file.size(final_output_loupe_path) / 1e6, 2), "MB\n\n")
  } else {
    cat("  4. final_output_object.cloupe - NOT CREATED (conversion failed or loupeR not installed)\n\n")
  }

  cat("Structure files:\n")
  cat("  - final_output_object_seurat_str.txt\n")
  if (h5ad_conversion_success) {
    cat("  - final_output_object_h5ad_structure.txt\n")
  }
  if (loupe_conversion_success) {
    cat("  - final_output_object_loupe_structure.txt\n")
  }
  cat("\n")
}

# ==============================================================================
# Output Directory Structure
# ==============================================================================
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("OUTPUT DIRECTORY STRUCTURE\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

cat("Output Root:", out_base, "\n\n")

# Count generated files
all_files <- list.files(out_base, recursive = TRUE, full.names = FALSE)
n_rds <- sum(grepl("\\.rds$", all_files, ignore.case = TRUE))
n_csv <- sum(grepl("\\.csv$", all_files, ignore.case = TRUE))
n_png <- sum(grepl("\\.png$", all_files, ignore.case = TRUE))
n_pdf <- sum(grepl("\\.pdf$", all_files, ignore.case = TRUE))
n_rdata <- sum(grepl("\\.RData$", all_files, ignore.case = TRUE))
n_h5ad <- sum(grepl("\\.h5ad$", all_files, ignore.case = TRUE))
n_cloupe <- sum(grepl("\\.cloupe$", all_files, ignore.case = TRUE))

cat("Generated files:\n")
cat("  Total files:", length(all_files), "\n")
cat("  RDS objects:", n_rds, "\n")
cat("  RData files:", n_rdata, "\n")
cat("  H5AD files:", n_h5ad, "\n")
cat("  Loupe files:", n_cloupe, "\n")
cat("  CSV tables:", n_csv, "\n")
cat("  PNG plots:", n_png, "\n")
cat("  PDF plots:", n_pdf, "\n")

# ==============================================================================
# Write Final README
# ==============================================================================

# Build dynamic README content based on detected metadata
readme_sex_info <- ""
if (!is.null(sex_col) && exists("clustered_obj") && !is.null(clustered_obj)) {
  sex_levels <- unique(clustered_obj@meta.data[[sex_col]])
  sex_levels <- sex_levels[!is.na(sex_levels)]
  for (sex_val in sex_levels) {
    count <- sum(clustered_obj@meta.data[[sex_col]] == sex_val, na.rm = TRUE)
    readme_sex_info <- paste0(readme_sex_info, sex_val, " cells: ", count, "\n")
  }
}

readme_content <- paste0(
  "================================================================================\n",
  "CHOROID PLEXUS scRNA-seq MULTI-SAMPLE ANALYSIS\n",
  "================================================================================\n\n",
  "Analysis Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
  "--------------------------------------------------------------------------------\n",
  "SAMPLES ANALYZED\n",
  "--------------------------------------------------------------------------------\n",
  samples_analyzed, "\n\n",
  "Batch variable: ", if (!is.null(params$batch_variable)) params$batch_variable else "batch", "\n",
  if (readme_sex_info != "") paste0("\n", readme_sex_info) else "",
  "\n--------------------------------------------------------------------------------\n",
  "METHODS SELECTED\n",
  "--------------------------------------------------------------------------------\n",
  "Normalization: ", if (exists("selected_normalization_method")) selected_normalization_method else "N/A", "\n",
  "Integration: ", if (exists("best_method")) best_method else "N/A", "\n",
  "Clustering: ", if (exists("clustering_method_used")) clustering_method_used else "N/A", "\n",
  "Batch Variable: ", if (!is.null(params$batch_variable)) params$batch_variable else "batch", "\n",
  "Batch Integration: ", as.character(isTRUE(params$run_batch_integration)), "\n\n",
  "--------------------------------------------------------------------------------\n",
  "FINAL OUTPUT FILES\n",
  "--------------------------------------------------------------------------------\n",
  "final_output_object.rds    - Seurat V5 object with all clustering/subclustering\n",
  "final_output_object.h5ad   - AnnData format for Python/scanpy compatibility\n",
  "final_output_object.cloupe - Loupe Browser format for 10x Genomics visualization\n",
  "final_output_object_seurat_str.txt - Full str() output of Seurat object\n",
  "final_output_object_h5ad_structure.txt - Full structure of H5AD object\n",
  "final_output_object_loupe_structure.txt - Summary of Loupe file contents\n\n",
  "--------------------------------------------------------------------------------\n",
  "DATA PROVENANCE\n",
  "--------------------------------------------------------------------------------\n",
  "Counts layer: scCDC-corrected counts (contamination-corrected from Step 7)\n",
  "Data layer: Normalized expression (", if (exists("selected_normalization_method")) selected_normalization_method else "N/A", ")\n",
  "Scale.data layer: Scaled expression for variable features only\n\n",
  "--------------------------------------------------------------------------------\n",
  "PIPELINE MODULES\n",
  "--------------------------------------------------------------------------------\n",
  "00. Environment Setup         - Package loading, directory creation\n",
  "01. Load Data                 - Import scCDC-corrected data\n",
  "02. QC Validation             - Quality control filtering\n",
  "03. Normalization             - SCT + scran + LogNormalize\n",
  "04. Integration               - Multi-method comparison\n",
  "05. CHOIR Clustering          - Hierarchical clustering\n",
  "06. scICE Subclustering       - Consistent subclustering\n",
  "07. Leiden Clustering         - Resolution testing\n",
  "08. Differential Expression   - MAST, edgeR, DESeq2\n",
  "09. Gene Visualization        - Expression plots\n",
  "10. Final Summary             - This summary\n",
  "11. HTML Report               - Interactive report\n\n",
  "--------------------------------------------------------------------------------\n",
  "FILE COUNTS\n",
  "--------------------------------------------------------------------------------\n",
  "Total files: ", length(all_files), "\n",
  "RDS objects: ", n_rds, "\n",
  "H5AD files: ", n_h5ad, "\n",
  "Loupe files: ", n_cloupe, "\n",
  "CSV tables: ", n_csv, "\n",
  "PNG plots: ", n_png, "\n",
  "PDF plots: ", n_pdf, "\n\n",
  "================================================================================\n"
)

writeLines(readme_content, file.path(out_base, "README.txt"))
cat("\nFinal README written to:", file.path(out_base, "README.txt"), "\n")

# ==============================================================================
# Session Information
# ==============================================================================
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("SESSION INFORMATION\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

session_file <- file.path(output_dirs$tables, "session_info.txt")
sink(session_file)
print(sessionInfo())
sink()
cat("Session info saved to:", session_file, "\n")

cat("\nKey package versions:\n")
cat("  Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("  SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
cat("  scCustomize:", as.character(packageVersion("scCustomize")), "\n")
cat("  ggplot2:", as.character(packageVersion("ggplot2")), "\n")
cat("  dplyr:", as.character(packageVersion("dplyr")), "\n")

# Check if loupeR is available and print version
if (requireNamespace("loupeR", quietly = TRUE)) {
  cat("  loupeR:", as.character(packageVersion("loupeR")), "\n")
} else {
  cat("  loupeR: not installed\n")
}

# ==============================================================================
# Final Message
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Output directory:", out_base, "\n")
cat("Total files generated:", length(all_files), "\n")

cat("\n>>> MODULE 10 COMPLETE <<<\n")
