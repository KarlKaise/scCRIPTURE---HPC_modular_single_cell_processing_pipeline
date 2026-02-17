#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 03: NORMALIZATION (MULTI-SAMPLE PIPELINE)
# ==============================================================================
#
# This module performs four normalization methods:
# 1. SCTransform (variance stabilization)
# 2. scran (deconvolution-based size factors)
# 3. LogNormalize (simple log normalization)
# 4. scKWARN (kernel weighted adjusted regularized normalization)
#
# All methods use Harmony integration for fair comparison when batch
# integration is enabled.
#
# INPUT: QC-filtered merged Seurat object from Module 02
# OUTPUT: Normalized objects with optional Harmony integration
#
# UPDATES:
# - 2026-01-15: Added scKWARN normalization method
# - 2026-01-15: scKWARN uses log1p() as recommended by authors
# - 2026-02-05: Fixed benchmark subsampling to use vars.to.regress in ScaleData
# - 2026-02-05: Strip pre-existing integration reductions before benchmarking
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 03: NORMALIZATION\n")
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
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scran)
  library(scuttle)
  library(reticulate) 
  library(SingleCellExperiment)
  library(tidyr)
})

# Check for scKWARN package
has_sckwarn <- requireNamespace("scKWARN", quietly = TRUE)
if (has_sckwarn) {
  suppressPackageStartupMessages(library(scKWARN))
  cat("scKWARN package: AVAILABLE\n")
} else {
  cat("scKWARN package: NOT AVAILABLE\n")
  cat("  Install with: devtools::install_github('cyhsuTN/scKWARN')\n")
}

out_base <- params$out_root
load(file.path(out_base, "objects", "pipeline_environment.RData"))
if (!exists("has_harmony")) has_harmony <- requireNamespace("harmony", quietly = TRUE)
if (!exists("has_glmGamPoi")) has_glmGamPoi <- requireNamespace("glmGamPoi", quietly = TRUE)
load(file.path(out_base, "objects", "02_qc_data.RData"))

if (has_harmony) library(harmony)

# Check for clustering metrics packages
has_mclust <- requireNamespace("mclust", quietly = TRUE)
has_aricode <- requireNamespace("aricode", quietly = TRUE)
cat("Clustering metrics packages:\n")
cat("  mclust (for ARI):", has_mclust, "\n")
cat("  aricode (for NMI):", has_aricode, "\n\n")

# --- Detect cell type annotations for biological conservation metrics ---
cat("--- Detecting cell type annotations for normalization benchmarking ---\n")

norm_bench_celltype_column <- NULL
norm_bench_has_celltype <- FALSE

# Priority 1: Explicit param
if (!is.null(params$celltype_column) && nchar(params$celltype_column) > 0) {
  norm_bench_celltype_column <- params$celltype_column
  cat("  Celltype column from params:", norm_bench_celltype_column, "\n")
}

# Priority 2: MapMyCells columns (preferred granularity order)
if (is.null(norm_bench_celltype_column)) {
  mapmycells_candidates <- c(
    "cluster_name_MapMyCells",       # ~22 types, best for silhouette/LISI
    "supercluster_name_MapMyCells",  # ~8 types, coarser
    "subcluster_name_MapMyCells"     # very fine, may fragment rare types
  )
  # We need to check against one of the sample objects since the merged
  # benchmarking objects inherit the same metadata columns
  check_obj <- sample_objects[[1]]
  for (col in mapmycells_candidates) {
    if (col %in% colnames(check_obj@meta.data)) {
      col_values <- check_obj@meta.data[[col]]
      n_unique <- length(unique(na.omit(col_values)))
      if (n_unique > 1) {
        norm_bench_celltype_column <- col
        cat("  Auto-detected MapMyCells column:", col, "(", n_unique, "types in first sample)\n")
        break
      }
    }
  }
}

# Priority 3: Common column names (same as Module 04)
if (is.null(norm_bench_celltype_column)) {
  check_obj <- sample_objects[[1]]
  common_cols <- c("cell_type", "celltype", "CellType", "cell_type_annotation",
                   "cluster_annotation", "annotation", "predicted_celltype",
                   "sctype_classification")
  for (col in common_cols) {
    if (col %in% colnames(check_obj@meta.data)) {
      n_unique <- length(unique(na.omit(check_obj@meta.data[[col]])))
      if (n_unique > 1) {
        norm_bench_celltype_column <- col
        cat("  Auto-detected celltype column:", col, "\n")
        break
      }
    }
  }
}

if (!is.null(norm_bench_celltype_column)) {
  norm_bench_has_celltype <- TRUE
  cat("  >>> Biological conservation metrics ENABLED for normalization benchmarking\n")
  cat("  >>> Using column:", norm_bench_celltype_column, "\n\n")
} else {
  cat("  No cell type annotations found.\n")
  cat("  Normalization benchmarking will use batch metrics only.\n")
  cat("  To enable: set params$celltype_column or add MapMyCells annotations.\n\n")
}

cat("Batch variable for integration:", params$batch_variable, "\n")
cat("Run batch integration:", params$run_batch_integration, "\n")
if (isTRUE(params$run_batch_integration)) {
  cat("All normalization methods will use Harmony for fair comparison.\n\n")
} else {
  cat("Batch integration disabled. Will merge without batch correction.\n\n")
}

sct_object <- NULL
scran_object <- NULL
lognorm_object <- NULL
sckwarn_object <- NULL

# Create directories for individual and merged normalized objects
individual_norm_dir <- file.path(output_dirs$objects, "individual_normalized")
dir.create(individual_norm_dir, showWarnings = FALSE, recursive = TRUE)

merged_norm_dir <- file.path(output_dirs$objects, "merged_normalized")
dir.create(merged_norm_dir, showWarnings = FALSE, recursive = TRUE)

# Create LogNormalize output directory
subdirs$norm_lognorm <- file.path(output_dirs$normalization, "LogNormalize")
dir.create(subdirs$norm_lognorm, showWarnings = FALSE, recursive = TRUE)

# Create scKWARN output directory
subdirs$norm_sckwarn <- file.path(output_dirs$normalization, "scKWARN")
dir.create(subdirs$norm_sckwarn, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Get sample names
# ==============================================================================
sample_names <- names(sample_objects)
cat("Samples to normalize:", paste(sample_names, collapse = ", "), "\n\n")

# ==============================================================================
# SCTRANSFORM NORMALIZATION
# ==============================================================================
if (isTRUE(params$run_sctransform)) {
  cat("\n--- SCTransform Normalization ---\n")

  sct_method_use <- if (has_glmGamPoi) "glmGamPoi" else "poisson"
  cat("Using method:", sct_method_use, "\n")

  # SCTransform each sample
  sct_samples <- list()

  for (samp in sample_names) {
    cat("\nSCTransform:", samp, "\n")
    sct_samples[[samp]] <- SCTransform(sample_objects[[samp]],
                                        vars.to.regress = params$sct_vars_to_regress,
                                        method = sct_method_use,
                                        verbose = TRUE)

    # Save individual normalized object
    sct_path <- file.path(individual_norm_dir, paste0(samp, "_SCTransform.rds"))
    saveRDS(sct_samples[[samp]], sct_path)
    cat("  Saved:", sct_path, "\n")
  }

  # Merge SCT objects
  cat("\nMerging SCTransform objects...\n")

  if (length(sct_samples) == 1) {
    merged_sct <- sct_samples[[1]]
  } else {
    merged_sct <- merge(sct_samples[[1]], sct_samples[2:length(sct_samples)],
                        add.cell.ids = sample_names)
  }

  cat("Merged SCT object:", ncol(merged_sct), "cells\n")

  tryCatch({
    merged_sct[["SCT"]] <- JoinLayers(merged_sct[["SCT"]])
  }, error = function(e) NULL)

  DefaultAssay(merged_sct) <- "SCT"
  merged_sct <- FindVariableFeatures(merged_sct, nfeatures = params$nfeatures_integration, verbose = FALSE)
  merged_sct <- ScaleData(merged_sct, verbose = FALSE)
  merged_sct <- RunPCA(merged_sct, npcs = 50, verbose = FALSE)
  merged_sct <- RunUMAP(merged_sct, reduction = "pca", dims = 1:params$dims_use,
                        reduction.name = "umap_unintegrated", verbose = FALSE)

  # Save merged normalized object BEFORE integration
  cat("\nSaving merged SCTransform object (before integration)...\n")
  merged_sct_unintegrated_path <- file.path(merged_norm_dir, "merged_SCTransform_unintegrated.rds")
  saveRDS(merged_sct, merged_sct_unintegrated_path)
  cat("  Saved:", merged_sct_unintegrated_path, "\n")

  p_sct_unintegrated <- DimPlot(merged_sct, reduction = "umap_unintegrated", group.by = "sex",
                                 cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
    ggtitle("SCT Merged (Unintegrated)")
  save_plot_multi(p_sct_unintegrated, "00_SCT_Merged_Unintegrated_UMAP",
                  output_dir = subdirs$norm_sct, width = 7, height = 6)

  # Conditional integration based on run_batch_integration parameter
  if (isTRUE(params$run_batch_integration) && has_harmony) {
    cat("\nRunning Harmony integration on SCTransform data...\n")
    merged_sct <- RunHarmony(merged_sct, group.by.vars = params$batch_variable,
                             dims.use = 1:params$dims_use, assay.use = "SCT")
    merged_sct <- RunUMAP(merged_sct, reduction = "harmony", dims = 1:params$dims_use,
                          reduction.name = "umap", verbose = FALSE)
    sct_integration_method <- "Harmony"
  } else {
    cat("\nSkipping batch integration (run_batch_integration = FALSE or Harmony not available)...\n")
    merged_sct[["umap"]] <- merged_sct[["umap_unintegrated"]]
    sct_integration_method <- "None (PCA only)"
  }

  cat("\n>>> SCT INTEGRATED OBJECT <<<\n")
  print_object_structure(merged_sct, paste0("SCT + ", sct_integration_method))

  # Visualization
  p_sct_sex <- DimPlot(merged_sct, reduction = "umap", group.by = "sex",
                       cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) + ggtitle("By Sex")
  p_sct_batch <- DimPlot(merged_sct, reduction = "umap", group.by = "batch") + ggtitle("By Batch")
  p_sct_combined <- p_sct_sex + p_sct_batch

  save_plot_multi(p_sct_combined, "01_SCT_Harmony_UMAP", output_dir = subdirs$norm_sct, width = 14, height = 6)

  sct_rds_path <- file.path(output_dirs$objects, "integrated_SCTransform_Harmony.rds")
  saveRDS(merged_sct, sct_rds_path)
  cat("\nSaved:", sct_rds_path, "\n")

  sct_object <- merged_sct

  write_readme(subdirs$norm_sct, "SCTransform Normalization",
               paste0("Method: ", sct_method_use, "\nVars regressed: ",
                      paste(params$sct_vars_to_regress, collapse = ", "), "\n",
                      "Features: ", params$nfeatures_integration, "\n",
                      "Integration: ", sct_integration_method, "\n",
                      "Batch variable: ", params$batch_variable),
               list("00_SCT_Merged_Unintegrated_UMAP.*" = "UMAP before integration",
                    "01_SCT_Harmony_UMAP.*" = "UMAP after integration"))
}

# ==============================================================================
# SCRAN NORMALIZATION
# ==============================================================================
if (isTRUE(params$run_scran)) {
  cat("\n--- scran Normalization ---\n")

  compute_scran_sf <- function(sce, sample_name) {
    cat("\nSize factors:", sample_name, "\n")
    set.seed(42)
    quick_clust <- quickCluster(sce, min.mean = params$scran_min_mean)
    sce <- computeSumFactors(sce, clusters = quick_clust, min.mean = params$scran_min_mean)
    sce <- logNormCounts(sce)
    cat("Size factor range:", round(min(sizeFactors(sce)), 3), "-", round(max(sizeFactors(sce)), 3), "\n")
    return(sce)
  }

  sce_to_seurat <- function(sce, original_seurat) {
    obj <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
    obj <- SetAssayData(obj, layer = "data", new.data = logcounts(sce))
    for (col in setdiff(colnames(original_seurat@meta.data), colnames(obj@meta.data))) {
      obj[[col]] <- original_seurat@meta.data[[col]]
    }
    obj$scran_size_factor <- sizeFactors(sce)
    return(obj)
  }

  scran_samples <- list()

  for (samp in sample_names) {
    sce <- seurat_to_sce(sample_objects[[samp]])
    sce <- compute_scran_sf(sce, samp)
    scran_samples[[samp]] <- sce_to_seurat(sce, sample_objects[[samp]])

    # Save individual normalized object
    scran_path <- file.path(individual_norm_dir, paste0(samp, "_scran.rds"))
    saveRDS(scran_samples[[samp]], scran_path)
    cat("  Saved:", scran_path, "\n")
  }

  # Merge scran objects
  cat("\nMerging scran objects...\n")

  if (length(scran_samples) == 1) {
    merged_scran <- scran_samples[[1]]
  } else {
    merged_scran <- merge(scran_samples[[1]], scran_samples[2:length(scran_samples)],
                          add.cell.ids = sample_names)
  }

  merged_scran <- FindVariableFeatures(merged_scran, nfeatures = params$nfeatures_integration)
  merged_scran <- ScaleData(merged_scran, vars.to.regress = params$sct_vars_to_regress)
  merged_scran <- RunPCA(merged_scran, npcs = 50, verbose = FALSE)
  merged_scran <- RunUMAP(merged_scran, reduction = "pca", dims = 1:params$dims_use,
                          reduction.name = "umap_unintegrated", verbose = FALSE)

  # Save merged normalized object BEFORE integration
  cat("\nSaving merged scran object (before integration)...\n")
  merged_scran_unintegrated_path <- file.path(merged_norm_dir, "merged_scran_unintegrated.rds")
  saveRDS(merged_scran, merged_scran_unintegrated_path)
  cat("  Saved:", merged_scran_unintegrated_path, "\n")

  p_scran_unintegrated <- DimPlot(merged_scran, reduction = "umap_unintegrated", group.by = "sex",
                                   cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
    ggtitle("scran Merged (Unintegrated)")
  save_plot_multi(p_scran_unintegrated, "00_scran_Merged_Unintegrated_UMAP",
                  output_dir = subdirs$norm_scran, width = 7, height = 6)

  # Conditional integration
  if (isTRUE(params$run_batch_integration) && has_harmony) {
    cat("\nRunning Harmony integration on scran data...\n")
    merged_scran <- RunHarmony(merged_scran, group.by.vars = params$batch_variable,
                               dims.use = 1:params$dims_use)
    merged_scran <- RunUMAP(merged_scran, reduction = "harmony", dims = 1:params$dims_use,
                            reduction.name = "umap", verbose = FALSE)
    scran_integration_method <- "Harmony"
  } else {
    cat("\nSkipping batch integration...\n")
    merged_scran[["umap"]] <- merged_scran[["umap_unintegrated"]]
    scran_integration_method <- "None (PCA only)"
  }

  cat("\n>>> SCRAN INTEGRATED OBJECT <<<\n")
  print_object_structure(merged_scran, paste0("scran + ", scran_integration_method))

  p_scran <- DimPlot(merged_scran, reduction = "umap", group.by = "sex",
                     cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
    ggtitle(paste0("scran + ", scran_integration_method))
  save_plot_multi(p_scran, "01_scran_Harmony_UMAP", output_dir = subdirs$norm_scran)

  scran_rds_path <- file.path(output_dirs$objects, "integrated_scran_Harmony.rds")
  saveRDS(merged_scran, scran_rds_path)
  cat("\nSaved:", scran_rds_path, "\n")

  scran_object <- merged_scran

  write_readme(subdirs$norm_scran, "scran Normalization",
               paste0("scran normalization with deconvolution size factors\n",
                      "Integration: ", scran_integration_method, "\n",
                      "Batch variable: ", params$batch_variable),
               list("00_scran_Merged_Unintegrated_UMAP.*" = "UMAP before integration",
                    "01_scran_Harmony_UMAP.*" = "UMAP after integration"))
}

# ==============================================================================
# LOGNORMALIZE NORMALIZATION
# ==============================================================================
if (isTRUE(params$run_lognorm)) {
  cat("\n--- LogNormalize Normalization ---\n")

  lognorm_samples <- list()

  for (samp in sample_names) {
    cat("\nLogNormalize:", samp, "\n")
    lognorm_samples[[samp]] <- NormalizeData(sample_objects[[samp]],
                                              normalization.method = "LogNormalize",
                                              scale.factor = 10000,
                                              verbose = TRUE)

    # Save individual normalized object
    lognorm_path <- file.path(individual_norm_dir, paste0(samp, "_LogNormalize.rds"))
    saveRDS(lognorm_samples[[samp]], lognorm_path)
    cat("  Saved:", lognorm_path, "\n")
  }

  # Merge LogNormalize objects
  cat("\nMerging LogNormalize objects...\n")

  if (length(lognorm_samples) == 1) {
    merged_lognorm <- lognorm_samples[[1]]
  } else {
    merged_lognorm <- merge(lognorm_samples[[1]], lognorm_samples[2:length(lognorm_samples)],
                            add.cell.ids = sample_names)
  }

  cat("Merged LogNormalize object:", ncol(merged_lognorm), "cells\n")

  tryCatch({
    merged_lognorm[["RNA"]] <- JoinLayers(merged_lognorm[["RNA"]])
  }, error = function(e) NULL)

  merged_lognorm <- FindVariableFeatures(merged_lognorm, nfeatures = params$nfeatures_integration, verbose = FALSE)
  merged_lognorm <- ScaleData(merged_lognorm, vars.to.regress = params$sct_vars_to_regress, verbose = FALSE)
  merged_lognorm <- RunPCA(merged_lognorm, npcs = 50, verbose = FALSE)
  merged_lognorm <- RunUMAP(merged_lognorm, reduction = "pca", dims = 1:params$dims_use,
                            reduction.name = "umap_unintegrated", verbose = FALSE)

  # Save merged normalized object BEFORE integration
  cat("\nSaving merged LogNormalize object (before integration)...\n")
  merged_lognorm_unintegrated_path <- file.path(merged_norm_dir, "merged_LogNormalize_unintegrated.rds")
  saveRDS(merged_lognorm, merged_lognorm_unintegrated_path)
  cat("  Saved:", merged_lognorm_unintegrated_path, "\n")

  p_lognorm_unintegrated <- DimPlot(merged_lognorm, reduction = "umap_unintegrated", group.by = "sex",
                                     cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
    ggtitle("LogNormalize Merged (Unintegrated)")
  save_plot_multi(p_lognorm_unintegrated, "00_LogNorm_Merged_Unintegrated_UMAP",
                  output_dir = subdirs$norm_lognorm, width = 7, height = 6)

  # Conditional integration
  if (isTRUE(params$run_batch_integration) && has_harmony) {
    cat("\nRunning Harmony integration on LogNormalize data...\n")
    merged_lognorm <- RunHarmony(merged_lognorm, group.by.vars = params$batch_variable,
                                  dims.use = 1:params$dims_use)
    merged_lognorm <- RunUMAP(merged_lognorm, reduction = "harmony", dims = 1:params$dims_use,
                              reduction.name = "umap", verbose = FALSE)
    lognorm_integration_method <- "Harmony"
  } else {
    cat("\nSkipping batch integration...\n")
    merged_lognorm[["umap"]] <- merged_lognorm[["umap_unintegrated"]]
    lognorm_integration_method <- "None (PCA only)"
  }

  cat("\n>>> LOGNORMALIZE INTEGRATED OBJECT <<<\n")
  print_object_structure(merged_lognorm, paste0("LogNormalize + ", lognorm_integration_method))

  # Visualization
  p_lognorm_sex <- DimPlot(merged_lognorm, reduction = "umap", group.by = "sex",
                            cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) + ggtitle("By Sex")
  p_lognorm_batch <- DimPlot(merged_lognorm, reduction = "umap", group.by = "batch") + ggtitle("By Batch")
  p_lognorm_combined <- p_lognorm_sex + p_lognorm_batch

  save_plot_multi(p_lognorm_combined, "01_LogNorm_Harmony_UMAP", output_dir = subdirs$norm_lognorm, width = 14, height = 6)

  lognorm_rds_path <- file.path(output_dirs$objects, "integrated_LogNormalize_Harmony.rds")
  saveRDS(merged_lognorm, lognorm_rds_path)
  cat("\nSaved:", lognorm_rds_path, "\n")

  lognorm_object <- merged_lognorm

  write_readme(subdirs$norm_lognorm, "LogNormalize Normalization",
               paste0("Standard log normalization (log1p(x/total * scale_factor))\n",
                      "Scale factor: 10000\n",
                      "Vars regressed in ScaleData: ", paste(params$sct_vars_to_regress, collapse = ", "), "\n",
                      "Features: ", params$nfeatures_integration, "\n",
                      "Integration: ", lognorm_integration_method, "\n",
                      "Batch variable: ", params$batch_variable),
               list("00_LogNorm_Merged_Unintegrated_UMAP.*" = "UMAP before integration",
                    "01_LogNorm_Harmony_UMAP.*" = "UMAP after integration"))
}

# ==============================================================================
# SCKWARN NORMALIZATION
# ==============================================================================
# scKWARN: Kernel Weighted Adjusted Regularized Normalization
# Reference: https://github.com/cyhsuTN/scKWARN
#
# The LocASN function computes adaptive size factors and returns:
#   - NormalizedData: counts / scaling_factor (linear scale)
#   - scalingFactor: per-cell normalization factors
#
# Per scKWARN documentation, log1p() should be applied for clustering:
#   "log1p(Result$NormalizedData) for clustering"
#
# This is consistent with Seurat's LogNormalize which also uses log1p().
# ==============================================================================

run_sckwarn_param <- if (!is.null(params$run_sckwarn)) params$run_sckwarn else TRUE

if (isTRUE(run_sckwarn_param) && has_sckwarn) {
  cat("\n--- scKWARN Normalization ---\n")
  cat("Using LocASN (Local Adaptive Size Normalization)\n")
  cat("Log transformation: log1p(x) as recommended by scKWARN authors\n")
  cat("This is consistent with Seurat's LogNormalize (natural log scale)\n\n")

  sckwarn_samples <- list()

  for (samp in sample_names) {
    cat("\nscKWARN:", samp, "\n")

    # Get counts matrix from sample object
    counts_mat <- tryCatch({
      GetAssayData(sample_objects[[samp]], assay = "RNA", layer = "counts")
    }, error = function(e) {
      tryCatch({
        GetAssayData(sample_objects[[samp]], assay = "RNA", slot = "counts")
      }, error = function(e2) {
        stop("Could not extract counts from sample ", samp, ": ", e2$message)
      })
    })

    cat("  Input counts matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "cells\n")

    # Run scKWARN LocASN normalization
    cat("  Running LocASN...\n")
    sckwarn_result <- tryCatch({
      scKWARN::LocASN(countmatrix = counts_mat)
    }, error = function(e) {
      cat("  ERROR in LocASN:", e$message, "\n")
      NULL
    })

    if (is.null(sckwarn_result)) {
      cat("  SKIPPING sample due to LocASN error\n")
      next
    }

    # Get normalized data and scaling factors
    norm_data <- sckwarn_result$NormalizedData
    scaling_factors <- sckwarn_result$scalingFactor

    cat("  Scaling factor range:", round(min(scaling_factors), 4), "-", round(max(scaling_factors), 4), "\n")
    cat("  Scaling factor mean:", round(mean(scaling_factors), 4), "\n")

    # Apply log1p() transformation as recommended by scKWARN authors
    # From documentation: "log1p(Result$NormalizedData) for clustering"
    # This is consistent with Seurat's LogNormalize (natural log scale)
    log_norm_data <- log1p(norm_data)

    cat("  Applied log1p(x) transformation (natural log, as per scKWARN documentation)\n")
    cat("  Log-normalized data range:", round(min(log_norm_data), 4), "-", round(max(log_norm_data), 4), "\n")

    # Create Seurat object with original counts
    sckwarn_obj <- CreateSeuratObject(
      counts = counts_mat,
      meta.data = sample_objects[[samp]]@meta.data
    )

    # Set the log-transformed normalized data as the "data" layer
    tryCatch({
      # Seurat v5 method
      sckwarn_obj[["RNA"]]$data <- log_norm_data
    }, error = function(e) {
      # Fallback for older Seurat versions
      tryCatch({
        sckwarn_obj <- SetAssayData(sckwarn_obj, layer = "data", new.data = log_norm_data)
      }, error = function(e2) {
        sckwarn_obj <- SetAssayData(sckwarn_obj, slot = "data", new.data = log_norm_data)
      })
    })

    # Store scaling factor in metadata
    sckwarn_obj$sckwarn_scaling_factor <- scaling_factors

    # Also store the linear normalized data as a separate assay for future use
    # This preserves the original scKWARN output without log transformation
    sckwarn_obj[["scKWARN_linear"]] <- CreateAssayObject(data = norm_data)

    sckwarn_samples[[samp]] <- sckwarn_obj

    # Save individual normalized object
    sckwarn_path <- file.path(individual_norm_dir, paste0(samp, "_scKWARN.rds"))
    saveRDS(sckwarn_obj, sckwarn_path)
    cat("  Saved:", sckwarn_path, "\n")
  }

  # Check if we have any successful samples
  if (length(sckwarn_samples) == 0) {
    cat("\nWARNING: No samples successfully processed with scKWARN. Skipping.\n")
  } else {
    # Merge scKWARN objects
    cat("\nMerging scKWARN objects...\n")

    if (length(sckwarn_samples) == 1) {
      merged_sckwarn <- sckwarn_samples[[1]]
    } else {
      merged_sckwarn <- merge(sckwarn_samples[[1]], sckwarn_samples[2:length(sckwarn_samples)],
                              add.cell.ids = names(sckwarn_samples))
    }

    cat("Merged scKWARN object:", ncol(merged_sckwarn), "cells\n")

    # Join layers if needed (Seurat v5)
    tryCatch({
      merged_sckwarn[["RNA"]] <- JoinLayers(merged_sckwarn[["RNA"]])
    }, error = function(e) NULL)

    # Standard preprocessing pipeline
    DefaultAssay(merged_sckwarn) <- "RNA"
    merged_sckwarn <- FindVariableFeatures(merged_sckwarn, nfeatures = params$nfeatures_integration, verbose = FALSE)
    merged_sckwarn <- ScaleData(merged_sckwarn, vars.to.regress = params$sct_vars_to_regress, verbose = FALSE)
    merged_sckwarn <- RunPCA(merged_sckwarn, npcs = 50, verbose = FALSE)
    merged_sckwarn <- RunUMAP(merged_sckwarn, reduction = "pca", dims = 1:params$dims_use,
                              reduction.name = "umap_unintegrated", verbose = FALSE)

    # Save merged normalized object BEFORE integration
    cat("\nSaving merged scKWARN object (before integration)...\n")
    merged_sckwarn_unintegrated_path <- file.path(merged_norm_dir, "merged_scKWARN_unintegrated.rds")
    saveRDS(merged_sckwarn, merged_sckwarn_unintegrated_path)
    cat("  Saved:", merged_sckwarn_unintegrated_path, "\n")

    p_sckwarn_unintegrated <- DimPlot(merged_sckwarn, reduction = "umap_unintegrated", group.by = "sex",
                                       cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
      ggtitle("scKWARN Merged (Unintegrated)")
    save_plot_multi(p_sckwarn_unintegrated, "00_scKWARN_Merged_Unintegrated_UMAP",
                    output_dir = subdirs$norm_sckwarn, width = 7, height = 6)

    # Conditional integration
    if (isTRUE(params$run_batch_integration) && has_harmony) {
      cat("\nRunning Harmony integration on scKWARN data...\n")
      merged_sckwarn <- RunHarmony(merged_sckwarn, group.by.vars = params$batch_variable,
                                    dims.use = 1:params$dims_use)
      merged_sckwarn <- RunUMAP(merged_sckwarn, reduction = "harmony", dims = 1:params$dims_use,
                                reduction.name = "umap", verbose = FALSE)
      sckwarn_integration_method <- "Harmony"
    } else {
      cat("\nSkipping batch integration...\n")
      merged_sckwarn[["umap"]] <- merged_sckwarn[["umap_unintegrated"]]
      sckwarn_integration_method <- "None (PCA only)"
    }

    cat("\n>>> SCKWARN INTEGRATED OBJECT <<<\n")
    print_object_structure(merged_sckwarn, paste0("scKWARN + ", sckwarn_integration_method))

    # Visualization
    p_sckwarn_sex <- DimPlot(merged_sckwarn, reduction = "umap", group.by = "sex",
                              cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) + ggtitle("By Sex")
    p_sckwarn_batch <- DimPlot(merged_sckwarn, reduction = "umap", group.by = "batch") + ggtitle("By Batch")
    p_sckwarn_combined <- p_sckwarn_sex + p_sckwarn_batch

    save_plot_multi(p_sckwarn_combined, "01_scKWARN_Harmony_UMAP", output_dir = subdirs$norm_sckwarn, width = 14, height = 6)

    sckwarn_rds_path <- file.path(output_dirs$objects, "integrated_scKWARN_Harmony.rds")
    saveRDS(merged_sckwarn, sckwarn_rds_path)
    cat("\nSaved:", sckwarn_rds_path, "\n")

    sckwarn_object <- merged_sckwarn

    write_readme(subdirs$norm_sckwarn, "scKWARN Normalization",
                 paste0("scKWARN: Kernel Weighted Adjusted Regularized Normalization\n",
                        "Method: LocASN (Local Adaptive Size Normalization)\n",
                        "Reference: https://github.com/cyhsuTN/scKWARN\n\n",
                        "Log transformation: log1p(x) = ln(x + 1)\n",
                        "This follows the scKWARN documentation recommendation:\n",
                        "  'log1p(Result$NormalizedData) for clustering'\n",
                        "And is consistent with Seurat's LogNormalize (natural log scale).\n\n",
                        "Assay contents:\n",
                        "  - 'RNA' data layer: log1p(NormalizedData) - for downstream analysis\n",
                        "  - 'scKWARN_linear' assay: linear normalized data (no log)\n",
                        "  - 'sckwarn_scaling_factor' in metadata: per-cell scaling factors\n\n",
                        "Vars regressed in ScaleData: ", paste(params$sct_vars_to_regress, collapse = ", "), "\n",
                        "Features: ", params$nfeatures_integration, "\n",
                        "Integration: ", sckwarn_integration_method, "\n",
                        "Batch variable: ", params$batch_variable),
                 list("00_scKWARN_Merged_Unintegrated_UMAP.*" = "UMAP before integration",
                      "01_scKWARN_Harmony_UMAP.*" = "UMAP after integration"))
  }
} else if (isTRUE(run_sckwarn_param) && !has_sckwarn) {
  cat("\n--- scKWARN Normalization ---\n")
  cat("SKIPPED: scKWARN package not available\n")
  cat("Install with: devtools::install_github('cyhsuTN/scKWARN')\n")
}

# ==============================================================================
# NORMALIZATION BENCHMARKING
# ==============================================================================
methods_available <- sum(c(!is.null(sct_object), !is.null(scran_object),
                           !is.null(lognorm_object), !is.null(sckwarn_object)))

if (isTRUE(params$run_normalization_benchmarking) && methods_available >= 2) {
  cat("\n--- Multi-Integration Normalization Benchmarking ---\n")
  cat("Using multiple integration methods to select best normalization via majority vote\n\n")

  # Determine which integration methods to use for benchmarking
  bench_int_methods <- params$norm_benchmark_integration_methods
  if (is.null(bench_int_methods) || length(bench_int_methods) == 0) {
    bench_int_methods <- c("harmony")
    cat("NOTE: norm_benchmark_integration_methods not set. Using Harmony only.\n")
  }

  # Check R-based method availability
  if (!exists("has_batchelor")) has_batchelor <- requireNamespace("batchelor", quietly = TRUE)
  if (!exists("has_SeuratWrappers")) has_SeuratWrappers <- requireNamespace("SeuratWrappers", quietly = TRUE)

  available_bench_methods <- c()
  if ("harmony" %in% bench_int_methods && has_harmony) {
    available_bench_methods <- c(available_bench_methods, "harmony")
    cat("  [OK] Harmony available for benchmarking\n")
  }
  if ("mnn" %in% bench_int_methods && has_batchelor && has_SeuratWrappers) {
    available_bench_methods <- c(available_bench_methods, "mnn")
    cat("  [OK] FastMNN available for benchmarking\n")
  }

  # RPCA and CCA use Seurat::IntegrateLayers() -- always available (Seurat built-in)
  # FIX 2026-02-13: These were specified in params but never added to
  # available_bench_methods, so they were silently skipped.
  if ("rpca" %in% bench_int_methods) {
    available_bench_methods <- c(available_bench_methods, "rpca")
    cat("  [OK] RPCA available for benchmarking (Seurat built-in)\n")
  }
  if ("cca" %in% bench_int_methods) {
    available_bench_methods <- c(available_bench_methods, "cca")
    cat("  [OK] CCA available for benchmarking (Seurat built-in)\n")
  }
  
  # Check Python-based method availability (only if Python methods requested)
  python_bench_methods <- intersect(bench_int_methods, c("scvi", "sccobra", "concord"))
  if (length(python_bench_methods) > 0) {
    cat("\nSetting up Python environment for normalization benchmarking...\n")
    suppressPackageStartupMessages(library(reticulate))
    Sys.setenv(RETICULATE_PYTHON = params$unified_python)
    use_python(params$unified_python, required = TRUE)

    tryCatch({
      np <- import("numpy")
      pd <- import("pandas")
      anndata <- import("anndata")
      sc <- import("scanpy")
      cat("  [OK] Core Python modules (numpy, pandas, anndata, scanpy)\n")
    }, error = function(e) {
      cat("  [WARNING] Failed to import core Python modules:", e$message, "\n")
      cat("  Python-based benchmarking methods will be skipped.\n")
      python_bench_methods <<- c()
    })

    if ("scvi" %in% python_bench_methods) {
      has_scvi_bench <- tryCatch({
        import("scvi"); cat("  [OK] scvi-tools\n"); TRUE
      }, error = function(e) { cat("  [WARNING] scvi not available\n"); FALSE })
      if (has_scvi_bench) available_bench_methods <- c(available_bench_methods, "scvi")
    }

    if ("sccobra" %in% python_bench_methods) {
      has_sccobra_bench <- tryCatch({
        sccobra_path <- "/scicore/home/doetsch/kaiser0001/GITHUB_repositories/scCobra"
        if (dir.exists(sccobra_path)) {
          py_run_string(sprintf("import sys; sys.path.insert(0, '%s')", sccobra_path))
          py_run_string("from scCobra import scCobra as _test_sccobra_bm")
          cat("  [OK] scCobra\n")
          TRUE
        } else {
          cat("  [WARNING] scCobra path not found:", sccobra_path, "\n")
          FALSE
        }
      }, error = function(e) { cat("  [WARNING] scCobra not available:", e$message, "\n"); FALSE })
      if (has_sccobra_bench) available_bench_methods <- c(available_bench_methods, "sccobra")
    }

    if ("concord" %in% python_bench_methods) {
      has_concord_bench <- tryCatch({
        py_run_string("import concord"); cat("  [OK] concord\n"); TRUE
      }, error = function(e) { cat("  [WARNING] concord not available:", e$message, "\n"); FALSE })
      if (has_concord_bench) available_bench_methods <- c(available_bench_methods, "concord")
    }
  }

  cat("\nAvailable benchmarking integration methods:", paste(available_bench_methods, collapse = ", "), "\n")

  # Build list of normalization objects to benchmark
  norm_objects <- list()
  if (!is.null(sct_object)) norm_objects[["SCTransform"]] <- sct_object
  if (!is.null(scran_object)) norm_objects[["scran"]] <- scran_object
  if (!is.null(lognorm_object)) norm_objects[["LogNormalize"]] <- lognorm_object
  if (!is.null(sckwarn_object)) norm_objects[["scKWARN"]] <- sckwarn_object

  cat("Normalization methods to benchmark:", paste(names(norm_objects), collapse = ", "), "\n")
  cat("Total combinations:", length(norm_objects), "normalizations x", length(available_bench_methods), "integrations =", length(norm_objects) * length(available_bench_methods), "\n\n")

  # Subsample objects for faster benchmarking
  bench_max_cells <- if (!is.null(params$norm_benchmark_max_cells)) params$norm_benchmark_max_cells else NULL

  if (!is.null(bench_max_cells) && bench_max_cells > 0) {
    cat("--- Subsampling for benchmarking (max", bench_max_cells, "cells) ---\n")
    norm_objects_bench <- list()
    for (norm_name in names(norm_objects)) {
      cat("\n", norm_name, ":\n")
      sub_obj <- subsample_for_benchmark(
        norm_objects[[norm_name]],
        max_cells = bench_max_cells,
        batch_var = params$batch_variable,
        seed = params$random_seed
      )
      DefaultAssay(sub_obj) <- DefaultAssay(norm_objects[[norm_name]])
      sub_obj <- FindVariableFeatures(sub_obj, nfeatures = params$nfeatures_integration, verbose = FALSE)
      # FIX 2026-02-05: Include vars.to.regress to match full pipeline preprocessing
      sub_obj <- ScaleData(sub_obj, vars.to.regress = params$sct_vars_to_regress, verbose = FALSE)
      sub_obj <- RunPCA(sub_obj, npcs = 50, verbose = FALSE)

      # FIX 2026-02-05: Strip pre-existing integration reductions so all methods
      # start from the same PCA-only state on the subsampled data.
      # Without this, Harmony (and any other method that checks for existing
      # reductions) would use stale full-data embeddings while other methods
      # recompute fresh on the subset.
      stale_reductions <- c("harmony", "integrated.mnn", "integrated.rpca", "integrated.cca", "umap", "umap_unintegrated")
      for (red_name in stale_reductions) {
        if (red_name %in% names(sub_obj@reductions)) {
          sub_obj[[red_name]] <- NULL
        }
      }

      norm_objects_bench[[norm_name]] <- sub_obj
    }
    cat("\n")
  } else {
    cat("--- No subsampling: benchmarking on full objects ---\n\n")
    norm_objects_bench <- norm_objects

    # FIX 2026-02-05: Even without subsampling, strip pre-existing integration
    # reductions to ensure fair recomputation across all methods.
    cat("--- Stripping pre-existing integration reductions for fair benchmarking ---\n")
    for (norm_name in names(norm_objects_bench)) {
      obj_tmp <- norm_objects_bench[[norm_name]]
      stale_reductions <- c("harmony", "integrated.mnn", "integrated.rpca", "integrated.cca", "umap", "umap_unintegrated")
      stripped <- c()
      for (red_name in stale_reductions) {
        if (red_name %in% names(obj_tmp@reductions)) {
          obj_tmp[[red_name]] <- NULL
          stripped <- c(stripped, red_name)
        }
      }
      if (length(stripped) > 0) {
        cat("  ", norm_name, ": removed", paste(stripped, collapse = ", "), "\n")
      }
      norm_objects_bench[[norm_name]] <- obj_tmp
    }
    cat("\n")
  }

  # ==========================================================================
  # Re-initialize Python bridge before benchmarking loop
  # ==========================================================================
  # After heavy R processing (4 normalizations, merging, Harmony, subsampling),
  # the reticulate-Python bridge may have become stale under memory pressure.
  # Re-warm the connection immediately before use.
  # ==========================================================================
  python_methods_in_bench <- intersect(available_bench_methods, c("scvi", "sccobra", "concord"))
  if (length(python_methods_in_bench) > 0) {
    cat("--- Re-initializing Python bridge for benchmarking ---\n")
    tryCatch({
      suppressPackageStartupMessages(library(reticulate))
      Sys.setenv(RETICULATE_PYTHON = params$unified_python)
      use_python(params$unified_python, required = TRUE)

      # Force a round-trip through the Python bridge to verify it's alive
      reticulate::py_run_string("import sys")
      cat("  [OK] Python bridge verified\n")

      # Re-import all needed modules
      reticulate::py_run_string("import numpy, pandas, anndata, scanpy")
      cat("  [OK] Core Python modules re-imported\n")

      if ("scvi" %in% python_methods_in_bench) {
        reticulate::py_run_string("import scvi")
        cat("  [OK] scvi re-imported\n")
      }
      if ("sccobra" %in% python_methods_in_bench) {
        sccobra_path <- "/scicore/home/doetsch/kaiser0001/GITHUB_repositories/scCobra"
        reticulate::py_run_string(sprintf("import sys; sys.path.insert(0, '%s')", sccobra_path))
        reticulate::py_run_string("from scCobra import scCobra")
        cat("  [OK] scCobra re-imported\n")
      }
      if ("concord" %in% python_methods_in_bench) {
        reticulate::py_run_string("import concord")
        cat("  [OK] concord re-imported\n")
      }
      cat("\n")
    }, error = function(e) {
      cat("  [WARNING] Python bridge failed:", e$message, "\n")
      cat("  Removing Python methods from benchmarking.\n\n")
      available_bench_methods <<- setdiff(available_bench_methods, c("scvi", "sccobra", "concord"))
    })
  }

  # Full benchmark matrix: normalization x integration
  # Python methods (scVI, scCobra, Concord) use raw counts internally and
  # normalize in Python, so they produce identical embeddings regardless of
  # which R-side normalization was applied. They are NOT valid voters for
  # selecting the best normalization method. They remain available for
  # Module 04 integration method comparison (where they compare integration
  # approaches on a single fixed normalization).
  python_methods <- c("scvi", "sccobra", "concord")
  norm_vote_methods <- setdiff(available_bench_methods, python_methods)

  if (length(norm_vote_methods) == 0) {
    cat("WARNING: No R-based integration methods available for normalization voting.\n")
    cat("  Python methods cannot distinguish between normalizations (they use raw counts).\n")
    cat("  Defaulting to SCTransform.\n")
    selected_normalization_method <- "SCTransform"
  } else {
    cat("Integration methods for normalization vote:", paste(norm_vote_methods, collapse = ", "), "\n")
    cat("Excluded from vote (use raw counts):", paste(intersect(available_bench_methods, python_methods), collapse = ", "), "\n\n")

    # Full benchmark matrix: normalization x integration (R-based methods only)
    norm_benchmark_full_df <- data.frame()
    for (norm_name in names(norm_objects_bench)) {
      cat("\n", paste(rep("-", 50), collapse = ""), "\n")
      cat("Normalization:", norm_name, "\n")
      cat(paste(rep("-", 50), collapse = ""), "\n")
      obj <- norm_objects_bench[[norm_name]]
      # Ensure batch column is set from batch_variable
      obj$batch <- obj@meta.data[[params$batch_variable]]
      for (int_method in norm_vote_methods) {
        emb <- run_integration_for_benchmark(obj, int_method, "batch", params$dims_use, params)
      if (!is.null(emb)) {
        batch_labels <- obj$batch

        # Batch correction metrics
        bv <- compute_batch_variance(emb, batch_labels)
        ba <- compute_batch_asw(emb, batch_labels)
        li <- compute_lisi_score(emb, batch_labels)

        # Biological conservation metrics (if annotations available)
        ct_asw <- NA_real_
        ct_lisi <- NA_real_
        nmi_val <- NA_real_
        ari_val <- NA_real_

        if (norm_bench_has_celltype && norm_bench_celltype_column %in% colnames(obj@meta.data)) {
          ct_labels <- obj@meta.data[[norm_bench_celltype_column]]
          ct_asw <- compute_celltype_asw(emb, ct_labels)
          ct_lisi <- compute_celltype_lisi(emb, ct_labels)
          nmi_val <- compute_nmi_score(emb, ct_labels)
          ari_val <- compute_ari_score(emb, ct_labels)
        }

        scores <- data.frame(
          normalization = norm_name,
          integration_method = int_method,
          batch_variance = bv,
          batch_asw = ba,
          lisi = li,
          celltype_asw = ct_asw,
          celltype_lisi = ct_lisi,
          nmi = nmi_val,
          ari = ari_val,
          stringsAsFactors = FALSE
        )

        norm_benchmark_full_df <- rbind(norm_benchmark_full_df, scores)

        cat("    batch_variance:", round(bv, 4),
            " batch_asw:", round(ba, 4),
            " lisi:", round(li, 4))
        if (norm_bench_has_celltype) {
          cat(" | ct_asw:", round(ct_asw, 4),
              " ct_lisi:", round(ct_lisi, 4),
              " nmi:", round(nmi_val, 4),
              " ari:", round(ari_val, 4))
        }
        cat("\n")
      }                         # close if (!is.null(emb))
      }                         # close for (int_method)
      gc()
    }                           # close for (norm_name)
  }                             # close else

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("FULL NORMALIZATION BENCHMARK MATRIX\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  print(norm_benchmark_full_df)
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  # Determine selection mode with safe fallback
  selection_mode <- if (!is.null(params$integration_selection_mode)) params$integration_selection_mode else "balanced"

  # Determine best normalization method
  norm_benchmark_vote_result <- NULL

  if (nrow(norm_benchmark_full_df) > 0 &&
      length(unique(norm_benchmark_full_df$normalization)) >= 2 &&
      length(unique(norm_benchmark_full_df$integration_method)) >= 2) {

    # Multiple integration methods available -> use majority vote
    norm_benchmark_vote_result <- select_normalization_by_majority_vote(
      norm_benchmark_full_df,
      selection_mode
    )
    selected_normalization_method <- norm_benchmark_vote_result$winner

  } else if (nrow(norm_benchmark_full_df) > 0) {

    # Only 1 integration method or 1 normalization -> fall back to composite score
    cat("Only 1 integration method available. Selecting by composite score (no majority vote).\n")

    norm_benchmark_full_df$bv_norm <- normalize_metric(norm_benchmark_full_df$batch_variance, lower_is_better = TRUE)
    norm_benchmark_full_df$asw_norm <- normalize_metric(norm_benchmark_full_df$batch_asw, lower_is_better = TRUE)
    norm_benchmark_full_df$lisi_norm <- normalize_metric(norm_benchmark_full_df$lisi, lower_is_better = FALSE)

    norm_benchmark_full_df$composite <- rowMeans(
      norm_benchmark_full_df[, c("bv_norm", "asw_norm", "lisi_norm")], na.rm = TRUE
    )

    best_idx <- which.max(norm_benchmark_full_df$composite)
    selected_normalization_method <- norm_benchmark_full_df$normalization[best_idx]

    cat("\n>>> BEST NORMALIZATION METHOD:", selected_normalization_method, "<<<\n\n")

  } else {
    cat("WARNING: No benchmark scores computed. Defaulting to SCTransform.\n")
    selected_normalization_method <- "SCTransform"
  }

  # Save detailed benchmark results
  write.csv(norm_benchmark_full_df,
            file.path(subdirs$norm_benchmark, "normalization_multi_integration_benchmark.csv"),
            row.names = FALSE)

  # Backward-compatible summary (average across integration methods per normalization)
  if (nrow(norm_benchmark_full_df) > 0) {
  if (norm_bench_has_celltype) {
      norm_benchmark_results <- aggregate(
        cbind(batch_variance, batch_asw, lisi, celltype_asw, celltype_lisi, nmi, ari) ~ normalization,
        data = norm_benchmark_full_df,
        FUN = mean, na.rm = TRUE,
        na.action = na.pass
      )
    } else {
      norm_benchmark_results <- aggregate(
        cbind(batch_variance, batch_asw, lisi) ~ normalization,
        data = norm_benchmark_full_df,
        FUN = mean, na.rm = TRUE,
        na.action = na.pass
      )
    }
    colnames(norm_benchmark_results)[1] <- "method"
    norm_benchmark_results$integration_method <- "multi_method_average"
    norm_benchmark_results$batch_variable <- params$batch_variable
  } else {
    norm_benchmark_results <- data.frame(
      method = character(), batch_variance = numeric(),
      batch_asw = numeric(), lisi = numeric(),
      integration_method = character(), batch_variable = character(),
      stringsAsFactors = FALSE
    )
  }

  write.csv(norm_benchmark_results,
            file.path(subdirs$norm_benchmark, "normalization_benchmark_scores.csv"),
            row.names = FALSE)

  # Save vote results if available
  if (!is.null(norm_benchmark_vote_result)) {
    write.csv(norm_benchmark_vote_result$per_integration_winners,
              file.path(subdirs$norm_benchmark, "normalization_majority_vote_details.csv"),
              row.names = FALSE)
    write.csv(norm_benchmark_vote_result$vote_table,
              file.path(subdirs$norm_benchmark, "normalization_vote_tally.csv"),
              row.names = FALSE)
  }

  cat("\n>>> NORMALIZATION BENCHMARKING SUMMARY <<<\n")
  print(norm_benchmark_results)
  cat("\n>>> BEST NORMALIZATION METHOD:", selected_normalization_method, "<<<\n")

  # Create comparison plots
  suppressPackageStartupMessages(library(tidyr))

  # Plot 1: Heatmap of composite scores (normalization x integration)
  if (nrow(norm_benchmark_full_df) > 0 &&
      length(unique(norm_benchmark_full_df$integration_method)) >= 1 &&
      length(unique(norm_benchmark_full_df$normalization)) >= 2) {

    plot_df <- norm_benchmark_full_df
    plot_df$bv_norm <- normalize_metric(plot_df$batch_variance, lower_is_better = TRUE)
    plot_df$asw_norm <- normalize_metric(plot_df$batch_asw, lower_is_better = TRUE)
    plot_df$lisi_norm <- normalize_metric(plot_df$lisi, lower_is_better = FALSE)
    plot_df$composite <- rowMeans(plot_df[, c("bv_norm", "asw_norm", "lisi_norm")], na.rm = TRUE)

    p_heatmap <- ggplot(plot_df, aes(x = integration_method, y = normalization, fill = composite)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", composite)), size = 3.5) +
      scale_fill_gradient2(low = "white", high = "#2166AC", midpoint = 0.5,
                           limits = c(0, 1), name = "Composite\nScore") +
      theme_minimal() +
      labs(title = "Normalization x Integration Benchmark (Composite Scores)",
           subtitle = paste0("Selection mode: ", selection_mode,
                            " | Winner: ", selected_normalization_method),
           x = "Integration Method", y = "Normalization Method") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_plot_multi(p_heatmap, "01_normalization_integration_heatmap",
                    output_dir = subdirs$norm_benchmark, width = 10, height = 6)
  }

  # Plot 2: Majority vote tally bar chart
  if (!is.null(norm_benchmark_vote_result)) {
    vote_df <- norm_benchmark_vote_result$vote_table

    p_votes <- ggplot(vote_df, aes(x = reorder(normalization, votes), y = votes, fill = normalization)) +
      geom_col() +
      geom_text(aes(label = votes), hjust = -0.2, size = 4) +
      coord_flip() +
      theme_minimal() +
      labs(title = "Normalization Method - Majority Vote Across Integration Methods",
           subtitle = paste0("Each of ", nrow(norm_benchmark_vote_result$per_integration_winners),
                            " integration methods votes for its best normalization | Winner: ",
                            selected_normalization_method),
           x = "Normalization Method", y = "Votes") +
      theme(legend.position = "none") +
      scale_fill_brewer(palette = "Set2") +
      ylim(0, max(vote_df$votes) * 1.3)

    save_plot_multi(p_votes, "02_normalization_majority_vote",
                    output_dir = subdirs$norm_benchmark, width = 8, height = 5)
  }

  # Plot 3: Backward-compatible bar chart (averaged metrics across integrations)
  norm_benchmark_long <- norm_benchmark_results %>%
    select(method, batch_variance, batch_asw, lisi) %>%
    pivot_longer(-method, names_to = "metric", values_to = "score") %>%
    filter(!is.na(score))

  if (nrow(norm_benchmark_long) > 0) {
    p_norm_bench <- ggplot(norm_benchmark_long, aes(x = method, y = score, fill = method)) +
      geom_col() +
      facet_wrap(~ metric, scales = "free_y", ncol = 3) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Normalization Methods Comparison (averaged across integration methods)",
           subtitle = paste0("Batch variable: ", params$batch_variable,
                            " | Winner: ", selected_normalization_method),
           x = "Method", y = "Score") +
      scale_fill_brewer(palette = "Set2")

    save_plot_multi(p_norm_bench, "03_normalization_benchmark_comparison",
                    output_dir = subdirs$norm_benchmark, width = 12, height = 6)
  }

} else {
  selected_normalization_method <- "SCTransform"
  norm_benchmark_full_df <- NULL
  norm_benchmark_vote_result <- NULL
  cat("\nSkipping benchmarking: fewer than 2 normalization methods available.\n")
}

# ==============================================================================
# Save for next module
# ==============================================================================
norm_file <- file.path(output_dirs$objects, "03_normalization_data.RData")
save(sct_object, scran_object, lognorm_object, sckwarn_object, selected_normalization_method,
     sample_objects, file = norm_file)
cat("\nNormalization data saved to:", norm_file, "\n")

# Save additional multi-integration benchmark details (separate file for backward compatibility)
if (exists("norm_benchmark_full_df") && !is.null(norm_benchmark_full_df)) {
  norm_bench_detail_file <- file.path(output_dirs$objects, "03_normalization_benchmark_details.RData")
  save(norm_benchmark_full_df, norm_benchmark_vote_result, file = norm_bench_detail_file)
  cat("Multi-integration benchmark details saved to:", norm_bench_detail_file, "\n")
}

cat("\n>>> MODULE 03 COMPLETE <<<\n")
cat("Selected normalization method:", selected_normalization_method, "\n")
cat("Individual normalized objects saved to:", individual_norm_dir, "\n")
cat("Merged normalized objects saved to:", merged_norm_dir, "\n")
if (!is.null(sckwarn_object)) {
  cat("scKWARN object includes:\n")
  cat("  - 'data' layer: log1p(NormalizedData) - natural log scale\n")
  cat("  - 'scKWARN_linear' assay: linear normalized data (no log)\n")
  cat("  - 'sckwarn_scaling_factor' in metadata\n")
}
