#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 04: MULTI-METHOD INTEGRATION BENCHMARKING (MULTI-SAMPLE PIPELINE)
# ==============================================================================
#
# This module performs batch correction using multiple methods and benchmarks
# their performance.
#
# Methods tested:
# - R-based: Harmony, CCA, RPCA, FastMNN
# - Python-based: scVI, Scanorama, BBKNN, scCobra
#
# INPUT: Normalized Seurat object from Module 03 (UNINTEGRATED version)
# OUTPUT: Integrated object with best method selected
#
# CRITICAL NOTES:
# 1. This module uses the UNINTEGRATED normalized object from Module 03
# 2. For CCA, RPCA, FastMNN: Must use RNA assay with LogNormalize
#    (SCTAssay has compatibility issues with these methods)
# 3. The object is prepared with proper layer structure before integration
#
# SELECTION MODES:
# - "batch_removal": Minimizes batch_variance (most aggressive)
#   Use for technical replicates only
# - "balanced": Maximizes composite score across all metrics (RECOMMENDED)
#   Use when batches have biological meaning (different conditions/sexes)
# - "conservative": Prioritizes LISI score (moderate mixing)
#   Use for preserving subtle biological differences
#
# METRICS:
# Batch Correction (higher normalized score = better correction):
# - batch_variance: Variance explained by batch (lower raw = better)
# - batch_asw: Batch silhouette width (closer to 0 = better)
# - lisi: Local Inverse Simpson Index (higher = better mixing)
#
# Biological Conservation (higher normalized score = better preservation):
# - celltype_asw: Cell type silhouette width (higher = better separation)
# - celltype_lisi: Cell type LISI (lower = better local purity)
# - nmi: Normalized Mutual Information (higher = better clustering match)
# - ari: Adjusted Rand Index (higher = better clustering match)
#
# UPDATES:
# - 2026-01-03: Improved Seurat v5 layer handling
# - 2026-01-03: Added better error handling for Python integration methods
# - 2026-01-03: Ensured JoinLayers before Python export
# - 2026-01-04: Fixed batch variable assignment - now properly uses params$batch_variable
# - 2026-01-06: Added integration_selection_mode for balanced vs aggressive selection
# - 2026-01-06: Uses normalize_metric from functions.R
# - 2026-01-10: CRITICAL FIX - Fixed batch_asw normalization to use absolute values
# - 2026-01-10: Added biological conservation metrics (celltype_asw, celltype_lisi, NMI, ARI)
# - 2026-01-10: Enhanced output reporting with detailed benchmark report
# - 2026-01-10: Added scIB-style composite scoring (40% bio + 60% batch)
# - 2026-02-09: Added per-label batch mixing validation after method selection
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 04: MULTI-METHOD INTEGRATION BENCHMARKING\n")
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
  library(future)
  library(cluster)
  library(reshape2)
  library(reticulate)
  library(Matrix)
  library(tidyr)
})

out_base <- params$out_root
load(file.path(out_base, "objects", "pipeline_environment.RData"))

# Defensive fallbacks if pipeline_environment.RData schema changes
if (!exists("has_harmony")) has_harmony <- requireNamespace("harmony", quietly = TRUE)
if (!exists("has_batchelor")) has_batchelor <- requireNamespace("batchelor", quietly = TRUE)
if (!exists("has_SeuratWrappers")) has_SeuratWrappers <- requireNamespace("SeuratWrappers", quietly = TRUE)
load(file.path(out_base, "objects", "03_normalization_data.RData"))

if (has_harmony) library(harmony)

# Check for batchelor/SeuratWrappers for FastMNN
cat("FastMNN dependencies:\n")
cat("  batchelor:", has_batchelor, "\n")
cat("  SeuratWrappers:", has_SeuratWrappers, "\n\n")

# Set Python environment
Sys.setenv(RETICULATE_PYTHON = params$unified_python)
use_python(params$unified_python, required = TRUE)

# ==============================================================================
# Skip if disabled
# ==============================================================================
if (!isTRUE(params$run_integration_benchmarking)) {
  cat("Integration benchmarking disabled. Skipping.\n")

  multi_integrated <- NULL
  best_method <- "none"
  best_reduction <- "pca"

  integration_file <- file.path(output_dirs$objects, "04_integration_data.RData")
  save(multi_integrated, best_method, best_reduction, file = integration_file)

  cat("\n>>> MODULE 04 SKIPPED <<<\n")
  quit(save = "no", status = 0)
}

# ==============================================================================
# Handle batch integration control
# ==============================================================================
cat("Batch variable:", params$batch_variable, "\n")
cat("Run batch integration:", params$run_batch_integration, "\n\n")

if (!isTRUE(params$run_batch_integration)) {
  cat("Batch integration disabled (run_batch_integration = FALSE).\n")
  cat("Will only compute unintegrated baseline.\n\n")

  params$integration_methods_r <- c()
  params$integration_methods_python <- c()
  params$run_python_integrations <- FALSE
  params$run_sccobra <- FALSE
}

# Future settings - increase for large objects
options(future.globals.maxSize = 30 * 1024^3)  # 30 GB
future::plan(sequential)

# ==============================================================================
# Import Python modules
# ==============================================================================
cat("--- Importing Python modules ---\n")

np <- import("numpy")
pd <- import("pandas")
anndata <- import("anndata")
sc <- import("scanpy")

cat("  numpy version:", np$`__version__`, "\n")
cat("  pandas version:", pd$`__version__`, "\n")
cat("  anndata version:", anndata$`__version__`, "\n")
cat("  scanpy version:", sc$`__version__`, "\n")

# Check for scvi
has_scvi_python <- tryCatch({
  scvi <- import("scvi")
  cat("  [OK] scvi-tools version:", scvi$`__version__`, "\n")
  TRUE
}, error = function(e) {
  cat("  [WARNING] scvi-tools not available:", e$message, "\n")
  FALSE
})

# Check for scanorama
has_scanorama_python <- tryCatch({
  scanorama <- import("scanorama")
  cat("  [OK] scanorama\n")
  TRUE
}, error = function(e) {
  cat("  [WARNING] scanorama not available\n")
  FALSE
})

# Check for bbknn
has_bbknn_python <- tryCatch({
  bbknn <- import("bbknn")
  cat("  [OK] bbknn\n")
  TRUE
}, error = function(e) {
  cat("  [WARNING] bbknn not available\n")
  FALSE
})

# Check for scCobra
has_sccobra_python <- tryCatch({
  sccobra_path <- "/scicore/home/doetsch/kaiser0001/GITHUB_repositories/scCobra"
  if (dir.exists(sccobra_path)) {
    py_run_string(sprintf("import sys; sys.path.insert(0, '%s')", sccobra_path))
    py_run_string("import scCobra")
    cat("  [OK] scCobra (from", sccobra_path, ")\n")
    TRUE
  } else {
    cat("  [WARNING] scCobra path not found:", sccobra_path, "\n")
    FALSE
  }
}, error = function(e) {
  cat("  [WARNING] scCobra not available:", e$message, "\n")
  FALSE
})

cat("\n")

# ==============================================================================
# Determine normalization method for integration
# ==============================================================================
normalization_method <- get_integration_normalization_method(params, output_dirs$objects)

# ==============================================================================
# Load the appropriate normalized object
# ==============================================================================
merged_raw <- load_normalized_object_for_integration(
  method = normalization_method,
  objects_dir = output_dirs$objects,
  use_unintegrated = TRUE
)

if (is.null(merged_raw)) {
  cat("\n!!! ERROR: Could not load normalized object for integration !!!\n")
  cat("Attempting fallback to 03_normalization_data.RData...\n")

  # Fallback: try to use objects from RData
  if (!is.null(lognorm_object)) {
    merged_raw <- lognorm_object
    normalization_method <- "LogNormalize"
    cat("Using LogNormalize object from RData\n")
  } else if (!is.null(sct_object)) {
    merged_raw <- sct_object
    normalization_method <- "SCTransform"
    cat("Using SCTransform object from RData\n")
  } else if (!is.null(scran_object)) {
    merged_raw <- scran_object
    normalization_method <- "scran"
    cat("Using scran object from RData\n")
  } else {
    stop("No normalized object available for integration")
  }
}

# ==============================================================================
# Print counts layer information
# ==============================================================================
print_counts_layer_info(merged_raw, assay = "RNA")

# ==============================================================================
# Prepare object for integration
# ==============================================================================
# This function:
# 1. Sets default assay to RNA (required for CCA/RPCA/FastMNN)
# 2. Joins layers, then splits by batch
# 3. Runs NormalizeData, FindVariableFeatures, ScaleData, RunPCA

merged_raw <- prepare_object_for_integration(
  obj = merged_raw,
  batch_var = params$batch_variable,
  nfeatures = params$nfeatures_integration,
  dims_use = params$dims_use
)

# ==============================================================================
# CRITICAL FIX: Ensure batch column is properly set from batch_variable
# ==============================================================================
# This ensures that obj$batch contains the actual batch labels (e.g., sample names)
# rather than a column that might not exist or have wrong values.
# This is essential for all integration methods and benchmarking to work correctly.
# ==============================================================================
merged_raw$batch <- merged_raw@meta.data[[params$batch_variable]]

cat("\n--- Batch variable assignment ---\n")
cat("  Using column '", params$batch_variable, "' as batch identity\n", sep = "")
cat("  Batch levels:", paste(unique(merged_raw$batch), collapse = ", "), "\n")
cat("  Cells per batch:\n")
print(table(merged_raw$batch))
cat("\n")

print_object_structure(merged_raw, "Merged Preprocessed for Integration")

# ==============================================================================
# Detect cell type column for biological conservation metrics
# ==============================================================================
cat("\n--- Detecting cell type annotations ---\n")

# Check for cell type column
celltype_column <- NULL
has_celltype_annotations <- FALSE

# First check params
if (!is.null(params$celltype_column) && params$celltype_column %in% colnames(merged_raw@meta.data)) {
  celltype_column <- params$celltype_column
  has_celltype_annotations <- TRUE
  cat("  Using celltype column from params:", celltype_column, "\n")
}

# If not in params, try common column names
if (!has_celltype_annotations) {
  common_celltype_cols <- c("cell_type", "celltype", "CellType", "cell_type_annotation",
                            "cluster_annotation", "annotation", "cell_annotation",
                            "predicted_celltype", "predicted.celltype", "sctype_classification")

  for (col in common_celltype_cols) {
    if (col %in% colnames(merged_raw@meta.data)) {
      # Check that it has meaningful values (not all NA, not all same value)
      col_values <- merged_raw@meta.data[[col]]
      n_unique <- length(unique(na.omit(col_values)))
      if (n_unique > 1 && n_unique < ncol(merged_raw) / 2) {
        celltype_column <- col
        has_celltype_annotations <- TRUE
        cat("  Auto-detected celltype column:", celltype_column, "\n")
        break
      }
    }
  }
}

if (has_celltype_annotations) {
  celltype_labels <- merged_raw@meta.data[[celltype_column]]
  n_celltypes <- length(unique(na.omit(celltype_labels)))
  cat("  Number of cell types:", n_celltypes, "\n")
  cat("  Cell type distribution:\n")
  print(table(celltype_labels, useNA = "ifany"))
} else {
  cat("  No cell type annotations found.\n")
  cat("  Biological conservation metrics will be skipped.\n")
  cat("  To enable: set params$celltype_column or add 'cell_type' column to metadata.\n")
}
cat("\n")

# ==============================================================================
# Define biological conservation metric functions
# ==============================================================================

# ==============================================================================
# Run integration methods
# ==============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("RUNNING INTEGRATION METHODS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

integration_results <- list()
method_names <- c("unintegrated")
integration_results[["unintegrated"]] <- merged_raw

get_integration_reduction <- function(method_name) {
  switch(method_name,
         "unintegrated" = "pca",
         "harmony" = "harmony",
         "cca" = "integrated.cca",
         "rpca" = "integrated.rpca",
         "mnn" = "integrated.mnn",
         "scvi" = "integrated.scvi",
         "scanorama" = "integrated.scanorama",
         "bbknn" = "integrated.bbknn",
         "sccobra" = "integrated.sccobra",
         "concord" = "integrated.concord",
         "pca")
}

# ==============================================================================
# R-based integration methods
# ==============================================================================

# Harmony
if ("harmony" %in% params$integration_methods_r && has_harmony) {
  cat("\n--- Running Harmony ---\n")
  cat("  Input: PCA embeddings <- ScaleData <- RNA$data (Module 03 normalization)\n")
  tryCatch({
    harmony_obj <- IntegrateLayers(
      merged_raw,
      method = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction = "harmony",
      verbose = FALSE
    )
    integration_results[["harmony"]] <- harmony_obj
    method_names <- c(method_names, "harmony")
    cat("  [SUCCESS] Harmony integration complete\n")
    cat("  New reduction: harmony with", ncol(Embeddings(harmony_obj, "harmony")), "dimensions\n")
  }, error = function(e) {
    print_integration_error("Harmony", e, merged_raw)
  })
}

# CCA
if ("cca" %in% params$integration_methods_r) {
  cat("\n--- Running CCA ---\n")
  cat("  Input: split RNA$data layers (Module 03 normalization)\n")
  cat("  Note: CCA requires RNA assay with proper layer structure\n")
  tryCatch({
    cca_obj <- IntegrateLayers(
      merged_raw,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.cca",
      verbose = FALSE
    )
    integration_results[["cca"]] <- cca_obj
    method_names <- c(method_names, "cca")
    cat("  [SUCCESS] CCA integration complete\n")
    cat("  New reduction: integrated.cca with", ncol(Embeddings(cca_obj, "integrated.cca")), "dimensions\n")
  }, error = function(e) {
    print_integration_error("CCA", e, merged_raw)
  })
}

# RPCA
if ("rpca" %in% params$integration_methods_r) {
  cat("\n--- Running RPCA ---\n")
  cat("  Input: split RNA$data layers (Module 03 normalization)\n")
  cat("  Note: RPCA requires RNA assay with proper layer structure\n")
  tryCatch({
    rpca_obj <- IntegrateLayers(
      merged_raw,
      method = RPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.rpca",
      verbose = FALSE
    )
    integration_results[["rpca"]] <- rpca_obj
    method_names <- c(method_names, "rpca")
    cat("  [SUCCESS] RPCA integration complete\n")
    cat("  New reduction: integrated.rpca with", ncol(Embeddings(rpca_obj, "integrated.rpca")), "dimensions\n")
  }, error = function(e) {
    print_integration_error("RPCA", e, merged_raw)
  })
}

# FastMNN
if ("mnn" %in% params$integration_methods_r && has_batchelor && has_SeuratWrappers) {
  cat("\n--- Running FastMNN ---\n")
  cat("  Input: split RNA$data layers (Module 03 normalization)\n")
  cat("  Note: FastMNN requires layers split by batch\n")
  tryCatch({
    library(SeuratWrappers)
    mnn_obj <- IntegrateLayers(
      object = merged_raw,
      method = FastMNNIntegration,
      new.reduction = "integrated.mnn",
      verbose = FALSE
    )
    integration_results[["mnn"]] <- mnn_obj
    method_names <- c(method_names, "mnn")
    cat("  [SUCCESS] FastMNN integration complete\n")
    cat("  New reduction: integrated.mnn with", ncol(Embeddings(mnn_obj, "integrated.mnn")), "dimensions\n")
  }, error = function(e) {
    print_integration_error("FastMNN", e, merged_raw)
  })
}

# ==============================================================================
# Python-based integration methods (scVI, Scanorama, BBKNN, scCobra)
# ==============================================================================
if (isTRUE(params$run_python_integrations)) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("PYTHON-BASED INTEGRATION METHODS\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  # scVI
  if ("scvi" %in% params$integration_methods_python && has_scvi_python) {
    cat("\n--- Running scVI ---\n")
    cat("  Input: RNA$counts (raw counts) -- scVI normalizes internally via negative binomial VAE\n")
    tryCatch({
      # Create a copy with joined layers for Python export
      temp_obj <- merged_raw

      # Ensure layers are joined before Python export
      tryCatch({
        temp_obj[["RNA"]] <- JoinLayers(temp_obj[["RNA"]])
      }, error = function(e) {
        cat("  Note: JoinLayers:", e$message, "\n")
      })

      var_features <- VariableFeatures(temp_obj)
      counts_data <- GetAssayData(temp_obj, assay = "RNA", layer = "counts")[var_features, ]
      counts_data <- as.matrix(counts_data)
      counts_data_t <- t(counts_data)

      cat("  Data dimensions:", nrow(counts_data_t), "cells x", ncol(counts_data_t), "genes\n")

      obs_df <- data.frame(
        batch = as.character(temp_obj$batch),
        row.names = colnames(temp_obj),
        stringsAsFactors = FALSE
      )
      var_df <- data.frame(
        gene = var_features,
        row.names = var_features,
        stringsAsFactors = FALSE
      )

      cat("  Creating AnnData object...\n")

      py_run_string("import numpy as np")
      py_run_string("import pandas as pd")
      py_run_string("import anndata")
      py_run_string("import scvi")

      py$counts_matrix <- r_to_py(counts_data_t)
      py$obs_data <- r_to_py(obs_df)
      py$var_data <- r_to_py(var_df)
      py$cell_names <- r_to_py(colnames(temp_obj))
      py$gene_names <- r_to_py(var_features)

      py_run_string("
adata_scvi = anndata.AnnData(
    X=np.array(counts_matrix, dtype='float32'),
    obs=pd.DataFrame(obs_data),
    var=pd.DataFrame(var_data)
)
adata_scvi.obs.index = list(cell_names)
adata_scvi.var.index = list(gene_names)
adata_scvi.obs['batch'] = adata_scvi.obs['batch'].astype('category')
")

      cat("  Setting up scVI model...\n")
      py_run_string("scvi.model.SCVI.setup_anndata(adata_scvi, batch_key='batch')")

      n_latent <- as.integer(params$dims_use)
      py_run_string(sprintf("model = scvi.model.SCVI(adata_scvi, n_latent=%d)", n_latent))

      cat("  Training scVI model (400 epochs)...\n")
      py_run_string("model.train(max_epochs=400, early_stopping=True)")

      cat("  Extracting latent representation...\n")
      py_run_string("latent = model.get_latent_representation()")

      latent <- py_to_r(py$latent)
      latent <- as.matrix(latent)
      rownames(latent) <- colnames(temp_obj)
      colnames(latent) <- paste0("scVI_", 1:ncol(latent))

      scvi_obj <- merged_raw
      scvi_obj[["integrated.scvi"]] <- CreateDimReducObject(
        embeddings = latent,
        key = "scvi_",
        assay = DefaultAssay(merged_raw)
      )

      integration_results[["scvi"]] <- scvi_obj
      method_names <- c(method_names, "scvi")
      cat("  [SUCCESS] scVI integration complete\n")
      cat("  New reduction: integrated.scvi with", ncol(latent), "dimensions\n")

      py_run_string("del adata_scvi, model, latent, counts_matrix, obs_data, var_data, cell_names, gene_names")

    }, error = function(e) {
      print_integration_error("scVI", e, merged_raw)
    })
  }

  # Scanorama
  if ("scanorama" %in% params$integration_methods_python && has_scanorama_python) {
    cat("\n--- Running Scanorama ---\n")
    cat("  Input: RNA$data (Module 03 normalization) -- Scanorama uses pre-normalized expression\n")
    tryCatch({
      scanorama <- import("scanorama")

      temp_obj <- merged_raw

      # Ensure layers are joined before Python export
      tryCatch({
        temp_obj[["RNA"]] <- JoinLayers(temp_obj[["RNA"]])
      }, error = function(e) {
        cat("  Note: JoinLayers:", e$message, "\n")
      })

      batch_labels <- temp_obj$batch
      batches <- unique(batch_labels)
      var_features <- VariableFeatures(temp_obj)

      cat("  Batches:", paste(batches, collapse = ", "), "\n")

      datasets <- list()
      genes_list <- list()
      cell_order <- c()

      for (batch in batches) {
        cells_batch <- colnames(temp_obj)[batch_labels == batch]
        cell_order <- c(cell_order, cells_batch)
        mat <- as.matrix(GetAssayData(temp_obj, assay = "RNA", layer = "data")[var_features, cells_batch])
        datasets <- c(datasets, list(t(mat)))
        genes_list <- c(genes_list, list(as.list(var_features)))
      }

      cat("  Running scanorama.integrate()...\n")
      result <- scanorama$integrate(datasets, genes_list)
      integrated_data <- result[[1]]

      corrected_list <- lapply(integrated_data, function(x) {
        tryCatch(as.matrix(x), error = function(e1) {
          tryCatch(as.matrix(x$toarray()), error = function(e2) as.matrix(np$array(x)))
        })
      })

      corrected_combined <- do.call(rbind, corrected_list)
      corrected_reordered <- corrected_combined[match(colnames(merged_raw), cell_order), , drop = FALSE]
      rownames(corrected_reordered) <- colnames(merged_raw)
      colnames(corrected_reordered) <- paste0("scanorama_", 1:ncol(corrected_reordered))

      scano_obj <- merged_raw
      scano_obj[["integrated.scanorama"]] <- CreateDimReducObject(
        embeddings = corrected_reordered, key = "scanorama_", assay = DefaultAssay(merged_raw))

      integration_results[["scanorama"]] <- scano_obj
      method_names <- c(method_names, "scanorama")
      cat("  [SUCCESS] Scanorama integration complete\n")
      cat("  New reduction: integrated.scanorama with", ncol(corrected_reordered), "dimensions\n")
    }, error = function(e) {
      print_integration_error("Scanorama", e, merged_raw)
    })
  }

  # BBKNN
  if ("bbknn" %in% params$integration_methods_python && has_bbknn_python) {
    cat("\n--- Running BBKNN ---\n")
    cat("  Input: RNA$data (Module 03 normalization) -- BBKNN uses pre-normalized expression\n")
    tryCatch({
      temp_obj <- merged_raw

      # Ensure layers are joined before Python export
      tryCatch({
        temp_obj[["RNA"]] <- JoinLayers(temp_obj[["RNA"]])
      }, error = function(e) {
        cat("  Note: JoinLayers:", e$message, "\n")
      })

      var_features <- VariableFeatures(temp_obj)
      norm_data <- as.matrix(GetAssayData(temp_obj, assay = "RNA", layer = "data")[var_features, ])

      obs_df <- data.frame(
        batch = as.character(temp_obj$batch),
        row.names = colnames(temp_obj),
        stringsAsFactors = FALSE
      )
      var_df <- data.frame(gene = var_features, row.names = var_features, stringsAsFactors = FALSE)

      py_run_string("import numpy as np")
      py_run_string("import pandas as pd")
      py_run_string("import anndata")
      py_run_string("import scanpy as sc")
      py_run_string("import bbknn")

      py$norm_matrix <- r_to_py(t(norm_data))
      py$obs_data <- r_to_py(obs_df)
      py$var_data <- r_to_py(var_df)
      py$cell_names <- r_to_py(colnames(temp_obj))
      py$gene_names <- r_to_py(var_features)

      py_run_string("
adata_bbknn = anndata.AnnData(
    X=np.array(norm_matrix, dtype='float32'),
    obs=pd.DataFrame(obs_data),
    var=pd.DataFrame(var_data)
)
adata_bbknn.obs.index = list(cell_names)
adata_bbknn.var.index = list(gene_names)
adata_bbknn.obs['batch'] = adata_bbknn.obs['batch'].astype('category')
")

      n_pcs <- as.integer(params$dims_use)
      py_run_string(sprintf("sc.pp.pca(adata_bbknn, n_comps=%d)", n_pcs))
      py_run_string(sprintf("bbknn.bbknn(adata_bbknn, batch_key='batch', n_pcs=%d)", n_pcs))
      py_run_string("sc.tl.umap(adata_bbknn)")

      py_run_string("umap_result = adata_bbknn.obsm['X_umap']")
      py_run_string("pca_result = adata_bbknn.obsm['X_pca']")

      umap_coords <- as.matrix(py_to_r(py$umap_result))
      pca_coords <- as.matrix(py_to_r(py$pca_result))

      rownames(umap_coords) <- colnames(temp_obj)
      colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
      rownames(pca_coords) <- colnames(temp_obj)
      colnames(pca_coords) <- paste0("PC_", 1:ncol(pca_coords))

      bbknn_obj <- merged_raw
      bbknn_obj[["umap_bbknn"]] <- CreateDimReducObject(embeddings = umap_coords, key = "UMAP_", assay = DefaultAssay(merged_raw))
      bbknn_obj[["integrated.bbknn"]] <- CreateDimReducObject(embeddings = pca_coords, key = "bbknn_", assay = DefaultAssay(merged_raw))

      integration_results[["bbknn"]] <- bbknn_obj
      method_names <- c(method_names, "bbknn")
      cat("  [SUCCESS] BBKNN integration complete\n")

      py_run_string("del adata_bbknn, umap_result, pca_result, norm_matrix, obs_data, var_data, cell_names, gene_names")

    }, error = function(e) {
      print_integration_error("BBKNN", e, merged_raw)
    })
  }

  # ==========================================================================
  # scCobra - Deep learning batch correction
  # ==========================================================================
  if (isTRUE(params$run_sccobra) && has_sccobra_python) {
    cat("\n--- Running scCobra (deep learning batch correction) ---\n")
    cat("  Input: RNA$counts (raw counts) -- scCobra normalizes internally (processed=False)\n")
    tryCatch({
      temp_obj <- merged_raw

      # Ensure layers are joined before Python export
      tryCatch({
        temp_obj[["RNA"]] <- JoinLayers(temp_obj[["RNA"]])
      }, error = function(e) {
        cat("  Note: JoinLayers:", e$message, "\n")
      })

      counts_mat <- GetAssayData(temp_obj, assay = "RNA", layer = "counts")
      counts_mat <- as.matrix(counts_mat)
      counts_mat_t <- t(counts_mat)

      cat("  Counts matrix:", nrow(counts_mat_t), "cells x", ncol(counts_mat_t), "genes\n")

      obs_df <- data.frame(
        batch = as.character(temp_obj$batch),
        sex = as.character(temp_obj$sex),
        row.names = colnames(temp_obj),
        stringsAsFactors = FALSE
      )
      var_df <- data.frame(
        gene = rownames(counts_mat),
        row.names = rownames(counts_mat),
        stringsAsFactors = FALSE
      )

      cat("  Batch distribution:\n")
      print(table(obs_df$batch))

      cat("  Creating AnnData object...\n")

      sccobra_outdir <- file.path(tempdir(), "sccobra_output")
      dir.create(sccobra_outdir, showWarnings = FALSE, recursive = TRUE)

      py_run_string("import numpy as np")
      py_run_string("import pandas as pd")
      py_run_string("import anndata")

      py$counts_matrix <- r_to_py(counts_mat_t)
      py$obs_data <- r_to_py(obs_df)
      py$var_data <- r_to_py(var_df)
      py$cell_names <- r_to_py(colnames(temp_obj))
      py$gene_names <- r_to_py(rownames(counts_mat))
      py$outdir <- r_to_py(sccobra_outdir)

      py_run_string("
adata_cobra = anndata.AnnData(
    X=np.array(counts_matrix, dtype='float32'),
    obs=pd.DataFrame(obs_data),
    var=pd.DataFrame(var_data)
)
adata_cobra.obs.index = list(cell_names)
adata_cobra.var.index = list(gene_names)
adata_cobra.obs['batch'] = adata_cobra.obs['batch'].astype('category')
print(f'  AnnData shape: {adata_cobra.shape}')
print(f'  batch dtype: {adata_cobra.obs[\"batch\"].dtype}')
")

      cat("  Running scCobra (this may take several minutes)...\n")

      py_run_string("
import sys
sys.path.insert(0, '/scicore/home/doetsch/kaiser0001/GITHUB_repositories/scCobra')
from scCobra import scCobra as run_scCobra

print('  scCobra function imported successfully')

adata_integrated = run_scCobra(
    data_list=adata_cobra,
    batch_name='batch',
    profile='RNA',
    processed=False,
    min_features=200,
    min_cells=3,
    n_top_features=2000,
    batch_size=64,
    lr=2e-4,
    max_iteration=30000,
    seed=42,
    gpu=0,
    outdir=outdir,
    ignore_umap=True,
    verbose=True
)
print(f'  scCobra completed')
print(f'  Result type: {type(adata_integrated)}')
print(f'  Result shape: {adata_integrated.shape}')
")

      py_run_string("
obsm_keys = list(adata_integrated.obsm.keys())
print(f'  Available obsm keys: {obsm_keys}')
for key in obsm_keys:
    shape = adata_integrated.obsm[key].shape
    print(f'    {key}: shape {shape}')
")

      obsm_keys <- py_to_r(py$obsm_keys)
      cat("  Available obsm keys:", paste(obsm_keys, collapse = ", "), "\n")

      latent_key <- NULL
      possible_keys <- c("X_scCobra", "scCobra", "latent", "X_latent", "X_emb", "emb")
      for (key in possible_keys) {
        if (key %in% obsm_keys) {
          latent_key <- key
          cat("  Found latent key:", latent_key, "\n")
          break
        }
      }

      if (is.null(latent_key) && length(obsm_keys) > 0) {
        latent_key <- obsm_keys[1]
        cat("  Using first available obsm key:", latent_key, "\n")
      }

      if (is.null(latent_key)) {
        stop("Could not find latent embeddings in scCobra output")
      }

      py_run_string(sprintf("latent_rep = adata_integrated.obsm['%s']", latent_key))
      py_run_string("cell_names_out = list(adata_integrated.obs_names)")
      py_run_string("print(f'  Latent shape: {latent_rep.shape}')")

      latent_rep <- py_to_r(py$latent_rep)
      cell_names_out <- py_to_r(py$cell_names_out)

      latent_rep <- as.matrix(latent_rep)

      cat("  Latent dimensions:", dim(latent_rep), "\n")

      common_cells <- intersect(cell_names_out, colnames(temp_obj))

      if (length(common_cells) == 0) {
        if (nrow(latent_rep) == ncol(temp_obj)) {
          cat("  Cell counts match. Assuming same order.\n")
          rownames(latent_rep) <- colnames(temp_obj)
        } else {
          stop("Cell count mismatch and no common cell names")
        }
      } else {
        rownames(latent_rep) <- cell_names_out
        latent_rep <- latent_rep[colnames(temp_obj), , drop = FALSE]
      }

      colnames(latent_rep) <- paste0("scCobra_", 1:ncol(latent_rep))

      sccobra_obj <- merged_raw
      sccobra_obj[["integrated.sccobra"]] <- CreateDimReducObject(
        embeddings = latent_rep,
        key = "sccobra_",
        assay = DefaultAssay(merged_raw)
      )

      integration_results[["sccobra"]] <- sccobra_obj
      method_names <- c(method_names, "sccobra")
      cat("  [SUCCESS] scCobra integration complete\n")
      cat("  New reduction: integrated.sccobra with", ncol(latent_rep), "dimensions\n")

      py_run_string("
try:
    del adata_cobra, adata_integrated, counts_matrix, obs_data, var_data
    del cell_names, gene_names, outdir, obsm_keys, latent_rep, cell_names_out
except:
    pass
")
      unlink(sccobra_outdir, recursive = TRUE)

    }, error = function(e) {
      print_integration_error("scCobra", e, merged_raw)
    })
  } else {
    if (!isTRUE(params$run_sccobra)) {
      cat("\nscCobra skipped: params$run_sccobra is not TRUE\n")
    }
    if (!has_sccobra_python) {
      cat("\nscCobra skipped: has_sccobra_python is FALSE (module not loaded)\n")
    }
  }

  # ==========================================================================
  # CONCORD - Contrastive Learning for Single-Cell Data Integration
  # (Nature Biotechnology 2025)
  # ==========================================================================
  # Key parameters (from Concord.default_params):
  #   n_epochs=15, batch_size=256, lr=1e-2, latent_dim=100
  #   domain_key: batch variable for integration
  #   input_feature: list of features to use
  #   normalize_total, log1p: preprocessing flags (default False = pre-normalized)
  # ==========================================================================

  # Check for CONCORD availability
  has_concord_python <- tryCatch({
    py_run_string("import concord")
    cat("  [OK] concord\n")
    TRUE
  }, error = function(e) {
    cat("  [WARNING] concord not available:", e$message, "\n")
    cat("  Install with: pip install concord-sc\n")
    FALSE
  })

  if (isTRUE(params$run_concord) && has_concord_python) {
    cat("\n--- Running CONCORD ---\n")
    cat("  Input: RNA$counts (raw counts) -- CONCORD applies normalize_total + log1p internally\n")
    cat("  Key insight: Dataset-aware sampling prevents batch effect learning\n")

    tryCatch({
      temp_obj <- merged_raw

      # Ensure layers are joined before Python export
      tryCatch({
        temp_obj[["RNA"]] <- JoinLayers(temp_obj[["RNA"]])
      }, error = function(e) {
        cat("  Note: JoinLayers:", e$message, "\n")
      })

      # Get counts matrix
      counts_data <- GetAssayData(temp_obj, assay = "RNA", layer = "counts")
      counts_data <- as.matrix(counts_data)
      counts_data_t <- t(counts_data)  # Cells x Genes for Python

      cat("  Data dimensions:", nrow(counts_data_t), "cells x", ncol(counts_data_t), "genes\n")

      # Prepare metadata
      obs_df <- data.frame(
        batch = as.character(temp_obj$batch),
        row.names = colnames(temp_obj),
        stringsAsFactors = FALSE
      )

      cat("  Batch distribution:\n")
      print(table(obs_df$batch))

      # Export to Python
      py_run_string("import numpy as np")
      py_run_string("import pandas as pd")
      py_run_string("import anndata")
      py_run_string("import scanpy as sc")
      py_run_string("import torch")
      py_run_string("import concord as ccd")

      py$counts_matrix <- r_to_py(counts_data_t)
      py$obs_data <- r_to_py(obs_df)
      py$cell_names <- r_to_py(colnames(temp_obj))
      py$gene_names <- r_to_py(rownames(counts_data))

      # Set parameters from params (with defaults matching CONCORD defaults)
      n_top_features <- if (!is.null(params$concord_n_top_features)) params$concord_n_top_features else 2000
      latent_dim <- if (!is.null(params$concord_n_latent)) params$concord_n_latent else 100  # CONCORD default is 100
      n_epochs <- if (!is.null(params$concord_max_epochs)) params$concord_max_epochs else 15  # CONCORD default is 15
      batch_size <- if (!is.null(params$concord_batch_size)) params$concord_batch_size else 256  # CONCORD default
      lr <- if (!is.null(params$concord_lr)) params$concord_lr else 1e-2  # CONCORD default is 0.01

      py$n_top_features <- as.integer(n_top_features)
      py$latent_dim <- as.integer(latent_dim)
      py$n_epochs <- as.integer(n_epochs)
      py$batch_size <- as.integer(batch_size)
      py$lr <- lr

      # Create temporary directory for CONCORD outputs
      concord_outdir <- file.path(tempdir(), "concord_output")
      dir.create(concord_outdir, showWarnings = FALSE, recursive = TRUE)
      py$save_dir <- r_to_py(concord_outdir)

      # Create AnnData object
      cat("  Creating AnnData object...\n")
      py_run_string("
# Create AnnData from counts
adata = anndata.AnnData(
    X=np.array(counts_matrix, dtype='float32'),
    obs=pd.DataFrame(obs_data)
)
adata.obs.index = list(cell_names)
adata.var.index = list(gene_names)
adata.obs['batch'] = adata.obs['batch'].astype('category')

print(f'  AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes')
print(f'  Batches: {adata.obs[\"batch\"].nunique()} unique values')
")

      # Select variable features BEFORE creating Concord object
      cat("  Selecting variable features (seurat_v3 method)...\n")
      py_run_string("
# Select top variable features using CONCORD's utility function
feature_list = ccd.ul.select_features(adata, n_top_features=n_top_features, flavor='seurat_v3')
print(f'  Selected {len(feature_list)} variable features')

# Normalize and log-transform (CONCORD expects this by default)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
print('  Data normalized and log-transformed')
")

      # Initialize and run CONCORD
      # All training params go to __init__ via kwargs
      cat("  Initializing CONCORD model...\n")
      cat("    n_epochs:", n_epochs, "\n")
      cat("    batch_size:", batch_size, "\n")
      cat("    lr:", lr, "\n")
      cat("    latent_dim:", latent_dim, "\n")

      py_run_string("
# Determine device
if torch.cuda.is_available():
    device = torch.device('cuda:0')
    print('  Using CUDA GPU')
elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
    device = torch.device('mps')
    print('  Using Apple MPS')
else:
    device = torch.device('cpu')
    print('  Using CPU')

# Initialize Concord with all parameters
# Key params: domain_key enables batch-aware sampling, input_feature specifies genes
cur_ccd = ccd.Concord(
    adata=adata,
    save_dir=save_dir,
    verbose=True,
    # Training parameters (override defaults)
    input_feature=feature_list,
    domain_key='batch',
    latent_dim=latent_dim,
    n_epochs=n_epochs,
    batch_size=batch_size,
    lr=lr,
    device=device
)
print('  CONCORD model initialized')
")

      # Train CONCORD model
      cat("  Training CONCORD (this may take several minutes)...\n")
      py_run_string("
# fit_transform only needs output_key - all training params were in __init__
cur_ccd.fit_transform(output_key='Concord')
print('  Training complete!')
print(f'  Latent dimensions: {adata.obsm[\"Concord\"].shape[1]}')
")

      # Retrieve latent representation
      py_run_string("latent = adata.obsm['Concord']")
      py_run_string("cell_names_out = list(adata.obs_names)")

      latent_rep <- py_to_r(py$latent)
      cell_names_out <- py_to_r(py$cell_names_out)

      latent_rep <- as.matrix(latent_rep)
      rownames(latent_rep) <- cell_names_out
      colnames(latent_rep) <- paste0("CONCORD_", 1:ncol(latent_rep))

      # Ensure cell order matches original
      latent_rep <- latent_rep[colnames(temp_obj), , drop = FALSE]

      cat("  Final latent dimensions:", ncol(latent_rep), "\n")

      # Create Seurat object with CONCORD reduction
      concord_obj <- merged_raw
      concord_obj[["integrated.concord"]] <- CreateDimReducObject(
        embeddings = latent_rep,
        key = "concord_",
        assay = DefaultAssay(merged_raw)
      )

      integration_results[["concord"]] <- concord_obj
      method_names <- c(method_names, "concord")
      cat("  [SUCCESS] CONCORD integration complete\n")
      cat("  New reduction: integrated.concord with", ncol(latent_rep), "dimensions\n")

      # Clean up Python memory
      py_run_string("
try:
    del adata, cur_ccd, latent, counts_matrix, obs_data, feature_list
    del cell_names, gene_names, cell_names_out
    import gc
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
except:
    pass
")
      # Clean up temp directory
      unlink(concord_outdir, recursive = TRUE)

    }, error = function(e) {
      print_integration_error("CONCORD", e, merged_raw)
    })
  } else {
    if (!isTRUE(params$run_concord)) {
      cat("\nCONCORD skipped: params$run_concord is not TRUE\n")
    }
    if (!has_concord_python) {
      cat("\nCONCORD skipped: has_concord_python is FALSE (module not loaded)\n")
    }
  }

}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("INTEGRATION SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Methods completed:", length(method_names), "\n")
cat("Methods:", paste(method_names, collapse = ", "), "\n")
cat("Normalization used:", normalization_method, "\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# ==============================================================================
# Compute UMAPs
# ==============================================================================
cat("\n--- Computing UMAPs ---\n")

for (method_name in method_names) {
  obj <- integration_results[[method_name]]
  red_name <- get_integration_reduction(method_name)
  umap_name <- paste0("umap_", method_name)

  if (method_name == "bbknn" && "umap_bbknn" %in% names(obj@reductions)) next
  if (!red_name %in% names(obj@reductions)) next

  tryCatch({
    nd <- ncol(Embeddings(obj, reduction = red_name))
    dims_to_use <- min(params$dims_use, nd)
    obj <- RunUMAP(obj, reduction = red_name, dims = 1:dims_to_use,
                   reduction.name = umap_name, verbose = FALSE)
    integration_results[[method_name]] <- obj
    cat("  Created:", umap_name, "\n")
  }, error = function(e) cat("  UMAP failed for", method_name, ":", e$message, "\n"))
}

# ==============================================================================
# Compute benchmark scores (BATCH CORRECTION + BIOLOGICAL CONSERVATION)
# ==============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("COMPUTING BENCHMARK SCORES\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

benchmark_scores <- list()

for (method_name in method_names) {
  obj <- integration_results[[method_name]]
  red_name <- get_integration_reduction(method_name)

  if (!red_name %in% names(obj@reductions)) next

  cat("Scoring:", method_name, "\n")
  emb <- Embeddings(obj, reduction = red_name)
  batch_labels <- obj$batch

  # --- Batch correction metrics ---
  batch_variance <- compute_batch_variance(emb, batch_labels)
  batch_asw <- compute_batch_asw(emb, batch_labels)
  lisi <- compute_lisi_score(emb, batch_labels)

  # --- Biological conservation metrics (if annotations available) ---
  celltype_asw <- NA_real_
  celltype_lisi <- NA_real_
  nmi <- NA_real_
  ari <- NA_real_

  if (has_celltype_annotations) {
    cat("  Computing biological conservation metrics...\n")
    celltype_labels <- obj@meta.data[[celltype_column]]

    celltype_asw <- compute_celltype_asw(emb, celltype_labels)
    celltype_lisi <- compute_celltype_lisi(emb, celltype_labels)
    nmi <- compute_nmi_score(emb, celltype_labels)
    ari <- compute_ari_score(emb, celltype_labels)

    cat("    celltype_asw:", round(celltype_asw, 4), "\n")
    cat("    celltype_lisi:", round(celltype_lisi, 4), "\n")
    cat("    NMI:", round(nmi, 4), "\n")
    cat("    ARI:", round(ari, 4), "\n")
  }

  benchmark_scores[[method_name]] <- data.frame(
    method = method_name,
    reduction = red_name,
    # Batch correction metrics
    batch_variance = batch_variance,
    batch_asw = batch_asw,
    lisi = lisi,
    # Biological conservation metrics
    celltype_asw = celltype_asw,
    celltype_lisi = celltype_lisi,
    nmi = nmi,
    ari = ari,
    # Metadata
    batch_variable = params$batch_variable,
    celltype_column = ifelse(has_celltype_annotations, celltype_column, NA_character_),
    batch_integration_enabled = params$run_batch_integration,
    normalization_method = normalization_method
  )
}

benchmark_df <- do.call(rbind, benchmark_scores)
rownames(benchmark_df) <- NULL

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("BENCHMARK RESULTS (RAW VALUES)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("--- Batch Correction Metrics ---\n")
cat("  batch_variance: Lower = better (less batch effect in embedding)\n")
cat("  batch_asw: Closer to 0 = better (batches are well-mixed)\n")
cat("  lisi: Higher = better (more mixing across batches)\n\n")

print(benchmark_df[, c("method", "batch_variance", "batch_asw", "lisi")])

if (has_celltype_annotations) {
  cat("\n--- Biological Conservation Metrics ---\n")
  cat("  celltype_asw: Higher = better (cell types well-separated)\n")
  cat("  celltype_lisi: Lower = better (cell types locally pure)\n")
  cat("  nmi: Higher = better (clustering matches annotations)\n")
  cat("  ari: Higher = better (clustering matches annotations)\n\n")

  print(benchmark_df[, c("method", "celltype_asw", "celltype_lisi", "nmi", "ari")])
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n\n")

# ==============================================================================
# Select best method based on selection mode
# ==============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("METHOD SELECTION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Get selection mode from params (default to balanced for new behavior)
selection_mode <- params$integration_selection_mode
if (is.null(selection_mode) || selection_mode == "") {
  selection_mode <- "balanced"
}

cat("Selection mode:", selection_mode, "\n")
cat("Cell type annotations available:", has_celltype_annotations, "\n\n")

# ------------------------------------------------------------------------------
# Compute normalized scores for all metrics
# ------------------------------------------------------------------------------
# Use normalize_metric from functions.R (or define locally if not available)
if (!exists("normalize_metric")) {
  normalize_metric <- function(x, lower_is_better = TRUE) {
    if (all(is.na(x)) || length(unique(na.omit(x))) <= 1) {
      return(rep(0.5, length(x)))
    }
    x_range <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    if (x_range == 0) {
      return(rep(0.5, length(x)))
    }
    normalized <- (x - min(x, na.rm = TRUE)) / x_range
    if (lower_is_better) {
      normalized <- 1 - normalized
    }
    return(normalized)
  }
}

# --- Batch correction metrics normalization ---
# batch_variance: lower is better (less batch effect in embedding)
# batch_asw: CLOSER TO 0 is better (batches are mixed) - USE ABSOLUTE VALUE!
# lisi: higher is better (more mixing across batches)

benchmark_df$batch_variance_norm <- normalize_metric(benchmark_df$batch_variance, lower_is_better = TRUE)
# CRITICAL FIX: Use absolute value for batch_asw (closer to 0 is better)
benchmark_df$batch_asw_norm <- normalize_metric(abs(benchmark_df$batch_asw), lower_is_better = TRUE)
benchmark_df$lisi_norm <- normalize_metric(benchmark_df$lisi, lower_is_better = FALSE)

# Batch correction composite (average of all batch metrics, higher = better correction)
benchmark_df$batch_score <- rowMeans(
  benchmark_df[, c("batch_variance_norm", "batch_asw_norm", "lisi_norm")],
  na.rm = TRUE
)

# --- Biological conservation metrics normalization ---
# celltype_asw: higher is better (cell types well-separated)
# celltype_lisi: lower is better (cell types locally pure)
# nmi: higher is better (clustering matches annotations)
# ari: higher is better (clustering matches annotations)

benchmark_df$celltype_asw_norm <- NA_real_
benchmark_df$celltype_lisi_norm <- NA_real_
benchmark_df$nmi_norm <- NA_real_
benchmark_df$ari_norm <- NA_real_
benchmark_df$bio_score <- NA_real_

if (has_celltype_annotations) {
  benchmark_df$celltype_asw_norm <- normalize_metric(benchmark_df$celltype_asw, lower_is_better = FALSE)
  benchmark_df$celltype_lisi_norm <- normalize_metric(benchmark_df$celltype_lisi, lower_is_better = TRUE)
  benchmark_df$nmi_norm <- normalize_metric(benchmark_df$nmi, lower_is_better = FALSE)
  benchmark_df$ari_norm <- normalize_metric(benchmark_df$ari, lower_is_better = FALSE)

  # Biological conservation composite
  benchmark_df$bio_score <- rowMeans(
    benchmark_df[, c("celltype_asw_norm", "celltype_lisi_norm", "nmi_norm", "ari_norm")],
    na.rm = TRUE
  )
}

# --- Get weights for composite score ---
bio_weight <- if (!is.null(params$bio_weight)) params$bio_weight else 0.4
batch_weight <- if (!is.null(params$batch_weight)) params$batch_weight else 0.6

# Ensure weights sum to 1
total_weight <- bio_weight + batch_weight
bio_weight <- bio_weight / total_weight
batch_weight <- batch_weight / total_weight

# ------------------------------------------------------------------------------
# Select based on mode
# ------------------------------------------------------------------------------
if (selection_mode == "batch_removal") {
  # Original behavior: minimize batch variance (most aggressive batch removal)
  cat("Strategy: Minimize batch_variance (aggressive batch removal)\n")
  cat("Best for: Technical replicates, same biology processed separately\n")
  cat("WARNING: May over-correct if batches have biological meaning!\n\n")

  best_idx <- which.min(benchmark_df$batch_variance)
  selection_criterion <- "min(batch_variance)"
  benchmark_df$composite_score <- benchmark_df$batch_variance_norm

} else if (selection_mode == "balanced") {
  # Balanced: use composite score across all metrics
  cat("Strategy: Maximize composite score across batch + biological metrics\n")
  cat("Best for: Batches with biological meaning (different conditions/sexes/timepoints)\n")
  cat("This preserves biological variation while correcting batch effects.\n\n")

  if (has_celltype_annotations) {
    # scIB-style composite: 40% bio + 60% batch (configurable)
    cat(sprintf("Using weighted composite: %.0f%% biological + %.0f%% batch\n\n",
                bio_weight * 100, batch_weight * 100))

    benchmark_df$composite_score <- bio_weight * benchmark_df$bio_score +
      batch_weight * benchmark_df$batch_score

    selection_criterion <- sprintf("%.2f * bio_score + %.2f * batch_score", bio_weight, batch_weight)
  } else {
    # No biological annotations: use batch score only
    cat("No cell type annotations - using batch metrics only.\n\n")
    benchmark_df$composite_score <- benchmark_df$batch_score
    selection_criterion <- "max(mean(batch_variance_norm, batch_asw_norm, lisi_norm))"
  }

  best_idx <- which.max(benchmark_df$composite_score)

} else if (selection_mode == "conservative") {
  # Conservative: prioritize methods with moderate batch correction
  cat("Strategy: Conservative batch correction (avoid over-integration)\n")
  cat("Best for: Exploratory analysis, preserving subtle biological differences\n\n")

  if (has_celltype_annotations) {
    # Prioritize biological conservation
    benchmark_df$composite_score <- 0.6 * benchmark_df$bio_score + 0.4 * benchmark_df$batch_score
    selection_criterion <- "0.6 * bio_score + 0.4 * batch_score (bio-prioritized)"
  } else {
    # Use LISI as primary metric (higher = more mixing, but not too aggressive)
    benchmark_df$composite_score <- benchmark_df$lisi_norm
    selection_criterion <- "max(lisi_norm)"
  }

  best_idx <- which.max(benchmark_df$composite_score)

} else {
  cat("WARNING: Unknown selection mode '", selection_mode, "'. Using balanced.\n", sep = "")
  benchmark_df$composite_score <- benchmark_df$batch_score
  best_idx <- which.max(benchmark_df$composite_score)
  selection_criterion <- "max(batch_score)"
}

best_method <- benchmark_df$method[best_idx]
best_reduction <- benchmark_df$reduction[best_idx]

# ------------------------------------------------------------------------------
# Print detailed selection summary
# ------------------------------------------------------------------------------
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("NORMALIZED SCORES (0-1 scale, higher = better)\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

# Create ranking dataframe
ranking_df <- benchmark_df[order(-benchmark_df$composite_score), ]
ranking_df$rank <- 1:nrow(ranking_df)

# Print batch metrics normalized
cat("--- Batch Correction (normalized) ---\n")
batch_norm_cols <- c("rank", "method", "batch_variance_norm", "batch_asw_norm", "lisi_norm", "batch_score")
print(ranking_df[, batch_norm_cols], row.names = FALSE, digits = 3)

if (has_celltype_annotations) {
  cat("\n--- Biological Conservation (normalized) ---\n")
  bio_norm_cols <- c("rank", "method", "celltype_asw_norm", "celltype_lisi_norm", "nmi_norm", "ari_norm", "bio_score")
  print(ranking_df[, bio_norm_cols], row.names = FALSE, digits = 3)
}

cat("\n--- Final Composite Scores ---\n")
final_cols <- c("rank", "method", "batch_score")
if (has_celltype_annotations) {
  final_cols <- c(final_cols, "bio_score")
}
final_cols <- c(final_cols, "composite_score")
print(ranking_df[, final_cols], row.names = FALSE, digits = 3)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SELECTION RESULT\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Selection mode:", selection_mode, "\n")
cat("Criterion:", selection_criterion, "\n")
cat("Cell type annotations:", ifelse(has_celltype_annotations, celltype_column, "not available"), "\n")
cat("\n>>> SELECTED METHOD:", best_method, "<<<\n")
cat(">>> Selected reduction:", best_reduction, "<<<\n")
cat(">>> Composite score:", round(benchmark_df$composite_score[best_idx], 4), "<<<\n\n")

# Provide interpretation
cat("--- Interpretation ---\n")
if (best_method == "unintegrated") {
  cat("The unintegrated baseline scored highest. This may indicate:\n")
  cat("  - Batch effects are minimal in your data\n")
  cat("  - Integration methods may be over-correcting\n")
  cat("  - Consider visual inspection of UMAPs\n\n")
} else {
  best_batch_score <- benchmark_df$batch_score[best_idx]
  best_bio_score <- benchmark_df$bio_score[best_idx]

  cat(sprintf("'%s' achieved:\n", best_method))
  cat(sprintf("  - Batch correction score: %.3f (1.0 = perfect batch removal)\n", best_batch_score))
  if (has_celltype_annotations) {
    cat(sprintf("  - Biological conservation score: %.3f (1.0 = perfect preservation)\n", best_bio_score))
  }

  # Compare to most aggressive method
  most_aggressive_idx <- which.min(benchmark_df$batch_variance)
  most_aggressive <- benchmark_df$method[most_aggressive_idx]
  if (best_method != most_aggressive) {
    cat(sprintf("\nNote: '%s' has lower batch_variance (%.4f vs %.4f) but\n",
                most_aggressive, benchmark_df$batch_variance[most_aggressive_idx],
                benchmark_df$batch_variance[best_idx]))
    cat(sprintf("'%s' was selected for better overall balance.\n", best_method))
  }
}

# Provide guidance based on selection
if (selection_mode == "batch_removal") {
  cat("\nNOTE: You selected 'batch_removal' mode which aggressively removes batch effects.\n")
  cat("      If your batches have biological meaning (e.g., different sexes, conditions),\n")
  cat("      consider using: params$integration_selection_mode <- 'balanced'\n\n")
}

cat(paste(rep("=", 60), collapse = ""), "\n\n")

# ==============================================================================
# PER-LABEL BATCH MIXING VALIDATION (added 2026-02-09)
# ==============================================================================
# Evaluates batch integration quality WITHIN each cell type label.
# Global metrics can be misleading — this checks whether batches overlap
# within each cell type, not just globally.
# ==============================================================================

run_perlabel_validation <- TRUE   # Set to FALSE to skip

if (run_perlabel_validation &&
    has_celltype_annotations &&
    !is.null(celltype_column) &&
    celltype_column %in% colnames(merged_raw@meta.data)) {

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("PER-LABEL BATCH MIXING VALIDATION\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Cell type column:", celltype_column, "\n")
  cat("Batch variable:", params$batch_variable, "\n\n")

  # ------------------------------------------------------------------
  # Step 1: Compute per-label batch mixing for all methods
  # ------------------------------------------------------------------
  # Loop through integration_results (each method has its own object)
  # to match the same pattern used in the main benchmark scoring loop.
  # ------------------------------------------------------------------

  perlabel_all_results <- list()

  for (method_name in method_names) {
    obj <- integration_results[[method_name]]
    red_name <- get_integration_reduction(method_name)

    if (!red_name %in% names(obj@reductions)) {
      cat("  [SKIP]", method_name, "- reduction not found:", red_name, "\n")
      next
    }

    cat("\n  Evaluating:", method_name, "(reduction:", red_name, ")\n")

    # Extract embedding
    emb <- Embeddings(obj, reduction = red_name)
    max_dims <- min(params$dims_use, ncol(emb))
    emb <- emb[, 1:max_dims, drop = FALSE]

    # Compute per-label metrics
    method_results <- compute_perlabel_batch_mixing(
      embedding = emb,
      batch_labels = obj$batch,
      celltype_labels = obj@meta.data[[celltype_column]],
      min_cells_per_label = 30,
      min_batches_per_label = 2
    )
    method_results$method <- method_name

    perlabel_all_results[[method_name]] <- method_results
  }

  # Combine all results
  perlabel_results <- do.call(rbind, perlabel_all_results)
  rownames(perlabel_results) <- NULL

  # ------------------------------------------------------------------
  # Step 2: Summarize results
  # ------------------------------------------------------------------
  perlabel_summary <- summarize_perlabel_batch_mixing(perlabel_results)

  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("PER-LABEL BATCH MIXING SUMMARY\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")

  cat("--- Method Summary (ordered by mean within-label LISI) ---\n")
  print(as.data.frame(perlabel_summary$summary), row.names = FALSE)

  if (nrow(perlabel_summary$flagged_labels) > 0) {
    cat("\n--- Flagged Cell Types (WARNING or FAIL) ---\n")
    print(as.data.frame(perlabel_summary$flagged_labels), row.names = FALSE)
  } else {
    cat("\n  No cell types flagged -- all labels show acceptable batch mixing.\n")
  }

  cat("\n--- Worst-Performing Cell Types Per Method (bottom 3 by LISI) ---\n")
  print(as.data.frame(perlabel_summary$worst_labels), row.names = FALSE)

  # ------------------------------------------------------------------
  # Step 3: Validate the selected method
  # ------------------------------------------------------------------
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("VALIDATION OF SELECTED METHOD:", best_method, "\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")

  if (best_method %in% perlabel_summary$summary$method) {
    sel_summary <- perlabel_summary$summary %>%
      dplyr::filter(method == best_method)

    cat("  Labels evaluated:", sel_summary$n_labels_evaluated, "\n")
    cat("  Mean within-label LISI:", round(sel_summary$mean_lisi, 4), "\n")
    cat("  Mean within-label batch ASW:", round(sel_summary$mean_batch_asw, 4), "\n")
    cat("  Labels OK:", sel_summary$n_ok, "/", sel_summary$n_labels_evaluated,
        "(", sel_summary$pct_ok, "%)\n")
    cat("  Labels WARNING:", sel_summary$n_warning, "\n")
    cat("  Labels FAIL:", sel_summary$n_fail, "\n")

    # Check if the selected method ranks well on per-label metrics
    lisi_rank <- which(perlabel_summary$summary$method == best_method)
    cat("  Per-label LISI rank:", lisi_rank, "of", nrow(perlabel_summary$summary), "\n")

    if (sel_summary$n_fail > 0) {
      cat("\n  *** ATTENTION: Selected method has", sel_summary$n_fail,
          "cell type(s) with FAILED batch mixing ***\n")
      failed_labels <- perlabel_summary$flagged_labels %>%
        dplyr::filter(method == best_method, flag == "FAIL")
      cat("  Failed labels:\n")
      for (i in seq_len(nrow(failed_labels))) {
        cat("    -", failed_labels$label[i],
            "(", failed_labels$n_cells[i], "cells,",
            "LISI:", round(failed_labels$lisi[i], 3),
            "batch_ASW:", round(failed_labels$batch_asw[i], 3), ")\n")
      }

      # Suggest alternatives that do better on failed labels
      cat("\n  Checking if alternative methods perform better on failed labels...\n")
      for (fl in failed_labels$label) {
        fl_scores <- perlabel_results %>%
          dplyr::filter(label == fl, !flag %in% c("SKIPPED_FEW_CELLS", "SKIPPED_SINGLE_BATCH")) %>%
          dplyr::arrange(dplyr::desc(lisi))
        if (nrow(fl_scores) > 0) {
          cat("    ", fl, "- best method:", fl_scores$method[1],
              "(LISI:", round(fl_scores$lisi[1], 3), ")\n")
        }
      }
    } else {
      cat("\n  VALIDATION PASSED: No cell types with failed batch mixing.\n")
    }
  }

  # ------------------------------------------------------------------
  # Step 4: Generate plots
  # ------------------------------------------------------------------
  if (!is.null(output_dirs$integration_benchmark)) {
    perlabel_plot_dir <- file.path(output_dirs$integration_benchmark, "per_label_batch_mixing")
  } else {
    perlabel_plot_dir <- file.path(output_dirs$benchmarking, "per_label_batch_mixing")
  }
  dir.create(perlabel_plot_dir, showWarnings = FALSE, recursive = TRUE)

  # Plot 1: Heatmap of within-label LISI (all methods x all labels)
  p_heatmap_lisi <- plot_perlabel_heatmap(
    perlabel_results,
    metric = "lisi",
    title = "Within-Label Batch LISI by Integration Method"
  )
  if (!is.null(p_heatmap_lisi)) {
    save_plot_multi(p_heatmap_lisi, "01_perlabel_lisi_heatmap",
                    output_dir = perlabel_plot_dir,
                    width = 12,
                    height = max(6, 0.4 * length(unique(perlabel_results$label))))
  }

  # Plot 2: Heatmap of within-label batch ASW
  p_heatmap_asw <- plot_perlabel_heatmap(
    perlabel_results,
    metric = "batch_asw",
    title = "Within-Label Batch ASW by Integration Method"
  )
  if (!is.null(p_heatmap_asw)) {
    save_plot_multi(p_heatmap_asw, "02_perlabel_batch_asw_heatmap",
                    output_dir = perlabel_plot_dir,
                    width = 12,
                    height = max(6, 0.4 * length(unique(perlabel_results$label))))
  }

  # Plot 3: Flag summary (stacked bar chart)
  p_flags <- plot_perlabel_flag_summary(perlabel_summary$summary)
  if (!is.null(p_flags)) {
    save_plot_multi(p_flags, "03_perlabel_flag_summary",
                    output_dir = perlabel_plot_dir, width = 10, height = 6)
  }

  # Plot 4: LISI comparison for top 3 methods
  top3_methods <- head(perlabel_summary$summary$method, 3)
  p_top3 <- plot_perlabel_lisi_comparison(perlabel_results, top3_methods)
  if (!is.null(p_top3)) {
    save_plot_multi(p_top3, "04_perlabel_lisi_top3_comparison",
                    output_dir = perlabel_plot_dir,
                    width = 10,
                    height = max(6, 0.3 * length(unique(perlabel_results$label))))
  }

  # ------------------------------------------------------------------
  # Step 5: Save results
  # ------------------------------------------------------------------
  write.csv(perlabel_results,
            file.path(perlabel_plot_dir, "perlabel_batch_mixing_all_methods.csv"),
            row.names = FALSE)

  write.csv(as.data.frame(perlabel_summary$summary),
            file.path(perlabel_plot_dir, "perlabel_method_summary.csv"),
            row.names = FALSE)

  if (nrow(perlabel_summary$flagged_labels) > 0) {
    write.csv(as.data.frame(perlabel_summary$flagged_labels),
              file.path(perlabel_plot_dir, "perlabel_flagged_labels.csv"),
              row.names = FALSE)
  }

  cat("\nPer-label batch mixing results saved to:", perlabel_plot_dir, "\n")

  # ------------------------------------------------------------------
  # Step 6: Save to environment for downstream use
  # ------------------------------------------------------------------
  perlabel_batch_mixing <- list(
    results = perlabel_results,
    summary = perlabel_summary,
    selected_method_validation = list(
      method = best_method,
      passed = ifelse(
        best_method %in% perlabel_summary$summary$method,
        perlabel_summary$summary$n_fail[perlabel_summary$summary$method == best_method] == 0,
        NA
      )
    )
  )

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("PER-LABEL VALIDATION COMPLETE\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
}

# ==============================================================================
# OPTIONAL - OVERRIDE SELECTION BASED ON PER-LABEL RESULTS
# ==============================================================================
#
# If the selected method fails per-label validation, you can optionally
# override the selection. Uncomment to enable automatic override:
# ==============================================================================

# if (exists("perlabel_batch_mixing") &&
#     !is.null(perlabel_batch_mixing$selected_method_validation$passed) &&
#     !perlabel_batch_mixing$selected_method_validation$passed) {
#
#   cat("\n*** OVERRIDE: Selected method failed per-label validation ***\n")
#
#   # Find best method that passes per-label validation
#   passing_methods <- perlabel_summary$summary %>%
#     dplyr::filter(n_fail == 0) %>%
#     dplyr::arrange(dplyr::desc(mean_lisi))
#
#   if (nrow(passing_methods) > 0) {
#     # Cross-reference with global composite scores
#     override_candidates <- passing_methods$method
#     cat("  Candidates passing per-label validation:", paste(override_candidates, collapse = ", "), "\n")
#
#     # Select the candidate with the best global composite score
#     global_scores <- benchmark_df %>%
#       dplyr::filter(method %in% override_candidates) %>%
#       dplyr::arrange(dplyr::desc(composite_score))
#
#     if (nrow(global_scores) > 0) {
#       new_selection <- global_scores$method[1]
#       cat("  Overriding selection:", best_method, "->", new_selection, "\n")
#       cat("  Reason: Better per-label batch mixing with acceptable global composite\n")
#       best_method <- new_selection
#       best_reduction <- global_scores$reduction[1]
#       best_idx <- which(benchmark_df$method == best_method)
#     }
#   } else {
#     cat("  No methods pass all per-label checks. Keeping original selection.\n")
#     cat("  Consider adjusting flag thresholds or inspecting problematic cell types.\n")
#   }
# }

# ==============================================================================
# END OF PER-LABEL BATCH MIXING VALIDATION
# ==============================================================================

multi_integrated <- integration_results[[best_method]]

# Join layers for downstream use
tryCatch({
  multi_integrated[["RNA"]] <- JoinLayers(multi_integrated[["RNA"]])
}, error = function(e) {
  cat("Note: JoinLayers for final object:", e$message, "\n")
})

umap_name <- paste0("umap_", best_method)
if (umap_name %in% names(multi_integrated@reductions)) {
  multi_integrated[["umap"]] <- multi_integrated[[umap_name]]
}

# ==============================================================================
# Save results
# ==============================================================================
cat("--- Saving results ---\n")

# Save benchmark scores with all computed metrics
write.csv(benchmark_df, file.path(output_dirs$benchmarking, "benchmark_scores.csv"), row.names = FALSE)

# Save selection summary
selection_summary <- data.frame(
  selection_mode = selection_mode,
  selection_criterion = selection_criterion,
  selected_method = best_method,
  selected_reduction = best_reduction,
  selected_composite_score = benchmark_df$composite_score[best_idx],
  selected_batch_score = benchmark_df$batch_score[best_idx],
  selected_bio_score = ifelse(has_celltype_annotations, benchmark_df$bio_score[best_idx], NA_real_),
  n_methods_tested = nrow(benchmark_df),
  batch_variable = params$batch_variable,
  celltype_column = ifelse(has_celltype_annotations, celltype_column, NA_character_),
  has_celltype_annotations = has_celltype_annotations,
  bio_weight = bio_weight,
  batch_weight = batch_weight,
  normalization_method = normalization_method,
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
write.csv(selection_summary, file.path(output_dirs$benchmarking, "selection_summary.csv"), row.names = FALSE)

# ==============================================================================
# Create detailed text report
# ==============================================================================
report_file <- file.path(output_dirs$benchmarking, "04_detailed_benchmark_report.txt")
sink(report_file)

cat("================================================================================\n")
cat("INTEGRATION BENCHMARK DETAILED REPORT\n")
cat("================================================================================\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("1. CONFIGURATION\n")
cat("--------------------------------------------------------------------------------\n")
cat("Batch variable:", params$batch_variable, "\n")
cat("Batch integration enabled:", params$run_batch_integration, "\n")
cat("Normalization method:", normalization_method, "\n")
cat("Selection mode:", selection_mode, "\n")
cat("Cell type annotations:", ifelse(has_celltype_annotations, celltype_column, "not available"), "\n")
cat("Bio weight:", bio_weight, "\n")
cat("Batch weight:", batch_weight, "\n")
cat("Methods tested:", paste(method_names, collapse = ", "), "\n\n")

cat("2. BATCH CORRECTION METRICS\n")
cat("--------------------------------------------------------------------------------\n")
cat("Metric interpretation:\n")
cat("  batch_variance: Variance explained by batch in embedding (lower = better)\n")
cat("  batch_asw: Batch silhouette width, closer to 0 = better mixing\n")
cat("  lisi: Local Inverse Simpson Index, higher = more batch mixing\n\n")

cat("Raw values:\n")
print(benchmark_df[order(benchmark_df$batch_variance), c("method", "batch_variance", "batch_asw", "lisi")])
cat("\n")

cat("Normalized scores (0-1, higher = better):\n")
print(benchmark_df[order(-benchmark_df$batch_score), c("method", "batch_variance_norm", "batch_asw_norm", "lisi_norm", "batch_score")])
cat("\n")

if (has_celltype_annotations) {
  cat("3. BIOLOGICAL CONSERVATION METRICS\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Metric interpretation:\n")
  cat("  celltype_asw: Cell type silhouette width, higher = better separation\n")
  cat("  celltype_lisi: Cell type LISI, lower = better local purity\n")
  cat("  nmi: Normalized Mutual Information with annotations, higher = better\n")
  cat("  ari: Adjusted Rand Index with annotations, higher = better\n\n")

  cat("Raw values:\n")
  print(benchmark_df[order(-benchmark_df$celltype_asw), c("method", "celltype_asw", "celltype_lisi", "nmi", "ari")])
  cat("\n")

  cat("Normalized scores (0-1, higher = better):\n")
  print(benchmark_df[order(-benchmark_df$bio_score), c("method", "celltype_asw_norm", "celltype_lisi_norm", "nmi_norm", "ari_norm", "bio_score")])
  cat("\n")
}

cat("4. FINAL RANKING\n")
cat("--------------------------------------------------------------------------------\n")
cat("Criterion:", selection_criterion, "\n\n")

ranking_print <- benchmark_df[order(-benchmark_df$composite_score),
                               c("method", "batch_score")]
if (has_celltype_annotations) {
  ranking_print$bio_score <- benchmark_df$bio_score[order(-benchmark_df$composite_score)]
}
ranking_print$composite_score <- benchmark_df$composite_score[order(-benchmark_df$composite_score)]
ranking_print$rank <- 1:nrow(ranking_print)
ranking_print <- ranking_print[, c("rank", setdiff(names(ranking_print), "rank"))]

print(ranking_print)
cat("\n")

cat("5. SELECTION RESULT\n")
cat("--------------------------------------------------------------------------------\n")
cat(">>> SELECTED METHOD:", best_method, "<<<\n")
cat(">>> Selected reduction:", best_reduction, "<<<\n")
cat(">>> Composite score:", round(benchmark_df$composite_score[best_idx], 4), "<<<\n\n")

cat("6. RECOMMENDATIONS\n")
cat("--------------------------------------------------------------------------------\n")
if (!has_celltype_annotations) {
  cat("- Cell type annotations were not available. Adding annotations would enable\n")
  cat("  biological conservation metrics for better method selection.\n")
  cat("- Set params$celltype_column or add 'cell_type' column to metadata.\n\n")
}

if (selection_mode == "batch_removal") {
  cat("- 'batch_removal' mode aggressively removes batch effects.\n")
  cat("- If batches have biological meaning, consider 'balanced' mode.\n\n")
}

cat("- Always visually inspect UMAPs to verify integration quality.\n")
cat("- Check that expected cell types cluster correctly after integration.\n")
cat("- Consider downstream analysis goals when evaluating integration.\n\n")

# Add per-label validation results to the report
if (exists("perlabel_batch_mixing")) {
  cat("7. PER-LABEL BATCH MIXING VALIDATION\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("This section evaluates batch mixing WITHIN each cell type label.\n\n")

  cat("Method Summary (ordered by mean within-label LISI):\n")
  print(as.data.frame(perlabel_summary$summary))
  cat("\n")

  if (nrow(perlabel_summary$flagged_labels) > 0) {
    cat("Flagged Cell Types (WARNING or FAIL):\n")
    print(as.data.frame(perlabel_summary$flagged_labels))
    cat("\n")
  } else {
    cat("No cell types flagged - all labels show acceptable batch mixing.\n\n")
  }

  cat("Selected method validation:\n")
  cat("  Method:", best_method, "\n")
  cat("  Passed:", perlabel_batch_mixing$selected_method_validation$passed, "\n\n")
}

cat("================================================================================\n")
cat("END OF REPORT\n")
cat("================================================================================\n")

sink()
cat("  Saved detailed report: 04_detailed_benchmark_report.txt\n")

# Save RDS files
saveRDS(multi_integrated, file.path(output_dirs$objects, "integrated_multi_methods.rds"))
saveRDS(integration_results, file.path(output_dirs$objects, "all_integration_results.rds"))

# Include per-label results in the RData save if they exist
if (exists("perlabel_batch_mixing")) {
  integration_file <- file.path(output_dirs$objects, "04_integration_data.RData")
  save(multi_integrated, best_method, best_reduction, benchmark_df,
       normalization_method, selection_mode, selection_criterion,
       has_celltype_annotations, celltype_column,
       perlabel_batch_mixing, file = integration_file)
} else {
  integration_file <- file.path(output_dirs$objects, "04_integration_data.RData")
  save(multi_integrated, best_method, best_reduction, benchmark_df,
       normalization_method, selection_mode, selection_criterion,
       has_celltype_annotations, celltype_column, file = integration_file)
}

cat("  Saved benchmark_scores.csv (with normalized metrics)\n")
cat("  Saved selection_summary.csv\n")
cat("  Saved integrated_multi_methods.rds\n")
cat("  Saved all_integration_results.rds\n")
cat("  Saved 04_integration_data.RData\n")

# ==============================================================================
# Visualizations
# ==============================================================================
cat("\n--- Creating visualizations ---\n")

umap_plots <- list()
for (method_name in method_names) {
  obj <- integration_results[[method_name]]
  umap_name <- paste0("umap_", method_name)
  if (!umap_name %in% names(obj@reductions)) next

  p <- DimPlot(obj, reduction = umap_name, group.by = "sex",
               cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
    ggtitle(method_name) + theme(legend.position = "bottom")
  umap_plots[[method_name]] <- p
}

if (length(umap_plots) > 0) {
  p_combined <- wrap_plots(umap_plots, ncol = min(3, length(umap_plots)))
  save_plot_multi(p_combined, "01_integration_comparison_UMAP",
                  output_dir = output_dirs$benchmarking, width = 14, height = 10)
}

# Raw metrics barplot
benchmark_long <- reshape2::melt(benchmark_df[, c("method", "batch_variance", "batch_asw", "lisi")],
                                  id.vars = "method", variable.name = "metric", value.name = "value")

p_benchmark <- ggplot(benchmark_long, aes(x = method, y = value, fill = method)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste0("Integration Benchmark (", normalization_method, ")",
                      if (isTRUE(params$run_batch_integration)) paste0(" - batch: ", params$batch_variable) else " - no batch correction"),
       subtitle = paste0("Selection mode: ", selection_mode, " | Selected: ", best_method))

save_plot_multi(p_benchmark, "02_benchmark_scores_barplot",
                output_dir = output_dirs$benchmarking, width = 12, height = 8)

# Composite score ranking plot
p_composite <- ggplot(benchmark_df, aes(x = reorder(method, composite_score), y = composite_score, fill = method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.3f", composite_score)), hjust = -0.1, size = 3) +
  coord_flip() +
  theme_minimal() +
  labs(title = paste0("Integration Method Ranking (", selection_mode, " mode)"),
       subtitle = paste0("Criterion: ", selection_criterion),
       x = "Method", y = "Composite Score (higher = better)") +
  theme(legend.position = "none") +
  ylim(0, max(benchmark_df$composite_score) * 1.15)

save_plot_multi(p_composite, "03_composite_score_ranking",
                output_dir = output_dirs$benchmarking, width = 10, height = 6)

# NEW: Composite score breakdown (stacked bar chart showing batch vs bio contributions)
if (has_celltype_annotations) {
  breakdown_df <- data.frame(
    method = rep(benchmark_df$method, 2),
    component = rep(c("Batch Correction", "Biological Conservation"), each = nrow(benchmark_df)),
    contribution = c(
      batch_weight * benchmark_df$batch_score,
      bio_weight * benchmark_df$bio_score
    )
  )

  # Order by composite score
  method_order <- benchmark_df$method[order(-benchmark_df$composite_score)]
  breakdown_df$method <- factor(breakdown_df$method, levels = method_order)

  p_breakdown <- ggplot(breakdown_df, aes(x = method, y = contribution, fill = component)) +
    geom_bar(stat = "identity", position = "stack") +
    coord_flip() +
    scale_fill_manual(values = c("Batch Correction" = "#4DAF4A", "Biological Conservation" = "#984EA3")) +
    theme_minimal() +
    labs(
      title = "Composite Score Breakdown by Component",
      subtitle = sprintf("Weights: %.0f%% Biological + %.0f%% Batch", bio_weight * 100, batch_weight * 100),
      x = "Method",
      y = "Weighted Score Contribution",
      fill = "Component"
    ) +
    theme(legend.position = "bottom")

  save_plot_multi(p_breakdown, "04_composite_score_breakdown",
                  output_dir = output_dirs$benchmarking, width = 10, height = 6)
}

# NEW: Biological metrics barplot (if available)
if (has_celltype_annotations) {
  bio_long <- reshape2::melt(benchmark_df[, c("method", "celltype_asw", "celltype_lisi", "nmi", "ari")],
                              id.vars = "method", variable.name = "metric", value.name = "value")

  p_bio <- ggplot(bio_long, aes(x = method, y = value, fill = method)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ metric, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Biological Conservation Metrics",
         subtitle = paste0("Cell type column: ", celltype_column))

  save_plot_multi(p_bio, "05_biological_conservation_metrics",
                  output_dir = output_dirs$benchmarking, width = 12, height = 8)
}

# NEW: Normalized scores heatmap
norm_cols <- c("batch_variance_norm", "batch_asw_norm", "lisi_norm")
if (has_celltype_annotations) {
  norm_cols <- c(norm_cols, "celltype_asw_norm", "celltype_lisi_norm", "nmi_norm", "ari_norm")
}

heatmap_df <- benchmark_df[, c("method", norm_cols)]
heatmap_long <- reshape2::melt(heatmap_df, id.vars = "method", variable.name = "metric", value.name = "score")

# Clean up metric names for display
heatmap_long$metric <- gsub("_norm$", "", heatmap_long$metric)
heatmap_long$metric <- gsub("_", " ", heatmap_long$metric)

# Order methods by composite score
method_order <- benchmark_df$method[order(-benchmark_df$composite_score)]
heatmap_long$method <- factor(heatmap_long$method, levels = rev(method_order))

p_heatmap <- ggplot(heatmap_long, aes(x = metric, y = method, fill = score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", score)), size = 3) +
  scale_fill_gradient2(low = "#D73027", mid = "#FFFFBF", high = "#1A9850",
                       midpoint = 0.5, limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Normalized Scores Heatmap (0-1 scale, higher = better)",
       x = "Metric", y = "Method", fill = "Score")

save_plot_multi(p_heatmap, "06_normalized_scores_heatmap",
                output_dir = output_dirs$benchmarking, width = 10, height = 8)

write_readme(output_dirs$benchmarking, "Integration Benchmarking",
             paste0("Multi-method integration comparison with comprehensive benchmarking.\n\n",
                    "=== CONFIGURATION ===\n",
                    "Normalization method: ", normalization_method, "\n",
                    "Batch variable: ", params$batch_variable, "\n",
                    "Batch integration enabled: ", params$run_batch_integration, "\n",
                    "Cell type annotations: ", ifelse(has_celltype_annotations, celltype_column, "not available"), "\n\n",
                    "=== SELECTION ===\n",
                    "Mode: ", selection_mode, "\n",
                    "Criterion: ", selection_criterion, "\n",
                    "Best method: ", best_method, "\n",
                    "Composite score: ", round(benchmark_df$composite_score[best_idx], 4), "\n\n",
                    "=== SELECTION MODES ===\n",
                    "- batch_removal: Minimizes batch_variance (most aggressive)\n",
                    "  Use when batches are technical replicates (same biology)\n\n",
                    "- balanced: Maximizes weighted composite of batch + bio metrics\n",
                    "  Default weights: 40% biological, 60% batch (configurable)\n",
                    "  Use when batches have biological meaning (recommended)\n\n",
                    "- conservative: Prioritizes biological conservation\n",
                    "  Use when preserving subtle biological differences is critical\n\n",
                    "=== METRICS ===\n",
                    "BATCH CORRECTION (higher normalized = better correction):\n",
                    "- batch_variance: Lower raw = less batch effect in embedding\n",
                    "- batch_asw: Closer to 0 = batches well-mixed\n",
                    "- lisi: Higher = more mixing across batches\n\n",
                    "BIOLOGICAL CONSERVATION (higher normalized = better preservation):\n",
                    "- celltype_asw: Higher = cell types well-separated\n",
                    "- celltype_lisi: Lower = cell types locally pure\n",
                    "- nmi: Higher = clustering matches annotations\n",
                    "- ari: Higher = clustering matches annotations\n\n",
                    "Methods tested: ", paste(method_names, collapse = ", "), "\n"),
             list("benchmark_scores.csv" = "All metrics (raw + normalized) for each method",
                  "selection_summary.csv" = "Selection configuration and result",
                  "04_detailed_benchmark_report.txt" = "Comprehensive text report",
                  "01_integration_comparison_UMAP.png" = "UMAP comparison across methods",
                  "02_benchmark_scores_barplot.png" = "Raw batch metrics barplot",
                  "03_composite_score_ranking.png" = "Methods ranked by composite score",
                  "04_composite_score_breakdown.png" = "Stacked bar showing batch vs bio contributions",
                  "05_biological_conservation_metrics.png" = "Raw biological metrics (if available)",
                  "06_normalized_scores_heatmap.png" = "Heatmap of all normalized scores",
                  "per_label_batch_mixing/" = "Per-cell-type batch mixing validation results"))

cat("\n>>> MODULE 04 COMPLETE <<<\n")