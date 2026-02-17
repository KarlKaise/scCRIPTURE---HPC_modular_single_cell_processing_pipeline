#!/usr/bin/env Rscript
# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
#
# This file contains all shared utility functions used across pipeline modules.
#
# UPDATES IN THIS VERSION:
# 1. Added layer detection functions for Seurat v5 compatibility
# 2. Added integration error helper functions
# 3. Added normalization method selection helper
# 4. Enhanced print statements for debugging
# 5. 2026-01-03: Added sce_to_seurat_v5() for scran normalization
# 6. 2026-01-03: Enhanced seurat_to_sce() with better v5 handling
# 7. 2026-01-06: Moved params.R functions here (load_sample_sheet, validate_params, etc.)
# 8. 2026-02-05: Fixed run_integration_for_benchmark() - removed early-return for
#    existing Harmony embeddings so benchmarking always recomputes on subsampled data
#
# ==============================================================================

# ==============================================================================
# SAVE PLOT IN MULTIPLE FORMATS
# ==============================================================================
save_plot_multi <- function(plot, filename, output_dir, width = 10, height = 8, dpi = 300) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  formats <- c("png", "tiff", "pdf", "svg")
  saved_files <- c()

  for (fmt in formats) {
    filepath <- file.path(output_dir, paste0(filename, ".", fmt))
    tryCatch({
      if (fmt == "tiff") {
        ggplot2::ggsave(filepath, plot, width = width, height = height, dpi = dpi, compression = "lzw")
      } else {
        ggplot2::ggsave(filepath, plot, width = width, height = height, dpi = dpi)
      }
      saved_files <- c(saved_files, filepath)
    }, error = function(e) {
      cat("Failed to save", filepath, ":", e$message, "\n")
    })
  }

  if (length(saved_files) > 0) {
    cat("Saved plots:\n")
    for (f in saved_files) {
      cat("  -", f, "\n")
    }
  }

  invisible(saved_files)
}

# ==============================================================================
# PRINT OBJECT STRUCTURE
# ==============================================================================
print_object_structure <- function(obj, name = "Object") {
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("OBJECT STRUCTURE:", name, "\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")

  cat("Class:", class(obj)[1], "\n")
  cat("Dimensions:", nrow(obj), "genes x", ncol(obj), "cells\n")

  # Assays
  assay_names <- names(obj@assays)
  if (length(assay_names) > 0) {
    cat("\nAssays:", paste(assay_names, collapse = ", "), "\n")
    cat("Default assay:", Seurat::DefaultAssay(obj), "\n")
    for (assay in assay_names) {
      assay_class <- class(obj[[assay]])[1]
      layers <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) character(0))
      if (length(layers) > 0) {
        cat("  ", assay, " (", assay_class, ") layers: ", paste(layers, collapse = ", "), "\n", sep = "")
      } else {
        cat("  ", assay, " (", assay_class, "): no layers\n", sep = "")
      }
    }
  } else {
    cat("\nAssays: <none>\n")
  }

  # Reductions
  red_names <- names(obj@reductions)
  if (length(red_names) > 0) {
    cat("\nReductions:", paste(red_names, collapse = ", "), "\n")
    for (red in red_names) {
      dims <- ncol(Seurat::Embeddings(obj, reduction = red))
      cat("  ", red, ":", dims, "dims\n")
    }
  } else {
    cat("\nReductions: <none>\n")
  }

  # Metadata
  meta_cols <- colnames(obj@meta.data)
  cat("\nMetadata columns (", length(meta_cols), "):\n", sep = "")
  cat("  ", paste(head(meta_cols, 15), collapse = ", "))
  if (length(meta_cols) > 15) {
    cat(" ... (+", length(meta_cols) - 15, " more)")
  }
  cat("\n")

  # Sex distribution
  if ("sex" %in% colnames(obj@meta.data)) {
    cat("\nSex distribution:\n")
    print(table(obj$sex))
  }

  # Sample distribution
  if ("sample_name" %in% colnames(obj@meta.data)) {
    cat("\nSample distribution:\n")
    print(table(obj$sample_name))
  }

  # Batch distribution
  if ("batch" %in% colnames(obj@meta.data)) {
    cat("\nBatch distribution:\n")
    print(table(obj$batch))
  }

  # Cluster distribution
  cluster_col <- NULL
  for (col in c("seurat_clusters", "CHOIR_clusters_0.05", "leiden_clusters")) {
    if (col %in% colnames(obj@meta.data)) {
      cluster_col <- col
      break
    }
  }
  if (!is.null(cluster_col)) {
    cat("\nCluster distribution (first 10):\n")
    cluster_table <- sort(table(obj@meta.data[[cluster_col]]), decreasing = TRUE)
    print(head(cluster_table, 10))
    if (length(cluster_table) > 10) {
      cat("  ... and", length(cluster_table) - 10, "more clusters\n")
    }
  }

  cat(paste(rep("-", 60), collapse = ""), "\n")
}

# ==============================================================================
# WRITE README FILES
# ==============================================================================
write_readme <- function(dir_path, title, description, files_info = NULL) {
  readme_path <- file.path(dir_path, "README.txt")

  content <- paste0(
    "=", paste(rep("=", 78), collapse = ""), "\n",
    title, "\n",
    "=", paste(rep("=", 78), collapse = ""), "\n\n",
    "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
    "Description:\n", description, "\n"
  )

  if (!is.null(files_info) && length(files_info) > 0) {
    content <- paste0(content, "\nFiles in this directory:\n")
    for (fname in names(files_info)) {
      content <- paste0(content, "  - ", fname, ": ", files_info[[fname]], "\n")
    }
  }

  writeLines(content, readme_path)
  cat("README written:", readme_path, "\n")
}

# ==============================================================================
# CONVERT SEURAT TO SINGLECELLEXPERIMENT (Seurat v5 compatible)
# ==============================================================================
#' Convert Seurat object to SingleCellExperiment
#'
#' This function handles Seurat v5's layer system by joining layers before
#' extracting counts.
#'
#' @param seurat_obj Seurat object
#' @param assay_use Assay to convert (default: "RNA")
#' @return SingleCellExperiment object
seurat_to_sce <- function(seurat_obj, assay_use = "RNA") {
  # Join layers if needed (Seurat v5 compatibility)
  tryCatch({
    layers <- SeuratObject::Layers(seurat_obj[[assay_use]])
    counts_layers <- grep("^counts", layers, value = TRUE)
    if (length(counts_layers) > 1) {
      cat("    Joining", length(counts_layers), "counts layers before SCE conversion...\n")
      seurat_obj[[assay_use]] <- SeuratObject::JoinLayers(seurat_obj[[assay_use]])
    }
  }, error = function(e) NULL)

  # Extract counts matrix
  counts_mat <- tryCatch({
    Seurat::GetAssayData(seurat_obj, assay = assay_use, layer = "counts")
  }, error = function(e) {
    tryCatch({
      Seurat::GetAssayData(seurat_obj, assay = assay_use, slot = "counts")
    }, error = function(e2) {
      stop("Could not extract counts matrix from Seurat object: ", e2$message)
    })
  })

  col_data <- seurat_obj@meta.data

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_mat),
    colData = col_data
  )

  return(sce)
}

# ==============================================================================
# CONVERT SINGLECELLEXPERIMENT TO SEURAT (Seurat v5 compatible) - NEW
# ==============================================================================
#' Convert SingleCellExperiment back to Seurat object
#'
#' This function creates a Seurat v5 object from an SCE, properly handling
#' the 'data' layer assignment.
#'
#' @param sce SingleCellExperiment object with counts and logcounts
#' @param original_seurat Original Seurat object (for metadata preservation)
#' @return Seurat object with counts and data layers
sce_to_seurat <- function(sce, original_seurat = NULL) {
  # Create new Seurat object with counts
  obj <- Seurat::CreateSeuratObject(
    counts = SingleCellExperiment::counts(sce),
    meta.data = as.data.frame(SummarizedExperiment::colData(sce))
  )

  # Add normalized data as 'data' layer (Seurat v5 style)
  # Use LayerData assignment for v5 compatibility
  if ("logcounts" %in% SummarizedExperiment::assayNames(sce)) {
    tryCatch({
      # Seurat v5 method
      obj[["RNA"]]$data <- SingleCellExperiment::logcounts(sce)
    }, error = function(e) {
      # Fallback for older Seurat versions
      tryCatch({
        obj <- Seurat::SetAssayData(obj, layer = "data",
                                    new.data = SingleCellExperiment::logcounts(sce))
      }, error = function(e2) {
        obj <- Seurat::SetAssayData(obj, slot = "data",
                                    new.data = SingleCellExperiment::logcounts(sce))
      })
    })
  }

  # Preserve additional metadata from original object if provided
  if (!is.null(original_seurat)) {
    original_meta_cols <- colnames(original_seurat@meta.data)
    current_meta_cols <- colnames(obj@meta.data)
    missing_cols <- setdiff(original_meta_cols, current_meta_cols)

    for (col in missing_cols) {
      # Only add if cells match
      if (all(rownames(obj@meta.data) %in% rownames(original_seurat@meta.data))) {
        obj[[col]] <- original_seurat@meta.data[rownames(obj@meta.data), col]
      }
    }
  }

  # Add scran size factor if available
  if (!is.null(SingleCellExperiment::sizeFactors(sce))) {
    obj$scran_size_factor <- SingleCellExperiment::sizeFactors(sce)
  }

  return(obj)
}

# ==============================================================================
# GET MITOCHONDRIAL GENE PATTERN
# ==============================================================================
get_mt_pattern <- function(gene_names) {
  if (any(grepl("^mt-", gene_names, ignore.case = FALSE))) {
    return("^mt-")
  }
  if (any(grepl("^MT-", gene_names, ignore.case = FALSE))) {
    return("^MT-")
  }
  return("^mt-|^MT-")
}

# ==============================================================================
# ENSURE QC COLUMNS EXIST
# ==============================================================================
ensure_qc_columns <- function(obj, sample_name) {
  mt_cols <- c("pct_counts_MT", "percent.mt", "percent_mt")
  mt_col_found <- intersect(mt_cols, colnames(obj@meta.data))

  if (length(mt_col_found) > 0) {
    obj$percent.mt <- obj@meta.data[[mt_col_found[1]]]
    cat(sample_name, "- Using MT column:", mt_col_found[1], "\n")
  } else {
    mt_pattern <- get_mt_pattern(rownames(obj))
    obj$percent.mt <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern)
    cat(sample_name, "- Calculated MT% using pattern:", mt_pattern, "\n")
  }

  if (!"nFeature_RNA" %in% colnames(obj@meta.data)) {
    if ("n_genes_by_counts" %in% colnames(obj@meta.data)) {
      obj$nFeature_RNA <- obj$n_genes_by_counts
    } else {
      obj$nFeature_RNA <- Matrix::colSums(Seurat::GetAssayData(obj, layer = "counts") > 0)
    }
  }

  if (!"nCount_RNA" %in% colnames(obj@meta.data)) {
    if ("total_counts" %in% colnames(obj@meta.data)) {
      obj$nCount_RNA <- obj$total_counts
    } else {
      obj$nCount_RNA <- Matrix::colSums(Seurat::GetAssayData(obj, layer = "counts"))
    }
  }

  return(obj)
}

# ==============================================================================
# DETECT AND VALIDATE COUNTS LAYER
# ==============================================================================
#' Detect available counts layers in a Seurat object
#'
#' @param obj Seurat object
#' @param assay Assay to check (default: "RNA")
#' @param preferred_layer Preferred layer name to use if available
#' @return List with layer name, validation status, and info message
detect_counts_layer <- function(obj, assay = "RNA", preferred_layer = "counts") {
  result <- list(
    layer = NULL,
    valid = FALSE,
    message = "",
    all_layers = character(0)
  )

  # Check if assay exists
  if (!assay %in% names(obj@assays)) {
    result$message <- paste("Assay", assay, "not found in object")
    return(result)
  }

  # Get all layers
  all_layers <- tryCatch({
    SeuratObject::Layers(obj[[assay]])
  }, error = function(e) {
    character(0)
  })

  result$all_layers <- all_layers

  if (length(all_layers) == 0) {
    result$message <- paste("No layers found in assay", assay)
    return(result)
  }

  # Look for counts layers
  counts_layers <- grep("^counts", all_layers, value = TRUE)

  if (length(counts_layers) == 0) {
    result$message <- paste("No counts layers found in assay", assay,
                            ". Available layers:", paste(all_layers, collapse = ", "))
    return(result)
  }

  # Check if preferred layer exists
  if (preferred_layer %in% all_layers) {
    result$layer <- preferred_layer
    result$valid <- TRUE
    result$message <- paste("Using preferred layer:", preferred_layer)
  } else if ("counts" %in% all_layers) {
    result$layer <- "counts"
    result$valid <- TRUE
    result$message <- "Using 'counts' layer"
  } else {
    # Use first counts layer found
    result$layer <- counts_layers[1]
    result$valid <- TRUE
    result$message <- paste("Using first available counts layer:", counts_layers[1])
  }

  return(result)
}

# ==============================================================================
# PRINT COUNTS LAYER INFO
# ==============================================================================
#' Print detailed information about counts layers
#'
#' @param obj Seurat object
#' @param assay Assay to check
print_counts_layer_info <- function(obj, assay = "RNA") {
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("COUNTS LAYER INFORMATION\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")

  if (!assay %in% names(obj@assays)) {
    cat("  ERROR: Assay", assay, "not found\n")
    return(invisible(NULL))
  }

  assay_class <- class(obj[[assay]])[1]
  cat("Assay:", assay, "(", assay_class, ")\n")

  all_layers <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) character(0))
  cat("Available layers:", paste(all_layers, collapse = ", "), "\n")

  # Check each counts layer
  counts_layers <- grep("^counts", all_layers, value = TRUE)
  if (length(counts_layers) > 0) {
    cat("\nCounts layers found:\n")
    for (layer in counts_layers) {
      mat <- tryCatch({
        Seurat::GetAssayData(obj, assay = assay, layer = layer)
      }, error = function(e) NULL)

      if (!is.null(mat)) {
        cat("  ", layer, ": ", nrow(mat), " genes x ", ncol(mat), " cells\n", sep = "")
        cat("    Sum:", format(sum(mat), big.mark = ","), "\n")
        cat("    Range:", min(mat), "-", max(mat), "\n")
      } else {
        cat("  ", layer, ": could not access\n", sep = "")
      }
    }
  } else {
    cat("  No counts layers found\n")
  }

  cat(paste(rep("-", 60), collapse = ""), "\n")
}

# ==============================================================================
# GET NORMALIZATION METHOD FOR INTEGRATION
# ==============================================================================
#' Determine which normalization method to use for integration
#'
#' @param params Pipeline parameters
#' @param objects_dir Directory containing Module 03 outputs
#' @return Character string with normalization method name
get_integration_normalization_method <- function(params, objects_dir) {
  method <- params$integration_normalization_method

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("DETERMINING NORMALIZATION METHOD FOR INTEGRATION\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Configured method:", method, "\n")

  if (method == "auto") {
    # Try to load the selected method from Module 03
    norm_data_file <- file.path(objects_dir, "03_normalization_data.RData")

    if (file.exists(norm_data_file)) {
      # Load into temporary environment
      temp_env <- new.env()
      load(norm_data_file, envir = temp_env)

      if ("selected_normalization_method" %in% ls(temp_env)) {
        method <- temp_env$selected_normalization_method
        cat("Auto-selected from Module 03 benchmarking:", method, "\n")
      } else {
        cat("WARNING: 'selected_normalization_method' not found in Module 03 output\n")
        method <- "LogNormalize"
        cat("Defaulting to:", method, "\n")
      }
    } else {
      cat("WARNING: Module 03 data not found at:", norm_data_file, "\n")
      method <- "LogNormalize"
      cat("Defaulting to:", method, "\n")
    }
  }

  # Validate method - now includes scKWARN
  valid_methods <- c("LogNormalize", "SCTransform", "scran", "scKWARN")
  if (!method %in% valid_methods) {
    cat("WARNING: Invalid method '", method, "'. Using LogNormalize.\n", sep = "")
    method <- "LogNormalize"
  }

  cat("\n>>> FINAL NORMALIZATION METHOD FOR INTEGRATION:", method, "<<<\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  return(method)
}

# ==============================================================================
# LOAD NORMALIZED OBJECT FOR INTEGRATION
# ==============================================================================
#' Load the appropriate normalized object for integration
#'
#' @param method Normalization method ("LogNormalize", "SCTransform", "scran", "scKWARN")
#' @param objects_dir Directory containing normalized objects
#' @param use_unintegrated Whether to use unintegrated version (TRUE) or integrated (FALSE)
#' @return Seurat object or NULL if not found
load_normalized_object_for_integration <- function(method, objects_dir, use_unintegrated = TRUE) {
  merged_norm_dir <- file.path(objects_dir, "merged_normalized")

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("LOADING NORMALIZED OBJECT FOR INTEGRATION\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Method:", method, "\n")
  cat("Use unintegrated:", use_unintegrated, "\n")
  cat("Looking in:", merged_norm_dir, "\n")

  # Define file paths based on method - now includes scKWARN
  if (use_unintegrated) {
    file_paths <- list(
      LogNormalize = file.path(merged_norm_dir, "merged_LogNormalize_unintegrated.rds"),
      SCTransform = file.path(merged_norm_dir, "merged_SCTransform_unintegrated.rds"),
      scran = file.path(merged_norm_dir, "merged_scran_unintegrated.rds"),
      scKWARN = file.path(merged_norm_dir, "merged_scKWARN_unintegrated.rds")
    )
  } else {
    file_paths <- list(
      LogNormalize = file.path(objects_dir, "integrated_LogNormalize_Harmony.rds"),
      SCTransform = file.path(objects_dir, "integrated_SCTransform_Harmony.rds"),
      scran = file.path(objects_dir, "integrated_scran_Harmony.rds"),
      scKWARN = file.path(objects_dir, "integrated_scKWARN_Harmony.rds")
    )
  }

  target_path <- file_paths[[method]]

  if (is.null(target_path)) {
    cat("ERROR: Unknown method:", method, "\n")
    return(NULL)
  }

  cat("Target file:", target_path, "\n")

  if (!file.exists(target_path)) {
    cat("ERROR: File not found:", target_path, "\n")

    # Try to find any available file
    cat("\nSearching for alternative files...\n")
    if (dir.exists(merged_norm_dir)) {
      available_files <- list.files(merged_norm_dir, pattern = "\\.rds$", full.names = TRUE)
      if (length(available_files) > 0) {
        cat("Available files in merged_normalized:\n")
        for (f in available_files) {
          cat("  -", basename(f), "\n")
        }
      } else {
        cat("  No RDS files found in merged_normalized directory\n")
      }
    } else {
      cat("  Directory does not exist:", merged_norm_dir, "\n")
    }

    return(NULL)
  }

  cat("\nLoading:", basename(target_path), "...\n")
  obj <- readRDS(target_path)

  cat("\n>>> LOADED OBJECT SUMMARY <<<\n")
  cat("Cells:", ncol(obj), "\n")
  cat("Genes:", nrow(obj), "\n")
  cat("Default assay:", Seurat::DefaultAssay(obj), "\n")

  # Print assay info
  for (assay_name in names(obj@assays)) {
    assay_class <- class(obj[[assay_name]])[1]
    layers <- tryCatch(SeuratObject::Layers(obj[[assay_name]]), error = function(e) character(0))
    cat("  ", assay_name, " (", assay_class, "): ", paste(layers, collapse = ", "), "\n", sep = "")
  }

  # Print reductions
  if (length(names(obj@reductions)) > 0) {
    cat("Reductions:", paste(names(obj@reductions), collapse = ", "), "\n")
  }

  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  return(obj)
}

# ==============================================================================
# INTEGRATION ERROR HELPER
# ==============================================================================
#' Print detailed error information for failed integration
#'
#' @param method_name Name of the integration method
#' @param error_obj Error object from tryCatch
#' @param obj Seurat object that was being processed
print_integration_error <- function(method_name, error_obj, obj = NULL) {
  cat("\n", paste(rep("!", 60), collapse = ""), "\n")
  cat("INTEGRATION FAILED:", method_name, "\n")
  cat(paste(rep("!", 60), collapse = ""), "\n")

  cat("\nError message:", conditionMessage(error_obj), "\n")

  # Print object info if available
  if (!is.null(obj)) {
    cat("\nObject details at time of failure:\n")
    cat("  Cells:", ncol(obj), "\n")
    cat("  Genes:", nrow(obj), "\n")
    cat("  Default assay:", Seurat::DefaultAssay(obj), "\n")

    # Check assay class
    default_assay <- Seurat::DefaultAssay(obj)
    assay_class <- class(obj[[default_assay]])[1]
    cat("  Assay class:", assay_class, "\n")

    # List layers
    layers <- tryCatch(SeuratObject::Layers(obj[[default_assay]]), error = function(e) character(0))
    cat("  Layers:", paste(layers, collapse = ", "), "\n")

    # Check batch variable
    if ("batch" %in% colnames(obj@meta.data)) {
      cat("  Batch levels:", paste(unique(obj$batch), collapse = ", "), "\n")
      cat("  Cells per batch:\n")
      print(table(obj$batch))
    }
  }

  # Provide suggestions based on error
  cat("\nPossible causes and solutions:\n")

  error_msg <- tolower(conditionMessage(error_obj))

  if (grepl("sctassay", error_msg) || grepl("joinlayers", error_msg) || grepl("no applicable method", error_msg)) {
    cat("  [!] SCTAssay doesn't support standard layer operations\n")
    cat("  [>] Solution: Use LogNormalize method for CCA/RPCA/FastMNN\n")
    cat("  [>] Set integration_normalization_method = 'LogNormalize' in params.R\n")
  }

  if (grepl("subscript out of bounds", error_msg)) {
    cat("  [!] Layer structure may be incompatible or data is missing\n")
    cat("  [>] Solution: Ensure layers are properly split by batch\n")
    cat("  [>] Check that the RNA assay has counts and data layers for each batch\n")
  }

  if (grepl("future.globals.maxsize", error_msg) || grepl("exceeds the maximum", error_msg)) {
    cat("  [!] Object is too large for parallel processing\n")
    cat("  [>] Solution: Increase future.globals.maxSize before integration:\n")
    cat("  [>]   options(future.globals.maxSize = 30 * 1024^3)  # 30 GB\n")
  }

  if (grepl("batch", error_msg) && grepl("must be specified", error_msg)) {
    cat("  [!] Batch information not properly set for FastMNN\n")
    cat("  [>] Solution: Ensure 'batch' column exists and layers are split by batch\n")
  }

  if (grepl("counts", error_msg) && grepl("null", error_msg)) {
    cat("  [!] Counts layer not found or not accessible\n")
    cat("  [>] Solution: Check that counts layer exists with detect_counts_layer()\n")
  }

  if (grepl("no variable features", error_msg) || grepl("variablefeatures", error_msg)) {
    cat("  [!] Variable features not found\n")
    cat("  [>] Solution: Run FindVariableFeatures() before integration\n")
  }

  cat(paste(rep("!", 60), collapse = ""), "\n\n")
}

# ==============================================================================
# PREPARE OBJECT FOR INTEGRATION
# ==============================================================================
#' Prepare a Seurat object for integration methods
#'
#' This function ensures the object is in the correct state for integration:
#' - Uses RNA assay (not SCT) for compatibility with all methods
#' - Joins layers first, then splits by batch
#' - Preserves existing normalization from Module 03 (does NOT re-normalize)
#' - Runs FindVariableFeatures, ScaleData, RunPCA
#'
#' CRITICAL (2026-02-09): The previous version called NormalizeData() here,
#' which overwrote scran/scKWARN/SCT normalization with LogNormalize.
#' This meant Module 03's normalization selection was silently ignored in
#' Module 04. The fix detects whether RNA$data already contains normalized
#' expression and skips NormalizeData() if so.
#'
#' @param obj Seurat object (with normalization already applied by Module 03)
#' @param batch_var Batch variable name in metadata
#' @param nfeatures Number of variable features
#' @param dims_use Number of PCA dimensions
#' @return Prepared Seurat object with layers split by batch and PCA computed
prepare_object_for_integration <- function(obj, batch_var = "batch",
                                            nfeatures = 3000, dims_use = 30) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("PREPARING OBJECT FOR INTEGRATION\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  # Always use RNA assay for integration (SCT has compatibility issues)
  if (!"RNA" %in% names(obj@assays)) {
    stop("RNA assay not found in object. Available assays: ",
         paste(names(obj@assays), collapse = ", "))
  }

  Seurat::DefaultAssay(obj) <- "RNA"
  cat("Default assay set to: RNA\n")

  # Print RNA assay info
  rna_class <- class(obj[["RNA"]])[1]
  cat("RNA assay class:", rna_class, "\n")

  # Ensure batch column exists
  if (!batch_var %in% colnames(obj@meta.data)) {
    stop("Batch variable '", batch_var, "' not found in metadata. ",
         "Available columns: ", paste(head(colnames(obj@meta.data), 10), collapse = ", "))
  }

  batch_levels <- unique(obj@meta.data[[batch_var]])
  cat("Batch variable:", batch_var, "\n")
  cat("Batch levels:", paste(batch_levels, collapse = ", "), "\n")
  cat("Cells per batch:\n")
  print(table(obj@meta.data[[batch_var]]))

  # Step 1: Join layers (in case they're already split)
  cat("\nStep 1: Joining RNA layers...\n")
  tryCatch({
    obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
    cat("  Layers joined successfully\n")
  }, error = function(e) {
    cat("  Note: JoinLayers returned:", e$message, "\n")
    cat("  Continuing...\n")
  })

  layers_after_join <- tryCatch(SeuratObject::Layers(obj[["RNA"]]), error = function(e) character(0))
  cat("  Layers after join:", paste(layers_after_join, collapse = ", "), "\n")

  # Step 2: Split by batch
  cat("\nStep 2: Splitting layers by batch...\n")
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[batch_var]])

  layers_after_split <- tryCatch(SeuratObject::Layers(obj[["RNA"]]), error = function(e) character(0))
  cat("  Layers after split:", paste(layers_after_split, collapse = ", "), "\n")

  # Step 3: Preprocessing -- PRESERVE EXISTING NORMALIZATION
  cat("\nStep 3: Running preprocessing...\n")

  # -------------------------------------------------------------------------
  # NORMALIZATION CHECK (FIX 2026-02-09)
  # -------------------------------------------------------------------------
  # Check if the object already has normalized data in the RNA data layer.
  # If a normalization method (scran, scKWARN, LogNormalize) was applied in
  # Module 03, the RNA$data layer already contains the correct normalized
  # expression. Re-running NormalizeData() would overwrite it with
  # LogNormalize, invalidating the normalization selection from Module 03.
  #
  # We only call NormalizeData() if the data layer appears to be missing.
  # -------------------------------------------------------------------------
  has_data_layer <- any(grepl("^data", layers_after_split))

  # Additional check: after splitting, data layers get names like "data.batch1"
  # If we only see "counts.batch1" layers but no "data.batch1", normalization
  # may not have been applied. Also check pre-split state.
  if (!has_data_layer) {
    # Try checking with joined layers (more reliable)
    has_data_layer <- any(grepl("^data", layers_after_join))
  }

  if (has_data_layer) {
    cat("  NormalizeData... SKIPPED\n")
    cat("    RNA$data layer already contains normalized expression from Module 03\n")
    cat("    Preserving selected normalization (no overwrite with LogNormalize)\n")

    # Report what normalization is being preserved
    # Try to detect from the loaded RData or metadata
    norm_source <- tryCatch({
      if (exists("normalization_method", envir = parent.frame())) {
        get("normalization_method", envir = parent.frame())
      } else if (exists("selected_normalization_method", envir = parent.frame())) {
        get("selected_normalization_method", envir = parent.frame())
      } else {
        "unknown (check Module 03 output)"
      }
    }, error = function(e) "unknown")
    cat("    Normalization source:", norm_source, "\n")
  } else {
    cat("  NormalizeData... APPLYING (no existing normalized data detected)\n")
    cat("    This typically means the object was loaded without prior normalization.\n")
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  }

  cat("  FindVariableFeatures (n =", nfeatures, ")...\n")
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE)

  cat("  ScaleData...\n")
  obj <- Seurat::ScaleData(obj, verbose = FALSE)

  cat("  RunPCA (npcs = 50)...\n")
  obj <- Seurat::RunPCA(obj, npcs = 50, verbose = FALSE)

  # -------------------------------------------------------------------------
  # DATA PROVENANCE SUMMARY
  # -------------------------------------------------------------------------
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("DATA PROVENANCE FOR INTEGRATION METHODS\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  cat("  RNA$data layer: ",
      ifelse(has_data_layer, "PRESERVED from Module 03", "freshly computed (LogNormalize)"), "\n")
  cat("  RNA$counts layer: raw counts (unchanged)\n")
  cat("\n")
  cat("  Methods using RNA$data (selected normalization):\n")
  cat("    Harmony  -- via PCA computed from scaled data layer\n")
  cat("    CCA      -- via split data layers directly\n")
  cat("    RPCA     -- via split data layers directly\n")
  cat("    FastMNN  -- via split data layers directly\n")
  cat("    Scanorama -- exports data layer to Python\n")
  cat("    BBKNN    -- exports data layer to Python\n")
  cat("\n")
  cat("  Methods using RNA$counts (raw counts, normalize internally):\n")
  cat("    scVI     -- negative binomial VAE, requires integer counts\n")
  cat("    scCobra  -- deep learning, requires raw counts\n")
  cat("    CONCORD  -- contrastive learning, applies own normalize_total + log1p\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")

  cat("\nObject prepared for integration\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  return(obj)
}

# ==============================================================================
# COMPUTE NORMALIZATION METRICS
# ==============================================================================
compute_normalization_metrics <- function(obj, reduction = "pca", dims_use = 30,
                                           method_name = "method", batch_var = "batch") {
  if (!reduction %in% names(obj@reductions)) {
    cat("  Warning: Reduction", reduction, "not found. Skipping metrics.\n")
    return(data.frame(
      method = method_name,
      batch_variance = NA,
      batch_asw = NA,
      lisi = NA
    ))
  }

  emb <- Seurat::Embeddings(obj, reduction = reduction)
  nd <- min(dims_use, ncol(emb))
  emb <- emb[, 1:nd, drop = FALSE]

  batch_labels <- obj@meta.data[[batch_var]]

  # Batch variance
  batch_var_score <- tryCatch({
    batches <- unique(batch_labels)
    batch_means <- lapply(batches, function(b) {
      colMeans(emb[batch_labels == b, , drop = FALSE])
    })
    overall_mean <- colMeans(emb)
    mean(sapply(batch_means, function(m) sum((m - overall_mean)^2)))
  }, error = function(e) NA)

  # Batch ASW
  batch_asw <- tryCatch({
    if (nrow(emb) > 5000) {
      set.seed(42)
      idx <- sample(nrow(emb), 5000)
      emb_sub <- emb[idx, ]
      batch_sub <- batch_labels[idx]
    } else {
      emb_sub <- emb
      batch_sub <- batch_labels
    }
    d <- dist(emb_sub)
    sil <- cluster::silhouette(as.numeric(factor(batch_sub)), d)
    mean(sil[, 3])
  }, error = function(e) NA)

  # LISI
  lisi_score <- tryCatch({
    if (requireNamespace("lisi", quietly = TRUE)) {
      if (nrow(emb) > 5000) {
        set.seed(42)
        idx <- sample(nrow(emb), 5000)
        emb_sub <- emb[idx, ]
        batch_sub <- batch_labels[idx]
      } else {
        emb_sub <- emb
        batch_sub <- batch_labels
      }
      meta_df <- data.frame(batch = batch_sub)
      lisi_result <- lisi::compute_lisi(emb_sub, meta_df, "batch")
      mean(lisi_result$batch)
    } else {
      NA
    }
  }, error = function(e) NA)

  return(data.frame(
    method = method_name,
    batch_variance = batch_var_score,
    batch_asw = batch_asw,
    lisi = lisi_score
  ))
}

# ==============================================================================
# COMPUTE BATCH VARIANCE
# ==============================================================================
compute_batch_variance <- function(embeddings, batch_labels) {
  batches <- unique(batch_labels)
  batch_means <- lapply(batches, function(b) {
    colMeans(embeddings[batch_labels == b, , drop = FALSE])
  })
  overall_mean <- colMeans(embeddings)
  batch_var <- mean(sapply(batch_means, function(m) sum((m - overall_mean)^2)))
  return(batch_var)
}

# ==============================================================================
# COMPUTE BATCH AVERAGE SILHOUETTE WIDTH
# ==============================================================================
compute_batch_asw <- function(embeddings, batch_labels) {
  tryCatch({
    if (nrow(embeddings) > 10000) {
      set.seed(42)
      idx <- sample(nrow(embeddings), 10000)
      embeddings <- embeddings[idx, ]
      batch_labels <- batch_labels[idx]
    }
    d <- dist(embeddings)
    sil <- cluster::silhouette(as.numeric(factor(batch_labels)), d)
    return(mean(sil[, 3]))
  }, error = function(e) {
    return(NA)
  })
}

# ==============================================================================
# COMPUTE LISI SCORE
# ==============================================================================
compute_lisi_score <- function(embeddings, batch_labels) {
  tryCatch({
    if (!requireNamespace("lisi", quietly = TRUE)) {
      return(NA)
    }
    if (nrow(embeddings) > 10000) {
      set.seed(42)
      idx <- sample(nrow(embeddings), 10000)
      embeddings <- embeddings[idx, ]
      batch_labels <- batch_labels[idx]
    }
    meta_df <- data.frame(batch = batch_labels)
    lisi_result <- lisi::compute_lisi(embeddings, meta_df, "batch")
    return(mean(lisi_result$batch))
  }, error = function(e) {
    return(NA)
  })
}

# ==============================================================================
# CREATE PSEUDOBULK WITH PSEUDO-REPLICATES
# ==============================================================================
create_pseudobulk_with_replicates <- function(obj, counts, group_var = "sex",
                                               n_splits = 3, seed = 42) {
  set.seed(seed)
  meta <- obj@meta.data

  group_labels <- meta[[group_var]]
  pseudo_sample <- character(nrow(meta))

  for (grp in unique(group_labels)) {
    cells_grp <- which(group_labels == grp)
    splits <- sample(1:n_splits, length(cells_grp), replace = TRUE)
    pseudo_sample[cells_grp] <- paste0(grp, "_rep", splits)
  }

  pb_list <- list()
  for (ps in unique(pseudo_sample)) {
    cells <- rownames(meta)[pseudo_sample == ps]
    if (length(cells) > 0) {
      pb_list[[ps]] <- Matrix::rowSums(counts[, cells, drop = FALSE])
    }
  }

  pb_counts <- do.call(cbind, pb_list)

  pb_meta <- data.frame(
    sample = colnames(pb_counts),
    group = factor(sub("_rep[0-9]+$", "", colnames(pb_counts))),
    replicate = as.integer(sub(".*_rep", "", colnames(pb_counts))),
    row.names = colnames(pb_counts)
  )

  return(list(counts = pb_counts, meta = pb_meta))
}

# ==============================================================================
# CASE-INSENSITIVE GENE MATCHING
# ==============================================================================
find_genes_in_object <- function(gene_list, obj) {
  available_genes <- c()
  gene_names <- rownames(obj)

  for (gene in gene_list) {
    # Strategy 1: Exact match
    if (gene %in% gene_names) {
      available_genes <- c(available_genes, gene)
      next
    }

    # Strategy 2: Case variations
    variations <- c(
      tolower(gene),
      toupper(gene),
      paste0(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene))))
    )

    found <- FALSE
    for (var in variations) {
      if (var %in% gene_names) {
        available_genes <- c(available_genes, var)
        found <- TRUE
        break
      }
    }

    # Strategy 3: Grep match
    if (!found) {
      matches <- grep(paste0("^", gene, "$"), gene_names, ignore.case = TRUE, value = TRUE)
      if (length(matches) > 0) {
        available_genes <- c(available_genes, matches[1])
      }
    }
  }

  return(unique(available_genes))
}

# ==============================================================================
# SAFE BOOLEAN CHECK
# ==============================================================================
safe_true <- function(x) {
  isTRUE(x)
}

# ==============================================================================
# PRINT MODULE HEADER
# ==============================================================================
print_module_header <- function(module_num, module_name) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("MODULE ", sprintf("%02d", module_num), ": ", module_name, "\n", sep = "")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
}

# ==============================================================================
# PRINT MODULE COMPLETE
# ==============================================================================
print_module_complete <- function(module_num) {
  cat("\n>>> MODULE ", sprintf("%02d", module_num), " COMPLETE <<<\n", sep = "")
}

# ==============================================================================
# METRIC NORMALIZATION FOR BENCHMARKING
# ==============================================================================
#' Min-max normalization for benchmark metrics
#'
#' Normalizes values to 0-1 scale. After transformation, higher values always
#' indicate "better" performance.
#'
#' @param x Numeric vector to normalize
#' @param lower_is_better If TRUE, inverts the scale (for metrics like batch_variance)
#' @return Numeric vector normalized to 0-1 scale (higher = better)
#'
#' @examples
#' # batch_variance: lower is better, so invert
#' normalize_metric(c(0.01, 0.5, 1.0), lower_is_better = TRUE)
#' # Returns: c(1.0, 0.5, 0.0)
#'
#' # LISI: higher is better, don't invert
#' normalize_metric(c(1.0, 1.5, 2.0), lower_is_better = FALSE)
#' # Returns: c(0.0, 0.5, 1.0)
normalize_metric <- function(x, lower_is_better = TRUE) {
  # Handle edge cases

if (all(is.na(x)) || length(unique(na.omit(x))) <= 1) {
    return(rep(0.5, length(x)))
  }

  x_range <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (x_range == 0) {
    return(rep(0.5, length(x)))
  }

  # Min-max normalization to 0-1
  normalized <- (x - min(x, na.rm = TRUE)) / x_range

  # Invert if lower is better (so higher = better after transformation)
  if (lower_is_better) {
    normalized <- 1 - normalized
  }

  return(normalized)
}

# ==============================================================================
# SUBSAMPLE SEURAT OBJECT FOR BENCHMARKING (STRATIFIED BY BATCH)
# ==============================================================================
#' Subsample a Seurat object while preserving batch proportions
#'
#' Used by Module 03 normalization benchmarking to speed up integration
#' method comparisons. Performs stratified sampling so each batch retains
#' its proportional representation. If the object already has fewer cells
#' than max_cells, it is returned unchanged.
#'
#' @param obj Seurat object
#' @param max_cells Maximum number of cells to retain (NULL or 0 = no subsampling)
#' @param batch_var Metadata column used for stratified sampling
#' @param seed Random seed for reproducibility
#' @return Subsampled Seurat object (or original if already small enough)
subsample_for_benchmark <- function(obj, max_cells = 5000, batch_var = "batch", seed = 42) {
  if (is.null(max_cells) || max_cells <= 0 || ncol(obj) <= max_cells) {
    cat("  Subsampling: not needed (", ncol(obj), " cells <= ",
        if (is.null(max_cells) || max_cells <= 0) "no limit" else max_cells, ")\n", sep = "")
    return(obj)
  }

  set.seed(seed)

  n_total <- ncol(obj)
  cat("  Subsampling:", n_total, "->", max_cells, "cells (stratified by", batch_var, ")\n")

  if (!batch_var %in% colnames(obj@meta.data)) {
    cat("    WARNING:", batch_var, "not found in metadata. Sampling uniformly.\n")
    keep_cells <- sample(colnames(obj), max_cells)
    sub_obj <- obj[, keep_cells]
    cat("    Subsampled:", ncol(sub_obj), "cells\n")
    return(sub_obj)
  }

  batch_labels <- obj@meta.data[[batch_var]]
  batches <- unique(batch_labels)

  batch_sizes <- table(batch_labels)
  batch_fracs <- batch_sizes / sum(batch_sizes)
  batch_alloc <- round(batch_fracs * max_cells)

  batch_alloc[batch_alloc < 1] <- 1
  while (sum(batch_alloc) > max_cells) {
    largest <- names(which.max(batch_alloc))
    batch_alloc[largest] <- batch_alloc[largest] - 1
  }

  keep_cells <- c()
  for (b in names(batch_alloc)) {
    cells_in_batch <- colnames(obj)[batch_labels == b]
    n_take <- min(batch_alloc[b], length(cells_in_batch))
    keep_cells <- c(keep_cells, sample(cells_in_batch, n_take))
  }

  sub_obj <- obj[, keep_cells]

  cat("    Per-batch allocation:\n")
  for (b in names(batch_alloc)) {
    original_n <- batch_sizes[b]
    sampled_n <- sum(obj@meta.data[keep_cells, batch_var] == b)
    cat("      ", b, ":", original_n, "->", sampled_n, "\n")
  }
  cat("    Total subsampled:", ncol(sub_obj), "cells\n")

  return(sub_obj)
}

# ==============================================================================
# RUN SINGLE INTEGRATION METHOD FOR NORMALIZATION BENCHMARKING
# ==============================================================================
#' Run a single integration method on a Seurat object and return embeddings
#'
#' Used by Module 03 to benchmark normalization methods across multiple
#' integration approaches. Returns the embedding matrix or NULL on failure.
#'
#' Supported methods:
#'   R-based (use existing normalization -- valid for normalization voting):
#'     - "harmony": Runs on existing PCA embeddings directly
#'     - "mnn": FastMNN on split RNA layers (skips SCT)
#'     - "rpca": Reciprocal PCA integration (skips SCT)
#'     - "cca": Canonical Correlation Analysis integration (skips SCT)
#'
#'   Python-based (use raw counts -- NOT valid for normalization voting):
#'     - "scvi", "sccobra", "concord"
#'
#' CRITICAL DESIGN PRINCIPLE:
#'   R-based methods NEVER call NormalizeData(). The input object already
#'   has normalization-specific data in RNA$data (for LogNorm/scran/scKWARN)
#'   or in the SCT assay (for SCTransform). Calling NormalizeData() would
#'   overwrite all normalizations with LogNormalize, invalidating the
#'   benchmark comparison.
#'
#' @param obj Seurat object with PCA computed and batch metadata
#' @param method Integration method name (see above)
#' @param batch_var Name of the metadata column containing batch labels
#' @param dims_use Number of dimensions for latent space
#' @param params Pipeline params list (for Concord hyperparameters)
#' @return Matrix of embeddings (cells x dims) or NULL on failure
run_integration_for_benchmark <- function(obj, method, batch_var, dims_use, params = NULL) {
  cat("    Running", method, "... ")

  tryCatch({

    # ==================================================================
    # HARMONY -- operates on PCA embeddings directly
    # ==================================================================
    if (method == "harmony") {
      if (!requireNamespace("harmony", quietly = TRUE)) {
        cat("[harmony not installed]\n")
        return(NULL)
      }

      temp <- obj

      if (!"pca" %in% names(temp@reductions)) {
        cat("[FAILED: PCA not found]\n")
        return(NULL)
      }

      nd <- min(dims_use, ncol(Seurat::Embeddings(temp, "pca")))
      if (nd < 2) {
        cat("[FAILED: too few PCA dims]\n")
        return(NULL)
      }

      # Remove stale reduction if present
      if ("harmony" %in% names(temp@reductions)) temp[["harmony"]] <- NULL

      # Run Harmony directly on existing PCA -- no NormalizeData, no layer ops
      temp <- harmony::RunHarmony(
        object        = temp,
        group.by.vars = batch_var,
        reduction.use = "pca",
        dims.use      = 1:nd,
        verbose       = FALSE
      )

      cat("[OK]\n")
      return(Seurat::Embeddings(temp, "harmony"))

    # ==================================================================
    # MNN (FastMNN) -- needs split RNA layers, preserves normalization
    # ==================================================================
    } else if (method == "mnn") {
      if (!requireNamespace("batchelor", quietly = TRUE) ||
          !requireNamespace("SeuratWrappers", quietly = TRUE)) {
        cat("[batchelor/SeuratWrappers not installed]\n")
        return(NULL)
      }

      library(SeuratWrappers)
      temp <- obj

      # SCT skip: RNA$data doesn't hold SCT normalization
      if (Seurat::DefaultAssay(temp) == "SCT") {
        cat("[SKIP: MNN cannot fairly evaluate SCT -- RNA$data is not SCT-normalized]\n")
        return(NULL)
      }

      if (!"RNA" %in% names(temp@assays)) {
        cat("[FAILED: RNA assay missing]\n")
        return(NULL)
      }

      Seurat::DefaultAssay(temp) <- "RNA"

      # Validate data layer exists
      rna_layers <- tryCatch(SeuratObject::Layers(temp[["RNA"]]), error = function(e) character(0))
      has_data_layer <- any(grepl("^data", rna_layers)) ||
        !is.null(tryCatch(temp[["RNA"]]$data, error = function(e) NULL))
      if (!has_data_layer) {
        cat("[SKIP: RNA data layer missing]\n")
        return(NULL)
      }

      # Join then split by batch
      tryCatch({
        temp[["RNA"]] <- SeuratObject::JoinLayers(temp[["RNA"]])
      }, error = function(e) NULL)
      temp[["RNA"]] <- split(temp[["RNA"]], f = temp@meta.data[[batch_var]])

      # NO NormalizeData() -- keep the benchmarked normalization
      temp <- Seurat::FindVariableFeatures(temp, nfeatures = 3000, verbose = FALSE)
      temp <- Seurat::ScaleData(temp, verbose = FALSE)
      temp <- Seurat::RunPCA(temp, npcs = 50, verbose = FALSE)

      if ("integrated.mnn" %in% names(temp@reductions)) temp[["integrated.mnn"]] <- NULL

      temp <- Seurat::IntegrateLayers(temp, method = SeuratWrappers::FastMNNIntegration,
                                       new.reduction = "integrated.mnn",
                                       verbose = FALSE)
      cat("[OK]\n")
      return(Seurat::Embeddings(temp, "integrated.mnn"))

    # ==================================================================
    # RPCA -- Reciprocal PCA integration, preserves normalization
    # ==================================================================
    } else if (method == "rpca") {
      temp <- obj

      # SCT skip: same reason as MNN
      if (Seurat::DefaultAssay(temp) == "SCT") {
        cat("[SKIP: RPCA cannot fairly evaluate SCT -- RNA$data is not SCT-normalized]\n")
        return(NULL)
      }

      if (!"RNA" %in% names(temp@assays)) {
        cat("[FAILED: RNA assay missing]\n")
        return(NULL)
      }

      Seurat::DefaultAssay(temp) <- "RNA"

      # Validate data layer exists
      rna_layers <- tryCatch(SeuratObject::Layers(temp[["RNA"]]), error = function(e) character(0))
      has_data_layer <- any(grepl("^data", rna_layers)) ||
        !is.null(tryCatch(temp[["RNA"]]$data, error = function(e) NULL))
      if (!has_data_layer) {
        cat("[SKIP: RNA data layer missing]\n")
        return(NULL)
      }

      # Join then split by batch
      tryCatch({
        temp[["RNA"]] <- SeuratObject::JoinLayers(temp[["RNA"]])
      }, error = function(e) NULL)
      temp[["RNA"]] <- split(temp[["RNA"]], f = temp@meta.data[[batch_var]])

      # NO NormalizeData() -- keep the benchmarked normalization
      temp <- Seurat::FindVariableFeatures(temp, nfeatures = 3000, verbose = FALSE)
      temp <- Seurat::ScaleData(temp, verbose = FALSE)
      temp <- Seurat::RunPCA(temp, npcs = 50, verbose = FALSE)

      if ("integrated.rpca" %in% names(temp@reductions)) temp[["integrated.rpca"]] <- NULL

      temp <- Seurat::IntegrateLayers(temp, method = Seurat::RPCAIntegration,
                                       new.reduction = "integrated.rpca",
                                       verbose = FALSE)
      cat("[OK]\n")
      return(Seurat::Embeddings(temp, "integrated.rpca"))

    # ==================================================================
    # CCA -- Canonical Correlation Analysis integration, preserves normalization
    # ==================================================================
    } else if (method == "cca") {
      temp <- obj

      # SCT skip: same reason as MNN
      if (Seurat::DefaultAssay(temp) == "SCT") {
        cat("[SKIP: CCA cannot fairly evaluate SCT -- RNA$data is not SCT-normalized]\n")
        return(NULL)
      }

      if (!"RNA" %in% names(temp@assays)) {
        cat("[FAILED: RNA assay missing]\n")
        return(NULL)
      }

      Seurat::DefaultAssay(temp) <- "RNA"

      # Validate data layer exists
      rna_layers <- tryCatch(SeuratObject::Layers(temp[["RNA"]]), error = function(e) character(0))
      has_data_layer <- any(grepl("^data", rna_layers)) ||
        !is.null(tryCatch(temp[["RNA"]]$data, error = function(e) NULL))
      if (!has_data_layer) {
        cat("[SKIP: RNA data layer missing]\n")
        return(NULL)
      }

      # Join then split by batch
      tryCatch({
        temp[["RNA"]] <- SeuratObject::JoinLayers(temp[["RNA"]])
      }, error = function(e) NULL)
      temp[["RNA"]] <- split(temp[["RNA"]], f = temp@meta.data[[batch_var]])

      # NO NormalizeData() -- keep the benchmarked normalization
      temp <- Seurat::FindVariableFeatures(temp, nfeatures = 3000, verbose = FALSE)
      temp <- Seurat::ScaleData(temp, verbose = FALSE)
      temp <- Seurat::RunPCA(temp, npcs = 50, verbose = FALSE)

      if ("integrated.cca" %in% names(temp@reductions)) temp[["integrated.cca"]] <- NULL

      temp <- Seurat::IntegrateLayers(temp, method = Seurat::CCAIntegration,
                                       new.reduction = "integrated.cca",
                                       verbose = FALSE)
      cat("[OK]\n")
      return(Seurat::Embeddings(temp, "integrated.cca"))

    # ==================================================================
    # scVI -- Python, uses raw counts (not valid for normalization vote)
    # ==================================================================
    } else if (method == "scvi") {
      library(reticulate)
      if (!py_module_available("scvi")) {
        cat("[scvi not available]\n")
        return(NULL)
      }

      temp <- obj
      Seurat::DefaultAssay(temp) <- "RNA"
      tryCatch({
        temp[["RNA"]] <- SeuratObject::JoinLayers(temp[["RNA"]])
      }, error = function(e) NULL)

      var_features <- Seurat::VariableFeatures(temp)
      if (length(var_features) == 0) {
        temp <- Seurat::FindVariableFeatures(temp, nfeatures = 3000, verbose = FALSE)
        var_features <- Seurat::VariableFeatures(temp)
      }

      if (length(var_features) == 0) {
        cat("[no variable features]\n")
        return(NULL)
      }

      counts_data_t <- t(as.matrix(Seurat::GetAssayData(temp, assay = "RNA", layer = "counts")[var_features, ]))
      obs_df <- data.frame(batch = as.character(temp@meta.data[[batch_var]]),
                           row.names = colnames(temp), stringsAsFactors = FALSE)
      var_df <- data.frame(gene = var_features, row.names = var_features, stringsAsFactors = FALSE)

      n_latent <- as.integer(dims_use)

      py_run_string("import numpy as np; import pandas as pd; import anndata; import scvi")
      py$bm_counts <- r_to_py(counts_data_t)
      py$bm_obs <- r_to_py(obs_df)
      py$bm_var <- r_to_py(var_df)
      py$bm_cells <- r_to_py(colnames(temp))
      py$bm_genes <- r_to_py(var_features)
      py$bm_n_latent <- r_to_py(n_latent)

      py_run_string("
bm_adata = anndata.AnnData(X=np.array(bm_counts, dtype='float32'),
                            obs=pd.DataFrame(bm_obs),
                            var=pd.DataFrame(bm_var))
bm_adata.obs.index = list(bm_cells)
bm_adata.var.index = list(bm_genes)
bm_adata.obs['batch'] = bm_adata.obs['batch'].astype('category')
scvi.model.SCVI.setup_anndata(bm_adata, batch_key='batch')
bm_model = scvi.model.SCVI(bm_adata, n_latent=int(bm_n_latent))
bm_model.train(max_epochs=400, early_stopping=True)
bm_latent = bm_model.get_latent_representation()
")

      latent <- as.matrix(py_to_r(py$bm_latent))
      rownames(latent) <- colnames(temp)
      colnames(latent) <- paste0("scVI_", 1:ncol(latent))

      py_run_string("
try:
    del bm_adata, bm_model, bm_latent, bm_counts, bm_obs, bm_var, bm_cells, bm_genes, bm_n_latent
    import gc; gc.collect()
except:
    pass
")
      gc()
      cat("[OK]\n")
      return(latent)

    # ==================================================================
    # scCobra -- Python, uses raw counts (not valid for normalization vote)
    # ==================================================================
    } else if (method == "sccobra") {
      sccobra_path <- "/scicore/home/doetsch/kaiser0001/GITHUB_repositories/scCobra"
      library(reticulate)
      sccobra_available <- tryCatch({
        if (!dir.exists(sccobra_path)) return(FALSE)
        py_run_string(sprintf("import sys; sys.path.insert(0, '%s')", sccobra_path))
        py_run_string("from scCobra import scCobra as _test_sccobra_bm")
        TRUE
      }, error = function(e) FALSE)

      if (!sccobra_available) {
        cat("[scCobra not available]\n")
        return(NULL)
      }

      temp <- obj
      Seurat::DefaultAssay(temp) <- "RNA"
      tryCatch({
        temp[["RNA"]] <- SeuratObject::JoinLayers(temp[["RNA"]])
      }, error = function(e) NULL)

      counts_mat <- as.matrix(Seurat::GetAssayData(temp, assay = "RNA", layer = "counts"))
      counts_mat_t <- t(counts_mat)
      obs_df <- data.frame(batch = as.character(temp@meta.data[[batch_var]]),
                           row.names = colnames(temp), stringsAsFactors = FALSE)
      var_df <- data.frame(gene = rownames(counts_mat), row.names = rownames(counts_mat),
                           stringsAsFactors = FALSE)

      bm_outdir <- file.path(tempdir(), "bm_sccobra_output")
      dir.create(bm_outdir, showWarnings = FALSE, recursive = TRUE)

      py_run_string("import numpy as np; import pandas as pd; import anndata")
      py$bm_counts <- r_to_py(counts_mat_t)
      py$bm_obs <- r_to_py(obs_df)
      py$bm_var <- r_to_py(var_df)
      py$bm_cells <- r_to_py(colnames(temp))
      py$bm_genes <- r_to_py(rownames(counts_mat))
      py$bm_outdir <- r_to_py(bm_outdir)

      py_run_string("
bm_adata = anndata.AnnData(X=np.array(bm_counts, dtype='float32'),
                            obs=pd.DataFrame(bm_obs),
                            var=pd.DataFrame(bm_var))
bm_adata.obs.index = list(bm_cells)
bm_adata.var.index = list(bm_genes)
bm_adata.obs['batch'] = bm_adata.obs['batch'].astype('category')

import sys
sys.path.insert(0, '/scicore/home/doetsch/kaiser0001/GITHUB_repositories/scCobra')
from scCobra import scCobra as run_scCobra_bm

bm_adata_int = run_scCobra_bm(
    data_list=bm_adata,
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
    outdir=bm_outdir,
    ignore_umap=True,
    verbose=False
)

bm_obsm_keys = list(bm_adata_int.obsm.keys())
bm_latent_key = None
for key in ['X_scCobra', 'scCobra', 'latent', 'X_latent', 'X_emb', 'emb']:
    if key in bm_obsm_keys:
        bm_latent_key = key
        break
if bm_latent_key is None and len(bm_obsm_keys) > 0:
    bm_latent_key = bm_obsm_keys[0]
if bm_latent_key is None:
    raise ValueError('No latent embedding found in scCobra output')

bm_latent = bm_adata_int.obsm[bm_latent_key]
bm_cell_names_out = list(bm_adata_int.obs_names)
")

      latent_rep <- as.matrix(py_to_r(py$bm_latent))
      cell_names_out <- py_to_r(py$bm_cell_names_out)

      common_cells <- intersect(cell_names_out, colnames(temp))
      if (length(common_cells) == 0 && nrow(latent_rep) == ncol(temp)) {
        rownames(latent_rep) <- colnames(temp)
      } else {
        rownames(latent_rep) <- cell_names_out
        latent_rep <- latent_rep[colnames(temp), , drop = FALSE]
      }
      colnames(latent_rep) <- paste0("scCobra_", 1:ncol(latent_rep))

      py_run_string("
try:
    del bm_adata, bm_adata_int, bm_latent, bm_counts, bm_obs, bm_var, bm_cells, bm_genes, bm_outdir, bm_obsm_keys, bm_latent_key, bm_cell_names_out
    import gc; gc.collect()
except:
    pass
")
      unlink(bm_outdir, recursive = TRUE)
      gc()
      cat("[OK]\n")
      return(latent_rep)

    # ==================================================================
    # Concord -- Python, uses raw counts (not valid for normalization vote)
    # ==================================================================
    } else if (method == "concord") {
      concord_available <- tryCatch({
        library(reticulate)
        py_run_string("import concord")
        TRUE
      }, error = function(e) FALSE)

      if (!concord_available) {
        cat("[concord not available]\n")
        return(NULL)
      }

      temp <- obj
      Seurat::DefaultAssay(temp) <- "RNA"
      tryCatch({
        temp[["RNA"]] <- SeuratObject::JoinLayers(temp[["RNA"]])
      }, error = function(e) NULL)

      counts_data <- as.matrix(Seurat::GetAssayData(temp, assay = "RNA", layer = "counts"))
      counts_data_t <- t(counts_data)
      obs_df <- data.frame(batch = as.character(temp@meta.data[[batch_var]]),
                           row.names = colnames(temp), stringsAsFactors = FALSE)

      bm_outdir <- file.path(tempdir(), "bm_concord_output")
      dir.create(bm_outdir, showWarnings = FALSE, recursive = TRUE)

      n_top <- if (!is.null(params$concord_n_top_features)) params$concord_n_top_features else 2000
      n_lat <- if (!is.null(params$concord_n_latent)) params$concord_n_latent else 100
      n_ep <- if (!is.null(params$concord_max_epochs)) params$concord_max_epochs else 15
      bsz <- if (!is.null(params$concord_batch_size)) params$concord_batch_size else 256
      lr_val <- if (!is.null(params$concord_lr)) params$concord_lr else 1e-2

      py_run_string("import numpy as np; import pandas as pd; import anndata; import scanpy as sc; import torch; import concord as ccd")
      py$bm_counts <- r_to_py(counts_data_t)
      py$bm_obs <- r_to_py(obs_df)
      py$bm_cells <- r_to_py(colnames(temp))
      py$bm_genes <- r_to_py(rownames(counts_data))
      py$bm_save_dir <- r_to_py(bm_outdir)
      py$bm_n_top <- r_to_py(as.integer(n_top))
      py$bm_latent_dim <- r_to_py(as.integer(n_lat))
      py$bm_n_epochs <- r_to_py(as.integer(n_ep))
      py$bm_batch_size <- r_to_py(as.integer(bsz))
      py$bm_lr <- r_to_py(lr_val)

      py_run_string("
bm_adata = anndata.AnnData(X=np.array(bm_counts, dtype='float32'),
                            obs=pd.DataFrame(bm_obs))
bm_adata.obs.index = list(bm_cells)
bm_adata.var.index = list(bm_genes)
bm_adata.obs['batch'] = bm_adata.obs['batch'].astype('category')

bm_feature_list = ccd.ul.select_features(bm_adata, n_top_features=int(bm_n_top), flavor='seurat_v3')
sc.pp.normalize_total(bm_adata)
sc.pp.log1p(bm_adata)

if torch.cuda.is_available():
    bm_device = torch.device('cuda:0')
elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
    bm_device = torch.device('mps')
else:
    bm_device = torch.device('cpu')

bm_ccd = ccd.Concord(
    adata=bm_adata,
    save_dir=bm_save_dir,
    verbose=False,
    input_feature=bm_feature_list,
    domain_key='batch',
    latent_dim=int(bm_latent_dim),
    n_epochs=int(bm_n_epochs),
    batch_size=int(bm_batch_size),
    lr=float(bm_lr),
    device=bm_device
)
bm_ccd.fit_transform(output_key='Concord')

bm_latent = bm_adata.obsm['Concord']
bm_cell_names_out = list(bm_adata.obs_names)
")

      latent_rep <- as.matrix(py_to_r(py$bm_latent))
      cell_names_out <- py_to_r(py$bm_cell_names_out)
      rownames(latent_rep) <- cell_names_out
      colnames(latent_rep) <- paste0("CONCORD_", 1:ncol(latent_rep))
      latent_rep <- latent_rep[colnames(temp), , drop = FALSE]

      py_run_string("
try:
    del bm_adata, bm_ccd, bm_latent, bm_counts, bm_obs, bm_cells, bm_genes, bm_feature_list, bm_cell_names_out, bm_save_dir
    import gc; gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
except:
    pass
")
      unlink(bm_outdir, recursive = TRUE)
      gc()
      cat("[OK]\n")
      return(latent_rep)
    }

    cat("[unknown method]\n")
    return(NULL)

  }, error = function(e) {
    cat("[FAILED:", conditionMessage(e), "]\n")
    return(NULL)
  })
}



# ==============================================================================
# ADDITIONS FOR MODULE 08 (DE HELPERS)
# ==============================================================================

# Safely get assay matrix across Seurat v4/v5 (layer vs slot)
safe_get_assay_mat <- function(seu, assay = NULL, layer = "data", slot_fallback = "data") {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(seu)

  m <- tryCatch({
    Seurat::GetAssayData(seu, assay = assay, layer = layer)
  }, error = function(e) NULL)

  if (is.null(m)) {
    m <- tryCatch({
      Seurat::GetAssayData(seu, assay = assay, slot = slot_fallback)
    }, error = function(e) NULL)
  }
  m
}

# Print first 5 rows of a CSV that was written (control check)
print_csv_head5 <- function(csv_path, label = NULL) {
  if (!file.exists(csv_path)) {
    cat("  [WARNING] CSV not found for head() check:", csv_path, "\n")
    return(invisible(NULL))
  }
  cat("\n--- CONTROL CHECK: First 5 rows of saved CSV",
      if (!is.null(label)) paste0(" (", label, ")") else "",
      " ---\n", sep = "")
  df <- tryCatch(utils::read.csv(csv_path, check.names = FALSE), error = function(e) NULL)
  if (is.null(df)) {
    cat("  [WARNING] Could not read CSV:", csv_path, "\n")
    return(invisible(NULL))
  }
  print(utils::head(df, 5))
  invisible(df)
}

# Add up/down direction column based on sign of LFC
add_direction_column <- function(df, lfc_col, group1_label = "Group1", group2_label = "Group2") {
  if (!lfc_col %in% colnames(df)) return(df)
  df$direction <- NA_character_
  df$direction[df[[lfc_col]] > 0] <- paste0("Up in ", group1_label)
  df$direction[df[[lfc_col]] < 0] <- paste0("Up in ", group2_label)
  df$direction[df[[lfc_col]] == 0] <- "No change"
  df
}

# Detect cluster column consistently
get_cluster_column <- function(obj) {
  for (col in c("leiden_clusters", "seurat_clusters")) {
    if (col %in% colnames(obj@meta.data)) return(col)
  }
  choir_cols <- grep("^CHOIR_clusters", colnames(obj@meta.data), value = TRUE)
  if (length(choir_cols) > 0) return(choir_cols[1])
  return(NULL)
}

# Load a Seurat object automatically from objects dir
load_seurat_object_auto <- function(objects_dir) {
  # Priority list (common pipeline outputs)
  preferred <- c(
    "07_final_object.rds",
    "07_leiden_final_object.rds",
    "07_leiden_clustered_object.rds",
    "07_leiden_object.rds",
    "scice_subclustered_object.rds",
    "choir_clustered_object.rds",
    "03_integrated_object.rds",
    "02_qc_filtered_object.rds"
  )
  preferred_paths <- file.path(objects_dir, preferred)
  preferred_paths <- preferred_paths[file.exists(preferred_paths)]

  # Any RDS
  rds_all <- list.files(objects_dir, pattern = "\\.rds$", full.names = TRUE)
  rds_paths <- unique(c(preferred_paths, rds_all))
  if (length(rds_paths) > 0) {
    # Use most recently modified among candidates
    info <- file.info(rds_paths)
    pick <- rds_paths[which.max(info$mtime)]
    obj <- readRDS(pick)
    if (!inherits(obj, "Seurat")) stop("File read but not a Seurat object: ", pick)
    return(list(obj = obj, path = pick))
  }

  # Else: RData containing Seurat
  rdata_paths <- list.files(objects_dir, pattern = "\\.RData$", full.names = TRUE)
  if (length(rdata_paths) == 0) return(NULL)

  info <- file.info(rdata_paths)
  pick <- rdata_paths[which.max(info$mtime)]

  e <- new.env(parent = emptyenv())
  load(pick, envir = e)
  seurat_names <- ls(e)[sapply(ls(e), function(n) inherits(e[[n]], "Seurat"))]

  if (length(seurat_names) == 0) {
    stop("Loaded RData but found no Seurat object inside: ", pick)
  }

  obj <- e[[seurat_names[1]]]
  return(list(obj = obj, path = pick, object_name = seurat_names[1]))
}

# MAST log2FC diagnostic
mast_log2fc_diagnostic <- function(seu, mast_df,
                                  group_var = "sex",
                                  group1 = "Male",
                                  group2 = "Female",
                                  lfc_col = "avg_log2FC",
                                  n_genes = 5) {
  if (is.null(mast_df) || nrow(mast_df) == 0) return(invisible(NULL))
  if (!lfc_col %in% colnames(mast_df)) {
    cat("\n[MAST log2FC DIAGNOSTIC] Column", lfc_col, "not present.\n")
    return(invisible(NULL))
  }

  cat("\n[MAST log2FC DIAGNOSTIC]\n")
  cat("  Seurat reports 'avg_log2FC' for MAST results.\n")
  cat("  This is derived from average expression in the normalized 'data' layer (log scale),\n")
  cat("  so values can look smaller than fold-changes computed directly from raw counts.\n")

  dat <- safe_get_assay_mat(seu, assay = Seurat::DefaultAssay(seu), layer = "data", slot_fallback = "data")
  if (is.null(dat)) {
    cat("  [WARNING] Cannot access assay 'data' layer/slot for diagnostic table.\n")
    return(invisible(NULL))
  }

  mast_df$gene <- if ("gene" %in% colnames(mast_df)) mast_df$gene else rownames(mast_df)
  ord <- order(abs(mast_df[[lfc_col]]), decreasing = TRUE)
  demo_genes <- mast_df$gene[ord]
  demo_genes <- demo_genes[!is.na(demo_genes)]
  demo_genes <- head(demo_genes, n_genes)

  cells1 <- colnames(seu)[seu@meta.data[[group_var]] == group1]
  cells2 <- colnames(seu)[seu@meta.data[[group_var]] == group2]

  demo <- data.frame(
    gene = demo_genes,
    mean_logdata_group1 = NA_real_,
    mean_logdata_group2 = NA_real_,
    diff_logdata = NA_real_,
    reported_avg_log2FC = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(demo_genes)) {
    g <- demo_genes[i]
    if (!g %in% rownames(dat)) next
    m1 <- mean(dat[g, cells1, drop = TRUE])
    m2 <- mean(dat[g, cells2, drop = TRUE])
    demo$mean_logdata_group1[i] <- m1
    demo$mean_logdata_group2[i] <- m2
    demo$diff_logdata[i] <- (m1 - m2)

    # Match row either by rowname or gene column
    idx <- which(mast_df$gene == g)
    if (length(idx) > 0) demo$reported_avg_log2FC[i] <- mast_df[[lfc_col]][idx[1]]
  }

  cat("  Example (means from assay 'data' log-normalized expression):\n")
  print(demo)
  invisible(demo)
}

# ==============================================================================
# FUNCTIONS MOVED FROM params.R
# ==============================================================================

# ==============================================================================
# LOAD SAMPLE SHEET
# ==============================================================================
#' Load and validate sample sheet for pipeline
#'
#' @param sample_sheet_path Path to sample sheet CSV
#' @param input_dir Directory containing input files
#' @param input_file_pattern Pattern for input file names with {sample_name} placeholder
#' @param files_in_subdirectories Whether files are in sample subdirectories
#' @param ventricle_filter Optional ventricle filter ("LV" or "4V")
#' @return Data frame with sample metadata and input file paths
load_sample_sheet <- function(sample_sheet_path, input_dir, input_file_pattern,
                               files_in_subdirectories, ventricle_filter = "") {

  if (!file.exists(sample_sheet_path)) {
    stop("Sample sheet not found: ", sample_sheet_path,
         "\nPlease create a CSV with columns: sample_name, sex, batch")
  }

  cat("Loading sample sheet from:", sample_sheet_path, "\n")

  sample_sheet <- read.csv(sample_sheet_path, stringsAsFactors = FALSE)

  # Handle different column naming conventions
  if ("sample_id" %in% colnames(sample_sheet) && !"sample_name" %in% colnames(sample_sheet)) {
    sample_sheet$sample_name <- sample_sheet$sample_id
    cat("  Mapped column: sample_id -> sample_name\n")
  }

  # Add include column if missing (default TRUE)
  if (!"include" %in% colnames(sample_sheet)) {
    sample_sheet$include <- TRUE
  }

  # Convert include to logical
  sample_sheet$include <- as.logical(sample_sheet$include)

  # Filter by ventricle if specified
  if (ventricle_filter != "" && "ventricle" %in% colnames(sample_sheet)) {
    cat("Filtering samples by ventricle:", ventricle_filter, "\n")
    sample_sheet <- sample_sheet[sample_sheet$ventricle == ventricle_filter, ]

    if (nrow(sample_sheet) == 0) {
      stop("No samples found for ventricle: ", ventricle_filter)
    }

    cat("Samples after filtering:", paste(sample_sheet$sample_name, collapse = ", "), "\n")
  }

  # Validate required columns
  required_cols <- c("sample_name", "sex", "batch")
  missing_cols <- setdiff(required_cols, colnames(sample_sheet))
  if (length(missing_cols) > 0) {
    stop("Sample sheet missing required columns: ", paste(missing_cols, collapse = ", "),
         "\nRequired columns: sample_name, sex, batch")
  }

  # Standardize sex values
  sample_sheet$sex <- tools::toTitleCase(tolower(sample_sheet$sex))
  valid_sex <- c("Male", "Female")
  invalid_sex <- setdiff(unique(sample_sheet$sex), valid_sex)
  if (length(invalid_sex) > 0) {
    stop("Invalid sex values: ", paste(invalid_sex, collapse = ", "),
         "\nAllowed values: Male, Female")
  }

  # Generate input file paths
  sample_sheet$input_file <- sapply(sample_sheet$sample_name, function(s) {
    filename <- gsub("\\{sample_name\\}", s, input_file_pattern)
    if (files_in_subdirectories) {
      file.path(input_dir, s, filename)
    } else {
      file.path(input_dir, filename)
    }
  })

  return(sample_sheet)
}

# ==============================================================================
# GET SAMPLES TO ANALYZE
# ==============================================================================
#' Get list of samples marked for inclusion in analysis
#'
#' @param params Pipeline parameters list containing sample_metadata
#' @return Character vector of sample names to analyze
get_samples_to_analyze <- function(params) {
  included_samples <- params$sample_metadata$sample_name[params$sample_metadata$include == TRUE]
  if (length(included_samples) == 0) {
    stop("No samples marked for inclusion in sample sheet")
  }
  return(included_samples)
}

# ==============================================================================
# GET INPUT PATHS
# ==============================================================================
#' Generate input file paths for all samples
#'
#' @param params Pipeline parameters list
#' @return Named character vector of input file paths
get_input_paths <- function(params) {
  samples <- params$samples_to_analyze
  meta <- params$sample_metadata

  paths <- sapply(samples, function(s) {
    idx <- which(meta$sample_name == s)
    if (length(idx) > 0 && "input_file" %in% colnames(meta)) {
      return(meta$input_file[idx])
    } else {
      filename <- gsub("\\{sample_name\\}", s, params$input_file_pattern)
      if (params$files_in_subdirectories) {
        return(file.path(params$input_dir, s, filename))
      } else {
        return(file.path(params$input_dir, filename))
      }
    }
  })
  names(paths) <- samples
  return(paths)
}

# ==============================================================================
# VALIDATE PARAMS
# ==============================================================================
#' Validate pipeline parameters and print configuration summary
#'
#' @param params Pipeline parameters list
#' @return Validated and potentially modified params list
validate_params <- function(params) {
  cat("\n")
  cat("================================================================================\n")
  cat("VALIDATING PIPELINE CONFIGURATION\n")
  cat("================================================================================\n\n")

  # Ensure logical parameters are properly set
  logical_params <- c(
    "run_mast", "run_pseudobulk_edger", "run_pseudobulk_deseq2",
    "run_integration_benchmarking", "run_choir_clustering",
    "run_python_integrations", "run_sctransform", "run_scran",
    "run_lognorm", "run_sckwarn", "run_clustering_quality", "run_normalization_benchmarking",
    "run_leiden_clustering", "run_sccobra", "run_batch_integration",
    "run_scice_subclustering", "run_idclust_subclustering", "filter_hemoglobin",
    "skip_all_filtering", "apply_qc_filtering",
    "filter_by_min_features", "filter_by_max_features",
    "filter_by_percent_mt", "filter_genes_by_min_cells",
    "filter_doublets", "use_doublet_vote_threshold",
    "files_in_subdirectories"
  )

  for (p in logical_params) {
    if (p %in% names(params)) {
      if (is.null(params[[p]])) {
        params[[p]] <- FALSE
      } else if (!is.logical(params[[p]])) {
        params[[p]] <- as.logical(params[[p]])
      }
    } else {
      params[[p]] <- FALSE
    }
  }

  # Validate integration_selection_mode
  valid_selection_modes <- c("batch_removal", "balanced", "conservative")
  if (!is.null(params$integration_selection_mode)) {
    if (!params$integration_selection_mode %in% valid_selection_modes) {
      cat("WARNING: Invalid integration_selection_mode '", params$integration_selection_mode,
          "'. Using 'balanced'.\n", sep = "")
      params$integration_selection_mode <- "balanced"
    }
  } else {
    params$integration_selection_mode <- "balanced"
  }

  # Print configuration
  cat("--- Project Configuration ---\n")
  cat("  Project root:", params$project_root, "\n")
  cat("  Dataset name:", params$dataset_name, "\n")
  cat("  Output directory:", params$out_root, "\n")
  cat("\n")

  cat("--- Input Configuration ---\n")
  cat("  Input directory:", params$input_dir, "\n")
  cat("  File pattern:", params$input_file_pattern, "\n")
  cat("  Files in subdirectories:", params$files_in_subdirectories, "\n")
  cat("  Counts layer:", params$counts_layer_to_use, "\n")
  cat("\n")

  # Validate input files
  cat("--- Validating Input Files ---\n")
  all_exist <- TRUE
  for (sample in names(params$input_paths)) {
    path <- params$input_paths[[sample]]
    if (!file.exists(path)) {
      cat("  [MISSING]", sample, ":", path, "\n")
      all_exist <- FALSE
    } else {
      size_mb <- round(file.info(path)$size / (1024^2), 1)
      cat("  [OK]", sample, "(", size_mb, "MB)\n")
    }
  }

  if (!all_exist) {
    cat("\n  NOTE: Some input files are missing.\n")
    cat("  This is expected if preprocessing has not completed yet.\n")
  }
  cat("\n")

  # Print filtering configuration
  cat("--- Filtering Configuration ---\n")
  if (isTRUE(params$skip_all_filtering)) {
    cat("  >>> ALL FILTERING DISABLED <<<\n")
  } else {
    cat("  apply_qc_filtering:", params$apply_qc_filtering, "\n")
    if (isTRUE(params$apply_qc_filtering)) {
      cat("    filter_by_min_features:", params$filter_by_min_features, "(", params$min_features, ")\n")
      cat("    filter_by_max_features:", params$filter_by_max_features, "(", params$max_features, ")\n")
      cat("    filter_by_percent_mt:", params$filter_by_percent_mt, "(", params$max_percent_mt, "%)\n")
    }
    cat("  filter_genes_by_min_cells:", params$filter_genes_by_min_cells, "(", params$min_cells_per_gene, ")\n")
    cat("  filter_doublets:", params$filter_doublets, "\n")
    cat("  filter_hemoglobin:", params$filter_hemoglobin, "\n")
  }
  cat("\n")

  # Print integration selection mode
  cat("--- Integration Selection Mode ---\n")
  cat("  Mode:", params$integration_selection_mode, "\n")
  if (params$integration_selection_mode == "batch_removal") {
    cat("  Strategy: Minimize batch_variance (aggressive batch removal)\n")
    cat("  Best for: Technical replicates (same biology, different processing)\n")
    cat("  WARNING: May over-correct if batches have biological meaning!\n")
  } else if (params$integration_selection_mode == "balanced") {
    cat("  Strategy: Maximize composite score across all batch metrics\n")
    cat("  Best for: Biological batches (different conditions/sexes/timepoints)\n")
    cat("  >>> RECOMMENDED for most studies <<<\n")
  } else if (params$integration_selection_mode == "conservative") {
    cat("  Strategy: Prioritize LISI (moderate mixing)\n")
    cat("  Best for: Preserving subtle biological differences\n")
  }
  cat("\n")

  return(params)
}

# ==============================================================================
# RUN scAURA CLUSTERING VIA PYTHON
# ==============================================================================
#' Run scAURA clustering on a Seurat object
#'
#' This function exports data to CSV, calls scaura_api.py via command line,
#' and loads results back into R.
#'
#' REPLACES the previous version that generated a wrapper script importing
#' nonexistent modules (model, utils, train) from the scAURA repo.
#'
#' REQUIRES: scaura_api.py placed in the scAURA repo directory alongside
#'           scAURA_gpu.py (or pointed to via params$scaura_api_script).
#'
#' @param obj Seurat object with normalized data and HVGs
#' @param params Pipeline parameters containing scAURA settings
#' @param output_dir Directory for scAURA outputs
#' @return List with embeddings (pre/post SSC), kmeans labels, ssc labels
run_scaura_clustering <- function(obj, params, output_dir) {

  cat("\n--- Running scAURA clustering ---\n")

  # Validate scAURA repository
  scaura_path <- params$scaura_repo_path
  if (!dir.exists(scaura_path)) {
    stop("scAURA repository not found at: ", scaura_path,
         "\nClone from: https://github.com/bozdaglab/scAURA")
  }

  # Locate the API script (scaura_api.py)
  api_script <- if (!is.null(params$scaura_api_script) && file.exists(params$scaura_api_script)) {
    params$scaura_api_script
  } else {
    file.path(scaura_path, "scaura_api.py")
  }

  if (!file.exists(api_script)) {
    stop("scaura_api.py not found at: ", api_script,
         "\nPlace scaura_api.py in the scAURA repo directory: ", scaura_path)
  }
  cat("  scAURA API script:", api_script, "\n")

  # Create output directory
  scaura_out <- file.path(output_dir, "scaura_output")
  dir.create(scaura_out, showWarnings = FALSE, recursive = TRUE)

  # Get HVGs
  hvgs <- Seurat::VariableFeatures(obj)
  if (length(hvgs) == 0) {
    cat("  No variable features found, running FindVariableFeatures...\n")
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = params$scaura_n_top_genes, verbose = FALSE)
    hvgs <- Seurat::VariableFeatures(obj)
  }

  n_genes <- min(length(hvgs), params$scaura_n_top_genes)
  hvgs <- hvgs[1:n_genes]
  cat("  Using", n_genes, "highly variable genes\n")

  # Extract normalized data for HVGs
  norm_data <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
  norm_data <- norm_data[hvgs, , drop = FALSE]
  norm_data <- as.matrix(norm_data)

  # Transpose: scAURA expects cells x genes
  data_for_scaura <- t(norm_data)
  cat("  Data matrix:", nrow(data_for_scaura), "cells x", ncol(data_for_scaura), "genes\n")

  # Save as CSV (cells as rows, first column = cell_id)
  input_csv <- file.path(scaura_out, "input_data.csv")
  input_df <- data.frame(cell_id = rownames(data_for_scaura), data_for_scaura, check.names = FALSE)
  write.csv(input_df, input_csv, row.names = FALSE)
  cat("  Saved input CSV:", input_csv, "\n")

  # Build CLI arguments
  python_cmd <- params$scaura_python

  cli_args <- c(
    api_script,
    "--input", input_csv,
    "--output", scaura_out,
    "--k", as.integer(params$scaura_k_clusters),
    "--seed", 42,
    "--hidden-dim", as.integer(params$scaura_hidden_dim),
    "--kmax", as.integer(params$scaura_kmax),
    "--tau", params$scaura_tau,
    "--tau-plus", params$scaura_tau_plus,
    "--pe", params$scaura_pe,
    "--pf", params$scaura_pf,
    "--lr", params$scaura_lr,
    "--epochs", as.integer(params$scaura_epochs),
    "--activation", if (!is.null(params$scaura_activation)) params$scaura_activation else "tanh",
    "--ssc-batchsize", as.integer(if (!is.null(params$scaura_ssc_batchsize)) params$scaura_ssc_batchsize else 64)
  )

  # GPU control
  if (isTRUE(params$scaura_use_gpu)) {
    gpu_id <- if (!is.null(params$scaura_gpu_id)) params$scaura_gpu_id else 0
    cli_args <- c(cli_args, "--gpu", as.integer(gpu_id))
  }

  # Self-training control
  if (!isTRUE(params$scaura_self_train)) {
    cli_args <- c(cli_args, "--no-self-train")
  }

  cat("  Command:", python_cmd, paste(basename(api_script), cli_args[-1], collapse = " "), "\n")

  # Execute scAURA
  cat("  Executing scAURA (this may take several minutes)...\n")

  result <- system2(
    python_cmd,
    args = cli_args,
    stdout = TRUE,
    stderr = TRUE
  )

  exit_code <- attr(result, "status")
  if (is.null(exit_code)) exit_code <- 0

  # Print Python output
  cat("  scAURA output:\n")
  for (line in result) {
    cat("    ", line, "\n")
  }

  if (exit_code != 0) {
    stop("scAURA failed with exit code: ", exit_code)
  }

  # Check for output files
  labels_file <- file.path(scaura_out, "cluster_labels.csv")
  embeddings_pre_file <- file.path(scaura_out, "embeddings_pre_ssc.npy")
  embeddings_post_file <- file.path(scaura_out, "embeddings_post_ssc.npy")

  if (!file.exists(labels_file)) {
    stop("scAURA failed: cluster_labels.csv not found in ", scaura_out)
  }

  # Load results using reticulate
  cat("  Loading scAURA results...\n")

  reticulate::use_python(python_cmd, required = TRUE)
  np <- reticulate::import("numpy")

  embeddings_pre <- np$load(embeddings_pre_file)
  embeddings_post <- np$load(embeddings_post_file)
  labels_df <- read.csv(labels_file, stringsAsFactors = FALSE)

  # Set row/col names
  rownames(embeddings_pre) <- labels_df$cell_id
  rownames(embeddings_post) <- labels_df$cell_id
  colnames(embeddings_pre) <- paste0("scAURA_", 1:ncol(embeddings_pre))
  colnames(embeddings_post) <- paste0("scAURA_", 1:ncol(embeddings_post))

  # Reorder to match Seurat object
  cell_order <- colnames(obj)
  embeddings_pre <- embeddings_pre[cell_order, , drop = FALSE]
  embeddings_post <- embeddings_post[cell_order, , drop = FALSE]
  labels_df <- labels_df[match(cell_order, labels_df$cell_id), ]

  cat("  scAURA results loaded successfully\n")
  cat("  Embedding dimensions:", ncol(embeddings_post), "\n")
  cat("  Clusters found:", length(unique(labels_df$ssc_cluster)), "\n")

  return(list(
    embeddings_pre_ssc = embeddings_pre,
    embeddings_post_ssc = embeddings_post,
    kmeans_labels = labels_df$kmeans_cluster,
    ssc_labels = labels_df$ssc_cluster,
    cell_ids = labels_df$cell_id,
    output_dir = scaura_out
  ))
}


# ==============================================================================
# CHANGE 1: Add these 4 functions to functions.R (copy from 04_integration.R)
# ==============================================================================
# INSERT LOCATION: After compute_lisi_score() in functions.R
#
# These are EXACT COPIES of the functions currently in 04_integration.R.
# After adding them to functions.R, REMOVE them from 04_integration.R to
# avoid double-definition.
# ==============================================================================

#' Compute cell type silhouette width (biological conservation)
#' Higher values = better separation of cell types = good
compute_celltype_asw <- function(embeddings, celltype_labels) {
  tryCatch({
    valid_idx <- !is.na(celltype_labels)
    if (sum(valid_idx) < 10) return(NA_real_)

    emb <- embeddings[valid_idx, , drop = FALSE]
    labels <- as.character(celltype_labels[valid_idx])

    if (length(unique(labels)) < 2) return(NA_real_)

    n_cells <- nrow(emb)
    if (n_cells > 5000) {
      set.seed(42)
      sample_idx <- sample(n_cells, 5000)
      emb <- emb[sample_idx, , drop = FALSE]
      labels <- labels[sample_idx]
    }

    dist_mat <- dist(emb)
    sil <- cluster::silhouette(as.integer(factor(labels)), dist_mat)
    mean(sil[, 3], na.rm = TRUE)
  }, error = function(e) {
    warning("celltype_asw computation failed: ", e$message)
    NA_real_
  })
}

#' Compute cell type LISI (local purity - biological conservation)
#' Lower values = better local purity of cell types = good
compute_celltype_lisi <- function(embeddings, celltype_labels, perplexity = 30) {
  tryCatch({
    valid_idx <- !is.na(celltype_labels)
    if (sum(valid_idx) < 10) return(NA_real_)

    emb <- embeddings[valid_idx, , drop = FALSE]
    labels <- as.character(celltype_labels[valid_idx])

    if (length(unique(labels)) < 2) return(NA_real_)

    n_cells <- nrow(emb)
    if (n_cells > 5000) {
      set.seed(42)
      sample_idx <- sample(n_cells, 5000)
      emb <- emb[sample_idx, , drop = FALSE]
      labels <- labels[sample_idx]
    }

    if (requireNamespace("lisi", quietly = TRUE)) {
      meta_df <- data.frame(celltype = labels, row.names = rownames(emb))
      lisi_result <- lisi::compute_lisi(emb, meta_df, "celltype", perplexity = perplexity)
      return(mean(lisi_result$celltype, na.rm = TRUE))
    }

    # Fallback: simplified LISI via k-NN
    k <- min(perplexity * 3, nrow(emb) - 1)
    dist_mat <- as.matrix(dist(emb))

    lisi_values <- sapply(1:nrow(emb), function(i) {
      neighbors <- order(dist_mat[i, ])[2:(k + 1)]
      neighbor_labels <- labels[neighbors]
      freqs <- table(neighbor_labels) / length(neighbor_labels)
      1 / sum(freqs^2)
    })

    mean(lisi_values, na.rm = TRUE)
  }, error = function(e) {
    warning("celltype_lisi computation failed: ", e$message)
    NA_real_
  })
}

#' Compute NMI between clustering and cell type annotations
#' Higher values = better agreement = good
compute_nmi_score <- function(embeddings, celltype_labels, resolution = 0.8) {
  tryCatch({
    valid_idx <- !is.na(celltype_labels)
    if (sum(valid_idx) < 10) return(NA_real_)

    emb <- embeddings[valid_idx, , drop = FALSE]
    labels <- as.character(celltype_labels[valid_idx])

    if (length(unique(labels)) < 2) return(NA_real_)

    # Build kNN graph directly on the embedding
    k <- min(20, nrow(emb) - 1)
    nn <- RANN::nn2(emb, k = k + 1)  # +1 because first neighbor is self
    nn_idx <- nn$nn.idx[, -1, drop = FALSE]  # remove self

    # Build SNN graph: edge weight = |shared neighbors| / (2k - |shared|)
    edges_from <- c()
    edges_to <- c()
    edges_weight <- c()

    for (i in 1:nrow(nn_idx)) {
      neighbors_i <- nn_idx[i, ]
      for (j in neighbors_i) {
        if (j > i) {
          shared <- length(intersect(nn_idx[i, ], nn_idx[j, ]))
          if (shared > 0) {
            weight <- shared / (2 * k - shared)
            edges_from <- c(edges_from, i)
            edges_to <- c(edges_to, j)
            edges_weight <- c(edges_weight, weight)
          }
        }
      }
    }

    if (length(edges_from) == 0) return(NA_real_)

    g <- igraph::make_empty_graph(n = nrow(emb), directed = FALSE)
    g <- igraph::add_edges(g, as.vector(rbind(edges_from, edges_to)))
    igraph::E(g)$weight <- edges_weight

    # Leiden clustering on the SNN graph
    clustering <- igraph::cluster_leiden(g, resolution_parameter = resolution,
                                          objective_function = "modularity")
    clusters <- as.character(igraph::membership(clustering))

    # Use aricode if available for fast NMI
    if (requireNamespace("aricode", quietly = TRUE)) {
      return(aricode::NMI(clusters, labels))
    }

    # Manual NMI computation
    compute_entropy <- function(x) {
      freqs <- table(x) / length(x)
      -sum(freqs * log(freqs + 1e-10))
    }

    compute_mutual_info <- function(x, y) {
      joint <- table(x, y) / length(x)
      px <- rowSums(joint)
      py <- colSums(joint)
      mi <- 0
      for (i in 1:nrow(joint)) {
        for (j in 1:ncol(joint)) {
          if (joint[i, j] > 0) {
            mi <- mi + joint[i, j] * log(joint[i, j] / (px[i] * py[j] + 1e-10) + 1e-10)
          }
        }
      }
      mi
    }

    h_clusters <- compute_entropy(clusters)
    h_labels <- compute_entropy(labels)
    mi <- compute_mutual_info(clusters, labels)

    nmi <- 2 * mi / (h_clusters + h_labels + 1e-10)
    return(nmi)
  }, error = function(e) {
    warning("NMI computation failed: ", e$message)
    NA_real_
  })
}

#' Compute ARI between clustering and cell type annotations
#' Higher values = better agreement = good
compute_ari_score <- function(embeddings, celltype_labels, resolution = 0.8) {
  tryCatch({
    valid_idx <- !is.na(celltype_labels)
    if (sum(valid_idx) < 10) return(NA_real_)

    emb <- embeddings[valid_idx, , drop = FALSE]
    labels <- as.character(celltype_labels[valid_idx])

    if (length(unique(labels)) < 2) return(NA_real_)

    # Build kNN graph directly on the embedding
    k <- min(20, nrow(emb) - 1)
    nn <- RANN::nn2(emb, k = k + 1)  # +1 because first neighbor is self
    nn_idx <- nn$nn.idx[, -1, drop = FALSE]  # remove self

    # Build SNN graph: edge weight = |shared neighbors| / (2k - |shared|)
    edges_from <- c()
    edges_to <- c()
    edges_weight <- c()

    for (i in 1:nrow(nn_idx)) {
      neighbors_i <- nn_idx[i, ]
      for (j in neighbors_i) {
        if (j > i) {
          shared <- length(intersect(nn_idx[i, ], nn_idx[j, ]))
          if (shared > 0) {
            weight <- shared / (2 * k - shared)
            edges_from <- c(edges_from, i)
            edges_to <- c(edges_to, j)
            edges_weight <- c(edges_weight, weight)
          }
        }
      }
    }

    if (length(edges_from) == 0) return(NA_real_)

    g <- igraph::make_empty_graph(n = nrow(emb), directed = FALSE)
    g <- igraph::add_edges(g, as.vector(rbind(edges_from, edges_to)))
    igraph::E(g)$weight <- edges_weight

    # Leiden clustering on the SNN graph
    clustering <- igraph::cluster_leiden(g, resolution_parameter = resolution,
                                          objective_function = "modularity")
    clusters <- as.character(igraph::membership(clustering))

    # Use mclust if available for ARI
    if (requireNamespace("mclust", quietly = TRUE)) {
      return(mclust::adjustedRandIndex(clusters, labels))
    }

    # Manual ARI computation
    contingency <- table(clusters, labels)
    n <- sum(contingency)
    sum_comb_rows <- sum(choose(rowSums(contingency), 2))
    sum_comb_cols <- sum(choose(colSums(contingency), 2))
    sum_comb_cells <- sum(choose(contingency, 2))
    expected <- sum_comb_rows * sum_comb_cols / choose(n, 2)
    max_index <- 0.5 * (sum_comb_rows + sum_comb_cols)
    ari <- (sum_comb_cells - expected) / (max_index - expected + 1e-10)
    return(ari)
  }, error = function(e) {
    warning("ARI computation failed: ", e$message)
    NA_real_
  })
}


# ==============================================================================
# CHANGE 2: Replace select_normalization_by_majority_vote() in functions.R
# ==============================================================================
#
# The updated version:
#  - Detects if bio columns (celltype_asw, celltype_lisi, nmi, ari) are present
#  - If present, computes weighted bio+batch composite (same weights as Module 04)
#  - If absent, falls back to batch-only composite (backward compatible)
#  - The function signature is UNCHANGED -- no breaking changes
#
# FIND AND REPLACE the entire select_normalization_by_majority_vote() function.
# ==============================================================================

select_normalization_by_majority_vote <- function(full_benchmark_df, selection_mode = "balanced") {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("MAJORITY VOTE FOR BEST NORMALIZATION\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Selection mode:", selection_mode, "\n")

  # Detect if biological conservation metrics are available
  bio_cols <- c("celltype_asw", "celltype_lisi", "nmi", "ari")
  has_bio <- all(bio_cols %in% colnames(full_benchmark_df)) &&
    any(!is.na(full_benchmark_df$celltype_asw))

  if (has_bio) {
    cat("Biological conservation metrics: AVAILABLE\n")
    cat("  Using scIB-style composite (40% bio + 60% batch) for voting\n\n")
  } else {
    cat("Biological conservation metrics: NOT AVAILABLE\n")
    cat("  Using batch metrics only for voting\n\n")
  }

  integration_methods <- unique(full_benchmark_df$integration_method)

  per_integration_winners <- data.frame(
    integration_method = character(),
    best_normalization = character(),
    score = numeric(),
    batch_score = numeric(),
    bio_score = numeric(),
    stringsAsFactors = FALSE
  )

  for (int_method in integration_methods) {
    subset_df <- full_benchmark_df[full_benchmark_df$integration_method == int_method, ]

    if (nrow(subset_df) < 2) {
      cat("  ", int_method, ": only", nrow(subset_df), "normalization(s), skipping vote\n")
      if (nrow(subset_df) == 1) {
        per_integration_winners <- rbind(per_integration_winners, data.frame(
          integration_method = int_method,
          best_normalization = subset_df$normalization[1],
          score = NA_real_,
          batch_score = NA_real_,
          bio_score = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      next
    }

    # Batch metrics (always available)
    subset_df$bv_norm <- normalize_metric(subset_df$batch_variance, lower_is_better = TRUE)
    subset_df$asw_norm <- normalize_metric(abs(subset_df$batch_asw), lower_is_better = TRUE)
    subset_df$lisi_norm <- normalize_metric(subset_df$lisi, lower_is_better = FALSE)
    subset_df$batch_composite <- rowMeans(
      subset_df[, c("bv_norm", "asw_norm", "lisi_norm")], na.rm = TRUE
    )

    # Bio metrics (if available)
    subset_df$bio_composite <- NA_real_
    if (has_bio) {
      subset_df$ct_asw_norm <- normalize_metric(subset_df$celltype_asw, lower_is_better = FALSE)
      subset_df$ct_lisi_norm <- normalize_metric(subset_df$celltype_lisi, lower_is_better = TRUE)
      subset_df$nmi_norm <- normalize_metric(subset_df$nmi, lower_is_better = FALSE)
      subset_df$ari_norm <- normalize_metric(subset_df$ari, lower_is_better = FALSE)
      subset_df$bio_composite <- rowMeans(
        subset_df[, c("ct_asw_norm", "ct_lisi_norm", "nmi_norm", "ari_norm")], na.rm = TRUE
      )
    }

    # Compute final composite based on selection mode
    if (selection_mode == "batch_removal") {
      subset_df$composite <- subset_df$bv_norm
    } else if (selection_mode == "conservative") {
      if (has_bio) {
        subset_df$composite <- 0.6 * subset_df$bio_composite + 0.4 * subset_df$batch_composite
      } else {
        subset_df$composite <- subset_df$lisi_norm
      }
    } else {
      # "balanced" (default)
      if (has_bio) {
        # scIB-style: 40% bio + 60% batch
        subset_df$composite <- 0.4 * subset_df$bio_composite + 0.6 * subset_df$batch_composite
      } else {
        subset_df$composite <- rowMeans(
          subset_df[, c("bv_norm", "asw_norm", "lisi_norm")], na.rm = TRUE
        )
      }
    }

    best_idx <- which.max(subset_df$composite)

    per_integration_winners <- rbind(per_integration_winners, data.frame(
      integration_method = int_method,
      best_normalization = subset_df$normalization[best_idx],
      score = subset_df$composite[best_idx],
      batch_score = subset_df$batch_composite[best_idx],
      bio_score = if (has_bio) subset_df$bio_composite[best_idx] else NA_real_,
      stringsAsFactors = FALSE
    ))

    cat("  ", int_method, "-> best:", subset_df$normalization[best_idx],
        "(composite:", round(subset_df$composite[best_idx], 4))
    if (has_bio) {
      cat(", batch:", round(subset_df$batch_composite[best_idx], 4),
          ", bio:", round(subset_df$bio_composite[best_idx], 4))
    }
    cat(")\n")
  }

  cat("\n--- Per-Integration Winners ---\n")
  print(per_integration_winners)
  cat("\n")

  # Tally votes
  vote_table <- as.data.frame(table(per_integration_winners$best_normalization),
                              stringsAsFactors = FALSE)
  colnames(vote_table) <- c("normalization", "votes")
  vote_table <- vote_table[order(-vote_table$votes), ]

  cat("--- Vote Tally ---\n")
  for (i in 1:nrow(vote_table)) {
    cat("  ", vote_table$normalization[i], ":", vote_table$votes[i], "vote(s)\n")
  }
  cat("\n")

  # Determine winner (break ties by average composite)
  max_votes <- max(vote_table$votes)
  tied <- vote_table$normalization[vote_table$votes == max_votes]

  if (length(tied) == 1) {
    winner <- tied
  } else {
    cat("Tie detected between:", paste(tied, collapse = ", "), "\n")
    cat("Breaking tie by average composite score across voting methods...\n")
    avg_scores <- sapply(tied, function(norm) {
      rows <- per_integration_winners[per_integration_winners$best_normalization == norm, ]
      mean(rows$score, na.rm = TRUE)
    })
    winner <- tied[which.max(avg_scores)]
    cat("Tie-break winner:", winner, "with avg score", round(max(avg_scores), 4), "\n")
  }

  cat("\n>>> MAJORITY VOTE WINNER:", winner, "<<<\n")
  cat("  Votes:", max_votes, "out of", nrow(per_integration_winners), "integration methods\n")
  if (has_bio) {
    cat("  Scoring used: scIB-style composite (40% bio + 60% batch)\n")
  } else {
    cat("  Scoring used: batch metrics only\n")
  }
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  return(list(
    winner = winner,
    vote_table = vote_table,
    per_integration_winners = per_integration_winners,
    full_benchmark = full_benchmark_df,
    has_bio_metrics = has_bio
  ))
}


#' Compute per-label batch mixing metrics
#'
#' For each cell type label, subsets the embedding and computes batch mixing
#' metrics only among cells of that type. This reveals whether batches are
#' well-mixed within each cell population, not just globally.
#'
#' @param embedding Matrix of cell embeddings (cells x dims). Can be PCA,
#'   Harmony, RPCA, or any integration output.
#' @param batch_labels Character/factor vector of batch assignments per cell.
#' @param celltype_labels Character/factor vector of cell type labels per cell.
#' @param min_cells_per_label Integer. Skip labels with fewer cells than this
#'   threshold (metrics become unreliable with very few cells). Default: 30.
#' @param min_batches_per_label Integer. Skip labels present in fewer than this
#'   many batches (batch mixing is undefined for single-batch labels). Default: 2.
#'
#' @return data.frame with columns:
#'   - label: cell type name
#'   - n_cells: number of cells with this label
#'   - n_batches: number of batches containing this label
#'   - batch_variance: variance of batch centroids in embedding space
#'   - batch_asw: average silhouette width by batch (closer to 0 = better mixing)
#'   - lisi: Local Inverse Simpson's Index for batch (higher = more mixing)
#'   - batch_entropy: Shannon entropy of batch proportions (higher = more balanced)
#'   - flag: "OK", "WARNING", or "FAIL" based on thresholds
#'
compute_perlabel_batch_mixing <- function(embedding,
                                           batch_labels,
                                           celltype_labels,
                                           min_cells_per_label = 30,
                                           min_batches_per_label = 2) {

  # Input validation
  stopifnot(nrow(embedding) == length(batch_labels))
  stopifnot(nrow(embedding) == length(celltype_labels))

  # Get unique labels, excluding NA
  unique_labels <- sort(unique(na.omit(celltype_labels)))
  cat("  Computing per-label batch mixing for", length(unique_labels), "cell types\n")

  results <- list()

  for (label in unique_labels) {
    # Subset to cells of this type
    idx <- which(celltype_labels == label)
    n_cells <- length(idx)

    # Get batch composition for this label
    label_batches <- batch_labels[idx]
    batch_counts <- table(label_batches)
    n_batches <- sum(batch_counts > 0)

    # Skip if too few cells or batches
    if (n_cells < min_cells_per_label) {
      results[[label]] <- data.frame(
        label = label,
        n_cells = n_cells,
        n_batches = n_batches,
        batch_variance = NA_real_,
        batch_asw = NA_real_,
        lisi = NA_real_,
        batch_entropy = NA_real_,
        flag = "SKIPPED_FEW_CELLS",
        stringsAsFactors = FALSE
      )
      next
    }

    if (n_batches < min_batches_per_label) {
      results[[label]] <- data.frame(
        label = label,
        n_cells = n_cells,
        n_batches = n_batches,
        batch_variance = NA_real_,
        batch_asw = NA_real_,
        lisi = NA_real_,
        batch_entropy = NA_real_,
        flag = "SKIPPED_SINGLE_BATCH",
        stringsAsFactors = FALSE
      )
      next
    }

    # Subset embedding
    emb_sub <- embedding[idx, , drop = FALSE]

    # --- Metric 1: Batch variance (variance of batch centroids) ---
    # Lower = better (batches occupy same region)
    bv <- tryCatch({
      compute_batch_variance(emb_sub, label_batches)
    }, error = function(e) NA_real_)

    # --- Metric 2: Batch ASW (silhouette width by batch) ---
    # Closer to 0 = better (batches are intermingled)
    ba <- tryCatch({
      compute_batch_asw(emb_sub, label_batches)
    }, error = function(e) NA_real_)

    # --- Metric 3: LISI (Local Inverse Simpson's Index) ---
    # Higher = better (more batch mixing in local neighborhoods)
    li <- tryCatch({
      compute_lisi_score(emb_sub, label_batches)
    }, error = function(e) NA_real_)

    # --- Metric 4: Batch entropy (composition balance) ---
    # Higher = more balanced batch representation
    # Normalized to [0, 1] where 1 = perfectly balanced
    be <- tryCatch({
      props <- as.numeric(batch_counts) / sum(batch_counts)
      props <- props[props > 0]
      raw_entropy <- -sum(props * log(props))
      max_entropy <- log(length(props))
      if (max_entropy > 0) raw_entropy / max_entropy else NA_real_
    }, error = function(e) NA_real_)

    # --- Flag assessment ---
    # Thresholds are heuristic; adjust based on your dataset
    flag <- "OK"
    if (!is.na(ba) && abs(ba) > 0.15) flag <- "WARNING"
    if (!is.na(li) && li < 1.2) flag <- "WARNING"
    if (!is.na(ba) && abs(ba) > 0.3) flag <- "FAIL"
    if (!is.na(li) && li < 1.05) flag <- "FAIL"

    results[[label]] <- data.frame(
      label = label,
      n_cells = n_cells,
      n_batches = n_batches,
      batch_variance = bv,
      batch_asw = ba,
      lisi = li,
      batch_entropy = be,
      flag = flag,
      stringsAsFactors = FALSE
    )
  }

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  return(result_df)
}


#' Compute per-label batch mixing for multiple integration methods
#'
#' Runs compute_perlabel_batch_mixing across all integration methods
#' stored in a Seurat object's reductions.
#'
#' @param seurat_obj Seurat object with multiple integration reductions.
#' @param methods Character vector of method names to evaluate.
#' @param reduction_map Named list mapping method names to reduction names
#'   in the Seurat object. E.g., list(harmony = "harmony", rpca = "integrated.rpca").
#' @param batch_var Character. Metadata column name for batch.
#' @param celltype_var Character. Metadata column name for cell type labels.
#' @param dims Integer. Number of embedding dimensions to use.
#' @param min_cells_per_label Integer. Minimum cells per label. Default: 30.
#' @param min_batches_per_label Integer. Minimum batches per label. Default: 2.
#'
#' @return data.frame with all per-label scores across all methods,
#'   with an additional 'method' column.
#'
compute_perlabel_batch_mixing_all_methods <- function(seurat_obj,
                                                       methods,
                                                       reduction_map,
                                                       batch_var,
                                                       celltype_var,
                                                       dims,
                                                       min_cells_per_label = 30,
                                                       min_batches_per_label = 2) {

  batch_labels <- seurat_obj@meta.data[[batch_var]]
  celltype_labels <- seurat_obj@meta.data[[celltype_var]]

  all_results <- list()

  for (method in methods) {
    reduction_name <- reduction_map[[method]]
    if (is.null(reduction_name) || !reduction_name %in% names(seurat_obj@reductions)) {
      cat("  [SKIP]", method, "- reduction not found:", reduction_name, "\n")
      next
    }

    cat("\n  Evaluating:", method, "(reduction:", reduction_name, ")\n")

    # Extract embedding
    emb <- Seurat::Embeddings(seurat_obj, reduction = reduction_name)
    max_dims <- min(dims, ncol(emb))
    emb <- emb[, 1:max_dims, drop = FALSE]

    # Compute per-label metrics
    method_results <- compute_perlabel_batch_mixing(
      embedding = emb,
      batch_labels = batch_labels,
      celltype_labels = celltype_labels,
      min_cells_per_label = min_cells_per_label,
      min_batches_per_label = min_batches_per_label
    )
    method_results$method <- method

    all_results[[method]] <- method_results
  }

  combined <- do.call(rbind, all_results)
  rownames(combined) <- NULL
  return(combined)
}


#' Summarize per-label batch mixing across methods
#'
#' Creates a summary comparing methods on average per-label metrics,
#' number of flagged labels, and worst-performing cell types.
#'
#' @param perlabel_df data.frame from compute_perlabel_batch_mixing_all_methods
#'   or from manual aggregation with a 'method' column.
#' @return list with:
#'   - summary: per-method averages and flag counts
#'   - flagged_labels: cell types with WARNING or FAIL for any method
#'   - worst_labels: per-method worst-performing cell types
#'
summarize_perlabel_batch_mixing <- function(perlabel_df) {

  # Filter to evaluable labels only
  eval_df <- perlabel_df[!perlabel_df$flag %in% c("SKIPPED_FEW_CELLS", "SKIPPED_SINGLE_BATCH"), ]

  # Per-method summary
  method_summary <- eval_df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      n_labels_evaluated = dplyr::n(),
      mean_batch_variance = mean(batch_variance, na.rm = TRUE),
      mean_batch_asw = mean(batch_asw, na.rm = TRUE),
      mean_lisi = mean(lisi, na.rm = TRUE),
      mean_batch_entropy = mean(batch_entropy, na.rm = TRUE),
      n_ok = sum(flag == "OK"),
      n_warning = sum(flag == "WARNING"),
      n_fail = sum(flag == "FAIL"),
      pct_ok = round(100 * sum(flag == "OK") / dplyr::n(), 1),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_lisi))

  # Identify problematic labels (FAIL in any method)
  flagged <- eval_df %>%
    dplyr::filter(flag %in% c("WARNING", "FAIL")) %>%
    dplyr::select(label, method, n_cells, n_batches, batch_asw, lisi, flag) %>%
    dplyr::arrange(label, method)

  # Per-method worst labels (lowest LISI)
  worst_per_method <- eval_df %>%
    dplyr::group_by(method) %>%
    dplyr::slice_min(order_by = lisi, n = 3, with_ties = FALSE) %>%
    dplyr::select(method, label, n_cells, n_batches, lisi, batch_asw, flag) %>%
    dplyr::arrange(method, lisi)

  return(list(
    summary = method_summary,
    flagged_labels = flagged,
    worst_labels = worst_per_method
  ))
}


# ==============================================================================
# PER-LABEL BATCH MIXING VISUALIZATION FUNCTIONS (added 2026-02-09)
# ==============================================================================

#' Plot per-label batch mixing heatmap
#'
#' Creates a heatmap showing LISI scores per cell type and integration method.
#' Labels are colored by flag status (green=OK, yellow=WARNING, red=FAIL).
#'
#' @param perlabel_df data.frame from compute_perlabel_batch_mixing_all_methods
#'   or from manual aggregation with a 'method' column.
#' @param metric Character. Which metric to display. Default: "lisi".
#' @param title Character. Plot title.
#' @return ggplot object
#'
plot_perlabel_heatmap <- function(perlabel_df,
                                  metric = "lisi",
                                  title = "Per-Label Batch Mixing (LISI)") {

  # Filter to evaluable labels
  plot_df <- perlabel_df[!perlabel_df$flag %in% c("SKIPPED_FEW_CELLS", "SKIPPED_SINGLE_BATCH"), ]

  if (nrow(plot_df) == 0) {
    cat("  No evaluable labels for heatmap\n")
    return(NULL)
  }

  # Order labels by average metric across methods (best at top)
  label_order <- plot_df %>%
    dplyr::group_by(label) %>%
    dplyr::summarise(mean_val = mean(.data[[metric]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_val)) %>%
    dplyr::pull(label)

  plot_df$label <- factor(plot_df$label, levels = rev(label_order))

  # Create cell size annotation (for label text)
  plot_df$cell_label <- paste0(round(plot_df[[metric]], 2))

  # Determine color scale based on metric
  if (metric == "lisi") {
    # Higher = better for LISI
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = method, y = label, fill = .data[[metric]])) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = cell_label), size = 2.8) +
      ggplot2::scale_fill_gradient2(
        low = "#D73027", mid = "#FEE08B", high = "#1A9850",
        midpoint = median(plot_df[[metric]], na.rm = TRUE),
        name = "LISI\n(higher=\nbetter)"
      )
  } else if (metric == "batch_asw") {
    # Closer to 0 = better for batch ASW
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = method, y = label, fill = abs(.data[[metric]]))) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = cell_label), size = 2.8) +
      ggplot2::scale_fill_gradient(
        low = "#1A9850", high = "#D73027",
        name = "|batch ASW|\n(lower=\nbetter)"
      )
  } else {
    # Generic: lower = better (batch_variance)
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = method, y = label, fill = .data[[metric]])) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = cell_label), size = 2.8) +
      ggplot2::scale_fill_gradient(
        low = "#1A9850", high = "#D73027",
        name = paste0(metric, "\n(lower=\nbetter)")
      )
  }

  # Add cell count annotation on the right
  label_counts <- plot_df %>%
    dplyr::group_by(label) %>%
    dplyr::summarise(n_cells = dplyr::first(n_cells), .groups = "drop")

  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Per cell-type evaluation of batch mixing | ",
                        length(unique(plot_df$label)), " cell types evaluated"),
      x = "Integration Method",
      y = "Cell Type Label"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank()
    )

  return(p)
}


#' Plot per-label flag summary (stacked bar chart)
#'
#' Shows proportion of OK / WARNING / FAIL labels for each method.
#'
#' @param summary_df data.frame from summarize_perlabel_batch_mixing()$summary
#' @return ggplot object
#'
plot_perlabel_flag_summary <- function(summary_df) {

  flag_df <- summary_df %>%
    dplyr::select(method, n_ok, n_warning, n_fail) %>%
    tidyr::pivot_longer(
      cols = c(n_ok, n_warning, n_fail),
      names_to = "flag",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      flag = factor(flag,
                    levels = c("n_fail", "n_warning", "n_ok"),
                    labels = c("FAIL", "WARNING", "OK"))
    )

  # Order methods by pct_ok
  method_order <- summary_df %>% dplyr::arrange(dplyr::desc(pct_ok)) %>% dplyr::pull(method)
  flag_df$method <- factor(flag_df$method, levels = method_order)

  p <- ggplot2::ggplot(flag_df, ggplot2::aes(x = method, y = count, fill = flag)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(
      values = c("OK" = "#4DAF4A", "WARNING" = "#FF7F00", "FAIL" = "#E41A1C"),
      name = "Status"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Per-Label Batch Mixing Status by Integration Method",
      subtitle = "How many cell types pass within-label batch mixing checks",
      x = "Integration Method",
      y = "Number of Cell Type Labels"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::coord_flip()

  return(p)
}


#' Plot per-label LISI comparison for top N methods
#'
#' Dot plot comparing within-label LISI across methods, one row per cell type.
#'
#' @param perlabel_df data.frame from compute_perlabel_batch_mixing_all_methods
#'   or from manual aggregation with a 'method' column.
#' @param top_methods Character vector of method names to compare
#' @return ggplot object
#'
plot_perlabel_lisi_comparison <- function(perlabel_df, top_methods) {

  plot_df <- perlabel_df %>%
    dplyr::filter(
      method %in% top_methods,
      !flag %in% c("SKIPPED_FEW_CELLS", "SKIPPED_SINGLE_BATCH")
    )

  if (nrow(plot_df) == 0) return(NULL)

  # Order labels by cell count
  label_order <- plot_df %>%
    dplyr::group_by(label) %>%
    dplyr::summarise(n = dplyr::first(n_cells), .groups = "drop") %>%
    dplyr::arrange(n) %>%
    dplyr::pull(label)

  plot_df$label <- factor(plot_df$label, levels = label_order)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = lisi, y = label, color = method, size = n_cells)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "Cells") +
    ggplot2::scale_color_brewer(palette = "Set1", name = "Method") +
    ggplot2::geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey50", alpha = 0.5) +
    ggplot2::annotate("text", x = 1.0, y = 0.5, label = "No mixing", hjust = 1.1,
             size = 3, color = "grey50") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Within-Label Batch LISI by Cell Type",
      subtitle = paste0("Comparing: ", paste(top_methods, collapse = ", "),
                        " | Higher LISI = better batch mixing within cell type"),
      x = "Batch LISI (within label)",
      y = "Cell Type"
    )

  return(p)
}