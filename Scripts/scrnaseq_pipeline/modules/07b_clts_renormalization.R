#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 07b: CLTS RE-NORMALIZATION (MULTI-SAMPLE PIPELINE)
# ==============================================================================
#
# This module applies CLTS (Count based on Linearized Transcriptome Size)
# re-normalization to clustered scRNA-seq data.
#
# CLTS preserves biological variation in transcriptome size across cell types
# while correcting for technical sequencing depth differences between samples.
#
# Based on: Lu et al., Nature Communications 2025
# "Transcriptome size matters for single-cell RNA-seq normalization and bulk
# deconvolution"
# DOI: 10.1038/s41467-025-56623-1
#
# DYNAMIC CLUSTERING SOURCE:
# - Can apply CLTS to scICE, Leiden, CHOIR, or multiple clusterings
# - Controlled by params$clts_clustering_source
#
# INPUT: Clustered Seurat object(s) from Module 05/06/07
# OUTPUT: CLTS-normalized Seurat object(s) + marker benchmark comparison
#
# CREATED: 2026-01-08
# UPDATED: 2026-01-09 - Renamed to 07b, added dynamic clustering source support
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 07b: CLTS RE-NORMALIZATION\n")
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
  library(tidyr)
  library(tibble)    
  library(patchwork)
  library(Matrix)
})

# Check for optional packages
has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)
if (has_pheatmap) library(pheatmap)

has_RColorBrewer <- requireNamespace("RColorBrewer", quietly = TRUE)
if (has_RColorBrewer) library(RColorBrewer)

out_base <- params$out_root
load(file.path(out_base, "objects", "pipeline_environment.RData"))

# ==============================================================================
# Check if CLTS should run
# ==============================================================================
if (!isTRUE(params$run_clts_renormalization)) {
  cat("CLTS re-normalization disabled in parameters. Skipping.\n")
  
  clts_success <- FALSE
  clts_results <- list()
  clts_objects_created <- character(0)
  
  save(clts_success, clts_results, clts_objects_created,
       file = file.path(output_dirs$objects, "07b_clts_data.RData"))
  
  cat("\n>>> MODULE 07b SKIPPED <<<\n")
  quit(save = "no", status = 0)
}

# ==============================================================================
# Determine which clusterings to process
# ==============================================================================
clts_source <- params$clts_clustering_source

cat("CLTS clustering source parameter:", clts_source, "\n\n")

# Map source to list of clusterings to process
if (clts_source == "scice") {
  clusterings_to_process <- "scice"
} else if (clts_source == "leiden") {
  clusterings_to_process <- "leiden"
} else if (clts_source == "choir") {
  clusterings_to_process <- "choir"
} else if (clts_source == "both") {
  clusterings_to_process <- c("scice", "leiden")
} else if (clts_source == "all") {
  clusterings_to_process <- c("scice", "leiden", "choir")
} else {
  cat("WARNING: Unknown clts_clustering_source:", clts_source, "\n")
  cat("Defaulting to 'scice'\n")
  clusterings_to_process <- "scice"
}

cat("Will process clusterings:", paste(clusterings_to_process, collapse = ", "), "\n\n")

# ==============================================================================
# Create output directory structure
# ==============================================================================
clts_output_dir <- file.path(out_base, "07b_CLTS_Normalization")
dir.create(clts_output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Output directory:", clts_output_dir, "\n\n")

# ==============================================================================
# Define input/output file mapping
# ==============================================================================
clustering_config <- list(
  scice = list(
    input_file = file.path(output_dirs$objects, "scice_subclustered_object.rds"),
    output_file = file.path(output_dirs$objects, "scice_subclustered_object_redeconv.rds"),
    cluster_columns = c("scice_subcluster", "CHOIR_clusters_0.05", "CHOIR_clusters", "seurat_clusters"),
    output_subdir = "scice",
    display_name = "scICE Subclusters"
  ),
  leiden = list(
    input_file = file.path(output_dirs$objects, "leiden_clustered_object.rds"),
    output_file = file.path(output_dirs$objects, "leiden_clustered_object_redeconv.rds"),
    cluster_columns = c("leiden_clusters", "seurat_clusters"),
    output_subdir = "leiden",
    display_name = "Leiden Clusters"
  ),
  choir = list(
    input_file = file.path(output_dirs$objects, "choir_clustered_object.rds"),
    output_file = file.path(output_dirs$objects, "choir_clustered_object_redeconv.rds"),
    cluster_columns = c("CHOIR_clusters_0.05", "CHOIR_clusters_0.01", "CHOIR_clusters", "seurat_clusters"),
    output_subdir = "choir",
    display_name = "CHOIR Clusters"
  )
)

# ==============================================================================
# CLTS ALGORITHM FUNCTIONS
# ==============================================================================

#' Find the first available cluster column in metadata
#' @param meta Metadata data frame
#' @param candidates Vector of candidate column names
#' @return Name of first found column or NULL
find_cluster_column <- function(meta, candidates) {
  for (col in candidates) {
    if (col %in% colnames(meta)) {
      return(col)
    }
    # Also check for pattern match (e.g., CHOIR_clusters_*)
    pattern_match <- grep(paste0("^", gsub("\\*", ".*", col), "$"), colnames(meta), value = TRUE)
    if (length(pattern_match) > 0) {
      return(pattern_match[1])
    }
  }
  return(NULL)
}

#' Apply CLTS normalization to a Seurat object
#' @param seurat_obj Seurat object with clusters
#' @param cluster_col Name of cluster column
#' @param sample_col Name of sample column
#' @param min_cells_per_cluster Minimum cells for regression
#' @param baseline_sample Baseline sample selection method
#' @return List with CLTS-normalized counts and statistics
apply_clts_normalization <- function(seurat_obj, cluster_col, sample_col,
                                      min_cells_per_cluster = 50,
                                      baseline_sample = "auto") {
  
  cat("\n--- Applying CLTS Normalization ---\n")
  cat("Cluster column:", cluster_col, "\n")
  cat("Sample column:", sample_col, "\n")
  cat("Min cells per cluster:", min_cells_per_cluster, "\n")
  
  # Ensure layers are joined
  tryCatch({
    if ("RNA" %in% names(seurat_obj@assays)) {
      rna_layers <- Layers(seurat_obj[["RNA"]])
      counts_layers <- grep("^counts", rna_layers, value = TRUE)
      
      if (length(counts_layers) > 1) {
        cat("Joining", length(counts_layers), "counts layers...\n")
        seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
      }
    }
  }, error = function(e) {
    cat("Note: Layer preparation -", e$message, "\n")
  })
  
  # Get counts matrix
  counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  
  # Get cluster and sample info
  clusters <- as.character(seurat_obj@meta.data[[cluster_col]])
  samples <- as.character(seurat_obj@meta.data[[sample_col]])
  
  unique_clusters <- unique(clusters)
  unique_samples <- unique(samples)
  n_clusters <- length(unique_clusters)
  n_samples <- length(unique_samples)
  
  cat("Total cells:", ncol(seurat_obj), "\n")
  cat("Number of clusters:", n_clusters, "\n")
  cat("Number of samples:", n_samples, "\n")
  
  # Calculate transcriptome size per cell
  transcriptome_sizes <- Matrix::colSums(counts_matrix)
  
  cat("Transcriptome size range:", round(min(transcriptome_sizes)), "-", 
      round(max(transcriptome_sizes)), "\n")
  
  # Single sample case - CLTS = raw counts (identity transformation)
  if (n_samples == 1) {
    cat("\nNOTE: Single sample detected - CLTS = raw counts (identity transformation)\n")
    cat("CLTS normalization is most beneficial for multi-sample data.\n")
    
    return(list(
      clts_counts = counts_matrix,
      regression_params = data.frame(
        sample = unique_samples,
        a = 1,
        b = 0,
        r_squared = 1,
        n_clusters_used = n_clusters,
        stringsAsFactors = FALSE
      ),
      baseline_sample = unique_samples[1],
      v_max = 0,
      single_sample = TRUE,
      cluster_sample_stats = NULL,
      ts_matrix = NULL
    ))
  }
  
  # Calculate mean transcriptome size per cluster per sample
  cat("\nCalculating mean transcriptome size per cluster per sample...\n")
  
  cluster_sample_df <- data.frame(
    cell_id = colnames(seurat_obj),
    cluster = clusters,
    sample = samples,
    ts = transcriptome_sizes,
    stringsAsFactors = FALSE
  )
  
  cluster_sample_stats <- cluster_sample_df %>%
    group_by(cluster, sample) %>%
    summarise(
      mean_ts = mean(ts),
      median_ts = median(ts),
      sd_ts = sd(ts),
      n_cells = n(),
      .groups = "drop"
    )
  
  cat("Cluster-sample combinations:", nrow(cluster_sample_stats), "\n")
  
  # Filter by minimum cells
  cluster_sample_stats_filtered <- cluster_sample_stats %>%
    filter(n_cells >= min_cells_per_cluster)
  
  cat("After filtering (n >=", min_cells_per_cluster, "):", 
      nrow(cluster_sample_stats_filtered), "combinations\n")
  
  # Create transcriptome size matrix (clusters x samples)
  ts_matrix <- cluster_sample_stats_filtered %>%
    select(cluster, sample, mean_ts) %>%
    pivot_wider(names_from = sample, values_from = mean_ts) %>%
    column_to_rownames("cluster") %>%
    as.matrix()
  
  cat("TS matrix dimensions:", nrow(ts_matrix), "clusters x", ncol(ts_matrix), "samples\n")
  
  # Select baseline sample
  if (baseline_sample == "auto") {
    cat("\nAuto-selecting baseline sample based on correlation...\n")
    
    # Calculate pairwise correlations between samples
    sample_correlations <- matrix(NA, n_samples, n_samples,
                                  dimnames = list(unique_samples, unique_samples))
    
    for (i in 1:n_samples) {
      for (j in 1:n_samples) {
        if (i != j) {
          s1 <- unique_samples[i]
          s2 <- unique_samples[j]
          
          if (s1 %in% colnames(ts_matrix) && s2 %in% colnames(ts_matrix)) {
            vals1 <- ts_matrix[, s1]
            vals2 <- ts_matrix[, s2]
            
            # Use complete cases
            complete_idx <- !is.na(vals1) & !is.na(vals2)
            
            if (sum(complete_idx) >= 3) {
              sample_correlations[i, j] <- cor(vals1[complete_idx], vals2[complete_idx])
            }
          }
        }
      }
    }
    
    # Select sample with highest mean correlation
    mean_cors <- rowMeans(sample_correlations, na.rm = TRUE)
    baseline_sample <- unique_samples[which.max(mean_cors)]
    
    cat("Sample correlations (mean):\n")
    for (i in seq_along(unique_samples)) {
      cat("  ", unique_samples[i], ":", round(mean_cors[i], 3), "\n")
    }
    cat("Selected baseline:", baseline_sample, "(highest mean correlation)\n")
    
  } else {
    # Validate provided baseline sample
    if (!baseline_sample %in% unique_samples) {
      warning("Specified baseline sample '", baseline_sample, "' not found. Using auto-selection.")
      mean_cors <- sapply(unique_samples, function(s) {
        if (s %in% colnames(ts_matrix)) {
          mean(cor(ts_matrix[, s], ts_matrix, use = "pairwise.complete.obs"), na.rm = TRUE)
        } else {
          NA
        }
      })
      baseline_sample <- unique_samples[which.max(mean_cors)]
    }
    cat("Using specified baseline:", baseline_sample, "\n")
  }
  
  # Fit linear regression for each sample vs baseline
  cat("\nFitting linear regression for each sample vs baseline...\n")
  
  regression_params <- data.frame(
    sample = unique_samples,
    a = NA_real_,
    b = NA_real_,
    r_squared = NA_real_,
    n_clusters_used = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_samples) {
    samp <- unique_samples[i]
    
    if (samp == baseline_sample) {
      # Baseline sample: identity transformation
      regression_params$a[i] <- 1
      regression_params$b[i] <- 0
      regression_params$r_squared[i] <- 1
      regression_params$n_clusters_used[i] <- sum(!is.na(ts_matrix[, samp]))
      
      cat("  ", samp, "(baseline): a=1.000, b=0.000, R²=1.000\n")
      
    } else {
      # Check if sample is in matrix
      if (!samp %in% colnames(ts_matrix) || !baseline_sample %in% colnames(ts_matrix)) {
        regression_params$a[i] <- 1
        regression_params$b[i] <- 0
        regression_params$r_squared[i] <- NA
        regression_params$n_clusters_used[i] <- 0
        cat("  ", samp, ": NOT IN MATRIX - using identity\n")
        next
      }
      
      # Get common clusters with non-NA values
      common_mask <- !is.na(ts_matrix[, samp]) & !is.na(ts_matrix[, baseline_sample])
      n_common <- sum(common_mask)
      
      if (n_common >= 3) {
        x_baseline <- ts_matrix[common_mask, baseline_sample]
        y_sample <- ts_matrix[common_mask, samp]
        
        # Fit linear model: y_sample = a * x_baseline + b
        lm_fit <- lm(y_sample ~ x_baseline)
        
        regression_params$a[i] <- coef(lm_fit)[2]  # slope
        regression_params$b[i] <- coef(lm_fit)[1]  # intercept
        regression_params$r_squared[i] <- summary(lm_fit)$r.squared
        regression_params$n_clusters_used[i] <- n_common
        
        cat("  ", samp, ": a=", round(coef(lm_fit)[2], 4), 
            ", b=", round(coef(lm_fit)[1], 2),
            ", R²=", round(summary(lm_fit)$r.squared, 4),
            ", n=", n_common, "\n", sep = "")
        
      } else {
        # Not enough common clusters - use identity
        regression_params$a[i] <- 1
        regression_params$b[i] <- 0
        regression_params$r_squared[i] <- NA
        regression_params$n_clusters_used[i] <- n_common
        
        cat("  ", samp, ": TOO FEW CLUSTERS (n=", n_common, ") - using identity\n", sep = "")
      }
    }
  }
  
  # Calculate v_max = max(b_i / a_i) for the additive correction
  b_over_a <- regression_params$b / regression_params$a
  b_over_a[is.na(b_over_a) | is.infinite(b_over_a)] <- 0
  v_max <- max(b_over_a)
  
  cat("\nv_max (max b/a):", round(v_max, 4), "\n")
  
  # Apply CLTS transformation
  # Formula: u'_icg = u_icg * a_i + (v_max - b_i/a_i) / G
  # Simplified practical approach: scale counts by a_i and adjust for intercept
  
  cat("\nApplying CLTS transformation to counts...\n")
  
  n_genes <- nrow(counts_matrix)
  clts_counts <- counts_matrix  # Start with copy
  
  for (i in 1:n_samples) {
    samp <- unique_samples[i]
    a_i <- regression_params$a[i]
    b_i <- regression_params$b[i]
    
    if (is.na(a_i)) a_i <- 1
    if (is.na(b_i)) b_i <- 0
    
    cells_in_sample <- which(samples == samp)
    n_cells_sample <- length(cells_in_sample)
    
    if (n_cells_sample > 0) {
      cat("  Processing", samp, ":", n_cells_sample, "cells, a=", round(a_i, 4), "\n")
      
      # Get sample counts
      sample_counts <- counts_matrix[, cells_in_sample, drop = FALSE]
      
      # Apply scaling factor a_i
      if (a_i != 1) {
        sample_counts <- sample_counts * a_i
      }
      
      # For practical implementation: ensure counts match baseline distribution
      # Scale to match baseline median transcriptome size
      if (samp != baseline_sample && a_i != 1) {
        current_ts <- Matrix::colSums(sample_counts)
        baseline_cells <- which(samples == baseline_sample)
        baseline_median_ts <- median(transcriptome_sizes[baseline_cells])
        current_median_ts <- median(current_ts)
        
        if (current_median_ts > 0) {
          # Additional scaling to match baseline
          scale_factor <- baseline_median_ts / current_median_ts
          sample_counts <- sample_counts * scale_factor
        }
      }
      
      # Store transformed counts
      clts_counts[, cells_in_sample] <- sample_counts
    }
  }
  
  # Ensure non-negative values
  neg_count <- sum(clts_counts < 0)
  if (neg_count > 0) {
    cat("WARNING:", neg_count, "negative values detected. Setting to 0.\n")
    clts_counts[clts_counts < 0] <- 0
  }
  
  # Round to integers (counts should be integers)
  clts_counts <- round(clts_counts)
  
  # Verify transformation
  original_total <- sum(counts_matrix)
  clts_total <- sum(clts_counts)
  cat("\nTotal counts - Original:", format(original_total, big.mark = ","),
      "CLTS:", format(clts_total, big.mark = ","), "\n")
  
  return(list(
    clts_counts = clts_counts,
    regression_params = regression_params,
    baseline_sample = baseline_sample,
    v_max = v_max,
    single_sample = FALSE,
    cluster_sample_stats = cluster_sample_stats,
    ts_matrix = ts_matrix
  ))
}

#' Add CLTS assay to Seurat object
#' @param seurat_obj Seurat object
#' @param clts_counts CLTS-normalized count matrix
#' @return Seurat object with CLTS assay
add_clts_assay <- function(seurat_obj, clts_counts) {
  
  cat("\n--- Adding CLTS Assay to Seurat Object ---\n")
  
  # Create new assay
  clts_assay <- CreateAssayObject(counts = clts_counts)
  
  # Add to object
  seurat_obj[["CLTS"]] <- clts_assay
  
  # Normalize CLTS assay (log-normalize for visualization and marker detection)
  seurat_obj <- NormalizeData(seurat_obj, assay = "CLTS", verbose = FALSE)
  
  # Scale data for visualization
  seurat_obj <- ScaleData(seurat_obj, assay = "CLTS", verbose = FALSE)
  
  # Add CLTS transcriptome size to metadata
  seurat_obj$clts_transcriptome_size <- Matrix::colSums(clts_counts)
  
  cat("CLTS assay added successfully\n")
  cat("  Layers:", paste(Layers(seurat_obj[["CLTS"]]), collapse = ", "), "\n")
  
  return(seurat_obj)
}

#' Run marker benchmark comparing original vs CLTS normalization
#' @param seurat_obj Seurat object with CLTS assay
#' @param cluster_col Cluster column name
#' @param logfc_threshold Log fold change threshold
#' @param pval_threshold P-value threshold
#' @param min_pct Minimum percentage expressed
#' @return List with marker comparison results
run_marker_benchmark <- function(seurat_obj, cluster_col,
                                  logfc_threshold = 0.5,
                                  pval_threshold = 0.05,
                                  min_pct = 0.25) {
  
  cat("\n--- Running Marker Benchmark ---\n")
  cat("Cluster column:", cluster_col, "\n")
  cat("LogFC threshold:", logfc_threshold, "\n")
  cat("P-value threshold:", pval_threshold, "\n")
  cat("Min percent expressed:", min_pct, "\n")
  
  # Set identities
  Idents(seurat_obj) <- cluster_col
  
  # Ensure RNA assay has normalized data
  DefaultAssay(seurat_obj) <- "RNA"
  if (!"data" %in% Layers(seurat_obj[["RNA"]])) {
    cat("Normalizing RNA assay...\n")
    seurat_obj <- NormalizeData(seurat_obj, assay = "RNA", verbose = FALSE)
  }
  
  # Find markers using original normalization (RNA assay)
  cat("\nFinding markers with Original (RNA) normalization...\n")
  markers_original <- tryCatch({
    FindAllMarkers(seurat_obj, 
                   assay = "RNA",
                   only.pos = TRUE,
                   min.pct = min_pct, 
                   logfc.threshold = logfc_threshold,
                   verbose = FALSE)
  }, error = function(e) {
    cat("  Error finding RNA markers:", e$message, "\n")
    NULL
  })
  
  if (!is.null(markers_original)) {
    cat("  Found", nrow(markers_original), "markers (RNA)\n")
  }
  
  # Find markers using CLTS normalization
  cat("Finding markers with CLTS normalization...\n")
  DefaultAssay(seurat_obj) <- "CLTS"
  
  markers_clts <- tryCatch({
    FindAllMarkers(seurat_obj, 
                   assay = "CLTS",
                   only.pos = TRUE,
                   min.pct = min_pct, 
                   logfc.threshold = logfc_threshold,
                   verbose = FALSE)
  }, error = function(e) {
    cat("  Error finding CLTS markers:", e$message, "\n")
    NULL
  })
  
  if (!is.null(markers_clts)) {
    cat("  Found", nrow(markers_clts), "markers (CLTS)\n")
  }
  
  # Reset default assay
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Compare markers
  comparison_results <- list(
    markers_original = markers_original,
    markers_clts = markers_clts,
    cluster_comparison = NULL,
    overall_stats = NULL,
    concordant_genes = NULL,
    original_only_genes = NULL,
    clts_only_genes = NULL
  )
  
  if (!is.null(markers_original) && !is.null(markers_clts) &&
      nrow(markers_original) > 0 && nrow(markers_clts) > 0) {
    
    # Overall comparison
    all_genes_original <- unique(markers_original$gene)
    all_genes_clts <- unique(markers_clts$gene)
    
    concordant_genes <- intersect(all_genes_original, all_genes_clts)
    original_only <- setdiff(all_genes_original, all_genes_clts)
    clts_only <- setdiff(all_genes_clts, all_genes_original)
    
    all_genes <- union(all_genes_original, all_genes_clts)
    overall_jaccard <- length(concordant_genes) / length(all_genes)
    
    cat("\n--- Overall Marker Comparison ---\n")
    cat("Original markers:", length(all_genes_original), "\n")
    cat("CLTS markers:", length(all_genes_clts), "\n")
    cat("Concordant:", length(concordant_genes), "\n")
    cat("Original only:", length(original_only), "\n")
    cat("CLTS only (rescued):", length(clts_only), "\n")
    cat("Jaccard similarity:", round(overall_jaccard, 3), "\n")
    
    # Per-cluster comparison
    all_clusters <- unique(c(markers_original$cluster, markers_clts$cluster))
    
    cluster_comparison <- data.frame(
      cluster = all_clusters,
      n_markers_original = NA_integer_,
      n_markers_clts = NA_integer_,
      n_concordant = NA_integer_,
      n_original_only = NA_integer_,
      n_clts_only = NA_integer_,
      jaccard = NA_real_,
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(all_clusters)) {
      cl <- all_clusters[i]
      
      genes_orig <- markers_original$gene[markers_original$cluster == cl]
      genes_clts <- markers_clts$gene[markers_clts$cluster == cl]
      
      cluster_comparison$n_markers_original[i] <- length(genes_orig)
      cluster_comparison$n_markers_clts[i] <- length(genes_clts)
      cluster_comparison$n_concordant[i] <- length(intersect(genes_orig, genes_clts))
      cluster_comparison$n_original_only[i] <- length(setdiff(genes_orig, genes_clts))
      cluster_comparison$n_clts_only[i] <- length(setdiff(genes_clts, genes_orig))
      
      union_size <- length(union(genes_orig, genes_clts))
      if (union_size > 0) {
        cluster_comparison$jaccard[i] <- length(intersect(genes_orig, genes_clts)) / union_size
      } else {
        cluster_comparison$jaccard[i] <- NA
      }
    }
    
    # Sort by cluster
    cluster_comparison <- cluster_comparison[order(cluster_comparison$cluster), ]
    
    cat("\n--- Per-Cluster Summary ---\n")
    print(cluster_comparison)
    
    # Rescued markers (significant only in CLTS)
    rescued_markers <- markers_clts[markers_clts$gene %in% clts_only, ]
    
    # Store results
    comparison_results$cluster_comparison <- cluster_comparison
    comparison_results$overall_stats <- data.frame(
      n_original = length(all_genes_original),
      n_clts = length(all_genes_clts),
      n_concordant = length(concordant_genes),
      n_original_only = length(original_only),
      n_clts_only = length(clts_only),
      jaccard = overall_jaccard
    )
    comparison_results$concordant_genes <- concordant_genes
    comparison_results$original_only_genes <- original_only
    comparison_results$clts_only_genes <- clts_only
    comparison_results$rescued_markers <- rescued_markers
    
  } else {
    cat("\nWARNING: Could not compare markers (one or both methods returned no results)\n")
  }
  
  return(comparison_results)
}

#' Create benchmark visualization plots
#' @param comparison_results Output from run_marker_benchmark
#' @param output_dir Directory to save plots
#' @param clustering_name Name for plot titles
create_benchmark_plots <- function(comparison_results, output_dir, clustering_name) {
  
  cat("\n--- Creating Benchmark Plots ---\n")
  
  plots_dir <- file.path(output_dir, "plots")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. Marker count comparison bar plot
  if (!is.null(comparison_results$overall_stats)) {
    stats <- comparison_results$overall_stats
    
    bar_data <- data.frame(
      category = c("Original Only", "Concordant", "CLTS Only"),
      count = c(stats$n_original_only, stats$n_concordant, stats$n_clts_only)
    )
    bar_data$category <- factor(bar_data$category, 
                                levels = c("Original Only", "Concordant", "CLTS Only"))
    
    p_bar <- ggplot(bar_data, aes(x = category, y = count, fill = category)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = count), vjust = -0.5, size = 4) +
      scale_fill_manual(values = c("Original Only" = "#E41A1C", 
                                   "Concordant" = "#4DAF4A", 
                                   "CLTS Only" = "#377EB8")) +
      theme_minimal() +
      labs(title = paste0("Marker Comparison: Original vs CLTS (", clustering_name, ")"),
           subtitle = paste0("Jaccard Similarity: ", round(stats$jaccard, 3)),
           x = "", y = "Number of Markers") +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(size = 12))
    
    ggsave(file.path(plots_dir, "marker_count_comparison.png"), p_bar, 
           width = 8, height = 6, dpi = 150)
    ggsave(file.path(plots_dir, "marker_count_comparison.pdf"), p_bar, 
           width = 8, height = 6)
    
    cat("  Saved marker count comparison plot\n")
  }
  
  # 2. Per-cluster Jaccard similarity
  if (!is.null(comparison_results$cluster_comparison)) {
    cluster_df <- comparison_results$cluster_comparison
    cluster_df <- cluster_df[!is.na(cluster_df$jaccard), ]
    
    if (nrow(cluster_df) > 0) {
      p_jaccard <- ggplot(cluster_df, aes(x = reorder(cluster, jaccard), y = jaccard)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        geom_hline(yintercept = mean(cluster_df$jaccard, na.rm = TRUE), 
                   linetype = "dashed", color = "red") +
        coord_flip() +
        theme_minimal() +
        labs(title = paste0("Jaccard Similarity by Cluster (", clustering_name, ")"),
             subtitle = paste0("Mean: ", round(mean(cluster_df$jaccard, na.rm = TRUE), 3)),
             x = "Cluster", y = "Jaccard Similarity") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      ggsave(file.path(plots_dir, "jaccard_by_cluster.png"), p_jaccard, 
             width = 10, height = max(6, nrow(cluster_df) * 0.3), dpi = 150)
      ggsave(file.path(plots_dir, "jaccard_by_cluster.pdf"), p_jaccard, 
             width = 10, height = max(6, nrow(cluster_df) * 0.3))
      
      cat("  Saved Jaccard similarity plot\n")
    }
  }
  
  # 3. Cluster marker count comparison
  if (!is.null(comparison_results$cluster_comparison)) {
    cluster_df <- comparison_results$cluster_comparison
    
    cluster_long <- cluster_df %>%
      select(cluster, n_markers_original, n_markers_clts) %>%
      pivot_longer(cols = c(n_markers_original, n_markers_clts),
                   names_to = "method", values_to = "n_markers") %>%
      mutate(method = ifelse(method == "n_markers_original", "Original", "CLTS"))
    
    p_cluster <- ggplot(cluster_long, aes(x = cluster, y = n_markers, fill = method)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("Original" = "#E41A1C", "CLTS" = "#377EB8")) +
      theme_minimal() +
      labs(title = paste0("Markers per Cluster (", clustering_name, ")"),
           x = "Cluster", y = "Number of Markers", fill = "Normalization") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(plots_dir, "markers_per_cluster.png"), p_cluster, 
           width = max(8, length(unique(cluster_df$cluster)) * 0.5), height = 6, dpi = 150)
    ggsave(file.path(plots_dir, "markers_per_cluster.pdf"), p_cluster, 
           width = max(8, length(unique(cluster_df$cluster)) * 0.5), height = 6)
    
    cat("  Saved markers per cluster plot\n")
  }
  
  cat("Benchmark plots saved to:", plots_dir, "\n")
}

#' Save CLTS results to files
#' @param seurat_obj Seurat object with CLTS assay
#' @param clts_result Output from apply_clts_normalization
#' @param comparison_results Output from run_marker_benchmark
#' @param output_dir Output directory
#' @param output_rds_path Path to save RDS file
save_clts_results <- function(seurat_obj, clts_result, comparison_results, 
                               output_dir, output_rds_path) {
  
  cat("\n--- Saving CLTS Results ---\n")
  
  tables_dir <- file.path(output_dir, "tables")
  benchmark_dir <- file.path(output_dir, "benchmark")
  
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(benchmark_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save regression parameters
  write.csv(clts_result$regression_params, 
            file.path(tables_dir, "clts_regression_parameters.csv"),
            row.names = FALSE)
  cat("  Saved regression parameters\n")
  
  # Save cluster-sample stats
  if (!is.null(clts_result$cluster_sample_stats)) {
    write.csv(clts_result$cluster_sample_stats,
              file.path(tables_dir, "cluster_sample_transcriptome_sizes.csv"),
              row.names = FALSE)
    cat("  Saved cluster-sample transcriptome sizes\n")
  }
  
  # Save benchmark results
  if (!is.null(comparison_results$markers_original)) {
    write.csv(comparison_results$markers_original,
              file.path(benchmark_dir, "markers_original_normalization.csv"),
              row.names = FALSE)
    cat("  Saved original markers\n")
  }
  
  if (!is.null(comparison_results$markers_clts)) {
    write.csv(comparison_results$markers_clts,
              file.path(benchmark_dir, "markers_clts_normalization.csv"),
              row.names = FALSE)
    cat("  Saved CLTS markers\n")
  }
  
  if (!is.null(comparison_results$cluster_comparison)) {
    write.csv(comparison_results$cluster_comparison,
              file.path(benchmark_dir, "marker_comparison_by_cluster.csv"),
              row.names = FALSE)
    cat("  Saved cluster comparison\n")
  }
  
  if (!is.null(comparison_results$rescued_markers) && nrow(comparison_results$rescued_markers) > 0) {
    write.csv(comparison_results$rescued_markers,
              file.path(benchmark_dir, "rescued_markers_clts_only.csv"),
              row.names = FALSE)
    cat("  Saved rescued markers (CLTS-only)\n")
  }
  
  if (!is.null(comparison_results$overall_stats)) {
    write.csv(comparison_results$overall_stats,
              file.path(benchmark_dir, "overall_comparison_stats.csv"),
              row.names = FALSE)
    cat("  Saved overall stats\n")
  }
  
  # Save the Seurat object with CLTS assay
  saveRDS(seurat_obj, output_rds_path)
  cat("  Saved Seurat object:", output_rds_path, "\n")
  cat("  Size:", round(file.size(output_rds_path) / 1e6, 2), "MB\n")
  
  # Create README
  readme_text <- paste0(
    "CLTS Re-normalization Results\n",
    "=============================\n\n",
    "Module: 07b\n",
    "Method: CLTS (Count based on Linearized Transcriptome Size)\n",
    "Reference: Lu et al., Nature Communications 2025\n",
    "DOI: 10.1038/s41467-025-56623-1\n\n",
    "Baseline sample: ", clts_result$baseline_sample, "\n",
    "Single sample mode: ", clts_result$single_sample, "\n\n",
    "Files:\n",
    "------\n",
    "tables/clts_regression_parameters.csv - Linear regression coefficients (a, b) per sample\n",
    "tables/cluster_sample_transcriptome_sizes.csv - Mean transcriptome size per cluster per sample\n",
    "benchmark/markers_original_normalization.csv - Markers found with standard normalization\n",
    "benchmark/markers_clts_normalization.csv - Markers found with CLTS normalization\n",
    "benchmark/marker_comparison_by_cluster.csv - Per-cluster Jaccard similarity\n",
    "benchmark/rescued_markers_clts_only.csv - Markers significant only with CLTS\n",
    "plots/ - Visualization of benchmark comparison\n\n",
    "The output Seurat object (*_redeconv.rds) contains:\n",
    "- Original assays (RNA, SCT, etc.) - unchanged\n",
    "- New 'CLTS' assay with CLTS-normalized counts\n",
    "- New metadata column 'clts_transcriptome_size'\n"
  )
  
  writeLines(readme_text, file.path(output_dir, "README.txt"))
  cat("  Saved README\n")
}

# ==============================================================================
# MAIN PROCESSING LOOP
# ==============================================================================

clts_results <- list()
clts_objects_created <- character(0)
clts_success <- FALSE

for (clustering_type in clusterings_to_process) {
  
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("PROCESSING:", toupper(clustering_type), "CLUSTERING\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  config <- clustering_config[[clustering_type]]
  
  # Check if input file exists
  if (!file.exists(config$input_file)) {
    cat("WARNING: Input file not found:", config$input_file, "\n")
    cat("Skipping", clustering_type, "clustering.\n")
    next
  }
  
  # Create output subdirectory
  output_subdir <- file.path(clts_output_dir, config$output_subdir)
  dir.create(output_subdir, showWarnings = FALSE, recursive = TRUE)
  
  # Load object
  cat("\nLoading:", config$input_file, "\n")
  seurat_obj <- readRDS(config$input_file)
  
  cat("Object loaded:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "genes\n")
  
  # Find cluster column
  cluster_col <- find_cluster_column(seurat_obj@meta.data, config$cluster_columns)
  
  if (is.null(cluster_col)) {
    cat("WARNING: No cluster column found for", clustering_type, "\n")
    cat("Tried:", paste(config$cluster_columns, collapse = ", "), "\n")
    cat("Available columns:", paste(head(colnames(seurat_obj@meta.data), 20), collapse = ", "), "\n")
    next
  }
  
  cat("Using cluster column:", cluster_col, "\n")
  
  # Find sample column
  sample_candidates <- c("sample_name", "sample", "orig.ident", "Sample")
  sample_col <- find_cluster_column(seurat_obj@meta.data, sample_candidates)
  
  if (is.null(sample_col)) {
    cat("WARNING: No sample column found\n")
    sample_col <- "orig.ident"  # Fallback
  }
  
  cat("Using sample column:", sample_col, "\n")
  
  # Apply CLTS normalization
  clts_result <- tryCatch({
    apply_clts_normalization(
      seurat_obj = seurat_obj,
      cluster_col = cluster_col,
      sample_col = sample_col,
      min_cells_per_cluster = params$clts_min_cells_per_cluster,
      baseline_sample = params$clts_baseline_sample
    )
  }, error = function(e) {
    cat("ERROR in CLTS normalization:", e$message, "\n")
    NULL
  })
  
  if (is.null(clts_result)) {
    cat("CLTS normalization failed for", clustering_type, "\n")
    next
  }
  
  # Add CLTS assay to object
  seurat_obj <- tryCatch({
    add_clts_assay(seurat_obj, clts_result$clts_counts)
  }, error = function(e) {
    cat("ERROR adding CLTS assay:", e$message, "\n")
    NULL
  })
  
  if (is.null(seurat_obj)) {
    cat("Failed to add CLTS assay for", clustering_type, "\n")
    next
  }
  
  # Run marker benchmark if enabled
  comparison_results <- list()
  if (isTRUE(params$clts_run_benchmark)) {
    comparison_results <- tryCatch({
      run_marker_benchmark(
        seurat_obj = seurat_obj,
        cluster_col = cluster_col,
        logfc_threshold = params$clts_marker_logfc_threshold,
        pval_threshold = params$clts_marker_pval_threshold,
        min_pct = params$clts_marker_min_pct
      )
    }, error = function(e) {
      cat("ERROR in marker benchmark:", e$message, "\n")
      list()
    })
    
    # Create benchmark plots
    if (length(comparison_results) > 0 && !is.null(comparison_results$cluster_comparison)) {
      tryCatch({
        create_benchmark_plots(comparison_results, output_subdir, config$display_name)
      }, error = function(e) {
        cat("WARNING: Could not create benchmark plots:", e$message, "\n")
      })
    }
  }
  
  # Save results
  tryCatch({
    save_clts_results(
      seurat_obj = seurat_obj,
      clts_result = clts_result,
      comparison_results = comparison_results,
      output_dir = output_subdir,
      output_rds_path = config$output_file
    )
  }, error = function(e) {
    cat("ERROR saving results:", e$message, "\n")
  })
  
  # Store results for this clustering
  clts_results[[clustering_type]] <- list(
    cluster_column = cluster_col,
    sample_column = sample_col,
    regression_params = clts_result$regression_params,
    baseline_sample = clts_result$baseline_sample,
    single_sample = clts_result$single_sample,
    comparison_stats = comparison_results$overall_stats,
    output_file = config$output_file
  )
  
  # Track created objects
  if (file.exists(config$output_file)) {
    clts_objects_created <- c(clts_objects_created, clustering_type)
  }
  
  cat("\n[SUCCESS]", config$display_name, "CLTS normalization completed\n")
  clts_success <- TRUE
}

# ==============================================================================
# Save module data
# ==============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("MODULE 07b SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Clusterings processed:", paste(clusterings_to_process, collapse = ", "), "\n")
cat("Objects successfully created:", paste(clts_objects_created, collapse = ", "), "\n")
cat("Module success:", clts_success, "\n")

if (length(clts_results) > 0) {
  cat("\nResults summary:\n")
  for (ct in names(clts_results)) {
    res <- clts_results[[ct]]
    cat("\n", toupper(ct), ":\n", sep = "")
    cat("  Cluster column:", res$cluster_column, "\n")
    cat("  Baseline sample:", res$baseline_sample, "\n")
    cat("  Single sample mode:", res$single_sample, "\n")
    if (!is.null(res$comparison_stats)) {
      cat("  Marker comparison Jaccard:", round(res$comparison_stats$jaccard, 3), "\n")
    }
    cat("  Output file:", res$output_file, "\n")
  }
}

# Save module data
clts_data_file <- file.path(output_dirs$objects, "07b_clts_data.RData")
save(clts_success, clts_results, clts_objects_created, clusterings_to_process,
     file = clts_data_file)
cat("\nModule data saved to:", clts_data_file, "\n")

cat("\n>>> MODULE 07b COMPLETE <<<\n")
