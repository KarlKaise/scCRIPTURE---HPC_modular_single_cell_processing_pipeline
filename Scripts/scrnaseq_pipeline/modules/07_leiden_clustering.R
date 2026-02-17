#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 07: LEIDEN CLUSTERING (MULTI-SAMPLE PIPELINE)
# ==============================================================================
#
# This module performs Leiden clustering with resolution testing.
# Can serve as:
# - Primary clustering method (if CHOIR failed/disabled)
# - Comparison method (if CHOIR succeeded)
#
# Features:
# - Multi-resolution testing
# - Clustree visualization
# - Quality metrics (silhouette, purity)
# - Optimal resolution selection
#
# INPUT: Integrated/normalized Seurat object
# OUTPUT: Leiden-clustered Seurat object
#
# UPDATES:
# - 2026-01-03: Added lognorm_object fallback in object selection
# - 2026-01-03: Improved JoinLayers handling for Seurat v5
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 07: LEIDEN CLUSTERING\n")
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
  library(reshape2)
  library(cluster)
  library(tidyr)
  library(Matrix)
})

out_base <- params$out_root
load(file.path(out_base, "objects", "pipeline_environment.RData"))

# Load CHOIR data (to check if it succeeded)
choir_file <- file.path(out_base, "objects", "05_choir_data.RData")
if (file.exists(choir_file)) {
  load(choir_file)
} else {
  choir_success <- FALSE
  clustering_method_used <- "none"
  clustered_obj <- NULL
}

# Load integration data
integration_file <- file.path(out_base, "objects", "04_integration_data.RData")
norm_file <- file.path(out_base, "objects", "03_normalization_data.RData")

if (file.exists(integration_file)) {
  load(integration_file)
  cat("Loaded integration data.\n")
} else if (file.exists(norm_file)) {
  load(norm_file)
  multi_integrated <- NULL
  best_method <- "none"
  best_reduction <- "pca"
  cat("Loaded normalization data.\n")
}

# ==============================================================================
# Check if Leiden should run
# ==============================================================================
if (!isTRUE(params$run_leiden_clustering)) {
  cat("Leiden clustering disabled in parameters.\n")

  if (choir_success) {
    cat("Using CHOIR clustering as final result.\n")
  } else {
    cat("WARNING: No clustering will be performed!\n")
  }

  leiden_file <- file.path(output_dirs$objects, "07_leiden_data.RData")
  save(choir_success, clustering_method_used, clustered_obj, file = leiden_file)

  cat("\n>>> MODULE 07 SKIPPED <<<\n")
  quit(save = "no", status = 0)
}

# ==============================================================================
# Determine purpose of Leiden clustering
# ==============================================================================
if (choir_success) {
  leiden_purpose <- "comparison"
  cat("Running Leiden as COMPARISON to CHOIR clustering.\n")
} else {
  leiden_purpose <- "primary"
  cat("Running Leiden as PRIMARY clustering method.\n")
}

# ==============================================================================
# Select input object
# ==============================================================================
cat("\n--- Selecting input object ---\n")

if (exists("multi_integrated") && !is.null(multi_integrated)) {
  leiden_input <- multi_integrated
  default_reduction <- best_reduction
  cat("Using integrated object with reduction:", default_reduction, "\n")
} else if (exists("sct_object") && !is.null(sct_object)) {
  leiden_input <- sct_object
  default_reduction <- if ("harmony" %in% names(sct_object@reductions)) "harmony" else "pca"
  cat("Using SCTransform object with reduction:", default_reduction, "\n")
} else if (exists("scran_object") && !is.null(scran_object)) {
  leiden_input <- scran_object
  default_reduction <- if ("harmony" %in% names(scran_object@reductions)) "harmony" else "pca"
  cat("Using scran object with reduction:", default_reduction, "\n")
} else if (exists("lognorm_object") && !is.null(lognorm_object)) {
  leiden_input <- lognorm_object
  default_reduction <- if ("harmony" %in% names(lognorm_object@reductions)) "harmony" else "pca"
  cat("Using LogNormalize object with reduction:", default_reduction, "\n")
} else {
  stop("No suitable object available for Leiden clustering")
}

# ==============================================================================
# Prepare object (Seurat v5 compatibility)
# ==============================================================================
cat("\n--- Preparing object ---\n")

# Join RNA layers if they are split
tryCatch({
  if ("RNA" %in% names(leiden_input@assays)) {
    rna_layers <- Layers(leiden_input[["RNA"]])
    counts_layers <- grep("^counts", rna_layers, value = TRUE)
    
    if (length(counts_layers) > 1) {
      cat("Found", length(counts_layers), "counts layers - joining...\n")
      leiden_input[["RNA"]] <- JoinLayers(leiden_input[["RNA"]])
      cat("Successfully joined RNA layers\n")
    } else {
      cat("RNA layers already unified\n")
    }
  }
}, error = function(e) {
  cat("Note: Layer preparation message -", e$message, "\n")
})

# Verify reduction exists
if (!default_reduction %in% names(leiden_input@reductions)) {
  cat("WARNING: Reduction", default_reduction, "not found. Falling back to pca.\n")
  default_reduction <- "pca"
}

if (!default_reduction %in% names(leiden_input@reductions)) {
  stop("No valid reduction found in object. Available: ", 
       paste(names(leiden_input@reductions), collapse = ", "))
}

nd <- ncol(Embeddings(leiden_input, reduction = default_reduction))
dims_to_use <- min(params$dims_use, nd)
cat("Using dims 1:", dims_to_use, "from", default_reduction, "\n")

print_object_structure(leiden_input, "Leiden Input")

# ==============================================================================
# Find Neighbors
# ==============================================================================
cat("\n--- Finding neighbors ---\n")

leiden_n_neighbors <- if (!is.null(params$leiden_n_neighbors)) params$leiden_n_neighbors else 20

leiden_obj <- FindNeighbors(leiden_input,
                            reduction = default_reduction,
                            dims = 1:dims_to_use,
                            k.param = leiden_n_neighbors,
                            verbose = FALSE)

cat("Neighbors computed with k =", leiden_n_neighbors, "\n")

# ==============================================================================
# Test multiple resolutions
# ==============================================================================
cat("\n--- Testing resolutions ---\n")

# Set defaults if not specified
leiden_resolutions <- if (!is.null(params$leiden_resolutions)) {
  params$leiden_resolutions
} else {
  c(0.3, 0.5, 0.8, 1.0, 1.2)
}

leiden_algorithm <- if (!is.null(params$leiden_algorithm)) params$leiden_algorithm else 4

resolution_results <- list()
quality_metrics <- list()

for (res in leiden_resolutions) {
  cat("  Resolution", res, "...")

  leiden_obj <- FindClusters(leiden_obj,
                             resolution = res,
                             algorithm = leiden_algorithm,
                             method = "igraph",
                             random.seed = 42,
                             verbose = FALSE)

  col_name <- paste0("leiden_res_", res)
  leiden_obj@meta.data[[col_name]] <- leiden_obj$seurat_clusters

  clusters <- leiden_obj@meta.data[[col_name]]
  n_clust <- length(unique(clusters))
  resolution_results[[as.character(res)]] <- n_clust

  emb <- Embeddings(leiden_obj, reduction = default_reduction)[, 1:dims_to_use]

  mean_sil <- NA
  pct_neg <- NA
  tryCatch({
    if (nrow(emb) > 5000) {
      set.seed(42)
      idx <- sample(nrow(emb), 5000)
      emb_sub <- emb[idx, ]
      clusters_sub <- clusters[idx]
    } else {
      emb_sub <- emb
      clusters_sub <- clusters
    }

    if (length(unique(clusters_sub)) > 1 && length(unique(clusters_sub)) < nrow(emb_sub)) {
      d <- dist(emb_sub)
      sil <- cluster::silhouette(as.numeric(factor(clusters_sub)), d)
      mean_sil <- mean(sil[, 3])
      pct_neg <- sum(sil[, 3] < 0) / length(sil[, 3]) * 100
    }
  }, error = function(e) NULL)

  mean_purity <- NA
  if (exists("has_bluster") && has_bluster) {
    tryCatch({
      library(bluster)
      purity_result <- neighborPurity(emb, clusters)
      mean_purity <- mean(purity_result$purity)
    }, error = function(e) NULL)
  }

  quality_metrics[[as.character(res)]] <- data.frame(
    resolution = res,
    n_clusters = n_clust,
    mean_silhouette = mean_sil,
    pct_negative_sil = pct_neg,
    mean_purity = mean_purity
  )

  cat(" ", n_clust, "clusters (silhouette:", round(mean_sil, 3), ")\n")
}

quality_df <- do.call(rbind, quality_metrics)
rownames(quality_df) <- NULL

res_summary <- data.frame(
  resolution = as.numeric(names(resolution_results)),
  n_clusters = unlist(resolution_results)
)

write.csv(res_summary,
          file.path(subdirs$resolution_testing, "resolution_summary.csv"),
          row.names = FALSE)

write.csv(quality_df,
          file.path(subdirs$quality_metrics, "clustering_quality.csv"),
          row.names = FALSE)

cat("\n>>> LEIDEN QUALITY SUMMARY <<<\n")
print(quality_df)

# ==============================================================================
# Create visualizations
# ==============================================================================
cat("\n--- Creating visualizations ---\n")

if (exists("has_clustree") && has_clustree) {
  library(clustree)
  tryCatch({
    p_clustree <- clustree(leiden_obj, prefix = "leiden_res_") +
      theme(legend.position = "bottom")
    save_plot_multi(p_clustree, "01_clustree_resolution",
                    output_dir = subdirs$plots_leiden, width = 12, height = 10)
    cat("  Created clustree plot\n")
  }, error = function(e) {
    cat("  Clustree plot failed:", e$message, "\n")
  })
}

# Set final resolution with default
final_resolution <- if (!is.null(params$final_resolution)) params$final_resolution else 0.8

p_res <- ggplot(res_summary, aes(x = resolution, y = n_clusters)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(xintercept = final_resolution, linetype = "dashed", color = "red") +
  annotate("text", x = final_resolution, y = max(res_summary$n_clusters),
           label = paste0("Final: ", final_resolution),
           hjust = -0.1, vjust = 1, color = "red") +
  theme_minimal() +
  labs(title = "Number of Clusters by Resolution",
       subtitle = paste0("Red line: selected resolution (", final_resolution, ")"),
       x = "Resolution", y = "Number of Clusters")

save_plot_multi(p_res, "02_resolution_vs_clusters",
                output_dir = subdirs$plots_leiden, width = 8, height = 6)

quality_long <- quality_df %>%
  select(resolution, mean_silhouette, mean_purity) %>%
  pivot_longer(-resolution, names_to = "metric", values_to = "value") %>%
  filter(!is.na(value))

if (nrow(quality_long) > 0) {
  p_quality <- ggplot(quality_long, aes(x = resolution, y = value, color = metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_vline(xintercept = final_resolution, linetype = "dashed", color = "grey50") +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    theme_minimal() +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Clustering Quality Across Resolutions",
         x = "Resolution", y = "Score")

  save_plot_multi(p_quality, "03_quality_metrics",
                  output_dir = subdirs$plots_leiden, width = 10, height = 12)
}

# ==============================================================================
# Set final clustering
# ==============================================================================
cat("\n--- Setting final clustering ---\n")

final_col <- paste0("leiden_res_", final_resolution)

# Verify the column exists
if (!final_col %in% colnames(leiden_obj@meta.data)) {
  cat("WARNING: Column", final_col, "not found. Using first available resolution.\n")
  available_cols <- grep("^leiden_res_", colnames(leiden_obj@meta.data), value = TRUE)
  if (length(available_cols) > 0) {
    final_col <- available_cols[1]
    cat("Using:", final_col, "\n")
  } else {
    stop("No leiden resolution columns found")
  }
}

leiden_obj$leiden_clusters <- leiden_obj@meta.data[[final_col]]
n_leiden_clusters <- length(unique(leiden_obj$leiden_clusters))

cat("Final resolution:", final_resolution, "\n")
cat("Final number of clusters:", n_leiden_clusters, "\n")

if (leiden_purpose == "primary" || !exists("clustered_obj") || is.null(clustered_obj)) {
  leiden_obj$seurat_clusters <- leiden_obj$leiden_clusters
  Idents(leiden_obj) <- "seurat_clusters"
  clustered_obj <- leiden_obj
  clustering_method_used <- "Leiden"
  cat("\nLeiden set as PRIMARY clustering method\n")
} else {
  cat("\nLeiden kept as COMPARISON (CHOIR is primary)\n")
}

# ==============================================================================
# Ensure UMAP exists
# ==============================================================================
if (!"umap" %in% names(leiden_obj@reductions)) {
  cat("\nRunning UMAP for visualization...\n")
  leiden_obj <- RunUMAP(leiden_obj, reduction = default_reduction,
                        dims = 1:dims_to_use, verbose = FALSE)
}

# ==============================================================================
# Create final UMAP plots
# ==============================================================================
p_leiden_final <- DimPlot(leiden_obj, reduction = "umap", group.by = "leiden_clusters",
                          label = TRUE, label.size = 4) +
  ggtitle(paste0("Leiden Clustering (res ", final_resolution, ", n=", n_leiden_clusters, ")"))

p_leiden_sex <- DimPlot(leiden_obj, reduction = "umap", group.by = "sex",
                        cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
  ggtitle("By Sex")

p_leiden_combined <- p_leiden_final + p_leiden_sex
save_plot_multi(p_leiden_combined, "04_final_leiden_UMAP",
                output_dir = subdirs$plots_leiden, width = 14, height = 6)

p_leiden_batch <- DimPlot(leiden_obj, reduction = "umap", group.by = "batch") +
  ggtitle("By Batch")
save_plot_multi(p_leiden_batch, "05_leiden_by_batch",
                output_dir = subdirs$plots_leiden, width = 8, height = 6)

p_leiden_sample <- DimPlot(leiden_obj, reduction = "umap", group.by = "sample_name") +
  ggtitle("By Sample")
save_plot_multi(p_leiden_sample, "06_leiden_by_sample",
                output_dir = subdirs$plots_leiden, width = 10, height = 6)

cluster_sizes <- as.data.frame(table(leiden_obj$leiden_clusters))
colnames(cluster_sizes) <- c("Cluster", "Cells")
cluster_sizes$Cluster <- factor(cluster_sizes$Cluster,
                                 levels = cluster_sizes$Cluster[order(cluster_sizes$Cells, decreasing = TRUE)])

p_sizes <- ggplot(cluster_sizes, aes(x = Cluster, y = Cells, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Leiden Cluster Sizes", x = "Cluster", y = "Number of Cells")

save_plot_multi(p_sizes, "07_leiden_cluster_sizes",
                output_dir = subdirs$plots_leiden, width = 10, height = 6)

# ==============================================================================
# Create cluster composition table
# ==============================================================================
cluster_comp <- table(Cluster = leiden_obj$leiden_clusters, Sex = leiden_obj$sex)
cluster_comp_df <- as.data.frame.matrix(cluster_comp)
cluster_comp_df$Cluster <- rownames(cluster_comp_df)

# Handle case where not all sex levels are present
sex_cols <- intersect(c("Female", "Male"), colnames(cluster_comp_df))
cluster_comp_df$Total <- rowSums(cluster_comp_df[, sex_cols, drop = FALSE])
cluster_comp_df <- cluster_comp_df[, c("Cluster", sex_cols, "Total")]

write.csv(cluster_comp_df,
          file.path(output_dirs$leiden_clustering, "leiden_cluster_composition.csv"),
          row.names = FALSE)

# Composition by sample
sample_comp <- table(Cluster = leiden_obj$leiden_clusters, Sample = leiden_obj$sample_name)
sample_comp_df <- as.data.frame.matrix(sample_comp)
sample_comp_df$Cluster <- rownames(sample_comp_df)
write.csv(sample_comp_df,
          file.path(output_dirs$leiden_clustering, "leiden_cluster_by_sample.csv"),
          row.names = FALSE)

# ==============================================================================
# Save objects
# ==============================================================================
leiden_rds_path <- file.path(output_dirs$objects, "leiden_clustered_object.rds")
saveRDS(leiden_obj, leiden_rds_path)
cat("\nSaved:", leiden_rds_path, "\n")
cat("Size:", round(file.size(leiden_rds_path) / 1e6, 2), "MB\n")

# Also save as final object
final_rds_path <- file.path(output_dirs$objects, "07_final_object.rds")
saveRDS(clustered_obj, final_rds_path)
cat("Saved final object:", final_rds_path, "\n")

if (clustering_method_used == "Leiden") {
  clustered_obj <- leiden_obj
}

cat("\n>>> FINAL CLUSTERED OBJECT <<<\n")
print_object_structure(clustered_obj, "Final Clustered")

leiden_file <- file.path(output_dirs$objects, "07_leiden_data.RData")
save(leiden_obj, clustered_obj, clustering_method_used, quality_df, file = leiden_file)
cat("\nLeiden data saved to:", leiden_file, "\n")

write_readme(output_dirs$leiden_clustering, "Leiden Clustering",
             paste0("Leiden clustering with resolution testing.\n\n",
                    "Purpose: ", leiden_purpose, "\n",
                    "Resolutions tested: ", paste(leiden_resolutions, collapse = ", "), "\n",
                    "Final resolution: ", final_resolution, "\n",
                    "Number of clusters: ", n_leiden_clusters, "\n",
                    "Algorithm: ", leiden_algorithm, " (Leiden)\n",
                    "k neighbors: ", leiden_n_neighbors, "\n"),
             list("resolution_summary.csv" = "Clusters per resolution",
                  "leiden_cluster_composition.csv" = "Cells per cluster by sex",
                  "leiden_cluster_by_sample.csv" = "Cells per cluster by sample",
                  "clustering_quality.csv" = "Quality metrics per resolution"))

cat("\n>>> MODULE 07 COMPLETE <<<\n")
cat("Final clustering method:", clustering_method_used, "\n")
cat("Number of clusters:", length(unique(clustered_obj$seurat_clusters)), "\n")
