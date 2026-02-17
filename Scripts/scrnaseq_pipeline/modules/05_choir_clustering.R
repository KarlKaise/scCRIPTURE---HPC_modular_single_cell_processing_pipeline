#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 05: CLUSTERING (CHOIR + scAURA)
# ==============================================================================
#
# This module performs clustering on the integrated data using one or both of:
# - CHOIR: Hierarchical clustering with random forest significance testing
# - scAURA: Graph debiased contrastive learning
#
# Both methods can run independently. The parameter scice_clustering_source
# controls which result is used as primary clustering for downstream analysis.
#
# INPUT: Integrated Seurat object from Module 04
# OUTPUT: Clustered Seurat object with cluster assignments
#
# UPDATES:
# - 2026-01-03: Improved JoinLayers handling for Seurat v5
# - 2026-01-03: Added explicit layer verification before CHOIR
# - 2026-02-04: Added scAURA clustering alongside CHOIR
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 05: CLUSTERING (CHOIR + scAURA)\n")
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
  library(patchwork)
  library(Matrix)
  library(reticulate)
})

out_base <- params$out_root
load(file.path(out_base, "objects", "pipeline_environment.RData"))

# Load integration data if available, otherwise use normalized data
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
  cat("Loaded normalization data (no integration).\n")
} else {
  stop("No input data found. Run previous modules first.")
}

# ==============================================================================
# Check which clustering methods to run
# ==============================================================================
run_choir <- isTRUE(params$run_choir_clustering) && has_CHOIR
run_scaura <- isTRUE(params$run_scaura_clustering)

cat("\n--- Clustering methods to run ---\n")
cat("  CHOIR:", run_choir, "\n")
if (isTRUE(params$run_choir_clustering) && !has_CHOIR) {
  cat("    (CHOIR requested but package not available)\n")
}
cat("  scAURA:", run_scaura, "\n")
cat("  Clustering source for scICE:", params$scice_clustering_source, "\n\n")

# ==============================================================================
# Skip if both methods disabled
# ==============================================================================
if (!run_choir && !run_scaura) {
  cat("Both CHOIR and scAURA clustering disabled. Skipping Module 05.\n")

  choir_success <- FALSE
  scaura_success <- FALSE
  clustering_method_used <- "none"
  clustered_obj <- NULL

  clustering_file <- file.path(output_dirs$objects, "05_clustering_data.RData")
  save(choir_success, scaura_success, clustering_method_used, clustered_obj, file = clustering_file)

  cat("\n>>> MODULE 05 SKIPPED <<<\n")
  quit(save = "no", status = 0)
}

if (run_choir) {
  library(CHOIR)
}

# ==============================================================================
# Initialize variables
# ==============================================================================
choir_success <- FALSE
scaura_success <- FALSE
clustering_method_used <- "none"
clustered_obj <- NULL
choir_clustered <- NULL
scaura_clustered <- NULL

# ==============================================================================
# Select input object
# ==============================================================================
cat("--- Selecting input object for clustering ---\n")

if (exists("multi_integrated") && !is.null(multi_integrated)) {
  clustering_input <- multi_integrated
  clustering_reduction <- best_reduction
  cat("Using integrated object with reduction:", clustering_reduction, "\n")
} else if (exists("sct_object") && !is.null(sct_object)) {
  clustering_input <- sct_object
  clustering_reduction <- "pca"
  cat("Using SCTransform object with reduction: pca\n")
} else if (exists("scran_object") && !is.null(scran_object)) {
  clustering_input <- scran_object
  clustering_reduction <- "pca"
  cat("Using scran object with reduction: pca\n")
} else if (exists("lognorm_object") && !is.null(lognorm_object)) {
  clustering_input <- lognorm_object
  clustering_reduction <- "pca"
  cat("Using LogNormalize object with reduction: pca\n")
} else {
  cat("WARNING: No suitable object available for clustering\n")

  clustering_file <- file.path(output_dirs$objects, "05_clustering_data.RData")
  save(choir_success, scaura_success, clustering_method_used, clustered_obj, file = clustering_file)

  cat("\n>>> MODULE 05 SKIPPED (no input) <<<\n")
  quit(save = "no", status = 0)
}

# ==============================================================================
# Prepare object for clustering (Seurat v5 compatibility)
# ==============================================================================
cat("\n--- Preparing object for clustering ---\n")

# Join RNA layers if they exist and are split
tryCatch({
  if ("RNA" %in% names(clustering_input@assays)) {
    rna_layers <- Layers(clustering_input[["RNA"]])
    counts_layers <- grep("^counts", rna_layers, value = TRUE)

    if (length(counts_layers) > 1) {
      cat("Found", length(counts_layers), "counts layers - joining...\n")
      clustering_input[["RNA"]] <- JoinLayers(clustering_input[["RNA"]])
      cat("Successfully joined RNA layers\n")
    } else {
      cat("RNA layers already unified\n")
    }

    # Verify counts layer is accessible
    counts_check <- tryCatch({
      GetAssayData(clustering_input, assay = "RNA", layer = "counts")
    }, error = function(e) NULL)

    if (is.null(counts_check)) {
      cat("WARNING: Could not access counts layer after JoinLayers\n")
    } else {
      cat("Verified: counts layer accessible (", nrow(counts_check), " genes x ",
          ncol(counts_check), " cells)\n", sep = "")
    }
  }
}, error = function(e) {
  cat("Note: Layer preparation message -", e$message, "\n")
})

# Also join SCT layers if present
tryCatch({
  if ("SCT" %in% names(clustering_input@assays)) {
    sct_layers <- Layers(clustering_input[["SCT"]])
    if (length(grep("^counts", sct_layers, value = TRUE)) > 1) {
      cat("Joining SCT layers...\n")
      clustering_input[["SCT"]] <- JoinLayers(clustering_input[["SCT"]])
    }
  }
}, error = function(e) {
  cat("Note: SCT layer join message -", e$message, "\n")
})

cat("\n>>> INPUT OBJECT FOR CLUSTERING <<<\n")
print_object_structure(clustering_input, "Clustering Input")

# Determine dimensions to use
if (!is.null(clustering_reduction) && clustering_reduction %in% names(clustering_input@reductions)) {
  nd <- ncol(Embeddings(clustering_input, reduction = clustering_reduction))
  dims_to_use <- min(params$dims_use, nd)
  cat("Using dimensions 1:", dims_to_use, "from reduction:", clustering_reduction, "\n")
} else {
  dims_to_use <- params$dims_use
  clustering_reduction <- "pca"
  cat("Using default PCA reduction\n")
}

# ==============================================================================
# Create output directories
# ==============================================================================
choir_plot_dir <- file.path(output_dirs$choir_clustering, "plots")
dir.create(choir_plot_dir, showWarnings = FALSE, recursive = TRUE)

scaura_output_dir <- file.path(output_dirs$choir_clustering, "scaura")
dir.create(scaura_output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# RUN CHOIR CLUSTERING
# ==============================================================================
if (run_choir) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("RUNNING CHOIR CLUSTERING\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  choir_input <- clustering_input

  # Set CHOIR parameters with defaults if not specified
  choir_alpha <- if (!is.null(params$choir_alpha)) params$choir_alpha else 0.05
  choir_use_assay <- if (!is.null(params$choir_use_assay)) params$choir_use_assay else "RNA"

  cat("Parameters:\n")
  cat("  alpha:", choir_alpha, "\n")
  cat("  p_adjust: bonferroni\n")
  cat("  use_assay:", choir_use_assay, "\n")
  cat("  n_cores: 4\n\n")

  tryCatch({
    choir_result <- CHOIR::CHOIR(
      object = choir_input,
      key = "CHOIR",
      alpha = choir_alpha,
      p_adjust = "bonferroni",
      use_assay = choir_use_assay,
      n_cores = 4,
      verbose = TRUE
    )

    cat("\n>>> CHOIR CLUSTERING COMPLETE <<<\n")

    if (is(choir_result, "Seurat")) {
      choir_clustered <- choir_result

      choir_cols <- grep("^CHOIR_clusters", colnames(choir_clustered@meta.data), value = TRUE)
      if (length(choir_cols) > 0) {
        cat("CHOIR cluster columns found:", paste(choir_cols, collapse = ", "), "\n")

        main_choir_col <- choir_cols[1]
        choir_clustered$choir_clusters <- choir_clustered@meta.data[[main_choir_col]]

        n_clusters <- length(unique(choir_clustered$choir_clusters))
        cat("Number of CHOIR clusters:", n_clusters, "\n")
        
        choir_success <- TRUE
      } else {
        cat("WARNING: No CHOIR cluster columns found in result\n")
      }
    } else {
      cat("CHOIR result type:", class(choir_result), "\n")
    }

    if (choir_success) {
      # Create cluster composition table
      cluster_comp <- table(Cluster = choir_clustered$choir_clusters, Sex = choir_clustered$sex)
      cluster_comp_df <- as.data.frame.matrix(cluster_comp)
      cluster_comp_df$Cluster <- rownames(cluster_comp_df)
      cluster_comp_df <- cluster_comp_df[, c("Cluster", names(cluster_comp_df)[names(cluster_comp_df) != "Cluster"])]

      cluster_comp_path <- file.path(output_dirs$choir_clustering, "choir_cluster_composition.csv")
      write.csv(cluster_comp_df, cluster_comp_path, row.names = FALSE)
      cat("Saved cluster composition:", cluster_comp_path, "\n")
    }

  }, error = function(e) {
    cat("\n[ERROR] CHOIR clustering failed:", conditionMessage(e), "\n")
    choir_success <<- FALSE
  })
}

# ==============================================================================
# RUN scAURA CLUSTERING
# ==============================================================================
if (run_scaura) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("RUNNING scAURA CLUSTERING\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  scaura_input <- clustering_input

  tryCatch({
    # Run scAURA using helper function
    scaura_results <- run_scaura_clustering(
      obj = scaura_input,
      params = params,
      output_dir = scaura_output_dir
    )

    cat("\n>>> scAURA CLUSTERING COMPLETE <<<\n")

    # Create scAURA clustered object
    scaura_clustered <- scaura_input

    # Add scAURA embeddings as a reduction
    scaura_clustered[["scaura"]] <- CreateDimReducObject(
      embeddings = scaura_results$embeddings_post_ssc,
      key = "scAURA_",
      assay = DefaultAssay(scaura_clustered)
    )

    # Add cluster labels
    scaura_clustered$scaura_clusters <- factor(scaura_results$ssc_labels)
    scaura_clustered$scaura_kmeans_clusters <- factor(scaura_results$kmeans_labels)

    n_clusters <- length(unique(scaura_clustered$scaura_clusters))
    cat("Number of scAURA clusters:", n_clusters, "\n")
    cat("scAURA embedding dimensions:", ncol(scaura_results$embeddings_post_ssc), "\n")

    scaura_success <- TRUE

    # Create cluster composition table
    cluster_comp <- table(Cluster = scaura_clustered$scaura_clusters, Sex = scaura_clustered$sex)
    cluster_comp_df <- as.data.frame.matrix(cluster_comp)
    cluster_comp_df$Cluster <- rownames(cluster_comp_df)
    cluster_comp_df <- cluster_comp_df[, c("Cluster", names(cluster_comp_df)[names(cluster_comp_df) != "Cluster"])]

    cluster_comp_path <- file.path(scaura_output_dir, "scaura_cluster_composition.csv")
    write.csv(cluster_comp_df, cluster_comp_path, row.names = FALSE)
    cat("Saved cluster composition:", cluster_comp_path, "\n")

  }, error = function(e) {
    cat("\n[ERROR] scAURA clustering failed:", conditionMessage(e), "\n")
    cat("Error details:", e$message, "\n")
    scaura_success <<- FALSE
  })
}

# ==============================================================================
# SELECT PRIMARY CLUSTERING METHOD
# ==============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SELECTING PRIMARY CLUSTERING METHOD\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Results:\n")
cat("  CHOIR success:", choir_success, "\n")
cat("  scAURA success:", scaura_success, "\n")
cat("  Selection mode:", params$scice_clustering_source, "\n\n")

# Determine which clustering to use as primary
selected_method <- "none"

if (params$scice_clustering_source == "choir") {
  if (choir_success) {
    selected_method <- "choir"
  } else {
    cat("WARNING: CHOIR requested but failed. ")
    if (scaura_success) {
      selected_method <- "scaura"
      cat("Falling back to scAURA.\n")
    } else {
      cat("No clustering available.\n")
    }
  }
} else if (params$scice_clustering_source == "scaura") {
  if (scaura_success) {
    selected_method <- "scaura"
  } else {
    cat("WARNING: scAURA requested but failed. ")
    if (choir_success) {
      selected_method <- "choir"
      cat("Falling back to CHOIR.\n")
    } else {
      cat("No clustering available.\n")
    }
  }
} else {
  # auto mode: CHOIR > scAURA
  if (choir_success) {
    selected_method <- "choir"
  } else if (scaura_success) {
    selected_method <- "scaura"
  }
}

cat(">>> SELECTED PRIMARY METHOD:", selected_method, "<<<\n\n")
clustering_method_used <- selected_method

# ==============================================================================
# CREATE FINAL CLUSTERED OBJECT
# ==============================================================================
if (selected_method != "none") {
  cat("--- Creating final clustered object ---\n")
  
  # Start with the appropriate base object
  if (selected_method == "choir" && choir_success) {
    clustered_obj <- choir_clustered
    clustered_obj$seurat_clusters <- clustered_obj$choir_clusters
    Idents(clustered_obj) <- "seurat_clusters"
    cat("  Set seurat_clusters from CHOIR\n")
    
    # Add scAURA results if available
    if (scaura_success) {
      clustered_obj$scaura_clusters <- scaura_clustered$scaura_clusters
      clustered_obj[["scaura"]] <- scaura_clustered[["scaura"]]
      cat("  Added scAURA clusters and embeddings\n")
    }
    
  } else if (selected_method == "scaura" && scaura_success) {
    clustered_obj <- scaura_clustered
    clustered_obj$seurat_clusters <- clustered_obj$scaura_clusters
    Idents(clustered_obj) <- "seurat_clusters"
    cat("  Set seurat_clusters from scAURA\n")
    
    # Add CHOIR results if available
    if (choir_success) {
      clustered_obj$choir_clusters <- choir_clustered$choir_clusters
      # Copy any CHOIR-specific columns
      choir_cols <- grep("^CHOIR_", colnames(choir_clustered@meta.data), value = TRUE)
      for (col in choir_cols) {
        clustered_obj[[col]] <- choir_clustered@meta.data[[col]]
      }
      cat("  Added CHOIR clusters\n")
    }
  }
  
  # Store which method was selected
  clustered_obj$primary_clustering_method <- selected_method
  
  # Run UMAP if not present
  if (!"umap" %in% names(clustered_obj@reductions)) {
    cat("\nRunning UMAP for visualization...\n")
    umap_reduction <- clustering_reduction
    
    # For scAURA, optionally use scAURA embeddings for UMAP
    if (selected_method == "scaura" && "scaura" %in% names(clustered_obj@reductions)) {
      umap_reduction <- "scaura"
      cat("  Using scAURA embeddings for UMAP\n")
    }
    
    if (umap_reduction %in% names(clustered_obj@reductions)) {
      n_dims <- ncol(Embeddings(clustered_obj, umap_reduction))
      clustered_obj <- RunUMAP(clustered_obj, reduction = umap_reduction,
                               dims = 1:min(30, n_dims), verbose = FALSE)
    } else if ("pca" %in% names(clustered_obj@reductions)) {
      clustered_obj <- RunUMAP(clustered_obj, reduction = "pca",
                               dims = 1:min(30, ncol(Embeddings(clustered_obj, "pca"))),
                               verbose = FALSE)
    }
  }
  
  # Create scAURA-specific UMAP if scAURA succeeded
  if (scaura_success && "scaura" %in% names(clustered_obj@reductions)) {
    cat("Creating scAURA-specific UMAP...\n")
    n_dims <- ncol(Embeddings(clustered_obj, "scaura"))
    clustered_obj <- RunUMAP(clustered_obj, reduction = "scaura",
                             dims = 1:min(30, n_dims),
                             reduction.name = "umap_scaura",
                             verbose = FALSE)
  }
}

# ==============================================================================
# CREATE VISUALIZATIONS
# ==============================================================================
if (!is.null(clustered_obj)) {
  cat("\n--- Creating visualizations ---\n")
  
  # Plot 1: Primary clusters
  p_clusters <- DimPlot(clustered_obj, reduction = "umap", group.by = "seurat_clusters",
                        label = TRUE, label.size = 4) +
    ggtitle(paste0("Primary Clusters (", selected_method, ", n=", 
                   length(unique(clustered_obj$seurat_clusters)), ")"))

  p_sex <- DimPlot(clustered_obj, reduction = "umap", group.by = "sex",
                   cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
    ggtitle("Clusters by Sex")

  p_combined <- p_clusters + p_sex
  save_plot_multi(p_combined, "01_primary_clusters",
                  output_dir = choir_plot_dir, width = 14, height = 6)
  
  # Plot 2: By Batch
  p_batch <- DimPlot(clustered_obj, reduction = "umap", group.by = "batch") +
    ggtitle("Clusters by Batch")
  save_plot_multi(p_batch, "02_clusters_by_batch",
                  output_dir = choir_plot_dir, width = 8, height = 6)
  
  # Plot 3: By Sample
  p_sample <- DimPlot(clustered_obj, reduction = "umap", group.by = "sample_name") +
    ggtitle("Clusters by Sample")
  save_plot_multi(p_sample, "03_clusters_by_sample",
                  output_dir = choir_plot_dir, width = 10, height = 6)
  
  # Plot 4: Cluster sizes
  cluster_sizes <- as.data.frame(table(clustered_obj$seurat_clusters))
  colnames(cluster_sizes) <- c("Cluster", "Cells")
  cluster_sizes$Cluster <- factor(cluster_sizes$Cluster,
                                   levels = cluster_sizes$Cluster[order(cluster_sizes$Cells, decreasing = TRUE)])

  p_sizes <- ggplot(cluster_sizes, aes(x = Cluster, y = Cells, fill = Cluster)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = paste0("Cluster Sizes (", selected_method, ")"), 
         x = "Cluster", y = "Number of Cells")

  save_plot_multi(p_sizes, "04_cluster_sizes",
                  output_dir = choir_plot_dir, width = 10, height = 6)
  
  # Plot 5: Method comparison (if both succeeded)
  if (choir_success && scaura_success) {
    cat("Creating method comparison plots...\n")
    
    p_choir <- DimPlot(clustered_obj, reduction = "umap", group.by = "choir_clusters",
                       label = TRUE, label.size = 3) +
      ggtitle(paste0("CHOIR (n=", length(unique(clustered_obj$choir_clusters)), ")"))
    
    p_scaura <- DimPlot(clustered_obj, reduction = "umap", group.by = "scaura_clusters",
                        label = TRUE, label.size = 3) +
      ggtitle(paste0("scAURA (n=", length(unique(clustered_obj$scaura_clusters)), ")"))
    
    p_comparison <- p_choir + p_scaura
    save_plot_multi(p_comparison, "05_method_comparison",
                    output_dir = choir_plot_dir, width = 14, height = 6)
    
    # If scAURA UMAP exists, show comparison
    if ("umap_scaura" %in% names(clustered_obj@reductions)) {
      p_umap_std <- DimPlot(clustered_obj, reduction = "umap", group.by = "scaura_clusters",
                            label = TRUE, label.size = 3) +
        ggtitle("scAURA clusters on standard UMAP")
      
      p_umap_scaura <- DimPlot(clustered_obj, reduction = "umap_scaura", group.by = "scaura_clusters",
                               label = TRUE, label.size = 3) +
        ggtitle("scAURA clusters on scAURA UMAP")
      
      p_umap_comp <- p_umap_std + p_umap_scaura
      save_plot_multi(p_umap_comp, "06_scaura_umap_comparison",
                      output_dir = choir_plot_dir, width = 14, height = 6)
    }
    
    # Cluster correspondence table
    correspondence <- table(CHOIR = clustered_obj$choir_clusters, 
                           scAURA = clustered_obj$scaura_clusters)
    write.csv(as.data.frame.matrix(correspondence), 
              file.path(output_dirs$choir_clustering, "cluster_correspondence_choir_scaura.csv"))
    cat("Saved cluster correspondence table\n")
  }
  
  # Save clustered object
  cluster_rds_path <- file.path(output_dirs$objects, "clustered_object.rds")
  saveRDS(clustered_obj, cluster_rds_path)
  cat("\nSaved:", cluster_rds_path, "\n")
  
  # Also save method-specific objects for reference
  if (choir_success) {
    choir_rds_path <- file.path(output_dirs$objects, "choir_clustered_object.rds")
    saveRDS(clustered_obj, choir_rds_path)
    cat("Saved:", choir_rds_path, "\n")
  }
  
  if (scaura_success) {
    scaura_rds_path <- file.path(output_dirs$objects, "scaura_clustered_object.rds")
    saveRDS(clustered_obj, scaura_rds_path)
    cat("Saved:", scaura_rds_path, "\n")
  }
  
  cat("\n>>> FINAL CLUSTERED OBJECT <<<\n")
  print_object_structure(clustered_obj, "Clustered Object")
}

# ==============================================================================
# Save data for next module
# ==============================================================================
clustering_file <- file.path(output_dirs$objects, "05_clustering_data.RData")
save(choir_success, scaura_success, clustering_method_used, clustered_obj, file = clustering_file)
cat("\nClustering data saved to:", clustering_file, "\n")

# Also save as 05_choir_data.RData for backward compatibility
choir_file <- file.path(output_dirs$objects, "05_choir_data.RData")
choir_success_compat <- choir_success || scaura_success
save(choir_success = choir_success_compat, 
     clustering_method_used, 
     clustered_obj, 
     file = choir_file)
cat("Backward-compatible file saved to:", choir_file, "\n")

# Write README
readme_content <- paste0(
  "Module 05: Clustering (CHOIR + scAURA)\n\n",
  "Methods run:\n",
  "  CHOIR: ", ifelse(run_choir, ifelse(choir_success, "SUCCESS", "FAILED"), "DISABLED"), "\n",
  "  scAURA: ", ifelse(run_scaura, ifelse(scaura_success, "SUCCESS", "FAILED"), "DISABLED"), "\n\n",
  "Primary clustering method: ", clustering_method_used, "\n",
  "Selection mode: ", params$scice_clustering_source, "\n\n",
  "Available cluster columns in output:\n",
  "  seurat_clusters: Primary clusters (from ", clustering_method_used, ")\n",
  if (choir_success) "  choir_clusters: CHOIR cluster assignments\n" else "",
  if (scaura_success) "  scaura_clusters: scAURA cluster assignments\n" else "",
  if (scaura_success) "\nscAURA reduction: 'scaura' (64-dimensional embeddings)\n" else ""
)

write_readme(output_dirs$choir_clustering, "Module 05: Clustering",
             readme_content,
             list("clustered_object.rds" = "Final clustered Seurat object",
                  "choir_cluster_composition.csv" = "CHOIR cluster composition by sex",
                  "scaura/" = "scAURA outputs (embeddings, labels, config)",
                  "plots/" = "Visualization plots",
                  "cluster_correspondence_choir_scaura.csv" = "Cluster correspondence between methods"))

cat("\n>>> MODULE 05 COMPLETE <<<\n")
cat("CHOIR success:", choir_success, "\n")
cat("scAURA success:", scaura_success, "\n")
cat("Primary method:", clustering_method_used, "\n")