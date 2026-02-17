#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 09: GENE EXPRESSION VISUALIZATION (MULTI-SAMPLE PIPELINE)
# ==============================================================================
#
# This module creates visualizations for genes of interest:
# - Feature plots (UMAP)
# - Violin plots by sex and sample
# - Dot plots by cluster
# - Heatmaps
# - Expression summary statistics
#
# INPUT: Clustered Seurat object from Module 07
# OUTPUT: Gene expression plots and summary tables
#
# UPDATES:
# - 2026-01-03: Improved object loading with multiple fallbacks
# - 2026-01-03: Enhanced JoinLayers handling for Seurat v5
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 09: GENE EXPRESSION VISUALIZATION\n")
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
  library(viridis)
  library(pheatmap)
  library(RColorBrewer)
  library(Matrix)
})

out_base <- params$out_root
load(file.path(out_base, "objects", "pipeline_environment.RData"))

# ==============================================================================
# Load clustered data (with multiple fallbacks)
# ==============================================================================
cat("--- Loading clustered data ---\n")

clustered_obj <- NULL

# Try multiple sources in order of preference
scice_rds <- file.path(out_base, "objects", "scice_subclustered_object.rds")
final_rds <- file.path(out_base, "objects", "07_final_object.rds")
leiden_file <- file.path(out_base, "objects", "07_leiden_data.RData")
choir_rds <- file.path(out_base, "objects", "choir_clustered_object.rds")

if (file.exists(scice_rds)) {
  clustered_obj <- readRDS(scice_rds)
  cat("Loaded scICE subclustered object.\n")
} else if (file.exists(final_rds)) {
  clustered_obj <- readRDS(final_rds)
  cat("Loaded final clustered object.\n")
} else if (file.exists(leiden_file)) {
  load(leiden_file)
  cat("Loaded Leiden clustered data.\n")
} else if (file.exists(choir_rds)) {
  clustered_obj <- readRDS(choir_rds)
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
viz_obj <- clustered_obj
DefaultAssay(viz_obj) <- "RNA"

# Join layers for visualization (Seurat v5)
tryCatch({
  rna_layers <- Layers(viz_obj[["RNA"]])
  counts_layers <- grep("^counts", rna_layers, value = TRUE)
  data_layers <- grep("^data", rna_layers, value = TRUE)
  
  if (length(counts_layers) > 1 || length(data_layers) > 1) {
    cat("Found split layers - joining for visualization...\n")
    viz_obj[["RNA"]] <- JoinLayers(viz_obj[["RNA"]])
    cat("Successfully joined RNA layers\n")
  } else {
    cat("RNA layers already unified\n")
  }
}, error = function(e) {
  cat("Note: Layer preparation message -", e$message, "\n")
})

cat("\n>>> OBJECT FOR VISUALIZATION <<<\n")
print_object_structure(viz_obj, "Visualization Input")

viz_dir <- output_dirs$gene_viz
plots_viz_dir <- subdirs$plots_gene_expression

# ==============================================================================
# Find genes of interest
# ==============================================================================
cat("\n--- Finding genes of interest ---\n")

# Combine all gene lists
sex_genes <- if (!is.null(params$sex_marker_genes)) params$sex_marker_genes else character(0)
cp_genes <- if (!is.null(params$cp_marker_genes)) params$cp_marker_genes else character(0)
goi_genes <- if (!is.null(params$genes_of_interest)) params$genes_of_interest else character(0)

all_genes <- unique(c(sex_genes, cp_genes, goi_genes))

if (length(all_genes) == 0) {
  cat("No genes of interest specified in params. Using example genes.\n")
  # Use some common marker genes as fallback
  all_genes <- c("Ttr", "Aqp1", "Folr1", "Slc12a2", "Clic6")
}

cat("Requested genes:", paste(all_genes, collapse = ", "), "\n\n")

available_genes <- find_genes_in_object(all_genes, viz_obj)

cat("Found", length(available_genes), "of", length(all_genes), "genes:\n")
if (length(available_genes) > 0) {
  cat("  ", paste(available_genes, collapse = ", "), "\n")
}

missing_genes <- setdiff(all_genes, available_genes)
if (length(missing_genes) > 0) {
  cat("\nMissing genes:\n")
  cat("  ", paste(missing_genes, collapse = ", "), "\n")
}

# ==============================================================================
# Create visualizations
# ==============================================================================

if (length(available_genes) == 0) {
  cat("\n*** No genes of interest found in dataset ***\n")
  cat("\nExample genes in dataset:\n")
  set.seed(42)
  cat("  ", paste(sample(rownames(viz_obj), min(20, nrow(viz_obj))), collapse = ", "), "\n")

} else {

  # --------------------------------------------------------------------------
  # Feature Plots (UMAP)
  # --------------------------------------------------------------------------
  cat("\n--- Creating Feature Plots ---\n")

  if ("umap" %in% names(viz_obj@reductions)) {
    genes_per_plot <- 4
    n_plots <- ceiling(length(available_genes) / genes_per_plot)

    for (i in 1:n_plots) {
      start_idx <- (i - 1) * genes_per_plot + 1
      end_idx <- min(i * genes_per_plot, length(available_genes))
      genes_subset <- available_genes[start_idx:end_idx]

      tryCatch({
        p <- FeaturePlot(viz_obj, features = genes_subset, reduction = "umap",
                         ncol = 2, cols = c("lightgrey", "darkblue"))
        save_plot_multi(p, paste0("01_feature_plot_", i),
                        output_dir = plots_viz_dir, width = 12, height = 10)
        cat("  Created feature plot", i, "for:", paste(genes_subset, collapse = ", "), "\n")
      }, error = function(e) {
        cat("  Feature plot", i, "failed:", e$message, "\n")
      })
    }
  } else {
    cat("  UMAP reduction not found. Skipping feature plots.\n")
  }

  # --------------------------------------------------------------------------
  # Violin Plots by Sex
  # --------------------------------------------------------------------------
  cat("\n--- Creating Violin Plots by Sex ---\n")

  for (gene in available_genes) {
    tryCatch({
      p <- VlnPlot(viz_obj, features = gene, group.by = "sex", pt.size = 0,
                   cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
        stat_summary(fun = median, geom = "point", size = 3, color = "black") +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(title = paste0(gene, " Expression by Sex"),
             x = "Sex", y = "Expression")

      save_plot_multi(p, paste0("02_violin_sex_", gene),
                      output_dir = plots_viz_dir, width = 6, height = 5)
      cat("  Created violin plot for:", gene, "\n")
    }, error = function(e) {
      cat("  Violin plot for", gene, "failed\n")
    })
  }

  # Combined violin plot
  tryCatch({
    p_combined <- VlnPlot(viz_obj, features = available_genes, group.by = "sex",
                          pt.size = 0, ncol = min(4, length(available_genes)),
                          cols = c("Female" = "#E41A1C", "Male" = "#377EB8"))
    save_plot_multi(p_combined, "02_violin_all_genes",
                    output_dir = plots_viz_dir,
                    width = min(4, length(available_genes)) * 4,
                    height = ceiling(length(available_genes) / 4) * 4)
    cat("  Created combined violin plot\n")
  }, error = function(e) {
    cat("  Combined violin plot failed\n")
  })

  # --------------------------------------------------------------------------
  # Violin Plots by Sample
  # --------------------------------------------------------------------------
  cat("\n--- Creating Violin Plots by Sample ---\n")

  for (gene in available_genes[1:min(3, length(available_genes))]) {
    tryCatch({
      p <- VlnPlot(viz_obj, features = gene, group.by = "sample_name",
                   pt.size = 0) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste0(gene, " by Sample"))

      save_plot_multi(p, paste0("02b_violin_sample_", gene),
                      output_dir = plots_viz_dir, width = 12, height = 6)
      cat("  Created sample violin for:", gene, "\n")
    }, error = function(e) {
      cat("  Sample violin for", gene, "failed\n")
    })
  }

  # --------------------------------------------------------------------------
  # Violin Plots by Cluster
  # --------------------------------------------------------------------------
  cat("\n--- Creating Violin Plots by Cluster ---\n")

  for (gene in available_genes[1:min(3, length(available_genes))]) {
    tryCatch({
      p <- VlnPlot(viz_obj, features = gene, group.by = "seurat_clusters",
                   pt.size = 0, split.by = "sex",
                   cols = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste0(gene, " by Cluster and Sex"))

      save_plot_multi(p, paste0("03_violin_cluster_", gene),
                      output_dir = plots_viz_dir, width = 12, height = 6)
      cat("  Created cluster violin for:", gene, "\n")
    }, error = function(e) {
      cat("  Cluster violin for", gene, "failed\n")
    })
  }

  # --------------------------------------------------------------------------
  # Dot Plot
  # --------------------------------------------------------------------------
  cat("\n--- Creating Dot Plot ---\n")

  tryCatch({
    p <- DotPlot(viz_obj, features = available_genes, group.by = "seurat_clusters") +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Gene Expression by Cluster")

    save_plot_multi(p, "04_dotplot_clusters",
                    output_dir = plots_viz_dir,
                    width = max(8, length(unique(viz_obj$seurat_clusters)) * 0.5),
                    height = max(6, length(available_genes) * 0.4))
    cat("  Created cluster dot plot\n")
  }, error = function(e) {
    cat("  Dot plot failed:", e$message, "\n")
  })

  # Dot plot by sex
  tryCatch({
    p <- DotPlot(viz_obj, features = available_genes, group.by = "sex") +
      coord_flip() +
      labs(title = "Gene Expression by Sex")

    save_plot_multi(p, "04_dotplot_sex",
                    output_dir = plots_viz_dir, width = 8, height = max(6, length(available_genes) * 0.4))
    cat("  Created sex dot plot\n")
  }, error = function(e) {
    cat("  Sex dot plot failed\n")
  })

  # Dot plot by sample
  tryCatch({
    p <- DotPlot(viz_obj, features = available_genes, group.by = "sample_name") +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Gene Expression by Sample")

    save_plot_multi(p, "04_dotplot_sample",
                    output_dir = plots_viz_dir, width = 12, height = max(6, length(available_genes) * 0.4))
    cat("  Created sample dot plot\n")
  }, error = function(e) {
    cat("  Sample dot plot failed\n")
  })

  # --------------------------------------------------------------------------
  # Heatmap
  # --------------------------------------------------------------------------
  cat("\n--- Creating Heatmap ---\n")

  tryCatch({
    avg_expr <- AverageExpression(viz_obj, features = available_genes,
                                   group.by = c("seurat_clusters", "sex"),
                                   assays = "RNA", return.seurat = FALSE)

    avg_mat <- avg_expr$RNA

    if (ncol(avg_mat) > 1 && nrow(avg_mat) > 1) {
      avg_mat_scaled <- t(scale(t(avg_mat)))

      p_heatmap <- pheatmap(avg_mat_scaled,
                            color = colorRampPalette(c("blue", "white", "red"))(100),
                            cluster_rows = TRUE,
                            cluster_cols = TRUE,
                            fontsize = 10,
                            main = "Gene Expression Heatmap\n(scaled by gene)",
                            silent = TRUE)

      for (fmt in c("png", "pdf")) {
        filepath <- file.path(plots_viz_dir, paste0("05_heatmap_expression.", fmt))
        if (fmt == "png") {
          png(filepath, width = 10, height = 8, units = "in", res = 300)
        } else {
          pdf(filepath, width = 10, height = 8)
        }
        print(p_heatmap)
        dev.off()
      }
      cat("  Created expression heatmap\n")
    }
  }, error = function(e) {
    cat("  Heatmap failed:", e$message, "\n")
  })

  # --------------------------------------------------------------------------
  # Expression Summary Table
  # --------------------------------------------------------------------------
  cat("\n--- Computing Expression Summary ---\n")

  tryCatch({
    expr_mat <- GetAssayData(viz_obj, layer = "data")
    sex_labels <- viz_obj$sex

    expr_summary <- data.frame(
      gene = available_genes,
      mean_female = sapply(available_genes, function(g) mean(expr_mat[g, sex_labels == "Female"])),
      mean_male = sapply(available_genes, function(g) mean(expr_mat[g, sex_labels == "Male"])),
      median_female = sapply(available_genes, function(g) median(expr_mat[g, sex_labels == "Female"])),
      median_male = sapply(available_genes, function(g) median(expr_mat[g, sex_labels == "Male"])),
      pct_female = sapply(available_genes, function(g) mean(expr_mat[g, sex_labels == "Female"] > 0) * 100),
      pct_male = sapply(available_genes, function(g) mean(expr_mat[g, sex_labels == "Male"] > 0) * 100)
    )

    expr_summary$log2FC <- log2((expr_summary$mean_male + 0.001) / (expr_summary$mean_female + 0.001))
    expr_summary$delta_pct <- expr_summary$pct_male - expr_summary$pct_female

    numeric_cols <- sapply(expr_summary, is.numeric)
    expr_summary[numeric_cols] <- lapply(expr_summary[numeric_cols], round, 4)

    write.csv(expr_summary, file.path(viz_dir, "gene_expression_summary.csv"), row.names = FALSE)

    cat("\n>>> EXPRESSION SUMMARY <<<\n")
    print(expr_summary)
    cat("\nSaved to:", file.path(viz_dir, "gene_expression_summary.csv"), "\n")

  }, error = function(e) {
    cat("  Expression summary failed:", e$message, "\n")
  })

  # --------------------------------------------------------------------------
  # Expression by Sample Summary
  # --------------------------------------------------------------------------
  cat("\n--- Computing Sample Expression Summary ---\n")

  tryCatch({
    sample_expr <- AverageExpression(viz_obj, features = available_genes,
                                      group.by = "sample_name",
                                      assays = "RNA", return.seurat = FALSE)

    sample_mat <- as.data.frame(t(sample_expr$RNA))
    sample_mat$sample <- rownames(sample_mat)
    sample_mat <- sample_mat[, c("sample", available_genes)]

    write.csv(sample_mat, file.path(viz_dir, "gene_expression_by_sample.csv"), row.names = FALSE)
    cat("  Saved sample expression summary\n")

  }, error = function(e) {
    cat("  Sample expression summary failed\n")
  })

  # --------------------------------------------------------------------------
  # Expression by Cluster Summary
  # --------------------------------------------------------------------------
  cat("\n--- Computing Cluster Expression Summary ---\n")

  tryCatch({
    cluster_expr <- AverageExpression(viz_obj, features = available_genes,
                                       group.by = "seurat_clusters",
                                       assays = "RNA", return.seurat = FALSE)

    cluster_mat <- as.data.frame(t(cluster_expr$RNA))
    cluster_mat$cluster <- rownames(cluster_mat)
    cluster_mat <- cluster_mat[, c("cluster", available_genes)]

    write.csv(cluster_mat, file.path(viz_dir, "gene_expression_by_cluster.csv"), row.names = FALSE)
    cat("  Saved cluster expression summary\n")

  }, error = function(e) {
    cat("  Cluster expression summary failed\n")
  })
}

# ==============================================================================
# Save visualization data
# ==============================================================================
viz_file <- file.path(output_dirs$objects, "09_visualization_data.RData")
save(available_genes, file = viz_file)
cat("\nVisualization data saved to:", viz_file, "\n")

# Write README
write_readme(viz_dir, "Gene Expression Visualization",
             paste0("Visualizations for genes of interest.\n\n",
                    "Requested: ", length(all_genes), " genes\n",
                    "Found: ", length(available_genes), " genes\n\n",
                    "Plots created:\n",
                    "- Feature plots (UMAP colored by expression)\n",
                    "- Violin plots (by sex, sample, and cluster)\n",
                    "- Dot plots (average expression and percent expressed)\n",
                    "- Heatmaps (scaled expression across clusters)"),
             list("gene_expression_summary.csv" = "Expression statistics by sex",
                  "gene_expression_by_sample.csv" = "Average expression per sample",
                  "gene_expression_by_cluster.csv" = "Average expression per cluster"))

cat("\n>>> MODULE 09 COMPLETE <<<\n")
cat("Genes visualized:", length(available_genes), "\n")
