#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 02: QC VALIDATION (MULTI-SAMPLE PIPELINE)
# ==============================================================================
#
# This module performs QC validation and filtering on the merged dataset.
#
# Steps:
# 1. Validate QC metrics exist for all samples
# 2. Apply QC thresholds (min_features, max_features, max_percent_mt) - OPTIONAL
# 3. Apply doublet filtering (if enabled) - OPTIONAL
# 4. Apply hemoglobin filtering (if enabled) - OPTIONAL
# 5. Filter genes (min_cells_per_gene) - OPTIONAL
# 6. Generate QC reports and visualizations
#
# INPUT: Merged Seurat object from Module 01
# OUTPUT: QC-filtered merged Seurat object
#
# FILTERING CONTROL FLAGS (set in params.R):
#   skip_all_filtering       - Master switch to bypass ALL filtering
#   apply_qc_filtering       - Enable/disable standard QC filtering
#   filter_by_min_features   - Filter by minimum gene count
#   filter_by_max_features   - Filter by maximum gene count
#   filter_by_percent_mt     - Filter by mitochondrial percentage
#   filter_doublets          - Filter cells marked as doublets
#   filter_hemoglobin        - Filter cells with high hemoglobin
#   filter_genes_by_min_cells - Filter genes by minimum cell count
#
# UPDATES:
# - 2026-01-05: Added comprehensive filtering control flags
# - 2026-01-05: Added optional doublet filtering
# - 2026-01-05: Made all filtering steps conditional
# - 2026-01-03: Added hemoglobin filtering support
# - 2026-01-03: Updated GetAssayData calls for Seurat v5 (layer= not slot=)
# - 2026-01-03: Improved JoinLayers handling for merged objects
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 02: QC VALIDATION\n")
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
  library(Matrix)
  library(tidyr)
})

out_base <- params$out_root
load(file.path(out_base, "objects", "pipeline_environment.RData"))
load(file.path(out_base, "objects", "01_loaded_data.RData"))

# ==============================================================================
# Seurat v5 FIX: Join multiple counts.* layers into a single "counts" layer
# ==============================================================================
if ("merged_obj" %in% ls()) {
  if ("RNA" %in% names(merged_obj@assays)) {
    counts_layers <- grep("^counts", Layers(merged_obj[["RNA"]]), value = TRUE)
    if (length(counts_layers) > 1) {
      cat("Seurat v5: detected multiple count layers in merged_obj:\n  ",
          paste(counts_layers, collapse = ", "), "\n", sep = "")
      cat("Joining layers to create unified 'counts' layer...\n")
      merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])
      cat("Layers after JoinLayers():\n  ",
          paste(Layers(merged_obj[["RNA"]]), collapse = ", "), "\n\n", sep = "")
    }
  }
}

cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ==============================================================================
# Filtering Configuration
# ==============================================================================
cat("--- Filtering Configuration ---\n")

# Get filtering flags with defaults
skip_all <- isTRUE(params$skip_all_filtering)
apply_qc <- isTRUE(params$apply_qc_filtering) %||% TRUE
filter_min_feat <- isTRUE(params$filter_by_min_features) %||% TRUE
filter_max_feat <- isTRUE(params$filter_by_max_features) %||% TRUE
filter_mt <- isTRUE(params$filter_by_percent_mt) %||% TRUE
filter_genes <- isTRUE(params$filter_genes_by_min_cells) %||% TRUE
filter_doublets <- isTRUE(params$filter_doublets)
filter_hb <- isTRUE(params$filter_hemoglobin)

# Print configuration
if (skip_all) {
  cat("\n>>> SKIP_ALL_FILTERING = TRUE <<<\n")
  cat(">>> ALL FILTERING STEPS WILL BE BYPASSED <<<\n\n")
} else {
  cat("Standard QC filtering (apply_qc_filtering):", apply_qc, "\n")
  if (apply_qc) {
    cat("  filter_by_min_features:", filter_min_feat, 
        if(filter_min_feat) paste0(" (threshold: ", params$min_features, ")") else "", "\n")
    cat("  filter_by_max_features:", filter_max_feat,
        if(filter_max_feat) paste0(" (threshold: ", params$max_features, ")") else "", "\n")
    cat("  filter_by_percent_mt:", filter_mt,
        if(filter_mt) paste0(" (threshold: ", params$max_percent_mt, "%)") else "", "\n")
  }
  cat("Doublet filtering (filter_doublets):", filter_doublets, "\n")
  if (filter_doublets) {
    cat("  doublet_column:", params$doublet_column %||% "doublet_consensus", "\n")
    cat("  doublet_value_to_remove:", params$doublet_value_to_remove %||% "Doublet", "\n")
  }
  cat("Hemoglobin filtering (filter_hemoglobin):", filter_hb, "\n")
  if (filter_hb) {
    cat("  max_percent_hb:", params$max_percent_hb %||% 2, "\n")
  }
  cat("Gene filtering (filter_genes_by_min_cells):", filter_genes,
      if(filter_genes) paste0(" (threshold: ", params$min_cells_per_gene, " cells)") else "", "\n")
}
cat("\n")

# Store hemoglobin threshold for later use
max_percent_hb <- params$max_percent_hb %||% 2

# ==============================================================================
# Pre-filtering summary
# ==============================================================================
cat("--- Pre-filtering Summary ---\n")

pre_filter_cells <- ncol(merged_obj)
pre_filter_genes <- nrow(merged_obj)

cat("Total cells before filtering:", pre_filter_cells, "\n")
cat("Total genes before filtering:", pre_filter_genes, "\n\n")

# Per-sample summary
pre_summary <- merged_obj@meta.data %>%
  group_by(sample_name) %>%
  summarise(
    n_cells = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_pct_mt = median(percent.mt, na.rm = TRUE),
    .groups = "drop"
  )

cat("Pre-filtering cells per sample:\n")
print(as.data.frame(pre_summary))
cat("\n")

# ==============================================================================
# Validate QC columns
# ==============================================================================
cat("--- Validating QC Columns ---\n")

required_cols <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
missing_cols <- setdiff(required_cols, colnames(merged_obj@meta.data))

if (length(missing_cols) > 0) {
  cat("WARNING: Missing QC columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("Attempting to compute missing columns...\n")

  counts_mat <- GetAssayData(merged_obj, assay = "RNA", layer = "counts")

  if ("nFeature_RNA" %in% missing_cols) {
    merged_obj$nFeature_RNA <- Matrix::colSums(counts_mat > 0)
    cat("  Computed nFeature_RNA\n")
  }
  if ("nCount_RNA" %in% missing_cols) {
    merged_obj$nCount_RNA <- Matrix::colSums(counts_mat)
    cat("  Computed nCount_RNA\n")
  }
  if ("percent.mt" %in% missing_cols) {
    mt_genes <- grep("^MT-|^mt-", rownames(merged_obj), value = TRUE)
    if (length(mt_genes) > 0) {
      mt_counts <- Matrix::colSums(counts_mat[mt_genes, , drop = FALSE])
      total_counts <- Matrix::colSums(counts_mat)
      merged_obj$percent.mt <- (mt_counts / total_counts) * 100
      cat("  Computed percent.mt from", length(mt_genes), "MT genes\n")
    } else {
      merged_obj$percent.mt <- 0
      cat("  WARNING: No MT genes found, setting percent.mt to 0\n")
    }
  }
}

cat("QC columns validated.\n\n")

# ==============================================================================
# Validate hemoglobin column (if filtering enabled)
# ==============================================================================
if (filter_hb && !skip_all) {
  cat("--- Validating Hemoglobin Column ---\n")

  if (!"percent.hb" %in% colnames(merged_obj@meta.data)) {
    cat("WARNING: percent.hb column not found. Computing now...\n")

    hb_pattern <- params$hemoglobin_pattern %||% "^HB[AB]-"
    counts_mat <- GetAssayData(merged_obj, assay = "RNA", layer = "counts")

    all_hb_matches <- grep(hb_pattern, rownames(merged_obj), value = TRUE, ignore.case = TRUE)
    hb_genes <- all_hb_matches[!grepl("HBEGF|Hbegf", all_hb_matches)]

    if (length(hb_genes) > 0) {
      hb_counts <- Matrix::colSums(counts_mat[hb_genes, , drop = FALSE])
      total_counts <- Matrix::colSums(counts_mat)
      merged_obj$percent.hb <- (hb_counts / total_counts) * 100
      merged_obj$percent.hb[is.na(merged_obj$percent.hb)] <- 0
      cat("  Computed percent.hb from", length(hb_genes), "genes:", paste(hb_genes, collapse = ", "), "\n")
    } else {
      merged_obj$percent.hb <- 0
      cat("  WARNING: No hemoglobin genes found with pattern:", hb_pattern, "\n")
    }
  } else {
    cat("  percent.hb column already exists\n")
  }

  cat("  Cells with >", max_percent_hb, "% Hb: ",
      sum(merged_obj$percent.hb > max_percent_hb, na.rm = TRUE), "\n", sep = "")
  cat("\n")
}

# ==============================================================================
# Validate doublet column (if filtering enabled)
# ==============================================================================
if (filter_doublets && !skip_all) {
  cat("--- Validating Doublet Column ---\n")
  
  doublet_col <- params$doublet_column %||% "doublet_consensus"
  doublet_val <- params$doublet_value_to_remove %||% "Doublet"
  
  if (!doublet_col %in% colnames(merged_obj@meta.data)) {
    cat("  WARNING: Doublet column '", doublet_col, "' not found in metadata.\n", sep = "")
    cat("  Available columns containing 'doublet':\n")
    doublet_cols <- grep("doublet", colnames(merged_obj@meta.data), value = TRUE, ignore.case = TRUE)
    if (length(doublet_cols) > 0) {
      cat("    ", paste(doublet_cols, collapse = ", "), "\n")
    } else {
      cat("    None found\n")
    }
    cat("  Doublet filtering will be SKIPPED.\n")
    filter_doublets <- FALSE
  } else {
    doublet_counts <- table(merged_obj@meta.data[[doublet_col]])
    cat("  Doublet column:", doublet_col, "\n")
    cat("  Values found:\n")
    print(doublet_counts)
    n_doublets <- sum(merged_obj@meta.data[[doublet_col]] == doublet_val, na.rm = TRUE)
    cat("  Cells to remove ('", doublet_val, "'): ", n_doublets, "\n", sep = "")
  }
  cat("\n")
}

# ==============================================================================
# Create pre-filtering QC plots
# ==============================================================================
cat("--- Creating Pre-filtering QC Plots ---\n")

# Histogram of nFeature
p_nfeature_hist <- ggplot(merged_obj@meta.data, aes(x = nFeature_RNA, fill = sample_name)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  geom_vline(xintercept = c(params$min_features, params$max_features),
             linetype = "dashed", color = "red") +
  facet_wrap(~ sample_name, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of Genes per Cell (Pre-filtering)",
       subtitle = paste("Red lines: min =", params$min_features, ", max =", params$max_features),
       x = "Number of Genes", y = "Count")

ggsave(file.path(output_dirs$qc, "02_nFeature_histogram_prefilter.png"),
       p_nfeature_hist, width = 12, height = 8, dpi = 300)

# Histogram of percent.mt
p_mt_hist <- ggplot(merged_obj@meta.data, aes(x = percent.mt, fill = sample_name)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  geom_vline(xintercept = params$max_percent_mt, linetype = "dashed", color = "red") +
  facet_wrap(~ sample_name, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of MT% per Cell (Pre-filtering)",
       subtitle = paste("Red line: max =", params$max_percent_mt),
       x = "Percent Mitochondrial", y = "Count")

ggsave(file.path(output_dirs$qc, "02_pctMT_histogram_prefilter.png"),
       p_mt_hist, width = 12, height = 8, dpi = 300)

# Histogram of percent.hb (if column exists)
if ("percent.hb" %in% colnames(merged_obj@meta.data)) {
  p_hb_hist <- ggplot(merged_obj@meta.data, aes(x = percent.hb, fill = sample_name)) +
    geom_histogram(bins = 50, alpha = 0.7) +
    geom_vline(xintercept = max_percent_hb, linetype = "dashed", color = "red") +
    facet_wrap(~ sample_name, scales = "free_y") +
    theme_minimal() +
    scale_x_continuous(limits = c(0, min(max(merged_obj$percent.hb, na.rm = TRUE) * 1.1, 50))) +
    labs(title = "Distribution of Hemoglobin % per Cell (Pre-filtering)",
         subtitle = paste("Red line: max =", max_percent_hb),
         x = "Percent Hemoglobin", y = "Count")

  ggsave(file.path(output_dirs$qc, "02_pctHb_histogram_prefilter.png"),
         p_hb_hist, width = 12, height = 8, dpi = 300)
}

cat("Pre-filtering plots saved.\n\n")

# ==============================================================================
# Apply Cell Filtering
# ==============================================================================
cat("--- Applying Cell Filtering ---\n")

# Initialize: all cells pass by default
cells_to_keep <- rep(TRUE, ncol(merged_obj))
names(cells_to_keep) <- colnames(merged_obj)

# Initialize tracking columns
merged_obj$pass_min_features <- TRUE
merged_obj$pass_max_features <- TRUE
merged_obj$pass_mt_filter <- TRUE
merged_obj$pass_doublet_filter <- TRUE
merged_obj$pass_hb_filter <- TRUE

# ==============================================================================
# MASTER SKIP SWITCH
# ==============================================================================
if (skip_all) {
  cat("\n>>> SKIPPING ALL CELL FILTERING (skip_all_filtering = TRUE) <<<\n\n")
  
} else {
  
  # ==========================================================================
  # DOUBLET FILTERING (applied first, before QC)
  # ==========================================================================
  if (filter_doublets) {
    cat("\nDoublet filtering:\n")
    
    doublet_col <- params$doublet_column %||% "doublet_consensus"
    
    if (isTRUE(params$use_doublet_vote_threshold)) {
      # Filter by vote count
      vote_col <- params$doublet_vote_column %||% "doublet_votes"
      vote_threshold <- params$doublet_vote_threshold %||% 2
      
      if (vote_col %in% colnames(merged_obj@meta.data)) {
        cells_pass_doublet <- merged_obj@meta.data[[vote_col]] < vote_threshold
        cells_pass_doublet[is.na(cells_pass_doublet)] <- TRUE
        cat("  Using vote threshold: <", vote_threshold, " votes\n")
        cat("  Cells passing:", sum(cells_pass_doublet), "\n")
        cat("  Cells failing (doublets):", sum(!cells_pass_doublet), "\n")
      } else {
        cat("  WARNING: Vote column '", vote_col, "' not found. Skipping.\n", sep = "")
        cells_pass_doublet <- TRUE
      }
    } else {
      # Filter by consensus label
      doublet_val <- params$doublet_value_to_remove %||% "Doublet"
      
      if (doublet_col %in% colnames(merged_obj@meta.data)) {
        cells_pass_doublet <- merged_obj@meta.data[[doublet_col]] != doublet_val
        cells_pass_doublet[is.na(cells_pass_doublet)] <- TRUE
        cat("  Using consensus column: '", doublet_col, "'\n", sep = "")
        cat("  Removing cells with value: '", doublet_val, "'\n", sep = "")
        cat("  Cells passing:", sum(cells_pass_doublet), "\n")
        cat("  Cells failing (doublets):", sum(!cells_pass_doublet), "\n")
      } else {
        cat("  WARNING: Doublet column '", doublet_col, "' not found. Skipping.\n", sep = "")
        cells_pass_doublet <- TRUE
      }
    }
    
    merged_obj$pass_doublet_filter <- cells_pass_doublet
    cells_to_keep <- cells_to_keep & cells_pass_doublet
  } else {
    cat("\nDoublet filtering: SKIPPED (filter_doublets = FALSE)\n")
  }
  
  # ==========================================================================
  # STANDARD QC FILTERING
  # ==========================================================================
  if (apply_qc) {
    cat("\nStandard QC filtering:\n")
    
    # Min features filter
    if (filter_min_feat) {
      cells_pass_min <- merged_obj$nFeature_RNA >= params$min_features
      merged_obj$pass_min_features <- cells_pass_min
      cells_to_keep <- cells_to_keep & cells_pass_min
      cat("  min_features (>=", params$min_features, "): ", 
          sum(cells_pass_min), " pass, ", sum(!cells_pass_min), " fail\n", sep = "")
    } else {
      cat("  min_features: SKIPPED\n")
    }
    
    # Max features filter
    if (filter_max_feat) {
      cells_pass_max <- merged_obj$nFeature_RNA <= params$max_features
      merged_obj$pass_max_features <- cells_pass_max
      cells_to_keep <- cells_to_keep & cells_pass_max
      cat("  max_features (<=", params$max_features, "): ", 
          sum(cells_pass_max), " pass, ", sum(!cells_pass_max), " fail\n", sep = "")
    } else {
      cat("  max_features: SKIPPED\n")
    }
    
    # MT percentage filter
    if (filter_mt) {
      cells_pass_mt <- merged_obj$percent.mt <= params$max_percent_mt
      cells_pass_mt[is.na(cells_pass_mt)] <- TRUE
      merged_obj$pass_mt_filter <- cells_pass_mt
      cells_to_keep <- cells_to_keep & cells_pass_mt
      cat("  max_percent_mt (<=", params$max_percent_mt, "): ", 
          sum(cells_pass_mt), " pass, ", sum(!cells_pass_mt), " fail\n", sep = "")
    } else {
      cat("  max_percent_mt: SKIPPED\n")
    }
    
  } else {
    cat("\nStandard QC filtering: SKIPPED (apply_qc_filtering = FALSE)\n")
  }
  
  # ==========================================================================
  # HEMOGLOBIN FILTERING
  # ==========================================================================
  if (filter_hb && "percent.hb" %in% colnames(merged_obj@meta.data)) {
    cat("\nHemoglobin filtering:\n")
    cells_pass_hb <- merged_obj$percent.hb <= max_percent_hb
    cells_pass_hb[is.na(cells_pass_hb)] <- TRUE
    merged_obj$pass_hb_filter <- cells_pass_hb
    cells_to_keep <- cells_to_keep & cells_pass_hb
    cat("  max_percent_hb (<=", max_percent_hb, "): ", 
        sum(cells_pass_hb), " pass, ", sum(!cells_pass_hb), " fail\n", sep = "")
  } else if (filter_hb) {
    cat("\nHemoglobin filtering: SKIPPED (percent.hb column not found)\n")
  } else {
    cat("\nHemoglobin filtering: SKIPPED (filter_hemoglobin = FALSE)\n")
  }
}

# Store combined filter result
merged_obj$pass_all_qc <- cells_to_keep

cat("\n  >>> TOTAL CELLS PASSING ALL FILTERS:", sum(cells_to_keep), "<<<\n\n")

# ==============================================================================
# Apply cell filter
# ==============================================================================
filtered_obj <- subset(merged_obj, cells = colnames(merged_obj)[cells_to_keep])

cat("Cells after filtering:", ncol(filtered_obj), "\n")
cat("Cells removed:", pre_filter_cells - ncol(filtered_obj),
    "(", round((pre_filter_cells - ncol(filtered_obj)) / pre_filter_cells * 100, 1), "%)\n\n")

# ==============================================================================
# Apply Gene Filtering
# ==============================================================================
cat("--- Applying Gene Filtering ---\n")

pre_gene_count <- nrow(filtered_obj)

if (!skip_all && filter_genes) {
  counts_mat <- GetAssayData(filtered_obj, assay = "RNA", layer = "counts")
  cells_per_gene <- Matrix::rowSums(counts_mat > 0)
  genes_to_keep <- names(cells_per_gene)[cells_per_gene >= params$min_cells_per_gene]

  cat("Genes before filtering:", pre_gene_count, "\n")
  cat("Genes with >=", params$min_cells_per_gene, "cells:", length(genes_to_keep), "\n")

  filtered_obj <- subset(filtered_obj, features = genes_to_keep)

  cat("Genes after filtering:", nrow(filtered_obj), "\n")
  cat("Genes removed:", pre_gene_count - nrow(filtered_obj), "\n\n")
} else {
  cat("Gene filtering: SKIPPED\n")
  if (skip_all) {
    cat("  Reason: skip_all_filtering = TRUE\n")
  } else {
    cat("  Reason: filter_genes_by_min_cells = FALSE\n")
  }
  cat("Genes retained:", nrow(filtered_obj), "\n\n")
}

# ==============================================================================
# Post-filtering summary
# ==============================================================================
cat("--- Post-filtering Summary ---\n")

post_filter_cells <- ncol(filtered_obj)
post_filter_genes <- nrow(filtered_obj)

# Per-sample filtering summary
filtering_summary <- merged_obj@meta.data %>%
  mutate(passed_filter = pass_all_qc) %>%
  group_by(sample_name, sex, batch) %>%
  summarise(
    cells_before = n(),
    cells_after = sum(passed_filter),
    cells_removed = sum(!passed_filter),
    pct_removed = round(sum(!passed_filter) / n() * 100, 1),
    .groups = "drop"
  )

cat("Filtering by sample:\n")
print(as.data.frame(filtering_summary))
cat("\n")

# Save filtering summary
write.csv(filtering_summary,
          file.path(output_dirs$tables, "02_filtering_summary.csv"),
          row.names = FALSE)

# Detailed filtering breakdown
filtering_breakdown <- merged_obj@meta.data %>%
  group_by(sample_name) %>%
  summarise(
    total_cells = n(),
    failed_doublet = sum(!pass_doublet_filter),
    failed_min_features = sum(!pass_min_features),
    failed_max_features = sum(!pass_max_features),
    failed_mt = sum(!pass_mt_filter),
    failed_hb = sum(!pass_hb_filter),
    passed_all = sum(pass_all_qc),
    .groups = "drop"
  )

cat("Detailed filtering breakdown:\n")
print(as.data.frame(filtering_breakdown))
cat("\n")

write.csv(filtering_breakdown,
          file.path(output_dirs$tables, "02_filtering_breakdown.csv"),
          row.names = FALSE)

# Post-filter per-sample summary
post_summary <- filtered_obj@meta.data %>%
  group_by(sample_name) %>%
  summarise(
    n_cells = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_pct_mt = median(percent.mt, na.rm = TRUE),
    .groups = "drop"
  )

cat("Post-filtering cells per sample:\n")
print(as.data.frame(post_summary))
cat("\n")

# ==============================================================================
# Create Post-filtering QC Plots
# ==============================================================================
cat("--- Creating Post-filtering QC Plots ---\n")

# Violin plots by sample
p_violin_post <- VlnPlot(
  filtered_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_name",
  pt.size = 0,
  ncol = 3
) +
  plot_annotation(title = "QC Metrics by Sample (Post-filtering)")

ggsave(file.path(output_dirs$qc, "02_qc_violin_postfilter_by_sample.png"),
       p_violin_post, width = 14, height = 5, dpi = 300)

# Violin plots by sex
p_violin_sex <- VlnPlot(
  filtered_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sex",
  pt.size = 0,
  ncol = 3
) +
  plot_annotation(title = "QC Metrics by Sex (Post-filtering)")

ggsave(file.path(output_dirs$qc, "02_qc_violin_postfilter_by_sex.png"),
       p_violin_sex, width = 10, height = 5, dpi = 300)

# Violin plots by batch
p_violin_batch <- VlnPlot(
  filtered_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "batch",
  pt.size = 0,
  ncol = 3
) +
  plot_annotation(title = "QC Metrics by Batch (Post-filtering)")

ggsave(file.path(output_dirs$qc, "02_qc_violin_postfilter_by_batch.png"),
       p_violin_batch, width = 10, height = 5, dpi = 300)

# Hemoglobin violin plot (post-filtering, if applicable)
if ("percent.hb" %in% colnames(filtered_obj@meta.data)) {
  p_hb_post <- VlnPlot(
    filtered_obj,
    features = "percent.hb",
    group.by = "sample_name",
    pt.size = 0
  ) +
    labs(title = "Hemoglobin Content by Sample (Post-filtering)")

  ggsave(file.path(output_dirs$qc, "02_hemoglobin_violin_postfilter.png"),
         p_hb_post, width = 10, height = 5, dpi = 300)
}

# Filtering comparison barplot
filtering_long <- filtering_summary %>%
  select(sample_name, cells_before, cells_after) %>%
  tidyr::pivot_longer(-sample_name, names_to = "status", values_to = "cells") %>%
  mutate(status = factor(status, levels = c("cells_before", "cells_after"),
                         labels = c("Before", "After")))

p_filter_comparison <- ggplot(filtering_long, aes(x = sample_name, y = cells, fill = status)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Before" = "#999999", "After" = "#2166AC")) +
  labs(title = "Cell Counts Before and After QC Filtering",
       x = "Sample", y = "Number of Cells", fill = "Status")

ggsave(file.path(output_dirs$qc, "02_filtering_comparison.png"),
       p_filter_comparison, width = 10, height = 6, dpi = 300)

# Scatter plot showing filtered cells
if (pre_filter_cells < 100000) {
  p_scatter_filtered <- ggplot(merged_obj@meta.data, aes(x = nFeature_RNA, y = percent.mt, color = pass_all_qc)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_hline(yintercept = params$max_percent_mt, linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(params$min_features, params$max_features), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                       labels = c("TRUE" = "Passed", "FALSE" = "Filtered")) +
    facet_wrap(~ sample_name) +
    theme_minimal() +
    labs(title = "Cells Filtered by QC Thresholds",
         x = "nFeature_RNA",
         y = "Percent MT",
         color = "QC Status")

  ggsave(file.path(output_dirs$qc, "02_scatter_filtered_cells.png"),
         p_scatter_filtered, width = 12, height = 10, dpi = 300)
}

cat("Post-filtering plots saved.\n\n")

# ==============================================================================
# Split object by sample for downstream processing
# ==============================================================================
cat("--- Preparing for Normalization ---\n")

sample_names <- unique(filtered_obj$sample_name)
sample_objects <- list()

for (samp in sample_names) {
  cells <- colnames(filtered_obj)[filtered_obj$sample_name == samp]
  sample_objects[[samp]] <- subset(filtered_obj, cells = cells)
  cat("  ", samp, ":", ncol(sample_objects[[samp]]), "cells\n")
}

cat("\nIndividual sample objects created.\n")

# ==============================================================================
# Save objects
# ==============================================================================
cat("\n--- Saving Objects ---\n")

qc_data_file <- file.path(output_dirs$objects, "02_qc_data.RData")
save(filtered_obj, sample_objects, filtering_summary, filtering_breakdown, file = qc_data_file)
cat("Saved:", qc_data_file, "\n")

qc_rds_file <- file.path(output_dirs$objects, "02_qc_filtered_object.rds")
saveRDS(filtered_obj, qc_rds_file)
cat("Saved:", qc_rds_file, "\n")

# ==============================================================================
# Write README
# ==============================================================================
readme_content <- paste0(
  "================================================================================\n",
  "MODULE 02: QC VALIDATION AND FILTERING\n",
  "================================================================================\n\n",
  "Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
  "FILTERING CONFIGURATION:\n",
  "  skip_all_filtering: ", skip_all, "\n",
  "  apply_qc_filtering: ", apply_qc, "\n",
  "  filter_by_min_features: ", filter_min_feat, "\n",
  "  filter_by_max_features: ", filter_max_feat, "\n",
  "  filter_by_percent_mt: ", filter_mt, "\n",
  "  filter_doublets: ", filter_doublets, "\n",
  "  filter_hemoglobin: ", filter_hb, "\n",
  "  filter_genes_by_min_cells: ", filter_genes, "\n\n",
  "QC THRESHOLDS:\n",
  "  min_features: ", params$min_features, "\n",
  "  max_features: ", params$max_features, "\n",
  "  max_percent_mt: ", params$max_percent_mt, "\n",
  "  min_cells_per_gene: ", params$min_cells_per_gene, "\n"
)

if (filter_hb) {
  readme_content <- paste0(readme_content,
    "  max_percent_hb: ", max_percent_hb, "\n",
    "  hemoglobin_pattern: ", params$hemoglobin_pattern %||% "^HB[AB]-", "\n"
  )
}

if (filter_doublets) {
  readme_content <- paste0(readme_content,
    "  doublet_column: ", params$doublet_column %||% "doublet_consensus", "\n",
    "  doublet_value_to_remove: ", params$doublet_value_to_remove %||% "Doublet", "\n"
  )
}

readme_content <- paste0(readme_content,
  "\nFILTERING SUMMARY:\n",
  "  Cells before: ", pre_filter_cells, "\n",
  "  Cells after: ", post_filter_cells, "\n",
  "  Cells removed: ", pre_filter_cells - post_filter_cells,
  " (", round((pre_filter_cells - post_filter_cells) / pre_filter_cells * 100, 1), "%)\n\n",
  "  Genes before: ", pre_filter_genes, "\n",
  "  Genes after: ", post_filter_genes, "\n",
  "  Genes removed: ", pre_filter_genes - post_filter_genes, "\n\n",
  "OUTPUT FILES:\n",
  "  02_qc_data.RData: Filtered object + sample objects + summaries\n",
  "  02_qc_filtered_object.rds: Filtered merged object\n",
  "  02_filtering_summary.csv: Per-sample filtering statistics\n",
  "  02_filtering_breakdown.csv: Detailed breakdown by filter type\n"
)

writeLines(readme_content, file.path(output_dirs$qc, "README_02_qc_validation.txt"))

# ==============================================================================
# Final Summary
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 02 COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Filtering Configuration:\n")
cat("  skip_all_filtering:", skip_all, "\n")
cat("  apply_qc_filtering:", apply_qc, "\n")
cat("  filter_doublets:", filter_doublets, "\n")
cat("  filter_hemoglobin:", filter_hb, "\n")
cat("  filter_genes_by_min_cells:", filter_genes, "\n")

cat("\nSummary:\n")
cat("  Samples:", length(sample_names), "\n")
cat("  Cells before filtering:", pre_filter_cells, "\n")
cat("  Cells after filtering:", post_filter_cells, "\n")
cat("  Genes before filtering:", pre_filter_genes, "\n")
cat("  Genes after filtering:", post_filter_genes, "\n")
cat("  Percent cells removed:", round((pre_filter_cells - post_filter_cells) / pre_filter_cells * 100, 1), "%\n")

if (filter_doublets && !skip_all) {
  cat("\nDoublet filtering:\n")
  cat("  Cells removed as doublets:", sum(!merged_obj$pass_doublet_filter), "\n")
}

if (filter_hb && !skip_all) {
  cat("\nHemoglobin filtering:\n")
  cat("  Threshold:", max_percent_hb, "%\n")
  cat("  Cells removed by Hb filter:", sum(!merged_obj$pass_hb_filter, na.rm = TRUE), "\n")
}

cat("\nCells per sex (post-filtering):\n")
print(table(filtered_obj$sex))

cat("\nCells per batch (post-filtering):\n")
print(table(filtered_obj$batch))

cat("\n>>> MODULE 02 COMPLETE <<<\n")
