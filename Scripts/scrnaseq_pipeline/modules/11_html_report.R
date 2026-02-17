#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 11: HTML REPORT GENERATOR (MULTI-SAMPLE PIPELINE)
# ==============================================================================
#
# This module generates a comprehensive interactive HTML report including:
# - Clickable navigation menu for all sections
# - Pipeline processing summary with step-by-step explanations
# - Benchmarking results with parameter explanations
# - Benchmarking visualization plots (scIB-style heatmaps)
# - All generated plots organized by analysis stage
# - Dark/light mode toggle
# - Print-friendly formatting
#
# INPUT: All previous module outputs from ${DOWNSTREAM_DIR}/10_Downstream_Analysis_${VENTRICLE}/
# OUTPUT: Interactive HTML report saved to same directory
#
# OUTPUT LOCATION:
#   ${DOWNSTREAM_DIR}/10_Downstream_Analysis_${VENTRICLE}/Analysis_Report_MultiSample.html
#
# UPDATES:
# - 2026-01-06: Fixed multibyte character encoding issues (replaced Unicode with HTML entities)
# - 2026-01-06: Updated for dataset-specific directory structure
# - 2026-01-06: Path resolution via environment variables (DOWNSTREAM_DIR)
# - 2026-01-04: Added scIB-style benchmarking visualization plots
# - 2026-01-03: Added clickable navigation sidebar
# - 2026-01-03: Added comprehensive pipeline processing summary
# - 2026-01-03: Added benchmarking parameter explanations
# - 2026-01-03: Enhanced plot gallery with all generated figures
# - 2026-01-03: Improved Seurat v5 compatibility
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 11: HTML REPORT GENERATOR\n")
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
  library(ggplot2)
  library(tidyr)
})

# Check for base64enc
if (!requireNamespace("base64enc", quietly = TRUE)) {
  stop("Package 'base64enc' is required. Install with: install.packages('base64enc')")
}
library(base64enc)

# Get output base from params (resolves from DOWNSTREAM_DIR + ventricle)
out_base <- params$out_root

cat("Output directory:", out_base, "\n")

# Load pipeline environment
tryCatch({
  load(file.path(out_base, "objects", "pipeline_environment.RData"))
  cat("Loaded pipeline environment\n")
}, error = function(e) {
  cat("Note: Could not load pipeline environment\n")
})

# Load all module data files
data_files <- c(
  "03_normalization_data.RData",
  "04_integration_data.RData",
  "05_choir_data.RData",
  "06_scice_data.RData",
  "07_leiden_data.RData",
  "08_de_data.RData",
  "09_visualization_data.RData"
)

for (f in data_files) {
  fp <- file.path(out_base, "objects", f)
  if (file.exists(fp)) {
    tryCatch({
      load(fp)
      cat("Loaded:", f, "\n")
    }, error = function(e) {
      cat("Could not load:", f, "\n")
    })
  }
}

# Try to load final object if not already loaded
if (!exists("clustered_obj") || is.null(clustered_obj)) {
  obj_files <- c(
    "final_clustered_object.rds",
    "07_final_object.rds",
    "scice_subclustered_object.rds",
    "leiden_clustered_object.rds",
    "choir_clustered_object.rds"
  )
  for (of in obj_files) {
    fp <- file.path(out_base, "objects", of)
    if (file.exists(fp)) {
      tryCatch({
        clustered_obj <- readRDS(fp)
        cat("Loaded object from:", of, "\n")
        break
      }, error = function(e) NULL)
    }
  }
}

# ==============================================================================
# Helper Functions
# ==============================================================================

# HTML escape special characters
html_escape <- function(text) {
  text <- gsub("&", "&amp;", as.character(text))
  text <- gsub("<", "&lt;", text)
  text <- gsub(">", "&gt;", text)
  text <- gsub('"', "&quot;", text)
  return(text)
}

# Convert image to base64
img_to_base64 <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch({
    base64encode(path)
  }, error = function(e) NULL)
}

# Embed image with base64
embed_image <- function(path, alt = "", max_width = "100%", clickable = TRUE) {
  b64 <- img_to_base64(path)
  if (is.null(b64)) {
    return(paste0('<div class="missing-img">Image not found: ', html_escape(basename(path)), '</div>'))
  }

  mime_type <- if (grepl("\\.png$", path, ignore.case = TRUE)) "image/png" else "image/jpeg"

  onclick <- if (clickable) ' onclick="openModal(this.src)"' else ""

  paste0('<img src="data:', mime_type, ';base64,', b64,
         '" alt="', html_escape(alt),
         '" style="max-width:', max_width, ';width:100%;cursor:',
         if (clickable) 'zoom-in' else 'default', ';"', onclick, '/>')
}

# Find PNG files in directory
find_pngs <- function(dir_path, recursive = FALSE) {
  if (!dir.exists(dir_path)) return(character(0))
  list.files(dir_path, pattern = "\\.png$", full.names = TRUE, recursive = recursive)
}

# Create HTML table from data frame
create_table <- function(df, max_rows = 50, class = "data-table") {
  if (is.null(df) || nrow(df) == 0) {
    return('<p class="no-data">No data available</p>')
  }

  df <- head(df, max_rows)

  # Header
  header <- paste0('<tr>', paste0('<th>', html_escape(colnames(df)), '</th>', collapse = ''), '</tr>')


  # Rows
  rows <- apply(df, 1, function(row) {
    cells <- sapply(row, function(cell) {
      val <- as.character(cell)
      # Format numbers
      if (!is.na(suppressWarnings(as.numeric(val)))) {
        num <- as.numeric(val)
        if (abs(num) < 0.001 && num != 0) {
          val <- sprintf("%.2e", num)
        } else if (abs(num) < 1) {
          val <- sprintf("%.4f", num)
        } else if (abs(num) > 1000) {
          val <- format(round(num, 1), big.mark = ",")
        }
      }
      paste0('<td>', html_escape(val), '</td>')
    })
    paste0('<tr>', paste(cells, collapse = ''), '</tr>')
  })

  paste0('<div class="table-wrapper"><table class="', class, '">',
         '<thead>', header, '</thead>',
         '<tbody>', paste(rows, collapse = ''), '</tbody>',
         '</table></div>',
         if (nrow(df) >= max_rows) paste0('<p class="table-note">Showing first ', max_rows, ' rows</p>') else '')
}

# Create image gallery
create_gallery <- function(paths, ncol = 2, titles = NULL) {
  if (length(paths) == 0) {
    return('<p class="no-data">No images available</p>')
  }

  if (is.null(titles)) {
    titles <- tools::file_path_sans_ext(basename(paths))
  }

  items <- mapply(function(path, title) {
    paste0('<div class="gallery-item">',
           '<div class="gallery-title">', html_escape(title), '</div>',
           '<div class="gallery-img">', embed_image(path), '</div>',
           '</div>')
  }, paths, titles, SIMPLIFY = FALSE)

  paste0('<div class="gallery" style="grid-template-columns:repeat(', ncol, ',1fr);">',
         paste(unlist(items), collapse = ''),
         '</div>')
}

# Create collapsible section
create_section <- function(id, title, content, icon = "") {
  paste0('<section id="', id, '" class="report-section">',
         '<div class="section-header" onclick="toggleSection(\'', id, '\')">',
         '<h2>', if (icon != "") paste0('<span class="icon">', icon, '</span> ') else '',
         html_escape(title), '</h2>',
         '<span class="toggle-icon">&#9660;</span>',
         '</div>',
         '<div class="section-content">', content, '</div>',
         '</section>')
}

# ==============================================================================
# Benchmarking Visualization Functions (scIB-style)
# ==============================================================================

#' Create scIB-style benchmarking heatmap plot
#'
#' @param bench_df Data frame with benchmarking results
#' @param title Plot title
#' @param output_path Path to save the PNG file
#' @param analysis_type "normalization" or "integration"
#' @return Path to saved plot
create_benchmark_plot <- function(bench_df, title = "Benchmarking Results",
                                  output_path = NULL, analysis_type = "integration") {

  if (is.null(bench_df) || nrow(bench_df) == 0) {
    cat("  No benchmarking data available for plotting\n")
    return(NULL)
  }

  # Identify method column
  method_col <- intersect(c("method", "Method", "normalization_method"), colnames(bench_df))[1]
  if (is.na(method_col)) {
    cat("  Could not find method column in benchmarking data\n")
    return(NULL)
  }

  # Identify available metric columns
  batch_correction_cols <- intersect(c("batch_variance", "batch_asw", "lisi", "kBET", "iLISI",
                                       "silhouette_batch", "graph_connectivity", "pcr_comparison"),
                                     colnames(bench_df))

  bio_conservation_cols <- intersect(c("cLISI", "silhouette_label", "isolated_labels",
                                       "kmeans_nmi", "kmeans_ari", "bio_conservation"),
                                     colnames(bench_df))

  # If no bio conservation metrics, use what's available
  if (length(bio_conservation_cols) == 0) {
    bio_conservation_cols <- c()
  }

  all_metric_cols <- c(batch_correction_cols, bio_conservation_cols)

  if (length(all_metric_cols) == 0) {
    cat("  No recognized metric columns found\n")
    return(NULL)
  }

  cat("  Found metric columns:", paste(all_metric_cols, collapse = ", "), "\n")

  # Prepare data for plotting
  plot_df <- bench_df[, c(method_col, all_metric_cols), drop = FALSE]
  colnames(plot_df)[1] <- "method"

  # Convert to numeric and handle NAs
  for (col in all_metric_cols) {
    plot_df[[col]] <- as.numeric(plot_df[[col]])
  }

  # Normalize metrics to 0-1 range
  normalize_metric <- function(x, invert = FALSE) {
    if (all(is.na(x))) return(rep(0.5, length(x)))
    x_min <- min(x, na.rm = TRUE)
    x_max <- max(x, na.rm = TRUE)
    if (x_max == x_min) return(rep(0.5, length(x)))
    normalized <- (x - x_min) / (x_max - x_min)
    if (invert) normalized <- 1 - normalized
    return(normalized)
  }

  # Define which metrics should be inverted (lower = better)
  invert_metrics <- c("batch_variance", "batch_asw", "silhouette_batch")

  # Create normalized data
  norm_df <- plot_df
  for (col in all_metric_cols) {
    should_invert <- col %in% invert_metrics
    norm_df[[paste0(col, "_norm")]] <- normalize_metric(plot_df[[col]], invert = should_invert)
  }

  # Calculate aggregate scores
  batch_norm_cols <- paste0(batch_correction_cols, "_norm")
  batch_norm_cols <- batch_norm_cols[batch_norm_cols %in% colnames(norm_df)]

  if (length(batch_norm_cols) > 0) {
    norm_df$batch_score <- rowMeans(norm_df[, batch_norm_cols, drop = FALSE], na.rm = TRUE)
  } else {
    norm_df$batch_score <- 0.5
  }

  if (length(bio_conservation_cols) > 0) {
    bio_norm_cols <- paste0(bio_conservation_cols, "_norm")
    bio_norm_cols <- bio_norm_cols[bio_norm_cols %in% colnames(norm_df)]
    if (length(bio_norm_cols) > 0) {
      norm_df$bio_score <- rowMeans(norm_df[, bio_norm_cols, drop = FALSE], na.rm = TRUE)
    } else {
      norm_df$bio_score <- 0.5
    }
  } else {
    norm_df$bio_score <- 0.5
  }

  # Overall score: weighted average (60% batch, 40% bio)
  norm_df$total_score <- 0.6 * norm_df$batch_score + 0.4 * norm_df$bio_score

  # Order by total score
  norm_df <- norm_df[order(norm_df$total_score, decreasing = TRUE), ]
  norm_df$method <- factor(norm_df$method, levels = rev(norm_df$method))

  # Pivot for ggplot
  metric_cols_for_pivot <- paste0(all_metric_cols, "_norm")
  metric_cols_for_pivot <- metric_cols_for_pivot[metric_cols_for_pivot %in% colnames(norm_df)]

  # Create the main heatmap data
  heat_data <- norm_df %>%
    select(method, all_of(metric_cols_for_pivot)) %>%
    pivot_longer(cols = -method, names_to = "metric", values_to = "value")

  # Clean metric names for display
  heat_data$metric <- gsub("_norm$", "", heat_data$metric)
  heat_data$metric <- gsub("_", " ", heat_data$metric)
  heat_data$metric <- tools::toTitleCase(heat_data$metric)

  # Determine metric category for coloring
  heat_data$category <- ifelse(
    tolower(gsub(" ", "_", heat_data$metric)) %in% tolower(batch_correction_cols),
    "Batch Correction",
    "Bio Conservation"
  )

  # Set factor levels for metrics
  batch_metrics_display <- tools::toTitleCase(gsub("_", " ", batch_correction_cols))
  bio_metrics_display <- tools::toTitleCase(gsub("_", " ", bio_conservation_cols))
  all_metrics_display <- c(batch_metrics_display, bio_metrics_display)
  heat_data$metric <- factor(heat_data$metric, levels = all_metrics_display)

  # Create aggregate score data
  agg_data <- norm_df %>%
    select(method, batch_score, bio_score, total_score) %>%
    pivot_longer(cols = -method, names_to = "score_type", values_to = "value")

  agg_data$score_type <- factor(agg_data$score_type,
                                levels = c("batch_score", "bio_score", "total_score"),
                                labels = c("Batch", "Bio", "Total"))

  # Create the heatmap plot
  p1 <- ggplot(heat_data, aes(x = metric, y = method)) +
    geom_point(aes(fill = value, color = category), shape = 21, size = 10, stroke = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 2.8, color = "white", fontface = "bold") +
    scale_fill_gradient2(low = "#f0f0f0", mid = "#7B68EE", high = "#2E0854",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Score\n(normalized)") +
    scale_color_manual(values = c("Batch Correction" = "#663399",
                                  "Bio Conservation" = "#228B22"),
                       name = "Category") +
    labs(title = title,
         subtitle = "Individual Metrics (higher = better after normalization)",
         x = "", y = "") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_text(size = 10, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))

  # Create aggregate score bar plot
  p2 <- ggplot(agg_data, aes(x = score_type, y = method, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3, color = "white", fontface = "bold") +
    scale_fill_gradient(low = "#87CEEB", high = "#00008B",
                        limits = c(0, 1), name = "Score") +
    labs(title = "Aggregate Scores",
         x = "", y = "") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 5)
    )

  # Combine plots if gridExtra is available
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined_plot <- gridExtra::grid.arrange(p1, p2, ncol = 2, widths = c(3, 1))
  } else {
    combined_plot <- p1
  }

  # Save plot
  if (!is.null(output_path)) {
    ggsave(output_path, plot = combined_plot, width = 14, height = max(6, nrow(norm_df) * 0.8),
           dpi = 150, bg = "white")
    cat("  Saved benchmarking plot to:", output_path, "\n")
    return(output_path)
  }

  return(combined_plot)
}

#' Create scIB-style benchmarking plot - COMPATIBLE WITH MODULE 04 OUTPUT
#'
#' Uses pre-computed normalized columns if available (from Module 04),
#' otherwise computes normalization on the fly.
#'
#' Purple = Batch correction, Green = Bio conservation
#'
#' @param bench_df Data frame with benchmarking results
#' @param title Plot title
#' @param output_path Path to save the PNG file
#' @param selection_mode Optional selection mode from Module 04
#' @param selection_criterion Optional selection criterion from Module 04
#' @return Path to saved plot
create_scib_style_plot <- function(bench_df, title = "Integration Benchmarking",
                                   output_path = NULL,
                                   selection_mode = NULL,
                                   selection_criterion = NULL) {

  if (is.null(bench_df) || nrow(bench_df) == 0) {
    cat("  No data available for scIB-style plot\n")
    return(NULL)
  }

  # Identify method column
 method_col <- intersect(c("method", "Method", "normalization_method"), colnames(bench_df))[1]
  if (is.na(method_col)) {
    cat("  Could not find method column\n")
    return(NULL)
  }

  # Separate batch correction from bio conservation metrics
  batch_metric_cols <- intersect(c("batch_variance", "batch_asw", "lisi", "iLISI", "kBET",
                                   "graph_connectivity", "pcr_comparison"), colnames(bench_df))
  bio_metric_cols <- intersect(c("cLISI", "silhouette_label", "isolated_labels",
                                 "kmeans_nmi", "kmeans_ari", "bio_conservation"), colnames(bench_df))

  has_bio_metrics <- length(bio_metric_cols) > 0
  all_metric_cols <- c(batch_metric_cols, bio_metric_cols)

  if (length(all_metric_cols) == 0) {
    cat("  No metric columns found\n")
    return(NULL)
  }

  cat("  Batch correction metrics:", paste(batch_metric_cols, collapse = ", "), "\n")
  if (has_bio_metrics) {
    cat("  Bio conservation metrics:", paste(bio_metric_cols, collapse = ", "), "\n")
  } else {
    cat("  NOTE: No bio conservation metrics - using batch metrics only\n")
  }

  # Prepare data
  df <- bench_df
  colnames(df)[colnames(df) == method_col] <- "method"
  n_methods <- nrow(df)

  # Convert to numeric
  for (col in all_metric_cols) {
    df[[col]] <- as.numeric(df[[col]])
  }

  # ============================================================================
  # CHECK FOR PRE-COMPUTED NORMALIZED COLUMNS (from Module 04)
  # ============================================================================
  precomputed_norm_cols <- c("batch_variance_norm", "batch_asw_norm", "lisi_norm")
  has_precomputed <- all(precomputed_norm_cols %in% colnames(df))

  if (has_precomputed) {
    cat("  Using pre-computed normalized values from Module 04\n")
    # Already have normalized columns, just ensure they're numeric
    for (col in precomputed_norm_cols) {
      df[[col]] <- as.numeric(df[[col]])
    }
  } else {
    cat("  Computing normalization on the fly\n")
    # Rank-based normalization to preserve granularity
    normalize_by_rank <- function(x, invert = FALSE) {
      if (all(is.na(x))) return(rep(0.5, length(x)))
      n <- length(x)
      if (n == 1) return(0.5)
      if (invert) {
        ranks <- rank(-x, na.last = "keep", ties.method = "average")
      } else {
        ranks <- rank(x, na.last = "keep", ties.method = "average")
      }
      normalized <- (ranks - 1) / (n - 1)
      return(normalized)
    }

    # Define which metrics to invert (lower original = better)
    invert_metrics <- c("batch_variance", "batch_asw")

    # Normalize all metrics
    for (col in all_metric_cols) {
      should_invert <- col %in% invert_metrics
      df[[paste0(col, "_norm")]] <- normalize_by_rank(df[[col]], invert = should_invert)
    }
  }

  # ============================================================================
  # CHECK FOR PRE-COMPUTED AGGREGATE SCORES (from Module 04)
  # ============================================================================
  has_batch_score <- "batch_score" %in% colnames(df)
  has_composite_score <- "composite_score" %in% colnames(df)

  if (has_batch_score) {
    cat("  Using pre-computed batch_score from Module 04\n")
    df$batch_score <- as.numeric(df$batch_score)
  } else {
    # Calculate aggregate scores
    batch_norm_cols <- paste0(batch_metric_cols, "_norm")
    batch_norm_cols <- batch_norm_cols[batch_norm_cols %in% colnames(df)]
    if (length(batch_norm_cols) > 0) {
      df$batch_score <- rowMeans(df[, batch_norm_cols, drop = FALSE], na.rm = TRUE)
    } else {
      df$batch_score <- 0.5
    }
  }

  # Bio score calculation (if metrics exist)
  if (has_bio_metrics) {
    bio_norm_cols <- paste0(bio_metric_cols, "_norm")
    bio_norm_cols <- bio_norm_cols[bio_norm_cols %in% colnames(df)]
    if (length(bio_norm_cols) > 0) {
      df$bio_score <- rowMeans(df[, bio_norm_cols, drop = FALSE], na.rm = TRUE)
    } else {
      df$bio_score <- 0.5
    }
  }

  # Total/composite score
  if (has_composite_score) {
    cat("  Using pre-computed composite_score from Module 04\n")
    df$total_score <- as.numeric(df$composite_score)
  } else if (has_bio_metrics) {
    df$total_score <- 0.6 * df$batch_score + 0.4 * df$bio_score
  } else {
    df$total_score <- df$batch_score
  }

  # Order by total score
  df <- df[order(df$total_score, decreasing = TRUE), ]
  df$method <- factor(df$method, levels = rev(df$method))

  # Create long format for metrics
  norm_cols_all <- paste0(all_metric_cols, "_norm")
  norm_cols_all <- norm_cols_all[norm_cols_all %in% colnames(df)]

  metric_data <- df %>%
    select(method, all_of(norm_cols_all)) %>%
    pivot_longer(cols = -method, names_to = "metric", values_to = "value") %>%
    mutate(
      metric_clean = gsub("_norm$", "", metric),
      metric_display = gsub("_", " ", metric_clean),
      metric_display = tools::toTitleCase(metric_display),
      category = ifelse(metric_clean %in% batch_metric_cols, "Batch", "Bio")
    )

  # Set metric order
  batch_display <- tools::toTitleCase(gsub("_", " ", batch_metric_cols))
  bio_display <- tools::toTitleCase(gsub("_", " ", bio_metric_cols))
  metric_data$metric_display <- factor(metric_data$metric_display,
                                       levels = c(batch_display, bio_display))

  # Score data - ONLY include Bio if we have bio metrics
  if (has_bio_metrics) {
    score_data <- df %>%
      select(method, batch_score, bio_score, total_score) %>%
      pivot_longer(cols = -method, names_to = "score_type", values_to = "value") %>%
      mutate(score_type = factor(score_type,
                                 levels = c("batch_score", "bio_score", "total_score"),
                                 labels = c("Batch", "Bio", "Total")))
  } else {
    score_data <- df %>%
      select(method, batch_score, total_score) %>%
      pivot_longer(cols = -method, names_to = "score_type", values_to = "value") %>%
      mutate(score_type = factor(score_type,
                                 levels = c("batch_score", "total_score"),
                                 labels = c("Batch", "Total")))
  }

  # ============================================================================
  # CREATE COLORS - Purple for batch, Green for bio
  # ============================================================================

  purple_palette <- colorRampPalette(c("#F5F0FA", "#D8BFD8", "#9370DB", "#663399", "#4B0082"))(101)
  green_palette <- colorRampPalette(c("#F0FFF0", "#90EE90", "#32CD32", "#228B22", "#006400"))(101)
  gray_palette <- colorRampPalette(c("#E8E8E8", "#1a1a2e"))(101)

  # Assign colors to metric data
  metric_data$fill_color <- mapply(function(val, cat) {
    val <- max(0, min(1, val))
    idx <- round(val * 100) + 1
    if (cat == "Batch") {
      purple_palette[idx]
    } else {
      green_palette[idx]
    }
  }, metric_data$value, metric_data$category)

  # Assign colors to score data
  score_data$bar_color <- mapply(function(val, type) {
    val <- max(0, min(1, val))
    idx <- round(val * 100) + 1
    type_char <- as.character(type)
    if (type_char == "Batch") {
      purple_palette[idx]
    } else if (type_char == "Bio") {
      green_palette[idx]
    } else {
      gray_palette[idx]
    }
  }, score_data$value, score_data$score_type)

  # Build subtitle with selection mode info if available
  if (!is.null(selection_mode) && selection_mode != "") {
    mode_label <- switch(selection_mode,
      "batch_removal" = "Aggressive (min batch_variance)",
      "balanced" = "Balanced (composite score)",
      "conservative" = "Conservative (max LISI)",
      selection_mode
    )
    plot_subtitle <- paste0("Selection: ", mode_label, " | Purple = Batch Correction (darker = better)")
  } else if (has_bio_metrics) {
    plot_subtitle <- "Purple = Batch Correction | Green = Bio Conservation (darker = better)"
  } else {
    plot_subtitle <- "Purple = Batch Correction (darker = better) | Bio metrics not available"
  }

  # Main metrics plot
  p_main <- ggplot(metric_data, aes(x = metric_display, y = method)) +
    geom_point(aes(fill = fill_color), shape = 21, size = 11, stroke = 0.8, color = "gray30") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 2.8, color = "white", fontface = "bold") +
    scale_fill_identity() +
    labs(x = "", y = "", title = title, subtitle = plot_subtitle) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )

  # Score bars plot
  n_score_cols <- length(unique(score_data$score_type))

  p_bars <- ggplot(score_data, aes(x = score_type, y = method)) +
    geom_tile(aes(fill = bar_color), width = 0.85, height = 0.85, color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3, color = "white", fontface = "bold") +
    scale_fill_identity() +
    labs(x = "", y = "", title = "Scores") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(10, 15, 10, 0)
    )

  # Combine plots
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    main_width <- length(all_metric_cols)
    score_width <- n_score_cols * 0.4

    # Build legend text
    if (!is.null(selection_mode) && selection_mode != "") {
      note_text <- paste0("Selection Mode: ", selection_mode,
                          " | Purple = Batch | ",
                          if (has_bio_metrics) "Green = Bio | " else "",
                          "Darker = Better")
    } else if (has_bio_metrics) {
      note_text <- "Purple = Batch correction | Green = Bio conservation | Darker = Better"
    } else {
      note_text <- "Note: Bio conservation not assessed | Ranking based on batch correction only"
    }

    legend_grob <- grid::textGrob(
      note_text,
      gp = grid::gpar(fontsize = 9, col = "gray40", fontface = "italic")
    )

    combined <- gridExtra::arrangeGrob(
      gridExtra::arrangeGrob(p_main, p_bars, ncol = 2, widths = c(main_width, score_width)),
      legend_grob,
      ncol = 1,
      heights = c(0.95, 0.05)
    )

    if (!is.null(output_path)) {
      ggsave(output_path, plot = combined,
             width = 11 + n_score_cols * 0.5,
             height = max(5, n_methods * 0.7 + 2),
             dpi = 150, bg = "white")
      cat("  Saved scIB-style plot to:", output_path, "\n")
      return(output_path)
    }
    return(combined)
  } else {
    if (!is.null(output_path)) {
      ggsave(output_path, plot = p_main,
             width = 10, height = max(5, n_methods * 0.7),
             dpi = 150, bg = "white")
      cat("  Saved plot to:", output_path, "\n")
      return(output_path)
    }
    return(p_main)
  }
}

# ==============================================================================
# Collect Statistics
# ==============================================================================
cat("\n--- Collecting statistics ---\n")

total_cells <- "N/A"
female_cells <- "N/A"
male_cells <- "N/A"
n_clusters <- "N/A"
n_samples <- "N/A"
sample_names <- "N/A"

if (exists("clustered_obj") && !is.null(clustered_obj)) {
  total_cells <- format(ncol(clustered_obj), big.mark = ",")
  female_cells <- format(sum(clustered_obj$sex == "Female", na.rm = TRUE), big.mark = ",")
  male_cells <- format(sum(clustered_obj$sex == "Male", na.rm = TRUE), big.mark = ",")
  n_clusters <- length(unique(clustered_obj$seurat_clusters))

  if ("sample_name" %in% colnames(clustered_obj@meta.data)) {
    n_samples <- length(unique(clustered_obj$sample_name))
    sample_names <- paste(unique(clustered_obj$sample_name), collapse = ", ")
  }
}

cat("Total cells:", total_cells, "\n")
cat("Female/Male:", female_cells, "/", male_cells, "\n")
cat("Clusters:", n_clusters, "\n")

# ==============================================================================
# Generate Benchmarking Plots
# ==============================================================================
cat("\n--- Generating benchmarking plots ---\n")

# Create directory for benchmarking plots
bench_plot_dir <- file.path(out_base, "plots", "benchmarking")
if (!dir.exists(bench_plot_dir)) {
  dir.create(bench_plot_dir, recursive = TRUE, showWarnings = FALSE)
}

# Normalization benchmarking plot
norm_bench_plot_path <- NULL
if (exists("normalization_benchmark_df") && !is.null(normalization_benchmark_df) && nrow(normalization_benchmark_df) > 0) {
  cat("  Creating normalization benchmarking plot...\n")
  norm_bench_plot_path <- file.path(bench_plot_dir, "normalization_benchmarking.png")
  create_scib_style_plot(normalization_benchmark_df,
                         title = "Normalization Method Comparison",
                         output_path = norm_bench_plot_path)
} else if (exists("norm_benchmark_df") && !is.null(norm_benchmark_df) && nrow(norm_benchmark_df) > 0) {
  cat("  Creating normalization benchmarking plot (alt name)...\n")
  norm_bench_plot_path <- file.path(bench_plot_dir, "normalization_benchmarking.png")
  create_scib_style_plot(norm_benchmark_df,
                         title = "Normalization Method Comparison",
                         output_path = norm_bench_plot_path)
} else {
  cat("  No normalization benchmarking data found\n")
}

# Integration benchmarking plot
int_bench_plot_path <- NULL
if (exists("benchmark_df") && !is.null(benchmark_df) && nrow(benchmark_df) > 0) {
  cat("  Creating integration benchmarking plot...\n")
  int_bench_plot_path <- file.path(bench_plot_dir, "integration_benchmarking.png")

  # Pass selection mode info if available from Module 04
  sel_mode <- if (exists("selection_mode")) selection_mode else NULL
  sel_crit <- if (exists("selection_criterion")) selection_criterion else NULL

  create_scib_style_plot(benchmark_df,
                         title = "Integration Method Comparison",
                         output_path = int_bench_plot_path,
                         selection_mode = sel_mode,
                         selection_criterion = sel_crit)
} else {
  cat("  No integration benchmarking data found\n")
}

# ==============================================================================
# Build HTML Report
# ==============================================================================
cat("\n--- Building HTML report ---\n")

# ============================================================================
# CSS Styles
# ============================================================================
css_styles <- '
<style>
/* ===== CSS Variables ===== */
:root {
  --bg-primary: #ffffff;
  --bg-secondary: #f8f9fa;
  --bg-tertiary: #e9ecef;
  --text-primary: #212529;
  --text-secondary: #6c757d;
  --accent: #0066cc;
  --accent-light: #e7f1ff;
  --border: #dee2e6;
  --success: #28a745;
  --warning: #ffc107;
  --danger: #dc3545;
  --shadow: 0 2px 8px rgba(0,0,0,0.1);
}

[data-theme="dark"] {
  --bg-primary: #1a1a2e;
  --bg-secondary: #16213e;
  --bg-tertiary: #0f3460;
  --text-primary: #e8e8e8;
  --text-secondary: #adb5bd;
  --accent: #4da6ff;
  --accent-light: #1a3a5c;
  --border: #3a3a5c;
  --shadow: 0 2px 8px rgba(0,0,0,0.3);
}

/* ===== Base Styles ===== */
* { box-sizing: border-box; }

body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  line-height: 1.6;
  color: var(--text-primary);
  background: var(--bg-primary);
  margin: 0;
  padding: 0;
}

/* ===== Layout ===== */
.layout {
  display: flex;
  min-height: 100vh;
}

/* ===== Sidebar Navigation ===== */
.sidebar {
  width: 280px;
  background: var(--bg-secondary);
  border-right: 1px solid var(--border);
  position: fixed;
  top: 0;
  left: 0;
  height: 100vh;
  overflow-y: auto;
  z-index: 1000;
  transition: transform 0.3s ease;
}

.sidebar-header {
  padding: 20px;
  background: linear-gradient(135deg, var(--accent), #004494);
  color: white;
}

.sidebar-header h1 {
  font-size: 1.2em;
  margin: 0 0 5px;
}

.sidebar-header .subtitle {
  font-size: 0.85em;
  opacity: 0.9;
}

.nav-menu {
  padding: 15px 0;
}

.nav-item {
  display: block;
  padding: 12px 20px;
  color: var(--text-primary);
  text-decoration: none;
  border-left: 3px solid transparent;
  transition: all 0.2s ease;
  font-size: 0.95em;
}

.nav-item:hover {
  background: var(--bg-tertiary);
  border-left-color: var(--accent);
}

.nav-item.active {
  background: var(--accent-light);
  border-left-color: var(--accent);
  font-weight: 600;
}

.nav-item .nav-icon {
  margin-right: 10px;
  width: 20px;
  display: inline-block;
  text-align: center;
}

.nav-item .nav-number {
  background: var(--bg-tertiary);
  color: var(--text-secondary);
  padding: 2px 8px;
  border-radius: 10px;
  font-size: 0.8em;
  margin-right: 8px;
}

/* ===== Main Content ===== */
.main-content {
  margin-left: 280px;
  flex: 1;
  padding: 30px 40px;
  max-width: 1200px;
}

/* ===== Header ===== */
.page-header {
  background: linear-gradient(135deg, var(--accent), #004494);
  color: white;
  padding: 40px;
  border-radius: 12px;
  margin-bottom: 30px;
  box-shadow: var(--shadow);
}

.page-header h1 {
  font-size: 2em;
  margin: 0 0 10px;
}

.page-header .meta {
  opacity: 0.9;
  font-size: 0.95em;
}

/* ===== Controls ===== */
.controls {
  display: flex;
  gap: 10px;
  margin-bottom: 25px;
  flex-wrap: wrap;
  position: sticky;
  top: 0;
  background: var(--bg-primary);
  padding: 15px 0;
  z-index: 100;
  border-bottom: 1px solid var(--border);
}

.btn {
  padding: 10px 18px;
  border: none;
  border-radius: 8px;
  cursor: pointer;
  font-size: 0.9em;
  font-weight: 500;
  transition: all 0.2s ease;
}

.btn-primary {
  background: var(--accent);
  color: white;
}

.btn-primary:hover {
  background: #0052a3;
}

.btn-secondary {
  background: var(--bg-tertiary);
  color: var(--text-primary);
}

.btn-secondary:hover {
  background: var(--border);
}

/* ===== Sections ===== */
.report-section {
  background: var(--bg-secondary);
  border: 1px solid var(--border);
  border-radius: 12px;
  margin-bottom: 20px;
  overflow: hidden;
  box-shadow: var(--shadow);
}

.section-header {
  padding: 18px 25px;
  background: var(--bg-tertiary);
  cursor: pointer;
  display: flex;
  justify-content: space-between;
  align-items: center;
  user-select: none;
}

.section-header:hover {
  background: var(--border);
}

.section-header h2 {
  margin: 0;
  font-size: 1.3em;
  color: var(--accent);
  display: flex;
  align-items: center;
  gap: 10px;
}

.section-header .icon {
  font-size: 1.2em;
}

.toggle-icon {
  color: var(--text-secondary);
  transition: transform 0.3s ease;
}

.section.collapsed .toggle-icon {
  transform: rotate(-90deg);
}

.section-content {
  padding: 25px;
  transition: max-height 0.3s ease;
}

.section.collapsed .section-content {
  display: none;
}

/* ===== Stats Grid ===== */
.stats-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
  gap: 20px;
  margin: 20px 0;
}

.stat-card {
  background: var(--bg-primary);
  border: 1px solid var(--border);
  border-radius: 10px;
  padding: 25px 20px;
  text-align: center;
  transition: transform 0.2s ease;
}

.stat-card:hover {
  transform: translateY(-3px);
  box-shadow: var(--shadow);
}

.stat-value {
  font-size: 2.2em;
  font-weight: 700;
  color: var(--accent);
  line-height: 1.2;
}

.stat-label {
  font-size: 0.9em;
  color: var(--text-secondary);
  margin-top: 8px;
  text-transform: uppercase;
  letter-spacing: 0.5px;
}

/* ===== Tables ===== */
.table-wrapper {
  overflow-x: auto;
  margin: 20px 0;
  border-radius: 8px;
  border: 1px solid var(--border);
}

.data-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.9em;
}

.data-table th,
.data-table td {
  padding: 12px 15px;
  text-align: left;
  border-bottom: 1px solid var(--border);
}

.data-table th {
  background: var(--bg-tertiary);
  font-weight: 600;
  color: var(--text-primary);
  position: sticky;
  top: 0;
}

.data-table tr:nth-child(even) {
  background: var(--bg-primary);
}

.data-table tr:hover {
  background: var(--accent-light);
}

.table-note {
  font-size: 0.85em;
  color: var(--text-secondary);
  font-style: italic;
  margin-top: 10px;
}

/* ===== Gallery ===== */
.gallery {
  display: grid;
  gap: 20px;
  margin: 20px 0;
}

.gallery-item {
  background: var(--bg-primary);
  border: 1px solid var(--border);
  border-radius: 10px;
  overflow: hidden;
}

.gallery-title {
  padding: 12px 15px;
  font-weight: 600;
  background: var(--bg-tertiary);
  font-size: 0.9em;
  border-bottom: 1px solid var(--border);
}

.gallery-img {
  padding: 15px;
  background: white;
}

.gallery-img img {
  border-radius: 6px;
}

/* ===== Info Boxes ===== */
.info-box {
  background: var(--accent-light);
  border-left: 4px solid var(--accent);
  padding: 15px 20px;
  border-radius: 0 8px 8px 0;
  margin: 20px 0;
}

.info-box h4 {
  margin: 0 0 8px;
  color: var(--accent);
}

.info-box p {
  margin: 0;
  color: var(--text-primary);
}

.warning-box {
  background: #fff3cd;
  border-left-color: var(--warning);
}

[data-theme="dark"] .warning-box {
  background: #3d3200;
}

/* ===== Processing Steps ===== */
.pipeline-flow {
  display: flex;
  flex-direction: column;
  gap: 15px;
  margin: 20px 0;
}

.pipeline-step {
  display: flex;
  align-items: flex-start;
  gap: 15px;
  padding: 20px;
  background: var(--bg-primary);
  border: 1px solid var(--border);
  border-radius: 10px;
  transition: all 0.2s ease;
}

.pipeline-step:hover {
  border-color: var(--accent);
  box-shadow: var(--shadow);
}

.step-number {
  width: 40px;
  height: 40px;
  background: var(--accent);
  color: white;
  border-radius: 50%;
  display: flex;
  align-items: center;
  justify-content: center;
  font-weight: 700;
  font-size: 1.1em;
  flex-shrink: 0;
}

.step-content {
  flex: 1;
}

.step-content h4 {
  margin: 0 0 5px;
  color: var(--accent);
}

.step-content p {
  margin: 0;
  color: var(--text-secondary);
  font-size: 0.95em;
}

.step-details {
  margin-top: 10px;
  padding-top: 10px;
  border-top: 1px dashed var(--border);
  font-size: 0.9em;
}

.step-details code {
  background: var(--bg-tertiary);
  padding: 2px 6px;
  border-radius: 4px;
  font-size: 0.9em;
}

/* ===== Benchmark Explanation ===== */
.metric-explanation {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 20px;
  margin: 20px 0;
}

.metric-card {
  background: var(--bg-primary);
  border: 1px solid var(--border);
  border-radius: 10px;
  padding: 20px;
}

.metric-card h4 {
  margin: 0 0 10px;
  color: var(--accent);
  display: flex;
  align-items: center;
  gap: 8px;
}

.metric-card .formula {
  background: var(--bg-tertiary);
  padding: 10px;
  border-radius: 6px;
  font-family: monospace;
  font-size: 0.9em;
  margin: 10px 0;
  overflow-x: auto;
}

.metric-card .interpretation {
  font-size: 0.9em;
  color: var(--text-secondary);
}

.metric-card .good-value {
  color: var(--success);
  font-weight: 600;
}

.metric-card .bad-value {
  color: var(--danger);
  font-weight: 600;
}

/* ===== Benchmark Plot ===== */
.benchmark-plot-container {
  background: white;
  border: 1px solid var(--border);
  border-radius: 10px;
  padding: 15px;
  margin: 20px 0;
  text-align: center;
}

.benchmark-plot-container img {
  max-width: 100%;
  height: auto;
  border-radius: 6px;
}

/* ===== Missing Content ===== */
.missing-img, .no-data {
  background: var(--bg-tertiary);
  padding: 30px;
  text-align: center;
  color: var(--text-secondary);
  border-radius: 8px;
  font-style: italic;
}

/* ===== Footer ===== */
.page-footer {
  margin-top: 40px;
  padding: 30px;
  text-align: center;
  background: var(--bg-secondary);
  border-radius: 12px;
  border: 1px solid var(--border);
}

.page-footer p {
  margin: 5px 0;
  color: var(--text-secondary);
}

/* ===== Modal ===== */
.modal {
  display: none;
  position: fixed;
  z-index: 2000;
  left: 0;
  top: 0;
  width: 100%;
  height: 100%;
  background: rgba(0,0,0,0.95);
  cursor: zoom-out;
}

.modal img {
  max-width: 95%;
  max-height: 95%;
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  border-radius: 8px;
  box-shadow: 0 0 50px rgba(0,0,0,0.5);
}

.modal-close {
  position: absolute;
  top: 20px;
  right: 30px;
  color: white;
  font-size: 40px;
  cursor: pointer;
  z-index: 2001;
}

/* ===== Mobile Responsive ===== */
@media (max-width: 900px) {
  .sidebar {
    transform: translateX(-100%);
  }

  .sidebar.open {
    transform: translateX(0);
  }

  .main-content {
    margin-left: 0;
    padding: 20px;
  }

  .mobile-menu-btn {
    display: block !important;
  }
}

.mobile-menu-btn {
  display: none;
  position: fixed;
  bottom: 20px;
  right: 20px;
  width: 60px;
  height: 60px;
  border-radius: 50%;
  background: var(--accent);
  color: white;
  border: none;
  font-size: 24px;
  cursor: pointer;
  z-index: 999;
  box-shadow: 0 4px 15px rgba(0,0,0,0.3);
}

/* ===== Print Styles ===== */
@media print {
  .sidebar, .controls, .modal, .mobile-menu-btn { display: none !important; }
  .main-content { margin-left: 0; padding: 0; }
  .report-section { break-inside: avoid; page-break-inside: avoid; }
  .section-content { display: block !important; }
}
</style>
'

# ============================================================================
# JavaScript
# ============================================================================
js_scripts <- '
<script>
// Toggle section collapse
function toggleSection(id) {
  const section = document.getElementById(id);
  section.classList.toggle("collapsed");
}

// Expand all sections
function expandAll() {
  document.querySelectorAll(".report-section").forEach(s => s.classList.remove("collapsed"));
}

// Collapse all sections
function collapseAll() {
  document.querySelectorAll(".report-section").forEach(s => s.classList.add("collapsed"));
}

// Toggle dark mode
function toggleDarkMode() {
  if (document.body.getAttribute("data-theme") === "dark") {
    document.body.removeAttribute("data-theme");
    localStorage.setItem("theme", "light");
  } else {
    document.body.setAttribute("data-theme", "dark");
    localStorage.setItem("theme", "dark");
  }
}

// Apply saved theme
if (localStorage.getItem("theme") === "dark") {
  document.body.setAttribute("data-theme", "dark");
}

// Modal functions
function openModal(src) {
  document.getElementById("modal").style.display = "block";
  document.getElementById("modalImg").src = src;
  document.body.style.overflow = "hidden";
}

function closeModal() {
  document.getElementById("modal").style.display = "none";
  document.body.style.overflow = "auto";
}

// Close modal on Escape key
document.addEventListener("keydown", function(e) {
  if (e.key === "Escape") closeModal();
});

// Mobile menu toggle
function toggleMobileMenu() {
  document.querySelector(".sidebar").classList.toggle("open");
}

// Smooth scroll to section
function scrollToSection(id) {
  const element = document.getElementById(id);
  if (element) {
    element.scrollIntoView({ behavior: "smooth", block: "start" });
    document.querySelector(".sidebar").classList.remove("open");
    element.classList.remove("collapsed");
  }
}

// Highlight active nav item on scroll
window.addEventListener("scroll", function() {
  const sections = document.querySelectorAll(".report-section");
  const navItems = document.querySelectorAll(".nav-item");

  let current = "";
  sections.forEach(section => {
    const sectionTop = section.offsetTop;
    if (window.pageYOffset >= sectionTop - 100) {
      current = section.getAttribute("id");
    }
  });

  navItems.forEach(item => {
    item.classList.remove("active");
    if (item.getAttribute("href") === "#" + current) {
      item.classList.add("active");
    }
  });
});
</script>
'

# ============================================================================
# HTML Header
# ============================================================================
html_header <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>scRNA-seq Analysis Report - Choroid Plexus</title>
', css_styles, '
</head>
<body>

<!-- Modal for image zoom -->
<div id="modal" class="modal" onclick="closeModal()">
  <span class="modal-close">&times;</span>
  <img id="modalImg" src="">
</div>

<!-- Mobile menu button -->
<button class="mobile-menu-btn" onclick="toggleMobileMenu()">&#9776;</button>
')

# ============================================================================
# Navigation Sidebar (using HTML entities instead of emojis)
# ============================================================================
nav_items <- c(
  "overview" = "&#128202; Overview",
  "pipeline" = "&#9881; Pipeline Summary",
  "qc" = "&#128269; Quality Control",
  "normalization" = "&#128200; Normalization",
  "integration" = "&#128279; Integration",
  "benchmarking" = "&#127942; Benchmarking",
  "clustering" = "&#127793; Clustering",
  "scice" = "&#128300; scICE Subclustering",
  "de" = "&#128201; Differential Expression",
  "visualization" = "&#127912; Gene Expression",
  "session" = "&#128196; Session Info"
)

nav_html <- paste0('
<nav class="sidebar">
  <div class="sidebar-header">
    <h1>&#128300; scRNA-seq Report</h1>
    <div class="subtitle">Choroid Plexus Analysis</div>
  </div>
  <div class="nav-menu">
', paste(sapply(seq_along(nav_items), function(i) {
  id <- names(nav_items)[i]
  label <- nav_items[i]
  paste0('<a href="#', id, '" class="nav-item" onclick="scrollToSection(\'', id, '\'); return false;">',
         '<span class="nav-number">', sprintf("%02d", i), '</span>', label, '</a>')
}), collapse = '\n'), '
  </div>
</nav>
')

# ============================================================================
# Get dataset and ventricle info for header
# ============================================================================
dataset_name_display <- if (!is.null(params$dataset_name) && params$dataset_name != "default") {
  params$dataset_name
} else {
  "Choroid Plexus"
}

ventricle_display <- if (!is.null(params$ventricle_filter) && params$ventricle_filter != "") {
  params$ventricle_filter
} else {
  "All"
}

# ============================================================================
# Main Content Start
# ============================================================================
main_start <- paste0('
<div class="layout">
', nav_html, '
<main class="main-content">

<!-- Page Header -->
<div class="page-header">
  <h1>&#128300; scRNA-seq Multi-Sample Analysis Report</h1>
  <div class="meta">
    <strong>Dataset:</strong> ', html_escape(dataset_name_display), '<br>
    <strong>Ventricle:</strong> ', html_escape(ventricle_display), '<br>
    <strong>Generated:</strong> ', format(Sys.time(), "%B %d, %Y at %H:%M:%S"), '<br>
    <strong>Pipeline:</strong> Seurat v5 Multi-Sample Pipeline
  </div>
</div>

<!-- Controls -->
<div class="controls">
  <button class="btn btn-primary" onclick="expandAll()">&#9660; Expand All</button>
  <button class="btn btn-primary" onclick="collapseAll()">&#9650; Collapse All</button>
  <button class="btn btn-secondary" onclick="toggleDarkMode()">&#9788; Toggle Dark Mode</button>
  <button class="btn btn-secondary" onclick="window.print()">&#128424; Print Report</button>
</div>
')

# ============================================================================
# SECTION 1: Overview
# ============================================================================
cat("  Building Section 1: Overview\n")

overview_content <- paste0('
<div class="stats-grid">
  <div class="stat-card">
    <div class="stat-value">', n_samples, '</div>
    <div class="stat-label">Samples</div>
  </div>
  <div class="stat-card">
    <div class="stat-value">', total_cells, '</div>
    <div class="stat-label">Total Cells</div>
  </div>
  <div class="stat-card">
    <div class="stat-value">', female_cells, '</div>
    <div class="stat-label">Female Cells</div>
  </div>
  <div class="stat-card">
    <div class="stat-value">', male_cells, '</div>
    <div class="stat-label">Male Cells</div>
  </div>
  <div class="stat-card">
    <div class="stat-value">', n_clusters, '</div>
    <div class="stat-label">Clusters</div>
  </div>
</div>

<h3>Analysis Configuration</h3>
')

# Determine normalization method display
norm_method_display <- if (exists("selected_normalization_method") && !is.null(selected_normalization_method)) {
  as.character(selected_normalization_method)
} else if (isTRUE(params$run_sctransform)) {
  "SCTransform"
} else if (isTRUE(params$run_scran)) {
  "scran"
} else if (isTRUE(params$run_lognorm)) {
  "LogNormalize"
} else {
  "N/A"
}

# Configuration table
config_df <- data.frame(
  Parameter = c(
    "Dataset",
    "Ventricle Filter",
    "Samples Analyzed",
    "Normalization Method",
    "Integration Method",
    "Clustering Method",
    "Batch Variable",
    "Batch Integration",
    "DE Comparison"
  ),
  Value = c(
    dataset_name_display,
    ventricle_display,
    sample_names,
    norm_method_display,
    if (exists("best_method")) as.character(best_method) else "N/A",
    if (exists("clustering_method_used")) as.character(clustering_method_used) else "N/A",
    if (!is.null(params$batch_variable)) params$batch_variable else "batch",
    as.character(isTRUE(params$run_batch_integration)),
    if (!is.null(params$de_comparison_scope)) params$de_comparison_scope else "global"
  ),
  stringsAsFactors = FALSE
)

overview_content <- paste0(overview_content, create_table(config_df))
section_overview <- create_section("overview", "1. Executive Overview", overview_content, "&#128202;")

# ============================================================================
# SECTION 2: Pipeline Processing Summary
# ============================================================================
cat("  Building Section 2: Pipeline Summary\n")

pipeline_content <- '
<div class="info-box">
  <h4>Pipeline Architecture</h4>
  <p>This analysis follows a modular pipeline design where each module performs specific tasks and passes results to subsequent modules. The pipeline is optimized for Seurat v5 with proper layer handling throughout.</p>
</div>

<h3>Processing Steps</h3>
<div class="pipeline-flow">

<div class="pipeline-step">
  <div class="step-number">00</div>
  <div class="step-content">
    <h4>Environment Setup</h4>
    <p>Initialize R environment, load required packages, create output directory structure, and set global parameters.</p>
    <div class="step-details">
      <strong>Key Actions:</strong> Load Seurat, create subdirectories for plots/objects/tables, initialize logging
    </div>
  </div>
</div>

<div class="pipeline-step">
  <div class="step-number">01</div>
  <div class="step-content">
    <h4>Data Loading</h4>
    <p>Import scCDC-corrected count matrices from preprocessing, extract corrected counts from the <code>scCDC_corrected</code> layer, and build initial Seurat objects.</p>
    <div class="step-details">
      <strong>Input:</strong> Step 7 preprocessed .rds files with scCDC contamination correction<br>
      <strong>Key Transformation:</strong> <code>pct_counts_MT</code> &rarr; <code>percent.mt</code>, calculate <code>percent.hb</code> from hemoglobin genes (^HB[AB]-)
    </div>
  </div>
</div>

<div class="pipeline-step">
  <div class="step-number">02</div>
  <div class="step-content">
    <h4>Quality Control Validation</h4>
    <p>Apply QC filters for gene counts, UMI counts, mitochondrial percentage, and hemoglobin percentage. Merge filtered samples.</p>
    <div class="step-details">
      <strong>Filters Applied:</strong><br>
      - nFeature_RNA: 200 - 8000<br>
      - percent.mt: &le; 15%<br>
      - percent.hb: &le; 2% (removes erythrocyte contamination)<br>
      <strong>Output:</strong> Merged Seurat object with all samples, JoinLayers applied after merge
    </div>
  </div>
</div>

<div class="pipeline-step">
  <div class="step-number">03</div>
  <div class="step-content">
    <h4>Normalization</h4>
    <p>Apply multiple normalization methods and compare results. Methods are benchmarked and the best performer is selected automatically.</p>
    <div class="step-details">
      <strong>Methods:</strong><br>
      - <strong>SCTransform:</strong> Variance-stabilizing transformation with regularized negative binomial regression<br>
      - <strong>scran:</strong> Deconvolution-based size factor estimation using cell pooling<br>
      - <strong>LogNormalize:</strong> Standard log-normalization with scale factor 10,000<br>
      <strong>Seurat v5:</strong> Uses <code>seurat_to_sce()</code> and <code>sce_to_seurat()</code> for SCE conversion, JoinLayers before operations
    </div>
  </div>
</div>

<div class="pipeline-step">
  <div class="step-number">04</div>
  <div class="step-content">
    <h4>Integration &amp; Batch Correction</h4>
    <p>Test multiple integration methods to find optimal batch correction while preserving biological signal. Benchmark all methods.</p>
    <div class="step-details">
      <strong>Methods Tested:</strong> Harmony, CCA, RPCA, FastMNN, scVI, Scanorama, BBKNN, scCobra, CONCORD<br>
      <strong>Seurat v5:</strong> JoinLayers before any Python method export, split layers for CCA/RPCA
    </div>
  </div>
</div>

<div class="pipeline-step">
  <div class="step-number">05</div>
  <div class="step-content">
    <h4>CHOIR Clustering</h4>
    <p>Apply CHOIR (Clustering Hierarchy Optimization by Iterative Random forests) for hierarchical cell type identification.</p>
    <div class="step-details">
      <strong>Algorithm:</strong> Uses random forests to identify clusters at multiple resolutions, builds hierarchy tree<br>
      <strong>Parameters:</strong> alpha = 0.05, uses Harmony/best reduction as input
    </div>
  </div>
</div>

<div class="pipeline-step">
  <div class="step-number">06</div>
  <div class="step-content">
    <h4>scICE Subclustering</h4>
    <p>Apply scICE (single-cell Iterative Clustering and Evaluation) for consistent subclustering of CHOIR clusters.</p>
    <div class="step-details">
      <strong>Purpose:</strong> Refine large clusters to identify subtypes while maintaining consistency<br>
      <strong>Seurat v5:</strong> JoinLayers before count matrix extraction for Julia export
    </div>
  </div>
</div>

<div class="pipeline-step">
  <div class="step-number">07</div>
  <div class="step-content">
    <h4>Leiden Clustering</h4>
    <p>Test multiple Leiden resolutions as comparison/fallback to CHOIR. Compute quality metrics for each resolution.</p>
    <div class="step-details">
      <strong>Resolutions Tested:</strong> 0.3, 0.5, 0.8, 1.0, 1.2, 1.5<br>
      <strong>Quality Metrics:</strong> Silhouette score, neighborhood purity, cluster size distribution
    </div>
  </div>
</div>

<div class="pipeline-step">
  <div class="step-number">08</div>
  <div class="step-content">
    <h4>Differential Expression</h4>
    <p>Compare gene expression between Male and Female cells using multiple DE methods for robust results.</p>
    <div class="step-details">
      <strong>Methods:</strong><br>
      - <strong>MAST:</strong> Hurdle model for single-cell data (accounts for dropout)<br>
      - <strong>edgeR:</strong> Pseudobulk with voom transformation<br>
      - <strong>DESeq2:</strong> Pseudobulk with negative binomial model<br>
      <strong>Note:</strong> MAST uses natural log scale; <code>log2FC_if_ln</code> column provides converted values
    </div>
  </div>
</div>

<div class="pipeline-step">
  <div class="step-number">09</div>
  <div class="step-content">
    <h4>Gene Visualization</h4>
    <p>Create comprehensive visualizations for genes of interest including UMAPs, violin plots, dot plots, and heatmaps.</p>
  </div>
</div>

<div class="pipeline-step">
  <div class="step-number">10</div>
  <div class="step-content">
    <h4>Final Summary</h4>
    <p>Compile analysis summary, per-sample statistics, and session information.</p>
  </div>
</div>

<div class="pipeline-step">
  <div class="step-number">11</div>
  <div class="step-content">
    <h4>HTML Report</h4>
    <p>Generate this interactive report with all results, plots, and documentation.</p>
  </div>
</div>

</div>
'

section_pipeline <- create_section("pipeline", "2. Pipeline Processing Summary", pipeline_content, "&#9881;")

# ============================================================================
# SECTION 3: Quality Control
# ============================================================================
cat("  Building Section 3: QC\n")

# Find QC plots
qc_dirs <- c(
  file.path(out_base, "plots", "qc"),
  file.path(out_base, "01_QC"),
  file.path(out_base, "02_QC_Validation")
)

qc_pngs <- character(0)
for (d in qc_dirs) {
  if (dir.exists(d)) {
    found <- find_pngs(d, recursive = TRUE)
    if (length(found) > 0) {
      qc_pngs <- c(qc_pngs, found)
    }
  }
}

qc_content <- '
<div class="info-box">
  <h4>Quality Control Filters</h4>
  <p>Cells are filtered based on multiple metrics to remove low-quality cells, doublets, and contamination.</p>
</div>

<h3>QC Thresholds Applied</h3>
'

qc_params_df <- data.frame(
  Metric = c("nFeature_RNA (genes)", "nCount_RNA (UMIs)", "percent.mt", "percent.hb"),
  Minimum = c(
    if (!is.null(params$min_features)) params$min_features else 200,
    "N/A",
    "N/A",
    "N/A"
  ),
  Maximum = c(
    if (!is.null(params$max_features)) params$max_features else 8000,
    "N/A",
    if (!is.null(params$max_percent_mt)) paste0(params$max_percent_mt, "%") else "15%",
    if (!is.null(params$max_percent_hb)) paste0(params$max_percent_hb, "%") else "2%"
  ),
  Rationale = c(
    "Remove empty droplets and doublets",
    "Remove low-quality cells",
    "Remove dying/stressed cells",
    "Remove erythrocyte contamination"
  ),
  stringsAsFactors = FALSE
)

qc_content <- paste0(qc_content, create_table(qc_params_df))

if (length(qc_pngs) > 0) {
  qc_content <- paste0(qc_content, '<h3>QC Visualizations</h3>', create_gallery(head(qc_pngs, 8), 2))
} else {
  qc_content <- paste0(qc_content, '<p class="no-data">No QC plots found</p>')
}

section_qc <- create_section("qc", "3. Quality Control", qc_content, "&#128269;")

# ============================================================================
# SECTION 4: Normalization
# ============================================================================
cat("  Building Section 4: Normalization\n")

norm_dirs <- c(
  file.path(out_base, "02_Normalization"),
  file.path(out_base, "03_Normalization"),
  file.path(out_base, "plots", "normalization")
)

norm_pngs <- character(0)
for (d in norm_dirs) {
  if (dir.exists(d)) {
    found <- find_pngs(d, recursive = TRUE)
    if (length(found) > 0) {
      norm_pngs <- c(norm_pngs, found)
    }
  }
}

norm_content <- paste0('
<div class="info-box">
  <h4>Selected Method: ', html_escape(norm_method_display), '</h4>
  <p>Normalization removes technical variation while preserving biological signal. The method was selected based on benchmarking metrics.</p>
</div>

<h3>Normalization Methods Comparison</h3>

<div class="metric-explanation">
  <div class="metric-card">
    <h4>SCTransform</h4>
    <p>Variance-stabilizing transformation using regularized negative binomial regression. Models each gene separately, accounting for sequencing depth.</p>
    <div class="formula">log(E[Y]) = &beta;<sub>0</sub> + &beta;<sub>1</sub> * log(UMI)</div>
    <p class="interpretation"><span class="good-value">&#10003; Recommended</span> for most datasets. Handles heteroscedasticity well.</p>
  </div>

  <div class="metric-card">
    <h4>scran</h4>
    <p>Deconvolution-based size factors computed by pooling cells and solving a linear system.</p>
    <div class="formula">SF<sub>i</sub> = median(counts<sub>i</sub> / pool_avg)</div>
    <p class="interpretation">Good for datasets with many zeros. Requires cell clustering first.</p>
  </div>

  <div class="metric-card">
    <h4>LogNormalize</h4>
    <p>Standard normalization: divide by total counts, multiply by scale factor, log-transform.</p>
    <div class="formula">log((count / total) * 10000 + 1)</div>
    <p class="interpretation">Simple baseline method. May not handle variance well.</p>
  </div>
</div>
')

# Add normalization benchmarking plot if available
if (!is.null(norm_bench_plot_path) && file.exists(norm_bench_plot_path)) {
  norm_content <- paste0(norm_content, '
<h3>Normalization Benchmarking</h3>
<div class="benchmark-plot-container">
', embed_image(norm_bench_plot_path, alt = "Normalization Benchmarking Plot"), '
</div>
<p class="table-note">Higher scores indicate better batch mixing after normalization. Metrics are normalized to 0-1 range with inversions applied where lower values originally indicated better performance.</p>
')
}

# Add normalization benchmarking table if available
if (exists("normalization_benchmark_df") && !is.null(normalization_benchmark_df) && nrow(normalization_benchmark_df) > 0) {
  norm_content <- paste0(norm_content, '<h3>Normalization Metrics (Raw Values)</h3>', create_table(normalization_benchmark_df))
} else if (exists("norm_benchmark_df") && !is.null(norm_benchmark_df) && nrow(norm_benchmark_df) > 0) {
  norm_content <- paste0(norm_content, '<h3>Normalization Metrics (Raw Values)</h3>', create_table(norm_benchmark_df))
}

if (length(norm_pngs) > 0) {
  norm_content <- paste0(norm_content, '<h3>Normalization Plots</h3>', create_gallery(head(norm_pngs, 8), 2))
}

section_norm <- create_section("normalization", "4. Normalization", norm_content, "&#128200;")

# ============================================================================
# SECTION 5: Integration
# ============================================================================
cat("  Building Section 5: Integration\n")

int_dirs <- c(
  file.path(out_base, "04_Benchmarking"),
  file.path(out_base, "03_Integration"),
  file.path(out_base, "plots", "integration")
)

int_pngs <- character(0)
for (d in int_dirs) {
  if (dir.exists(d)) {
    found <- find_pngs(d, recursive = TRUE)
    if (length(found) > 0) {
      int_pngs <- c(int_pngs, found)
    }
  }
}

int_content <- paste0('
<div class="info-box">
  <h4>Selected Method: ',
  if (exists("best_method")) html_escape(best_method) else "N/A",
  '</h4>
  <p>Integration corrects for batch effects while preserving biological variation between conditions.</p>
</div>

<h3>Integration Methods Tested</h3>

<div class="metric-explanation">
  <div class="metric-card">
    <h4>Harmony</h4>
    <p>Iteratively adjusts PCA embeddings to remove batch effects using soft k-means clustering.</p>
    <p class="interpretation"><span class="good-value">&#10003; Fast</span>, works well for most datasets. Operates on reduced dimensions.</p>
  </div>

  <div class="metric-card">
    <h4>CCA (Seurat)</h4>
    <p>Canonical Correlation Analysis identifies shared variation between batches.</p>
    <p class="interpretation">Good when batches share cell types. Can over-correct if batches are very different.</p>
  </div>

  <div class="metric-card">
    <h4>RPCA (Seurat)</h4>
    <p>Reciprocal PCA - more conservative than CCA, projects each batch onto others.</p>
    <p class="interpretation">Better when batches have unique populations. Preserves more biological signal.</p>
  </div>

  <div class="metric-card">
    <h4>FastMNN</h4>
    <p>Mutual Nearest Neighbors - finds pairs of similar cells across batches.</p>
    <p class="interpretation">Robust to differences in cell type composition between batches.</p>
  </div>

  <div class="metric-card">
    <h4>scVI</h4>
    <p>Deep learning variational autoencoder that learns batch-corrected latent space.</p>
    <p class="interpretation"><span class="good-value">&#10003; Powerful</span> for complex batch effects. Requires GPU for speed.</p>
  </div>

  <div class="metric-card">
    <h4>Scanorama</h4>
    <p>Panoramic stitching approach - aligns datasets like image panoramas.</p>
    <p class="interpretation">Works well when datasets have partial overlap.</p>
  </div>
  
  <div class="metric-card">
    <h4>Scanorama</h4>
    <p>Panoramic stitching approach - aligns datasets like image panoramas.</p>
    <p class="interpretation">Works well when datasets have partial overlap.</p>
  </div>

  <div class="metric-card">
    <h4>CONCORD</h4>
    <p>Contrastive learning with dataset-aware sampling. Each minibatch contains cells from only one batch, preventing the model from learning batch effects.</p>
    <p class="interpretation"><span class="good-value">&#10003; State-of-the-art</span> (Nature Biotechnology 2025). Preserves biological signal while removing batch effects.</p>
  </div>
</div>
')

if (length(int_pngs) > 0) {
  int_content <- paste0(int_content, '<h3>Integration Plots</h3>', create_gallery(head(int_pngs, 8), 2))
}

section_int <- create_section("integration", "5. Integration", int_content, "&#128279;")

# ============================================================================
# SECTION 6: Benchmarking
# ============================================================================
cat("  Building Section 6: Benchmarking\n")

bench_content <- '
<div class="info-box">
  <h4>Benchmarking Integration Quality</h4>
  <p>Multiple metrics are computed to evaluate how well each integration method corrects batch effects while preserving biological signal. The visualization below uses the scIB (single-cell Integration Benchmarking) framework style.</p>
</div>

<h3>Benchmarking Metrics Explained</h3>

<div class="metric-explanation">
  <div class="metric-card">
    <h4>&#9632; Batch Variance (kBET)</h4>
    <p>Measures how well batches are mixed locally using k-nearest neighbor acceptance testing.</p>
    <div class="formula">kBET = 1 - rejection_rate</div>
    <p class="interpretation">
      <span class="good-value">Higher is better</span> (close to 1.0)<br>
      Values near 1.0 indicate batches are well-mixed locally.<br>
      <span class="bad-value">Low values</span> suggest batch separation persists.
    </p>
  </div>

  <div class="metric-card">
    <h4>&#9632; Batch ASW (Silhouette Width)</h4>
    <p>Average Silhouette Width for batch labels - measures batch separation in embedding space.</p>
    <div class="formula">ASW = (b - a) / max(a, b)</div>
    <p class="interpretation">
      <span class="good-value">Lower is better</span> (close to 0 or negative)<br>
      0 = batches perfectly mixed<br>
      <span class="bad-value">Positive values</span> indicate batch clustering
    </p>
  </div>

  <div class="metric-card">
    <h4>&#9632; LISI (Local Inverse Simpson Index)</h4>
    <p>Measures diversity of batch labels in local neighborhoods.</p>
    <div class="formula">LISI = 1 / &Sigma;(p<sub>i</sub>&sup2;)</div>
    <p class="interpretation">
      <span class="good-value">Higher is better</span><br>
      LISI = n_batches means perfect mixing<br>
      LISI = 1 means no mixing (cells only near same batch)
    </p>
  </div>

  <div class="metric-card">
    <h4>&#9632; Bio Conservation (cLISI)</h4>
    <p>Cell-type LISI - ensures biological populations are preserved after integration.</p>
    <div class="formula">cLISI for cell type labels</div>
    <p class="interpretation">
      <span class="good-value">Lower is better</span> (close to 1.0)<br>
      Cell types should remain clustered, not artificially mixed.
    </p>
  </div>

  <div class="metric-card">
    <h4>&#9632; Overall Score</h4>
    <p>Weighted combination of batch mixing and bio conservation.</p>
    <div class="formula">Score = 0.6 * batch_mixing + 0.4 * bio_conservation</div>
    <p class="interpretation">
      <span class="good-value">Higher is better</span><br>
      Balances batch correction with biological signal preservation.
    </p>
  </div>
</div>
'

# Add integration benchmarking plot if available
if (!is.null(int_bench_plot_path) && file.exists(int_bench_plot_path)) {
  bench_content <- paste0(bench_content, '
<h3>Integration Benchmarking Visualization</h3>
<div class="benchmark-plot-container">
', embed_image(int_bench_plot_path, alt = "Integration Benchmarking Plot"), '
</div>
<p class="table-note">
<strong>Reading the plot:</strong> Each row represents an integration method. Circle colors show normalized metric values (darker = better performance).
The "Total" column shows the aggregate score combining all metrics. Methods are ordered by total score (best at top).
<br><br>
<strong>Metric normalization:</strong> batch_variance and batch_asw are inverted (original: lower = better &rarr; display: higher = better).
LISI is shown as-is (higher = better mixing).
</p>
')
}

# Add benchmark results table if available
if (exists("benchmark_df") && !is.null(benchmark_df) && nrow(benchmark_df) > 0) {
  bench_content <- paste0(bench_content, '<h3>Integration Benchmarking Results (Raw Values)</h3>', create_table(benchmark_df))

  # Highlight best method
  if (exists("best_method") && !is.null(best_method)) {
    bench_content <- paste0(bench_content, '
    <div class="info-box">
      <h4>&#10003; Selected Method: ', html_escape(best_method), '</h4>
      <p>This method achieved the best overall score balancing batch correction and biological signal preservation.</p>
    </div>')
  }
} else {
  bench_content <- paste0(bench_content, '<p class="no-data">No integration benchmarking results table available</p>')
}

section_bench <- create_section("benchmarking", "6. Integration Benchmarking", bench_content, "&#127942;")

# ============================================================================
# SECTION 7: Clustering
# ============================================================================
cat("  Building Section 7: Clustering\n")

# Find clustering plots
choir_dirs <- c(
  file.path(out_base, "05_CHOIR_Clustering", "plots"),
  file.path(out_base, "05_CHOIR_Clustering"),
  file.path(out_base, "plots", "choir")
)

leiden_dirs <- c(
  file.path(out_base, "plots", "leiden_clustering"),
  file.path(out_base, "07_Leiden_Clustering"),
  file.path(out_base, "06_Leiden_Clustering")
)

choir_pngs <- character(0)
for (d in choir_dirs) {
  if (dir.exists(d)) {
    found <- find_pngs(d, recursive = TRUE)
    if (length(found) > 0) choir_pngs <- c(choir_pngs, found)
  }
}

leiden_pngs <- character(0)
for (d in leiden_dirs) {
  if (dir.exists(d)) {
    found <- find_pngs(d, recursive = TRUE)
    if (length(found) > 0) leiden_pngs <- c(leiden_pngs, found)
  }
}

clust_content <- paste0('
<div class="info-box">
  <h4>Primary Method: ',
  if (exists("clustering_method_used")) html_escape(clustering_method_used) else "N/A",
  '</h4>
  <p>Final number of clusters: <strong>', n_clusters, '</strong></p>
</div>

<h3>CHOIR Clustering</h3>
<p>CHOIR (Clustering Hierarchy Optimization by Iterative Random forests) builds a hierarchical clustering tree using random forests to identify cell populations at multiple resolutions.</p>
')

if (length(choir_pngs) > 0) {
  clust_content <- paste0(clust_content, create_gallery(head(choir_pngs, 6), 2))
} else {
  clust_content <- paste0(clust_content, '<p class="no-data">No CHOIR plots found</p>')
}

clust_content <- paste0(clust_content, '
<h3>Leiden Clustering</h3>
<p>Leiden algorithm tested at multiple resolutions for comparison. Quality metrics computed for each resolution.</p>
')

if (length(leiden_pngs) > 0) {
  clust_content <- paste0(clust_content, create_gallery(head(leiden_pngs, 6), 2))
} else {
  clust_content <- paste0(clust_content, '<p class="no-data">No Leiden plots found</p>')
}

section_clust <- create_section("clustering", "7. Clustering", clust_content, "&#127793;")

# ============================================================================
# SECTION 8: scICE Subclustering - COMPLETE UPDATED VERSION
# ============================================================================
cat("  Building Section 8: scICE\n")

scice_content <- '
<div class="info-box">
  <h4>scICE Subclustering</h4>
  <p>scICE (single-cell Iterative Clustering and Evaluation) provides consistent subclustering to identify cell subtypes within major CHOIR clusters. The method uses scLENS dimensionality reduction followed by iterative clustering to find stable subclusters.</p>
</div>
'

# Add scICE results table if available
if (exists("scice_results") && length(scice_results) > 0) {

  scice_df <- data.frame(
    Cluster = names(scice_results),
    Success = sapply(scice_results, function(x) if (isTRUE(x$success)) "Yes" else "No"),
    Method = sapply(scice_results, function(x) if (!is.null(x$method)) x$method else "N/A"),
    Input_Cells = sapply(scice_results, function(x) if (!is.null(x$input_cells)) format(x$input_cells, big.mark=",") else "N/A"),
    Output_Cells = sapply(scice_results, function(x) if (!is.null(x$n_cells)) format(x$n_cells, big.mark=",") else "N/A"),
    stringsAsFactors = FALSE
  )

  n_successful <- sum(sapply(scice_results, function(x) isTRUE(x$success)))
  scice_content <- paste0(scice_content,
    '<p>Processed <strong>', length(scice_results), '</strong> clusters, <strong>',
    n_successful, '</strong> successful.</p>',
    create_table(scice_df))
} else {
  scice_content <- paste0(scice_content, '<p class="no-data">scICE subclustering was not performed or no results available.</p>')
}

# ============================================================================
# IMPROVED: Find scICE plots with COMPLETE directory and file pattern search
# ============================================================================

# Define all possible scICE directory names (case variations)
scice_dirs <- c(
  file.path(out_base, "06_scICE_subclustering"),      # lowercase 's' - actual directory name
  file.path(out_base, "06_scICE_Subclustering"),      # uppercase 'S' - fallback
  file.path(out_base, "06_scice_subclustering"),      # all lowercase
  file.path(out_base, "plots", "scice"),
  file.path(out_base, "plots", "scICE")
)

# Initialize collectors
scice_main_plot <- NULL
scice_comparison_plots <- character(0)
scice_umap_dist_plots <- character(0)
scice_umap_cluster_plots <- character(0)  # NEW: umap_l_*.png files
scice_ic_plots <- character(0)
scice_other_plots <- character(0)

cat("  Searching for scICE plots...\n")

for (d in scice_dirs) {
  if (!dir.exists(d)) next
  
  cat("    Found scICE directory:", d, "\n")
  
  # -------------------------------------------------------------------------
  # 1. Look for main comparison plot at root level
  # -------------------------------------------------------------------------
  main_plot_patterns <- c(
    "scICE_comparison_plot.png",
    "scice_comparison_plot.png",
    "comparison_plot.png",
    "scICE_overview.png"
  )
  
  for (pattern in main_plot_patterns) {
    main_plot <- file.path(d, pattern)
    if (file.exists(main_plot) && is.null(scice_main_plot)) {
      scice_main_plot <- main_plot
      cat("      Found main comparison plot:", basename(main_plot), "\n")
      break
    }
  }
  
  # Also check for any comparison plots
  root_pngs <- list.files(d, pattern = "\\.png$", full.names = TRUE, recursive = FALSE)
  comparison_matches <- root_pngs[grepl("comparison|overview|summary", basename(root_pngs), ignore.case = TRUE)]
  if (length(comparison_matches) > 0) {
    scice_comparison_plots <- c(scice_comparison_plots, comparison_matches)
  }
  
  # -------------------------------------------------------------------------
  # 2. Look for per-cluster plots in output/cluster_*/ subdirectories
  # -------------------------------------------------------------------------
  output_dir <- file.path(d, "output")
  if (dir.exists(output_dir)) {
    cluster_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
    cluster_dirs <- cluster_dirs[grepl("cluster_", basename(cluster_dirs))]
    
    cat("      Found", length(cluster_dirs), "cluster output directories\n")
    
    for (cdir in cluster_dirs) {
      cluster_name <- basename(cdir)
      
      # Get all PNG files in this cluster directory
      cluster_pngs <- list.files(cdir, pattern = "\\.png$", full.names = TRUE, 
                                  recursive = FALSE, ignore.case = TRUE)
      
      if (length(cluster_pngs) == 0) next
      
      # Categorize by type
      for (png_file in cluster_pngs) {
        fname <- basename(png_file)
        
        if (grepl("^umap_dist", fname, ignore.case = TRUE)) {
          # UMAP distribution plot
          scice_umap_dist_plots <- c(scice_umap_dist_plots, png_file)
          
        } else if (grepl("^umap_l_", fname, ignore.case = TRUE)) {
          # NEW: Cluster label UMAP visualization (umap_l_2.png, umap_l_3.png, etc.)
          scice_umap_cluster_plots <- c(scice_umap_cluster_plots, png_file)
          
        } else if (grepl("^ic_plot", fname, ignore.case = TRUE)) {
          # IC (Information Criterion) plot
          scice_ic_plots <- c(scice_ic_plots, png_file)
          
        } else if (grepl("umap", fname, ignore.case = TRUE)) {
          # Other UMAP-related plots
          scice_umap_cluster_plots <- c(scice_umap_cluster_plots, png_file)
          
        } else {
          # Other plots
          scice_other_plots <- c(scice_other_plots, png_file)
        }
      }
    }
  }
  
  # -------------------------------------------------------------------------
  # 3. Also check for plots directly in the scICE directory
  # -------------------------------------------------------------------------
  direct_pngs <- list.files(d, pattern = "\\.png$", full.names = TRUE, recursive = FALSE)
  for (png_file in direct_pngs) {
    fname <- basename(png_file)
    # Skip if already captured as main plot
    if (!is.null(scice_main_plot) && png_file == scice_main_plot) next
    if (png_file %in% scice_comparison_plots) next
    
    if (grepl("umap_dist", fname, ignore.case = TRUE)) {
      scice_umap_dist_plots <- c(scice_umap_dist_plots, png_file)
    } else if (grepl("umap_l_", fname, ignore.case = TRUE)) {
      scice_umap_cluster_plots <- c(scice_umap_cluster_plots, png_file)
    } else if (grepl("ic_plot", fname, ignore.case = TRUE)) {
      scice_ic_plots <- c(scice_ic_plots, png_file)
    }
  }
}

# Remove duplicates
scice_comparison_plots <- unique(scice_comparison_plots)
scice_umap_dist_plots <- unique(scice_umap_dist_plots)
scice_umap_cluster_plots <- unique(scice_umap_cluster_plots)
scice_ic_plots <- unique(scice_ic_plots)
scice_other_plots <- unique(scice_other_plots)

# ============================================================================
# Sort plots by cluster number for consistent ordering
# ============================================================================
sort_by_cluster_number <- function(paths) {
  if (length(paths) == 0) return(paths)
  
  # Extract cluster numbers from paths
  cluster_nums <- sapply(paths, function(p) {
    # Try to extract from directory name (cluster_1, cluster_10, etc.)
    dir_match <- regmatches(dirname(p), regexpr("cluster_([0-9]+)", dirname(p)))
    if (length(dir_match) > 0 && nchar(dir_match) > 0) {
      return(as.numeric(gsub("cluster_", "", dir_match)))
    }
    # Try to extract from filename
    file_match <- regmatches(basename(p), regexpr("_([0-9]+)", basename(p)))
    if (length(file_match) > 0 && nchar(file_match) > 0) {
      return(as.numeric(gsub("_", "", file_match)))
    }
    return(999)  # Default for unmatched
  })
  
  paths[order(cluster_nums)]
}

scice_umap_dist_plots <- sort_by_cluster_number(scice_umap_dist_plots)
scice_umap_cluster_plots <- sort_by_cluster_number(scice_umap_cluster_plots)
scice_ic_plots <- sort_by_cluster_number(scice_ic_plots)

# Report what was found
cat("  scICE plot summary:\n")
cat("    Main comparison plot:", if (!is.null(scice_main_plot)) "Found" else "Not found", "\n")
cat("    UMAP distribution plots:", length(scice_umap_dist_plots), "\n")
cat("    Cluster UMAP visualizations:", length(scice_umap_cluster_plots), "\n")
cat("    IC plots:", length(scice_ic_plots), "\n")
cat("    Other plots:", length(scice_other_plots), "\n")

# ============================================================================
# Build HTML content for scICE plots
# ============================================================================

# Helper function to create titles from paths
create_cluster_title <- function(path, suffix = "") {
  # Extract cluster number
  dir_name <- basename(dirname(path))
  if (grepl("cluster_", dir_name)) {
    cluster_num <- gsub("cluster_", "Cluster ", dir_name)
  } else {
    # Try from filename
    fname <- tools::file_path_sans_ext(basename(path))
    cluster_num <- gsub("_", " ", fname)
    cluster_num <- tools::toTitleCase(cluster_num)
  }
  
  if (suffix != "") {
    paste(cluster_num, "-", suffix)
  } else {
    cluster_num
  }
}

# -------------------------------------------------------------------------
# Add main comparison plot prominently
# -------------------------------------------------------------------------
if (!is.null(scice_main_plot) && file.exists(scice_main_plot)) {
  scice_content <- paste0(scice_content, 
    '<h3>scICE Subclustering Overview</h3>',
    '<div class="benchmark-plot-container">',
    embed_image(scice_main_plot, alt = "scICE Comparison Plot"),
    '</div>',
    '<p class="table-note">Comparison of CHOIR clusters with scICE subclustering results showing cluster assignments and consistency.</p>'
  )
}

# -------------------------------------------------------------------------
# Add Cluster UMAP Visualizations (umap_l_*.png) - NEW PRIORITY SECTION
# -------------------------------------------------------------------------
if (length(scice_umap_cluster_plots) > 0) {
  # Create descriptive titles
  umap_cluster_titles <- sapply(scice_umap_cluster_plots, function(p) {
    fname <- tools::file_path_sans_ext(basename(p))
    dir_name <- basename(dirname(p))
    
    # Extract k value from filename (e.g., umap_l_2.png -> k=2)
    k_match <- regmatches(fname, regexpr("[0-9]+", fname))
    k_val <- if (length(k_match) > 0) k_match else "?"
    
    # Extract cluster number from directory
    cluster_match <- regmatches(dir_name, regexpr("[0-9]+", dir_name))
    cluster_num <- if (length(cluster_match) > 0) cluster_match else "?"
    
    paste0("Cluster ", cluster_num, " (k=", k_val, ")")
  })
  
  # Determine how many to show
  n_show <- min(20, length(scice_umap_cluster_plots))
  
  scice_content <- paste0(scice_content, 
    '<h3>Cluster Label UMAPs (', length(scice_umap_cluster_plots), ' visualizations)</h3>',
    '<p>UMAP visualizations showing the consistent subcluster assignments identified by scICE for each CHOIR cluster. The k value indicates the number of subclusters found.</p>',
    create_gallery(head(scice_umap_cluster_plots, n_show), ncol = 4, 
                   titles = head(umap_cluster_titles, n_show))
  )
  
  if (length(scice_umap_cluster_plots) > n_show) {
    scice_content <- paste0(scice_content, 
      '<p class="table-note">Showing first ', n_show, ' of ', 
      length(scice_umap_cluster_plots), ' cluster UMAP visualizations</p>')
  }
}

# -------------------------------------------------------------------------
# Add UMAP Distribution plots
# -------------------------------------------------------------------------
if (length(scice_umap_dist_plots) > 0) {
  umap_dist_titles <- sapply(scice_umap_dist_plots, function(p) {
    dir_name <- basename(dirname(p))
    cluster_match <- regmatches(dir_name, regexpr("[0-9]+", dir_name))
    cluster_num <- if (length(cluster_match) > 0) cluster_match else "?"
    paste0("Cluster ", cluster_num)
  })
  
  n_show <- min(12, length(scice_umap_dist_plots))
  
  scice_content <- paste0(scice_content, 
    '<h3>UMAP Distributions (', length(scice_umap_dist_plots), ' clusters)</h3>',
    '<p>scLENS-derived UMAP embeddings for each cluster before subcluster assignment.</p>',
    create_gallery(head(scice_umap_dist_plots, n_show), ncol = 3, 
                   titles = head(umap_dist_titles, n_show))
  )
  
  if (length(scice_umap_dist_plots) > n_show) {
    scice_content <- paste0(scice_content, 
      '<p class="table-note">Showing first ', n_show, ' of ', 
      length(scice_umap_dist_plots), ' UMAP distribution plots</p>')
  }
}

# -------------------------------------------------------------------------
# Add IC plots
# -------------------------------------------------------------------------
if (length(scice_ic_plots) > 0) {
  ic_titles <- sapply(scice_ic_plots, function(p) {
    dir_name <- basename(dirname(p))
    cluster_match <- regmatches(dir_name, regexpr("[0-9]+", dir_name))
    cluster_num <- if (length(cluster_match) > 0) cluster_match else "?"
    paste0("Cluster ", cluster_num, " IC")
  })
  
  n_show <- min(12, length(scice_ic_plots))
  
  scice_content <- paste0(scice_content, 
    '<h3>Information Criterion (IC) Plots (', length(scice_ic_plots), ' clusters)</h3>',
    '<p>IC plots showing the clustering quality across different k values. The optimal k is selected where IC approaches the threshold (1.005).</p>',
    create_gallery(head(scice_ic_plots, n_show), ncol = 3, 
                   titles = head(ic_titles, n_show))
  )
  
  if (length(scice_ic_plots) > n_show) {
    scice_content <- paste0(scice_content, 
      '<p class="table-note">Showing first ', n_show, ' of ', 
      length(scice_ic_plots), ' IC plots</p>')
  }
}

# -------------------------------------------------------------------------
# Add any other plots
# -------------------------------------------------------------------------
if (length(scice_other_plots) > 0) {
  other_titles <- tools::file_path_sans_ext(basename(scice_other_plots))
  
  n_show <- min(6, length(scice_other_plots))
  
  scice_content <- paste0(scice_content, 
    '<h3>Additional scICE Visualizations</h3>',
    create_gallery(head(scice_other_plots, n_show), ncol = 2, 
                   titles = head(other_titles, n_show))
  )
}

# -------------------------------------------------------------------------
# Add note if no plots were found
# -------------------------------------------------------------------------
total_plots <- sum(
  !is.null(scice_main_plot),
  length(scice_umap_dist_plots),
  length(scice_umap_cluster_plots),
  length(scice_ic_plots),
  length(scice_other_plots)
)

if (total_plots == 0) {
  scice_content <- paste0(scice_content, 
    '<div class="info-box warning-box">',
    '<h4>No scICE Plots Found</h4>',
    '<p>scICE subclustering plots were not found in the expected locations. ',
    'Expected directory: <code>06_scICE_subclustering/output/cluster_*/</code></p>',
    '<p>Check that Module 06 completed successfully.</p>',
    '</div>'
  )
} else {
  cat("  Total scICE plots found:", total_plots, "\n")
}

# Create the section
section_scice <- create_section("scice", "8. scICE Subclustering", scice_content, "&#128300;")

# ============================================================================
# SECTION 9: Differential Expression
# ============================================================================
cat("  Building Section 9: DE\n")

de_dirs <- c(
  file.path(out_base, "plots", "de"),
  file.path(out_base, "08_Differential_Expression"),
  file.path(out_base, "07_DE_Analysis")
)

de_pngs <- character(0)
for (d in de_dirs) {
  if (dir.exists(d)) {
    found <- find_pngs(d, recursive = TRUE)
    if (length(found) > 0) de_pngs <- c(de_pngs, found)
  }
}

# Also check subdirs if they exist
if (exists("subdirs") && !is.null(subdirs$plots_de)) {
  found <- find_pngs(subdirs$plots_de, recursive = TRUE)
  if (length(found) > 0) de_pngs <- c(de_pngs, found)
}

de_content <- '
<div class="info-box">
  <h4>Differential Expression Analysis</h4>
  <p>Comparing gene expression between Male and Female cells using multiple methods for robust results.</p>
</div>

<h3>DE Methods</h3>
<div class="metric-explanation">
  <div class="metric-card">
    <h4>MAST</h4>
    <p>Hurdle model designed for single-cell data. Models both the rate of expression (logistic) and the level of expression (Gaussian).</p>
    <p class="interpretation"><strong>Note:</strong> Reports logFC in natural log scale. Use <code>log2FC_if_ln</code> column for log2 scale.</p>
  </div>

  <div class="metric-card">
    <h4>edgeR (Pseudobulk)</h4>
    <p>Aggregates cells into pseudobulk samples, then applies negative binomial GLM with voom transformation.</p>
    <p class="interpretation">More robust for small cell numbers per group. Requires biological replicates.</p>
  </div>

  <div class="metric-card">
    <h4>DESeq2 (Pseudobulk)</h4>
    <p>Pseudobulk aggregation with negative binomial model and shrinkage estimators for dispersion.</p>
    <p class="interpretation">Conservative but reliable. Good for datasets with few replicates.</p>
  </div>
</div>
'

# Add DE summary if available
if (exists("de_summary") && !is.null(de_summary) && nrow(de_summary) > 0) {
  de_content <- paste0(de_content, '<h3>DE Summary</h3>', create_table(de_summary))
}

# Add top DE genes for each method
if (exists("all_de_results") && length(all_de_results) > 0) {
  # Define mapping of method names to their significance column names
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
    "findmarkers" = "p_val_adj"
  )

  for (method in names(all_de_results)) {
    res <- all_de_results[[method]]

    # Skip if result is NULL or not a data frame
    if (is.null(res) || !is.data.frame(res)) {
      next
    }

    # Try to get sig_col from mapping first
    sig_col <- sig_col_mapping[[method]]

    # If not in mapping, try to auto-detect
    if (is.null(sig_col)) {
      possible_sig_cols <- c("FDR", "adj.P.Val", "padj", "p_val_adj", "q_value",
                             "qvalue", "p.adjust", "BH", "fdr", "adjusted_pvalue")
      detected_cols <- intersect(possible_sig_cols, colnames(res))
      if (length(detected_cols) > 0) {
        sig_col <- detected_cols[1]
      } else {
        next
      }
    }

    if (sig_col %in% colnames(res)) {
      sig_genes <- res[!is.na(res[[sig_col]]) & res[[sig_col]] < 0.05, ]
      if (nrow(sig_genes) > 0) {
        top_genes <- head(sig_genes[order(sig_genes[[sig_col]]), ], 15)
        # Select relevant columns
        display_cols <- intersect(c("gene", "logFC", "log2FoldChange", "FDR", "adj.P.Val", "padj", "direction"), colnames(top_genes))
        top_genes_df <- as.data.frame(top_genes)[, display_cols, drop = FALSE]
        de_content <- paste0(de_content,
          '<h3>Top ', method, ' Significant Genes (n=', nrow(sig_genes), ' total)</h3>',
          create_table(top_genes_df))
      }
    }
  }
}

if (length(de_pngs) > 0) {
  de_content <- paste0(de_content, '<h3>Volcano Plots</h3>', create_gallery(de_pngs, 2))
}

section_de <- create_section("de", "9. Differential Expression", de_content, "&#128201;")

# ============================================================================
# SECTION 10: Gene Expression Visualization
# ============================================================================
cat("  Building Section 10: Gene Visualization\n")

viz_dirs <- c(
  file.path(out_base, "plots", "gene_expression"),
  file.path(out_base, "09_Gene_Visualization"),
  file.path(out_base, "08_Gene_Visualization")
)

viz_pngs <- character(0)
for (d in viz_dirs) {
  if (dir.exists(d)) {
    found <- find_pngs(d, recursive = TRUE)
    if (length(found) > 0) viz_pngs <- c(viz_pngs, found)
  }
}

viz_content <- '
<div class="info-box">
  <h4>Gene Expression Visualization</h4>
  <p>Visualizations for genes of interest including feature plots, violin plots, dot plots, and heatmaps.</p>
</div>
'

# Add available genes if known
if (exists("available_genes") && length(available_genes) > 0) {
  viz_content <- paste0(viz_content,
    '<p><strong>Genes visualized:</strong> ', html_escape(paste(available_genes, collapse = ", ")), '</p>')
}

if (length(viz_pngs) > 0) {
  viz_content <- paste0(viz_content, '<h3>Expression Plots</h3>', create_gallery(head(viz_pngs, 12), 2))
} else {
  viz_content <- paste0(viz_content, '<p class="no-data">No gene expression plots found</p>')
}

section_viz <- create_section("visualization", "10. Gene Expression", viz_content, "&#127912;")

# ============================================================================
# SECTION 11: Session Info
# ============================================================================
cat("  Building Section 11: Session Info\n")

session_content <- paste0('
<h3>R Environment</h3>
<p><strong>R Version:</strong> ', R.version.string, '</p>
<p><strong>Platform:</strong> ', R.version$platform, '</p>

<h3>Output Location</h3>
<p><code>', html_escape(out_base), '</code></p>

<h3>Key Package Versions</h3>
')

# Get package versions
pkg_list <- c("Seurat", "SeuratObject", "ggplot2", "dplyr", "Matrix", "MAST", "edgeR", "DESeq2", "harmony", "scran")
pkg_versions <- sapply(pkg_list, function(p) {
  tryCatch(as.character(packageVersion(p)), error = function(e) "Not installed")
})

pkg_df <- data.frame(
  Package = pkg_list,
  Version = pkg_versions,
  stringsAsFactors = FALSE
)

session_content <- paste0(session_content, create_table(pkg_df))

session_content <- paste0(session_content, '
<h3>Environment Variables</h3>
')

# Show environment variables if available
env_vars_df <- data.frame(
  Variable = c("PIPELINE_DIR", "PREPROCESS_DIR", "DOWNSTREAM_DIR", "DATASET_NAME", "VENTRICLE_FILTER"),
  Value = c(
    Sys.getenv("PIPELINE_DIR", unset = "(not set)"),
    Sys.getenv("PREPROCESS_DIR", unset = "(not set)"),
    Sys.getenv("DOWNSTREAM_DIR", unset = "(not set)"),
    Sys.getenv("DATASET_NAME", unset = "(not set)"),
    Sys.getenv("VENTRICLE_FILTER", unset = "(not set)")
  ),
  stringsAsFactors = FALSE
)

session_content <- paste0(session_content, create_table(env_vars_df))

session_content <- paste0(session_content, '
<h3>Report Generated</h3>
<p>', format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), '</p>
')

section_session <- create_section("session", "11. Session Information", session_content, "&#128196;")

# ============================================================================
# Footer
# ============================================================================
html_footer <- paste0('
</main>
</div>

<footer class="page-footer">
  <p><strong>scRNA-seq Multi-Sample Analysis Pipeline</strong></p>
  <p>Doetsch Lab | University of Basel</p>
  <p>Report generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>
</footer>

', js_scripts, '
</body>
</html>
')

# ============================================================================
# Assemble Full Report
# ============================================================================
cat("\n--- Assembling final report ---\n")

full_html <- paste0(
  html_header,
  main_start,
  section_overview,
  section_pipeline,
  section_qc,
  section_norm,
  section_int,
  section_bench,
  section_clust,
  section_scice,
  section_de,
  section_viz,
  section_session,
  html_footer
)

# ============================================================================
# Write Report
# ============================================================================
html_path <- file.path(out_base, "Analysis_Report_MultiSample.html")
writeLines(full_html, html_path)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat(">>> HTML REPORT GENERATED <<<\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Output file:", html_path, "\n")
if (file.exists(html_path)) {
  cat("File size:", round(file.size(html_path) / 1e6, 2), "MB\n")
}

cat("\nBenchmarking plots generated:\n")
if (!is.null(norm_bench_plot_path) && file.exists(norm_bench_plot_path)) {
  cat("  - Normalization:", norm_bench_plot_path, "\n")
}
if (!is.null(int_bench_plot_path) && file.exists(int_bench_plot_path)) {
  cat("  - Integration:", int_bench_plot_path, "\n")
}

cat("\nReport sections:\n")
cat("  1. Executive Overview\n")
cat("  2. Pipeline Processing Summary\n")
cat("  3. Quality Control\n")
cat("  4. Normalization (with benchmarking plot)\n")
cat("  5. Integration\n")
cat("  6. Integration Benchmarking (with visualization plot)\n")
cat("  7. Clustering\n")
cat("  8. scICE Subclustering\n")
cat("  9. Differential Expression\n")
cat("  10. Gene Expression\n")
cat("  11. Session Information\n")

cat("\n>>> MODULE 11 COMPLETE <<<\n")