# ==============================================================================
# MODULE 03b: ALRA IMPUTATION (NORMALIZED-DATA-BASED)
# ==============================================================================
#
# This module performs dropout imputation using ALRA (Adaptively-thresholded
# Low Rank Approximation) on normalized data AFTER Module 03 normalization.
#
# WHEN TO RUN:
#   - When params$imputation_method = "alra" or "both"
#   - Runs AFTER Module 03 (Normalization)
#   - ONLY compatible with LogNormalize and scran (NOT SCTransform)
#   - Automatically SKIPPED if SCTransform wins normalization benchmarking
#
# INPUT:
#   - 03_normalization_data.RData from Module 03
#   - merged_normalized/merged_{method}_unintegrated.rds
#
# OUTPUT:
#   - Updated 03_normalization_data.RData with ALRA assay
#   - Updated merged_normalized/merged_{method}_unintegrated.rds
#   - normalized_object_with_alra.rds
#   - QC plots and statistics
#
# REFERENCE:
#   https://github.com/KlugerLab/ALRA
#   Linderman et al., Nature Methods 2022
#
# UPDATES:
# - 2026-01-12: Fixed cell/gene name preservation after ALRA processing
# - 2026-01-12: Added tryCatch for cleaner error messages
# - 2026-01-10: Fixed Seurat v5 compatibility - Assays() returns S4, use names()
# - 2026-01-10: Fixed normalization method detection priority
# - 2026-01-09: Created module for ALRA imputation post-normalization
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 03b: ALRA IMPUTATION (NORMALIZED-DATA-BASED)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")


# ==============================================================================
# REQUIRED LIBRARIES
# ==============================================================================
# Load these libraries when running module independently
# ==============================================================================

required_packages <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "ALRA")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package '", pkg, "' is not installed.\n",
         "Install with: install.packages('", pkg, "')", 
         if (pkg == "ALRA") " or devtools::install_github('KlugerLab/ALRA')" else "",
         call. = FALSE)
  }
}

library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)
# Note: ALRA is loaded later after availability check

cat("Required libraries loaded successfully\n")


# ==============================================================================
# SETUP
# ==============================================================================

# Get pipeline directory
pipeline_dir <- Sys.getenv("PIPELINE_DIR", unset = "")
if (pipeline_dir == "") {
  script_dir <- tryCatch({
    dirname(sys.frame(1)$ofile)
  }, error = function(e) ".")
  pipeline_dir <- normalizePath(file.path(script_dir, ".."))
}

# Load environment
env_file <- file.path(pipeline_dir, "config", "params.R")
source(env_file)

# Load pipeline environment from Module 00
env_data_file <- file.path(params$out_root, "objects", "pipeline_environment.RData")
if (file.exists(env_data_file)) {
  load(env_data_file)
  cat("Loaded pipeline environment from Module 00\n")
} else {
  stop("Pipeline environment not found. Run Module 00 first.")
}

cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ==============================================================================
# CHECK IF ALRA IMPUTATION SHOULD RUN
# ==============================================================================

# Check imputation_method parameter
run_alra_imputation <- params$imputation_method %in% c("alra", "both")

if (!run_alra_imputation) {
  cat(">>> ALRA IMPUTATION SKIPPED <<<\n")
  cat("    imputation_method =", params$imputation_method, "\n")
  cat("    Set imputation_method = 'alra' or 'both' to enable ALRA imputation\n")
  cat("\n>>> MODULE 03b COMPLETE (skipped) <<<\n")
  quit(save = "no", status = 0)
}

cat("ALRA imputation enabled (imputation_method =", params$imputation_method, ")\n\n")

# ==============================================================================
# CHECK ALRA AVAILABILITY
# ==============================================================================

if (!has_alra) {
  cat(">>> ALRA IMPUTATION SKIPPED <<<\n")
  cat("    ALRA R package not available\n")
  cat("    To install ALRA:\n")
  cat("      install.packages('ALRA')  # If on CRAN\n")
  cat("      # or\n")
  cat("      devtools::install_github('KlugerLab/ALRA')\n")
  cat("\n>>> MODULE 03b COMPLETE (skipped due to missing ALRA) <<<\n")
  quit(save = "no", status = 0)
}

# Load ALRA package
library(ALRA)
cat("ALRA package loaded\n\n")

# ==============================================================================
# DETERMINE SELECTED NORMALIZATION METHOD
# ==============================================================================
# Priority cascade for finding the normalization method:
#   1. Load from 03_normalization_data.RData (MOST RELIABLE)
#   2. Check integration_winner.rds (if Module 04 ran)
#   3. Check selected_normalization_method.txt
#   4. Check benchmark results CSV
#   5. Infer from filenames (fallback)
#   6. Default from params
# ==============================================================================

cat("--- Checking selected normalization method ---\n")

selected_normalization_method <- NULL

# PRIORITY 1: Load from Module 03 normalization data (MOST RELIABLE)
norm_data_file <- file.path(output_dirs$objects, "03_normalization_data.RData")
if (file.exists(norm_data_file)) {
  cat("Loading normalization data from Module 03...\n")
  norm_env <- new.env()
  load(norm_data_file, envir = norm_env)
  if ("selected_normalization_method" %in% ls(norm_env)) {
    selected_normalization_method <- norm_env$selected_normalization_method
    cat(">>> Found normalization method from Module 03:", selected_normalization_method, "<<<\n")
  }
}

# PRIORITY 2: Check integration_winner.rds (if Module 04 ran)
if (is.null(selected_normalization_method)) {
  winner_file <- file.path(output_dirs$objects, "integration_winner.rds")
  if (file.exists(winner_file)) {
    winner <- readRDS(winner_file)
    if (!is.null(winner$normalization_method)) {
      selected_normalization_method <- winner$normalization_method
      cat("Found normalization method from integration winner:", selected_normalization_method, "\n")
    }
  }
}

# PRIORITY 3: Check selected_normalization_method.txt
if (is.null(selected_normalization_method)) {
  method_file <- file.path(output_dirs$objects, "selected_normalization_method.txt")
  if (file.exists(method_file)) {
    selected_normalization_method <- trimws(readLines(method_file, n = 1))
    cat("Found normalization method from text file:", selected_normalization_method, "\n")
  }
}

# PRIORITY 4: Check benchmark results CSV
if (is.null(selected_normalization_method)) {
  benchmark_file <- file.path(output_dirs$normalization, "benchmarking", "normalization_benchmark_results.csv")
  if (file.exists(benchmark_file)) {
    bench_results <- read.csv(benchmark_file)
    if ("method" %in% colnames(bench_results) && "composite_score" %in% colnames(bench_results)) {
      selected_normalization_method <- bench_results$method[which.max(bench_results$composite_score)]
      cat("Found normalization method from benchmark CSV:", selected_normalization_method, "\n")
    }
  }
}

# PRIORITY 5: Infer from available files (fallback)
if (is.null(selected_normalization_method)) {
  cat("No explicit method selection found. Checking available normalized objects...\n")
  merged_norm_dir <- file.path(output_dirs$objects, "merged_normalized")
  if (dir.exists(merged_norm_dir)) {
    available_files <- list.files(merged_norm_dir, pattern = "merged_.*_unintegrated\\.rds$")
    # Extract method names from filenames
    available_methods <- gsub("merged_(.*)_unintegrated\\.rds", "\\1", available_files)
    cat("  Available methods:", paste(available_methods, collapse = ", "), "\n")

    # Preference order: scran > LogNormalize > SCTransform
    # (since SCTransform is incompatible with ALRA anyway)
    if ("scran" %in% available_methods) {
      selected_normalization_method <- "scran"
    } else if ("LogNormalize" %in% available_methods) {
      selected_normalization_method <- "LogNormalize"
    } else if ("SCTransform" %in% available_methods) {
      selected_normalization_method <- "SCTransform"
    }

    if (!is.null(selected_normalization_method)) {
      cat("Inferred normalization method from available objects:", selected_normalization_method, "\n")
    }
  }
}

# PRIORITY 6: Default from params
if (is.null(selected_normalization_method)) {
  if (!is.null(params$integration_normalization_method) &&
      params$integration_normalization_method != "auto") {
    selected_normalization_method <- params$integration_normalization_method
    cat("Using default normalization method from params:", selected_normalization_method, "\n")
  } else {
    selected_normalization_method <- "scran"
    cat("Using fallback default normalization method: scran\n")
  }
}

cat("\n>>> Selected normalization method:", selected_normalization_method, "<<<\n\n")

# ==============================================================================
# CHECK ALRA COMPATIBILITY
# ==============================================================================
# ALRA only works with LogNormalize and scran (NOT SCTransform)
# ==============================================================================

alra_compatible <- params$alra_compatible_methods
if (is.null(alra_compatible)) {
  alra_compatible <- c("LogNormalize", "scran")
}

if (!selected_normalization_method %in% alra_compatible) {
  cat(">>> ALRA IMPUTATION SKIPPED <<<\n")
  cat("    Selected normalization method:", selected_normalization_method, "\n")
  cat("    ALRA is ONLY compatible with:", paste(alra_compatible, collapse = ", "), "\n")
  cat("    SCTransform already handles technical noise via negative binomial regression.\n")
  cat("\n>>> MODULE 03b COMPLETE (skipped due to incompatible normalization) <<<\n")
  quit(save = "no", status = 0)
}

cat("Normalization method", selected_normalization_method, "is compatible with ALRA\n\n")

# ==============================================================================
# LOAD NORMALIZED OBJECT
# ==============================================================================

cat("--- Loading normalized object ---\n")

# Construct path to normalized object
merged_norm_dir <- file.path(output_dirs$objects, "merged_normalized")
input_file <- file.path(merged_norm_dir,
                        paste0("merged_", selected_normalization_method, "_unintegrated.rds"))

if (!file.exists(input_file)) {
  # Try alternative locations
  alt_paths <- c(
    file.path(output_dirs$objects, paste0("integrated_", selected_normalization_method, "_Harmony.rds")),
    file.path(output_dirs$objects, paste0("normalized_", selected_normalization_method, ".rds"))
  )
  for (alt in alt_paths) {
    if (file.exists(alt)) {
      input_file <- alt
      break
    }
  }
}

if (!file.exists(input_file)) {
  stop("Normalized object not found: ", input_file, "\nRun Module 03 first.")
}

seurat_obj <- readRDS(input_file)
cat("Loaded:", input_file, "\n")
cat("Cells:", ncol(seurat_obj), "\n")
cat("Genes:", nrow(seurat_obj), "\n")

# ==============================================================================
# FIX FOR SEURAT V5: Assays() returns S4 object, not character vector
# ==============================================================================
# In Seurat v5, Assays(seurat_obj) returns a SimpleAssays S4 object.
# To get assay names as a character vector, use names(seurat_obj@assays)
# ==============================================================================

available_assays <- names(seurat_obj@assays)
cat("Assays:", paste(available_assays, collapse = ", "), "\n\n")

# ==============================================================================
# SETUP OUTPUT DIRECTORIES
# ==============================================================================

alra_out_dir <- output_dirs$imputation_alra
alra_plots_dir <- subdirs$alra_plots
alra_tables_dir <- subdirs$alra_tables

dir.create(alra_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(alra_plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(alra_tables_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# JOIN LAYERS (Required for Seurat v5 with multiple samples)
# ==============================================================================
# In Seurat v5, merged objects have separate layers per sample (data.1, data.2, etc.)
# GetAssayData cannot extract data when multiple layers exist.
# JoinLayers combines them into a single layer.
# Reference: https://satijalab.org/seurat/articles/seurat5_essential_commands
# ==============================================================================

cat("--- Joining layers (Seurat v5 requirement) ---\n")
cat("Layers before joining:", paste(Layers(seurat_obj), collapse = ", "), "\n")

seurat_obj <- JoinLayers(seurat_obj)

cat("Layers after joining:", paste(Layers(seurat_obj), collapse = ", "), "\n\n")

# ==============================================================================
# PREPARE DATA FOR ALRA
# ==============================================================================

cat("--- Preparing data for ALRA ---\n")

# Get normalized data layer
# For scran/LogNormalize, data is in "data" layer of RNA assay
if ("RNA" %in% available_assays) {
  norm_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
} else {
  norm_matrix <- GetAssayData(seurat_obj, layer = "data")
}

cat("Normalized matrix:", nrow(norm_matrix), "genes x", ncol(norm_matrix), "cells\n")

# ==============================================================================
# CRITICAL FIX: Store original dimension names before ALRA processing
# ==============================================================================
# ALRA does NOT preserve row/column names during SVD reconstruction.
# We must store them and restore them after processing.
# ==============================================================================

original_genes <- rownames(norm_matrix)
original_cells <- colnames(norm_matrix)

cat("Stored original dimension names:\n")
cat("  Genes:", length(original_genes), "(first 5:", paste(head(original_genes, 5), collapse = ", "), "...)\n")
cat("  Cells:", length(original_cells), "(first 5:", paste(head(original_cells, 5), collapse = ", "), "...)\n")

# Convert to dense matrix if needed (ALRA works with dense matrices)
if (inherits(norm_matrix, "sparseMatrix")) {
  cat("Converting sparse matrix to dense for ALRA...\n")
  norm_dense <- as.matrix(norm_matrix)
} else {
  norm_dense <- norm_matrix
}

# Calculate pre-imputation statistics
pre_zeros <- sum(norm_dense == 0)
pre_total <- length(norm_dense)
pre_sparsity <- pre_zeros / pre_total * 100

cat("Pre-ALRA sparsity:", round(pre_sparsity, 2), "%\n")
cat("Non-zero entries:", pre_total - pre_zeros, "\n\n")

# ==============================================================================
# DETERMINE OPTIMAL K (SVD RANK)
# ==============================================================================

cat("--- Determining optimal SVD rank (k) ---\n")

k_value <- params$alra_k

if (is.null(k_value)) {
  cat("Auto-detecting k using choose_k()...\n")

  # ALRA's choose_k function finds the optimal rank
  # Need to transpose: ALRA expects cells x genes
  k_result <- tryCatch({
    choose_k(t(norm_dense), noise_start = 80)
  }, error = function(e) {
    cat("Warning: choose_k failed, using default k=20\n")
    list(k = 20)
  })

  k_value <- k_result$k
  cat("Auto-detected k:", k_value, "\n")

  # Save k selection plot if available
  if (!is.null(k_result$plot)) {
    tryCatch({
      png(file.path(alra_plots_dir, "alra_k_selection.png"), width = 800, height = 600)
      print(k_result$plot)
      dev.off()
      cat("Saved k selection plot\n")
    }, error = function(e) {
      cat("Could not save k selection plot:", e$message, "\n")
    })
  }
} else {
  cat("Using user-specified k:", k_value, "\n")
}

cat("\n")

# ==============================================================================
# RUN ALRA IMPUTATION
# ==============================================================================

cat("--- Running ALRA imputation ---\n")
cat("Parameters:\n")
cat("  k (SVD rank):", k_value, "\n")
cat("  q (power iterations):", params$alra_q, "\n")
cat("  quantile_prob:", params$alra_quantile_prob, "\n")

alra_start <- Sys.time()

# ALRA expects cells x genes matrix
alra_result <- tryCatch({
  alra(
    A_norm = t(norm_dense),
    k = k_value,
    q = params$alra_q,
    quantile.prob = params$alra_quantile_prob
  )
}, error = function(e) {
  stop("ALRA imputation failed: ", e$message, call. = FALSE)
})

alra_duration <- difftime(Sys.time(), alra_start, units = "mins")
cat("ALRA completed in", round(as.numeric(alra_duration), 2), "minutes\n\n")

# ==============================================================================
# PROCESS ALRA OUTPUT
# ==============================================================================

cat("--- Processing ALRA output ---\n")

# ALRA returns list with:
#   - A_norm_rank_k_cor_res: the imputed matrix (cells x genes)
#   - A_norm_rank_k: the low-rank approximation
#   - rand_svd: the randomized SVD result

# Get imputed matrix and transpose back to genes x cells
imputed_matrix <- t(alra_result[[1]])

cat("Imputed matrix dimensions:", nrow(imputed_matrix), "genes x", ncol(imputed_matrix), "cells\n")

# ==============================================================================
# CRITICAL FIX: Restore original dimension names
# ==============================================================================
# ALRA does not preserve dimension names - we must restore them explicitly
# ==============================================================================

cat("Restoring original dimension names...\n")

# Verify dimensions match
if (nrow(imputed_matrix) != length(original_genes)) {
  stop("Gene count mismatch after ALRA: expected ", length(original_genes), 
       ", got ", nrow(imputed_matrix), call. = FALSE)
}
if (ncol(imputed_matrix) != length(original_cells)) {
  stop("Cell count mismatch after ALRA: expected ", length(original_cells), 
       ", got ", ncol(imputed_matrix), call. = FALSE)
}

# Restore dimension names
rownames(imputed_matrix) <- original_genes
colnames(imputed_matrix) <- original_cells

cat("  Restored gene names:", length(rownames(imputed_matrix)), "\n")
cat("  Restored cell names:", length(colnames(imputed_matrix)), "\n")

# Verify names match Seurat object
seurat_cells <- colnames(seurat_obj)
if (!all(colnames(imputed_matrix) == seurat_cells)) {
  cat("  WARNING: Cell order verification...\n")
  if (setequal(colnames(imputed_matrix), seurat_cells)) {
    cat("  Cell names match but order differs - reordering...\n")
    imputed_matrix <- imputed_matrix[, seurat_cells]
  } else {
    stop("Cell name mismatch between imputed matrix and Seurat object", call. = FALSE)
  }
}

cat("  Dimension names verified and restored successfully\n\n")

# Calculate post-imputation statistics
post_zeros <- sum(imputed_matrix == 0)
post_total <- length(imputed_matrix)
post_sparsity <- post_zeros / post_total * 100

cat("Post-ALRA sparsity:", round(post_sparsity, 2), "%\n")
cat("Sparsity reduction:", round(pre_sparsity - post_sparsity, 2), "percentage points\n\n")

# ==============================================================================
# ADD ALRA ASSAY TO SEURAT OBJECT
# ==============================================================================

cat("--- Adding ALRA assay to Seurat object ---\n")

# Convert to sparse matrix for storage efficiency
imputed_sparse <- Matrix::Matrix(imputed_matrix, sparse = TRUE)

# Verify sparse matrix has correct names
cat("Sparse matrix dimensions:", nrow(imputed_sparse), "x", ncol(imputed_sparse), "\n")
cat("Sparse matrix has rownames:", !is.null(rownames(imputed_sparse)), "\n")
cat("Sparse matrix has colnames:", !is.null(colnames(imputed_sparse)), "\n")

# Create new assay with ALRA-imputed data
# Note: ALRA returns normalized log-space values, not counts
alra_assay <- CreateAssay5Object(data = imputed_sparse)

cat("Created ALRA assay with", nrow(alra_assay), "features\n")

# Add to Seurat object with error handling to suppress verbose output
tryCatch({
  seurat_obj[["ALRA"]] <- alra_assay
  cat("Successfully added 'ALRA' assay to Seurat object\n")
}, error = function(e) {
  # Provide clean error message without verbose object dump
  stop("Failed to add ALRA assay to Seurat object: ", conditionMessage(e),
       "\n\nDebug info:",
       "\n  Seurat object cells: ", ncol(seurat_obj),
       "\n  ALRA assay cells: ", ncol(alra_assay),
       "\n  Seurat cell names (first 3): ", paste(head(colnames(seurat_obj), 3), collapse = ", "),
       "\n  ALRA cell names (first 3): ", paste(head(colnames(imputed_sparse), 3), collapse = ", "),
       call. = FALSE)
})

cat("Available assays:", paste(names(seurat_obj@assays), collapse = ", "), "\n\n")

# ==============================================================================
# GENERATE QC PLOTS
# ==============================================================================

cat("--- Generating QC plots ---\n")

# 1. Sparsity comparison
sparsity_data <- data.frame(
  Stage = c("Original (normalized)", "After ALRA"),
  Sparsity = c(pre_sparsity, post_sparsity)
)

p_sparsity <- ggplot(sparsity_data, aes(x = Stage, y = Sparsity, fill = Stage)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Sparsity)), vjust = -0.5) +
  scale_fill_manual(values = c("Original (normalized)" = "#E74C3C", "After ALRA" = "#3498DB")) +
  labs(title = "Matrix Sparsity Before/After ALRA Imputation",
       subtitle = paste("Normalization method:", selected_normalization_method),
       y = "Sparsity (% zeros)") +
  theme_minimal() +
  theme(legend.position = "none") +
  ylim(0, 100)

ggsave(file.path(alra_plots_dir, "alra_sparsity_comparison.png"), p_sparsity,
       width = 6, height = 5, dpi = 300)

# 2. Gene-level zero comparison
gene_zeros_before <- rowSums(norm_dense == 0)
gene_zeros_after <- rowSums(imputed_matrix == 0)

zeros_comparison <- data.frame(
  Gene = rownames(norm_dense),
  Before = gene_zeros_before,
  After = gene_zeros_after,
  Reduction = gene_zeros_before - gene_zeros_after
)

p_zeros <- ggplot(zeros_comparison, aes(x = Before, y = After)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Zero Counts per Gene: Before vs After ALRA",
       subtitle = paste("Normalization:", selected_normalization_method),
       x = "Zeros per gene (before ALRA)",
       y = "Zeros per gene (after ALRA)") +
  theme_minimal()

ggsave(file.path(alra_plots_dir, "alra_gene_zeros_comparison.png"), p_zeros,
       width = 8, height = 6, dpi = 300)

# 3. Distribution of imputation effect
p_reduction <- ggplot(zeros_comparison, aes(x = Reduction)) +
  geom_histogram(bins = 50, fill = "#3498DB", alpha = 0.7) +
  geom_vline(xintercept = mean(zeros_comparison$Reduction), color = "red", linetype = "dashed") +
  labs(title = "Distribution of Zero-Count Reduction per Gene",
       subtitle = paste("Mean reduction:", round(mean(zeros_comparison$Reduction), 1), "zeros/gene"),
       x = "Reduction in zero counts",
       y = "Number of genes") +
  theme_minimal()

ggsave(file.path(alra_plots_dir, "alra_zero_reduction_distribution.png"), p_reduction,
       width = 8, height = 5, dpi = 300)

# 4. Value distribution comparison (sample of genes)
set.seed(42)
sample_genes <- sample(rownames(norm_dense), min(100, nrow(norm_dense)))

value_data <- data.frame(
  Value = c(as.vector(norm_dense[sample_genes, ]),
            as.vector(imputed_matrix[sample_genes, ])),
  Stage = rep(c("Original", "ALRA"), each = length(sample_genes) * ncol(norm_dense))
)

p_values <- ggplot(value_data, aes(x = Value, fill = Stage)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Expression Values",
       subtitle = paste("Sample of", length(sample_genes), "genes"),
       x = "Normalized expression value",
       y = "Density") +
  theme_minimal()

ggsave(file.path(alra_plots_dir, "alra_value_distribution.png"), p_values,
       width = 8, height = 5, dpi = 300)

cat("QC plots saved to:", alra_plots_dir, "\n\n")

# ==============================================================================
# SAVE STATISTICS
# ==============================================================================

cat("--- Saving statistics ---\n")

# Overall statistics
stats_overall <- data.frame(
  Metric = c(
    "Normalization method",
    "Total genes",
    "Total cells",
    "SVD rank (k)",
    "Power iterations (q)",
    "Quantile probability",
    "Pre-ALRA sparsity (%)",
    "Post-ALRA sparsity (%)",
    "Sparsity reduction (pp)",
    "Mean zeros reduced per gene",
    "Runtime (minutes)"
  ),
  Value = c(
    selected_normalization_method,
    nrow(seurat_obj),
    ncol(seurat_obj),
    k_value,
    params$alra_q,
    params$alra_quantile_prob,
    round(pre_sparsity, 2),
    round(post_sparsity, 2),
    round(pre_sparsity - post_sparsity, 2),
    round(mean(zeros_comparison$Reduction), 2),
    round(as.numeric(alra_duration), 2)
  )
)

write.csv(stats_overall, file.path(alra_tables_dir, "alra_summary_statistics.csv"), row.names = FALSE)

# Per-gene statistics
write.csv(zeros_comparison, file.path(alra_tables_dir, "alra_per_gene_statistics.csv"), row.names = FALSE)

cat("Statistics saved to:", alra_tables_dir, "\n\n")

# ==============================================================================
# SAVE OUTPUT OBJECTS
# ==============================================================================

cat("--- Saving output objects ---\n")

# Save standard output file
output_file <- file.path(alra_out_dir, "normalized_object_with_alra.rds")
saveRDS(seurat_obj, output_file)
cat("Saved:", output_file, "\n")

# ==============================================================================
# UPDATE MODULE 03 OUTPUT FOR MODULE 04 COMPATIBILITY
# ==============================================================================
# Module 04 loads:
#   1. 03_normalization_data.RData (contains normalized objects)
#   2. merged_normalized/merged_{method}_unintegrated.rds
# We need to update both with the ALRA assay.
# ==============================================================================

cat("\n--- Updating Module 03 outputs for Module 04 compatibility ---\n")

# Update the merged_normalized file
unintegrated_file <- file.path(merged_norm_dir,
                               paste0("merged_", selected_normalization_method, "_unintegrated.rds"))
if (file.exists(unintegrated_file)) {
  saveRDS(seurat_obj, unintegrated_file)
  cat("Updated:", unintegrated_file, "\n")
}

# Update 03_normalization_data.RData
if (file.exists(norm_data_file)) {
  # Load existing data
  load(norm_data_file)

  # Update the appropriate object based on normalization method
  if (selected_normalization_method == "scran" && exists("scran_object")) {
    scran_object <- seurat_obj
    cat("Updated scran_object with ALRA assay\n")
  } else if (selected_normalization_method == "LogNormalize" && exists("lognorm_object")) {
    lognorm_object <- seurat_obj
    cat("Updated lognorm_object with ALRA assay\n")
  }

  # Save updated normalization data
  save_list <- c("selected_normalization_method")
  if (exists("sct_object")) save_list <- c(save_list, "sct_object")
  if (exists("scran_object")) save_list <- c(save_list, "scran_object")
  if (exists("lognorm_object")) save_list <- c(save_list, "lognorm_object")
  if (exists("sample_objects")) save_list <- c(save_list, "sample_objects")
  if (exists("normalization_benchmark")) save_list <- c(save_list, "normalization_benchmark")

  save(list = save_list, file = norm_data_file)
  cat("Updated:", norm_data_file, "\n")
}

# ==============================================================================
# UPDATE ENVIRONMENT
# ==============================================================================

# Save updated environment with ALRA info
alra_completed <- TRUE
alra_stats <- stats_overall
alra_k_used <- k_value

save(
  params,
  output_dirs,
  subdirs,
  out_base,
  pipeline_dir,
  input_validation,
  has_SeuratIntegrate,
  has_CHOIR,
  has_glmGamPoi,
  has_DESeq2,
  has_harmony,
  has_lisi,
  has_kBET,
  has_batchelor,
  has_SeuratWrappers,
  has_MAST,
  has_clustree,
  has_bluster,
  has_variancePartition,
  has_alra,
  alra_completed,
  alra_stats,
  alra_k_used,
  python_packages_ok,
  has_afmf,
  afmf_python,
  has_julia,
  julia_version,
  julia_bin,
  has_scice_env,
  scice_env,
  scice_pkg_dir,
  has_scice_source,
  scice_source,
  has_scice_packages,
  has_scice,
  file = env_data_file
)

cat("Updated pipeline environment\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("MODULE 03b SUMMARY: ALRA IMPUTATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Input:\n")
cat("  Normalized object:", input_file, "\n")
cat("  Normalization method:", selected_normalization_method, "\n")
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")

cat("\nALRA parameters:\n")
cat("  k (SVD rank):", k_value, "\n")
cat("  q (power iterations):", params$alra_q, "\n")
cat("  quantile_prob:", params$alra_quantile_prob, "\n")

cat("\nResults:\n")
cat("  Pre-ALRA sparsity:", round(pre_sparsity, 2), "%\n")
cat("  Post-ALRA sparsity:", round(post_sparsity, 2), "%\n")
cat("  Sparsity reduction:", round(pre_sparsity - post_sparsity, 2), "pp\n")
cat("  Mean zeros reduced per gene:", round(mean(zeros_comparison$Reduction), 1), "\n")
cat("  Runtime:", round(as.numeric(alra_duration), 2), "minutes\n")

cat("\nOutput:\n")
cat("  ALRA object:", output_file, "\n")
cat("  New assay: 'ALRA'\n")
cat("  Plots:", alra_plots_dir, "\n")
cat("  Statistics:", alra_tables_dir, "\n")

cat("\nModule 04 compatibility:\n")
cat("  Updated:", unintegrated_file, "\n")
cat("  Updated:", norm_data_file, "\n")

cat("\nDownstream usage (Module 04+):\n")
cat("  use_alra_for_downstream:", isTRUE(params$use_alra_for_downstream), "\n")
if (isTRUE(params$use_alra_for_downstream)) {
  cat("  -> Module 04+ will use 'ALRA' assay for downstream analysis\n")
} else {
  cat("  -> Module 04+ will use original normalized data (ALRA available for comparison)\n")
}

cat("\n>>> MODULE 03b COMPLETE <<<\n")