# ==============================================================================
# MODULE 02b: afMF IMPUTATION (COUNTS-BASED)
# ==============================================================================
#
# This module performs dropout imputation using afMF (adaptive low-rank Full
# Matrix Factorization) on raw counts BEFORE normalization.
#
# WHEN TO RUN:
#   - When params$imputation_method = "afmf" or "both"
#   - Runs AFTER Module 02 (QC/Merge) and BEFORE Module 03 (Normalization)
#
# INPUT:
#   - 02_qc_data.RData from Module 02 (contains filtered_obj and sample_objects)
#
# OUTPUT:
#   - Updated 02_qc_data.RData with imputed assay added to sample_objects
#   - merged_object_imputed.rds with "imputed" assay containing afMF counts
#   - QC plots comparing original vs imputed
#   - Statistics on imputation effects
#
# IMPUTATION STRATEGY:
#   - Imputation is performed PER-SAMPLE to avoid information leakage
#   - This is critical for downstream DE analysis (e.g., Male vs Female)
#   - Imputing on merged data could blur biological differences between groups
#
# afMF API (from https://github.com/GO3295/SCImputation):
#   - Input: pandas DataFrame (genes x cells)
#   - Output: pandas DataFrame (genes x cells)
#   - NO keyword arguments (max_iter, tol are NOT parameters of afMF function)
#   - Uses helper script: utils/python/run_afmf_imputation.py
#
# DOWNSTREAM USAGE:
#   - If use_afmf_for_normalization = TRUE: Module 03 uses imputed counts
#   - If use_afmf_for_normalization = FALSE: Module 03 uses original counts
#     (imputed assay still available for comparison)
#
# UPDATES:
# - 2026-01-10: Fixed to use helper script run_afmf_imputation.py with correct API
# - 2026-01-10: Fixed input to load from 02_qc_data.RData (not merged_object.rds)
# - 2026-01-10: Changed to per-sample imputation to avoid information leakage
# - 2026-01-09: Renamed parameters (run_imputation -> imputation_method)
# - 2026-01-09: Added imputation_method check ("afmf" or "both")
# - 2026-01-09: Updated output directory to 02b_Imputation_afMF
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 02b: afMF IMPUTATION (COUNTS-BASED)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

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
# CHECK IF afMF IMPUTATION SHOULD RUN
# ==============================================================================

# Check imputation_method parameter
run_afmf_imputation <- params$imputation_method %in% c("afmf", "both")

if (!run_afmf_imputation) {
  cat(">>> afMF IMPUTATION SKIPPED <<<\n")
  cat("    imputation_method =", params$imputation_method, "\n")
  cat("    Set imputation_method = 'afmf' or 'both' to enable afMF imputation\n")
  cat("\n>>> MODULE 02b COMPLETE (skipped) <<<\n")
  quit(save = "no", status = 0)
}

cat("afMF imputation enabled (imputation_method =", params$imputation_method, ")\n\n")

# ==============================================================================
# CHECK afMF AVAILABILITY
# ==============================================================================

if (!has_afmf) {
  cat(">>> afMF IMPUTATION SKIPPED <<<\n")
  cat("    afMF Python environment not available\n")
  cat("    To install afMF:\n")
  cat("      conda create -n afMF_SCImputation_env python=3.10 -y\n")
  cat("      conda activate afMF_SCImputation_env\n")
  cat("      pip install numpy scipy pandas scikit-learn\n")
  cat("      cd ~/GITHUB_repositories/SCImputation/afMF && pip install .\n")
  cat("\n>>> MODULE 02b COMPLETE (skipped due to missing afMF) <<<\n")
  quit(save = "no", status = 0)
}

# ==============================================================================
# LOCATE HELPER SCRIPT
# ==============================================================================
# The helper script run_afmf_imputation.py correctly calls afMF(dat) with
# a pandas DataFrame and NO keyword arguments.
# ==============================================================================

afmf_script <- file.path(pipeline_dir, "utils", "python", "run_afmf_imputation.py")

if (!file.exists(afmf_script)) {
  # Try alternative locations
  alt_paths <- c(
    file.path(params$pipeline_scripts_dir, "utils", "python", "run_afmf_imputation.py"),
    file.path(dirname(pipeline_dir), "utils", "python", "run_afmf_imputation.py")
  )
  for (alt in alt_paths) {
    if (file.exists(alt)) {
      afmf_script <- alt
      break
    }
  }
}

if (!file.exists(afmf_script)) {
  stop("afMF helper script not found at: ", afmf_script,
       "\nExpected location: ", pipeline_dir, "/utils/python/run_afmf_imputation.py")
}

cat("Using afMF helper script:", afmf_script, "\n\n")

# ==============================================================================
# LOAD INPUT DATA
# ==============================================================================

cat("--- Loading QC-filtered data from Module 02 ---\n")

# Load 02_qc_data.RData which contains filtered_obj and sample_objects
qc_data_file <- file.path(output_dirs$objects, "02_qc_data.RData")
if (!file.exists(qc_data_file)) {
  stop("Input file not found: ", qc_data_file, "\nRun Module 02 first.")
}

load(qc_data_file)
cat("Loaded:", qc_data_file, "\n")

# Verify we have sample_objects
if (!exists("sample_objects") || length(sample_objects) == 0) {
  stop("sample_objects not found in 02_qc_data.RData. Run Module 02 first.")
}

sample_names <- names(sample_objects)
cat("Samples to impute:", paste(sample_names, collapse = ", "), "\n")
cat("Total cells in merged object:", ncol(filtered_obj), "\n")
for (samp in sample_names) {
  cat("  ", samp, ":", ncol(sample_objects[[samp]]), "cells\n")
}
cat("\n")

# ==============================================================================
# SETUP OUTPUT DIRECTORIES
# ==============================================================================

afmf_out_dir <- output_dirs$imputation_afmf
afmf_plots_dir <- subdirs$afmf_plots
afmf_tables_dir <- subdirs$afmf_tables

dir.create(afmf_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(afmf_plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(afmf_tables_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# PER-SAMPLE afMF IMPUTATION
# ==============================================================================
# CRITICAL: Imputation is performed on each sample separately to avoid
# information leakage between biological groups (e.g., Male vs Female).
# This ensures that downstream DE analysis is not biased.
# ==============================================================================

cat("--- Per-sample afMF imputation ---\n")
cat("NOTE: Imputing each sample separately to preserve biological differences.\n")
cat("      This is critical for accurate downstream DE analysis.\n\n")

# Store imputed sample objects
imputed_sample_objects <- list()

# Store per-sample statistics
all_sample_stats <- list()

# Overall statistics
total_pre_zeros <- 0
total_post_zeros <- 0
total_elements <- 0

for (samp in sample_names) {
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("Processing sample:", samp, "\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  seurat_obj <- sample_objects[[samp]]
  cat("Cells:", ncol(seurat_obj), "\n")
  cat("Genes:", nrow(seurat_obj), "\n")
  
  # ==========================================================================
  # PREPARE COUNT MATRIX FOR afMF
  # ==========================================================================
  
  cat("\n  Preparing count matrix for afMF...\n")
  
  # Get counts layer
  counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  
  # Convert to dense matrix if sparse
  if (inherits(counts_matrix, "sparseMatrix")) {
    cat("  Converting sparse matrix to dense for afMF...\n")
    counts_dense <- as.matrix(counts_matrix)
  } else {
    counts_dense <- counts_matrix
  }
  
  cat("  Matrix dimensions:", nrow(counts_dense), "genes x", ncol(counts_dense), "cells\n")
  
  # Calculate pre-imputation statistics
  pre_zeros <- sum(counts_dense == 0)
  pre_total <- length(counts_dense)
  pre_sparsity <- pre_zeros / pre_total * 100
  
  cat("  Pre-imputation sparsity:", round(pre_sparsity, 2), "%\n")
  cat("  Non-zero entries:", pre_total - pre_zeros, "\n")
  
  # ==========================================================================
  # FILTER GENES FOR IMPUTATION
  # ==========================================================================
  
  cat("\n  Filtering genes for imputation...\n")
  
  min_cells <- params$afmf_min_cells_expressing
  cells_expressing <- rowSums(counts_dense > 0)
  genes_to_impute <- cells_expressing >= min_cells
  
  cat("  Minimum cells expressing:", min_cells, "\n")
  cat("  Genes passing filter:", sum(genes_to_impute), "/", length(genes_to_impute), "\n")
  
  # Subset matrix for imputation
  counts_for_imputation <- counts_dense[genes_to_impute, , drop = FALSE]
  cat("  Matrix for imputation:", nrow(counts_for_imputation), "genes x",
      ncol(counts_for_imputation), "cells\n")
  
  # ==========================================================================
  # SAVE TEMPORARY FILES FOR PYTHON
  # ==========================================================================
  
  cat("\n  Preparing files for Python afMF...\n")
  
  temp_dir <- file.path(afmf_out_dir, "temp", samp)
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save count matrix as CSV (genes x cells, with row/col names)
  temp_counts_file <- file.path(temp_dir, "counts_matrix.csv")
  write.csv(counts_for_imputation, temp_counts_file, row.names = TRUE)
  cat("  Saved counts to:", temp_counts_file, "\n")
  
  # Output file for imputed counts
  temp_output_file <- file.path(temp_dir, "imputed_counts.csv")
  
  # Log file
  log_file <- file.path(temp_dir, "afmf_log.txt")
  
  # ==========================================================================
  # RUN afMF IMPUTATION VIA HELPER SCRIPT
  # ==========================================================================
  # The helper script run_afmf_imputation.py correctly calls:
  #   imputed_dat = afMF(dat)
  # where dat is a pandas DataFrame. NO keyword arguments.
  # ==========================================================================
  
  cat("\n  Running afMF imputation...\n")
  cat("  Python:", afmf_python, "\n")
  cat("  Helper script:", afmf_script, "\n")
  
  afmf_start <- Sys.time()
  
  # Call the helper script with proper arguments
  # The script handles: loading, filtering, imputation, saving
  result <- system2(
    command = afmf_python,
    args = c(
      afmf_script,
      temp_counts_file,
      temp_output_file,
      "--min-genes", "0",  # Already filtered in R
      "--log-file", log_file
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  afmf_duration <- difftime(Sys.time(), afmf_start, units = "mins")
  
  # Print Python output
  cat("\n  --- Python output ---\n")
  cat(paste(result, collapse = "\n"), "\n")
  cat("  ---------------------\n")
  
  # Check exit status
  exit_status <- attr(result, "status")
  if (!is.null(exit_status) && exit_status != 0) {
    cat("  ERROR: afMF imputation failed for sample", samp, "\n")
    cat("  Exit status:", exit_status, "\n")
    if (file.exists(log_file)) {
      cat("  Log file contents:\n")
      cat(readLines(log_file), sep = "\n")
    }
    cat("  Skipping this sample and continuing with others...\n")
    imputed_sample_objects[[samp]] <- seurat_obj  # Keep original
    next
  }
  
  # Check if output file was created
  if (!file.exists(temp_output_file)) {
    cat("  ERROR: afMF imputation failed for sample", samp, "- output file not created\n")
    cat("  Skipping this sample and continuing with others...\n")
    imputed_sample_objects[[samp]] <- seurat_obj  # Keep original
    next
  }
  
  cat("  afMF completed in", round(as.numeric(afmf_duration), 2), "minutes\n")
  
  # ==========================================================================
  # LOAD IMPUTED COUNTS AND CREATE ASSAY
  # ==========================================================================
  
  cat("\n  Loading imputed counts...\n")
  
  imputed_df <- read.csv(temp_output_file, row.names = 1, check.names = FALSE)
  imputed_matrix <- as.matrix(imputed_df)
  
  cat("  Imputed matrix dimensions:", nrow(imputed_matrix), "x", ncol(imputed_matrix), "\n")
  
  # Ensure non-negative values
  imputed_matrix[imputed_matrix < 0] <- 0
  
  # Round to integers for count data
  imputed_matrix <- round(imputed_matrix)
  
  # Calculate post-imputation statistics for imputed genes
  post_zeros_imputed <- sum(imputed_matrix == 0)
  post_total_imputed <- length(imputed_matrix)
  post_sparsity_imputed <- post_zeros_imputed / post_total_imputed * 100
  
  cat("  Post-imputation sparsity (imputed genes):", round(post_sparsity_imputed, 2), "%\n")
  
  # ==========================================================================
  # CREATE FULL IMPUTED MATRIX (INCLUDING NON-IMPUTED GENES)
  # ==========================================================================
  
  cat("\n  Creating full imputed matrix...\n")
  
  # Start with original counts
  full_imputed_matrix <- counts_dense
  
  # Replace imputed genes
  full_imputed_matrix[genes_to_impute, ] <- imputed_matrix
  
  # Calculate overall post-imputation statistics
  post_zeros <- sum(full_imputed_matrix == 0)
  post_total <- length(full_imputed_matrix)
  post_sparsity <- post_zeros / post_total * 100
  
  cat("  Full matrix post-imputation sparsity:", round(post_sparsity, 2), "%\n")
  cat("  Sparsity reduction:", round(pre_sparsity - post_sparsity, 2), "percentage points\n")
  
  # Update totals for overall statistics
  total_pre_zeros <- total_pre_zeros + pre_zeros
  total_post_zeros <- total_post_zeros + post_zeros
  total_elements <- total_elements + pre_total
  
  # Convert back to sparse matrix for storage efficiency
  full_imputed_sparse <- Matrix::Matrix(full_imputed_matrix, sparse = TRUE)
  
  cat("  Full imputed matrix:", nrow(full_imputed_sparse), "genes x",
      ncol(full_imputed_sparse), "cells\n")
  
  # ==========================================================================
  # ADD IMPUTED ASSAY TO SEURAT OBJECT
  # ==========================================================================
  
  cat("\n  Adding imputed assay to Seurat object...\n")
  
  # Create new assay with imputed counts
  imputed_assay <- CreateAssay5Object(counts = full_imputed_sparse)
  
  # Add to Seurat object
  seurat_obj[["imputed"]] <- imputed_assay
  
  cat("  Added 'imputed' assay to Seurat object\n")
  cat("  Available assays:", paste(names(seurat_obj@assays), collapse = ", "), "\n")
  
  # Store in list
  imputed_sample_objects[[samp]] <- seurat_obj
  
  # ==========================================================================
  # STORE PER-SAMPLE STATISTICS
  # ==========================================================================
  
  # Gene-level statistics for this sample
  gene_zeros_before <- rowSums(counts_for_imputation == 0)
  gene_zeros_after <- rowSums(imputed_matrix == 0)
  
  zeros_comparison <- data.frame(
    Gene = rownames(counts_for_imputation),
    Before = gene_zeros_before,
    After = gene_zeros_after,
    Reduction = gene_zeros_before - gene_zeros_after
  )
  
  all_sample_stats[[samp]] <- list(
    n_cells = ncol(seurat_obj),
    n_genes = nrow(seurat_obj),
    genes_imputed = sum(genes_to_impute),
    pre_sparsity = pre_sparsity,
    post_sparsity = post_sparsity,
    sparsity_reduction = pre_sparsity - post_sparsity,
    mean_zeros_reduced = mean(zeros_comparison$Reduction),
    runtime_minutes = as.numeric(afmf_duration),
    zeros_comparison = zeros_comparison
  )
  
  # ==========================================================================
  # CLEANUP TEMPORARY FILES FOR THIS SAMPLE
  # ==========================================================================
  
  unlink(temp_dir, recursive = TRUE)
  cat("  Removed temporary directory:", temp_dir, "\n")
}

# ==============================================================================
# MERGE IMPUTED SAMPLE OBJECTS
# ==============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("MERGING IMPUTED SAMPLES\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

if (length(imputed_sample_objects) == 1) {
  merged_imputed_obj <- imputed_sample_objects[[1]]
} else {
  merged_imputed_obj <- merge(
    imputed_sample_objects[[1]], 
    imputed_sample_objects[2:length(imputed_sample_objects)],
    add.cell.ids = sample_names
  )
}

cat("Merged imputed object:", ncol(merged_imputed_obj), "cells\n")
cat("Available assays:", paste(names(merged_imputed_obj@assays), collapse = ", "), "\n\n")

# ==============================================================================
# GENERATE QC PLOTS
# ==============================================================================

cat("--- Generating QC plots ---\n")

# Calculate overall statistics
total_pre_sparsity <- total_pre_zeros / total_elements * 100
total_post_sparsity <- total_post_zeros / total_elements * 100

# 1. Overall sparsity comparison
sparsity_data <- data.frame(
  Stage = c("Original", "After afMF"),
  Sparsity = c(total_pre_sparsity, total_post_sparsity)
)

p_sparsity <- ggplot(sparsity_data, aes(x = Stage, y = Sparsity, fill = Stage)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Sparsity)), vjust = -0.5) +
  scale_fill_manual(values = c("Original" = "#E74C3C", "After afMF" = "#27AE60")) +
  labs(title = "Matrix Sparsity Before/After afMF Imputation",
       subtitle = paste("Per-sample imputation across", length(sample_names), "samples"),
       y = "Sparsity (% zeros)") +
  theme_minimal() +
  theme(legend.position = "none") +
  ylim(0, 100)

ggsave(file.path(afmf_plots_dir, "sparsity_comparison.png"), p_sparsity,
       width = 6, height = 5, dpi = 300)

# 2. Per-sample sparsity comparison
sample_sparsity_data <- do.call(rbind, lapply(names(all_sample_stats), function(samp) {
  data.frame(
    Sample = samp,
    Stage = c("Original", "After afMF"),
    Sparsity = c(all_sample_stats[[samp]]$pre_sparsity, 
                 all_sample_stats[[samp]]$post_sparsity)
  )
}))

p_sample_sparsity <- ggplot(sample_sparsity_data, aes(x = Sample, y = Sparsity, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Original" = "#E74C3C", "After afMF" = "#27AE60")) +
  labs(title = "Per-Sample Sparsity Before/After afMF Imputation",
       y = "Sparsity (% zeros)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 100)

ggsave(file.path(afmf_plots_dir, "per_sample_sparsity_comparison.png"), p_sample_sparsity,
       width = 8, height = 5, dpi = 300)

# 3. Combined gene-level zeros comparison (across all samples)
all_zeros_comparison <- do.call(rbind, lapply(names(all_sample_stats), function(samp) {
  df <- all_sample_stats[[samp]]$zeros_comparison
  df$Sample <- samp
  df
}))

p_zeros <- ggplot(all_zeros_comparison, aes(x = Before, y = After)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  facet_wrap(~ Sample) +
  labs(title = "Zero Counts per Gene: Before vs After afMF",
       subtitle = "Per-sample imputation",
       x = "Zeros per gene (before)",
       y = "Zeros per gene (after)") +
  theme_minimal()

ggsave(file.path(afmf_plots_dir, "gene_zeros_comparison.png"), p_zeros,
       width = 10, height = 8, dpi = 300)

# 4. Distribution of imputation effect
p_reduction <- ggplot(all_zeros_comparison, aes(x = Reduction, fill = Sample)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  labs(title = "Distribution of Zero-Count Reduction per Gene",
       subtitle = paste("Mean reduction:", round(mean(all_zeros_comparison$Reduction), 1), "zeros/gene"),
       x = "Reduction in zero counts",
       y = "Number of genes") +
  theme_minimal()

ggsave(file.path(afmf_plots_dir, "zero_reduction_distribution.png"), p_reduction,
       width = 8, height = 5, dpi = 300)

cat("QC plots saved to:", afmf_plots_dir, "\n\n")

# ==============================================================================
# SAVE STATISTICS
# ==============================================================================

cat("--- Saving statistics ---\n")

# Overall statistics
stats_overall <- data.frame(
  Metric = c(
    "Total samples",
    "Total genes",
    "Total cells",
    "Pre-imputation sparsity (%)",
    "Post-imputation sparsity (%)",
    "Sparsity reduction (pp)",
    "Mean zeros reduced per gene",
    "Total runtime (minutes)"
  ),
  Value = c(
    length(sample_names),
    nrow(merged_imputed_obj),
    ncol(merged_imputed_obj),
    round(total_pre_sparsity, 2),
    round(total_post_sparsity, 2),
    round(total_pre_sparsity - total_post_sparsity, 2),
    round(mean(all_zeros_comparison$Reduction), 2),
    round(sum(sapply(all_sample_stats, function(x) x$runtime_minutes)), 2)
  )
)

write.csv(stats_overall, file.path(afmf_tables_dir, "afmf_summary_statistics.csv"), row.names = FALSE)

# Per-sample statistics
per_sample_stats <- do.call(rbind, lapply(names(all_sample_stats), function(samp) {
  stats <- all_sample_stats[[samp]]
  data.frame(
    Sample = samp,
    Cells = stats$n_cells,
    Genes = stats$n_genes,
    Genes_Imputed = stats$genes_imputed,
    Pre_Sparsity = round(stats$pre_sparsity, 2),
    Post_Sparsity = round(stats$post_sparsity, 2),
    Sparsity_Reduction = round(stats$sparsity_reduction, 2),
    Mean_Zeros_Reduced = round(stats$mean_zeros_reduced, 2),
    Runtime_Minutes = round(stats$runtime_minutes, 2)
  )
}))

write.csv(per_sample_stats, file.path(afmf_tables_dir, "afmf_per_sample_statistics.csv"), row.names = FALSE)

# Per-gene statistics
write.csv(all_zeros_comparison, file.path(afmf_tables_dir, "afmf_per_gene_statistics.csv"), row.names = FALSE)

cat("Statistics saved to:", afmf_tables_dir, "\n\n")

# ==============================================================================
# SAVE OUTPUT OBJECTS
# ==============================================================================

cat("--- Saving output objects ---\n")

# Save merged imputed object
output_file <- file.path(output_dirs$objects, "merged_object_imputed.rds")
saveRDS(merged_imputed_obj, output_file)
cat("Saved:", output_file, "\n")
cat("File size:", round(file.info(output_file)$size / (1024^2), 2), "MB\n")

# Update sample_objects with imputed versions for Module 03 compatibility
sample_objects <- imputed_sample_objects

# Update filtered_obj with the merged imputed object
filtered_obj <- merged_imputed_obj

# Save updated 02_qc_data.RData for Module 03 to use
qc_data_file_updated <- file.path(output_dirs$objects, "02_qc_data.RData")
save(filtered_obj, sample_objects, filtering_summary, filtering_breakdown, file = qc_data_file_updated)
cat("Updated:", qc_data_file_updated, "\n")
cat("  - sample_objects now contain 'imputed' assay\n")
cat("  - filtered_obj now contains merged imputed data\n\n")

# ==============================================================================
# UPDATE ENVIRONMENT
# ==============================================================================

# Save updated environment with imputation info
afmf_completed <- TRUE
afmf_stats <- stats_overall
afmf_per_sample_stats <- per_sample_stats

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
  python_packages_ok,
  has_afmf,
  afmf_python,
  afmf_completed,
  afmf_stats,
  afmf_per_sample_stats,
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
cat("MODULE 02b SUMMARY: afMF IMPUTATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Input:\n")
cat("  QC data file:", qc_data_file, "\n")
cat("  Samples:", paste(sample_names, collapse = ", "), "\n")
cat("  Total cells:", ncol(merged_imputed_obj), "\n")
cat("  Total genes:", nrow(merged_imputed_obj), "\n")

cat("\nImputation strategy:\n")
cat("  Method: Per-sample imputation (avoids information leakage)\n")
cat("  Min cells expressing threshold:", params$afmf_min_cells_expressing, "\n")
cat("  afMF API: afMF(DataFrame) via helper script (no keyword args)\n")

cat("\nPer-sample results:\n")
print(per_sample_stats[, c("Sample", "Cells", "Pre_Sparsity", "Post_Sparsity", "Sparsity_Reduction")])

cat("\nOverall results:\n")
cat("  Pre-imputation sparsity:", round(total_pre_sparsity, 2), "%\n")
cat("  Post-imputation sparsity:", round(total_post_sparsity, 2), "%\n")
cat("  Sparsity reduction:", round(total_pre_sparsity - total_post_sparsity, 2), "pp\n")
cat("  Total runtime:", round(sum(sapply(all_sample_stats, function(x) x$runtime_minutes)), 2), "minutes\n")

cat("\nOutput:\n")
cat("  Merged imputed object:", output_file, "\n")
cat("  Updated QC data:", qc_data_file_updated, "\n")
cat("  New assay: 'imputed'\n")
cat("  Plots:", afmf_plots_dir, "\n")
cat("  Statistics:", afmf_tables_dir, "\n")

cat("\nDownstream usage (Module 03):\n")
cat("  use_afmf_for_normalization:", isTRUE(params$use_afmf_for_normalization), "\n")
if (isTRUE(params$use_afmf_for_normalization)) {
  cat("  -> Module 03 will use 'imputed' assay for normalization\n")
} else {
  cat("  -> Module 03 will use original 'RNA' assay (imputed available for comparison)\n")
}

cat("\n>>> MODULE 02b COMPLETE <<<\n")
