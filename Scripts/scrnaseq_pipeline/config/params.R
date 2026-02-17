
# ==============================================================================
# CONFIGURATION PARAMETERS - SCRNASEQ DOWNSTREAM ANALYSIS PIPELINE
# ==============================================================================
#
# This configuration file controls all aspects of the downstream analysis.
# Parameters are organized following the pipeline module structure.
#
# =============================================================================
# PIPELINE MODULE STRUCTURE
# =============================================================================
#
# Module 00: Environment Setup
# Module 01: Load & QC (filtering, doublet removal)
# Module 02: Merge Samples
# Module 02b: Imputation - afMF (optional, counts-based)
# Module 03: Normalization (SCTransform, scran, LogNormalize)
# Module 03b: Imputation - ALRA (optional, normalized-data-based)
# Module 04: Integration & Benchmarking
# Module 05: CHOIR Clustering
# Module 05b: scCobra
# Module 06: scICE Subclustering
# Module 07: Leiden Clustering
# Module 07b: CLTS Re-normalization
# Module 08: Differential Expression
# Module 09: Visualization
#
# =============================================================================
# KEY ANALYSIS PARAMETERS - QUICK REFERENCE
# =============================================================================
#
# IMPUTATION (Modules 02b and 03b):
# ------------------------------------------------------------------------------
# Two imputation methods are available at different pipeline stages:
#
# 1. afMF (Module 02b) - COUNTS-BASED IMPUTATION:
#    - Runs BEFORE normalization on raw counts
#    - Uses Low-rank Full Matrix Factorization
#    - Output: "imputed" assay with imputed counts
#    - Always runs if imputation_method = "afmf" or "both"
#
# 2. ALRA (Module 03b) - NORMALIZED-DATA IMPUTATION:
#    - Runs AFTER normalization on log-normalized data
#    - Uses Adaptively-thresholded Low Rank Approximation
#    - Output: "ALRA" assay with imputed normalized data
#    - ONLY works with LogNormalize and scran (NOT SCTransform)
#    - Automatically SKIPPED if SCTransform wins benchmarking
#
# IMPUTATION METHOD SELECTION:
#   params$imputation_method <- "both"   # Default: run both methods
#
#   Options:
#     "none"  - No imputation (original data only)
#     "afmf"  - Only afMF (counts-based, Module 02b)
#     "alra"  - Only ALRA (normalized-data, Module 03b)*
#     "both"  - Both methods (afMF in 02b, ALRA in 03b)*
#
#   * ALRA automatically skipped if SCTransform is selected as best
#     normalization method from benchmarking.
#
# DOWNSTREAM USAGE:
#   params$use_afmf_for_normalization <- FALSE  # Use imputed counts in Module 03
#   params$use_alra_for_downstream <- FALSE     # Use ALRA data in Module 04+
#
# INTEGRATION METHOD SELECTION (Module 04):
# ------------------------------------------------------------------------------
#   params$integration_selection_mode <- "balanced"  # RECOMMENDED for most studies
#
#   Options:
#     "batch_removal"  - Minimizes batch_variance (most aggressive correction)
#                        Winner criterion: which.min(batch_variance)
#                        Use when: Batches are TECHNICAL replicates (same biology,
#                                  different processing times/labs/technicians)
#                        Risk: May over-correct and remove real biological signal
#                        Example: Same cell line processed on different days
#
#     "balanced"       - Maximizes composite score across ALL metrics
#                        Winner criterion: which.max(mean(batch_var_norm, asw_norm, lisi_norm))
#                        Use when: Batches have BIOLOGICAL meaning (different conditions,
#                                  sexes, timepoints, treatments, genotypes)
#                        Benefit: Removes technical noise while preserving biology
#                        Example: Male vs Female samples, Treatment vs Control
#                        >>> RECOMMENDED FOR MOST BIOLOGICAL COMPARISONS <<<
#
#     "conservative"   - Prioritizes LISI score (moderate mixing)
#                        Winner criterion: which.max(lisi_norm)
#                        Use when: Preserving subtle biological differences is critical
#                        Example: Rare cell populations, subtle disease states
#
# BATCH VARIABLE:
# ------------------------------------------------------------------------------
#   params$batch_variable <- "sample_name"
#
#   This defines what constitutes a "batch" for integration:
#     "sample_name"  - Each sample is a batch (most common)
#     "batch"        - Use explicit batch column from sample sheet
#     "sex"          - NOT recommended (removes sex differences!)
#
# CLTS RE-NORMALIZATION (Module 07b):
# ------------------------------------------------------------------------------
#   params$clts_clustering_source <- "scice"
#
#   Options:
#     "scice"  - Apply CLTS using scICE subclusters (DEFAULT)
#     "leiden" - Apply CLTS using Leiden clusters
#     "choir"  - Apply CLTS using CHOIR clusters
#     "both"   - Apply CLTS to BOTH scICE and Leiden
#     "all"    - Apply CLTS to scICE, Leiden, AND CHOIR
#
# DIFFERENTIAL EXPRESSION SETTINGS (Module 08):
# ------------------------------------------------------------------------------
#   params$de_comparison_variable <- "sex"     # What to compare
#   params$de_group1 <- "Male"                 # Numerator (positive fold change)
#   params$de_group2 <- "Female"               # Denominator (reference)
#   params$de_comparison_scope <- "whole_dataset"  # or "per_cluster"
#   params$de_object_sources <- c("scice", "scice_redeconv")  # Objects for DE
#
# CLUSTERING METHODS:
# ------------------------------------------------------------------------------
#   params$run_choir_clustering <- TRUE        # Hierarchical clustering
#   params$run_scice_subclustering <- TRUE     # Julia-based subclustering
#   params$run_leiden_clustering <- TRUE       # Graph-based clustering
#
# =============================================================================
#
# QUICK START - CONFIGURE THESE SECTIONS:
#   1. INPUT FILES CONFIGURATION (lines ~150-200) - Where are your input files?
#   2. SAMPLE SHEET (auto-loaded from project root or config folder)
#   3. FILTERING OPTIONS (lines ~350-400) - What filtering to apply?
#   4. ANALYSIS OPTIONS (lines ~500+) - What analyses to run?
#
# OUTPUT LOCATION:
#   Output is automatically placed based on environment variables:
#   - PREPROCESS_DIR -> for reading preprocessing outputs
#   - DOWNSTREAM_DIR -> for writing downstream analysis outputs
#
# See Define_input_paths.txt for detailed documentation.
#
# CONFIGURATOR TAGS:
#   Lines marked with # <<CFG:param_name>> are managed by the pipeline
#   configurator. The configurator only modifies values on tagged lines.
#   All infrastructure, logic, and validation code is preserved.
#
# UPDATES:
# - 2026-02-06: Added <<CFG>> configurator tags to all adjustable parameters
# - 2026-02-05: Added norm_benchmark_integration_methods parameter (was referenced
#               by 03_normalization.R but never defined, defaulting to Harmony-only)
# - 2026-01-09: Added unified imputation_method parameter ("none"/"afmf"/"alra"/"both")
# - 2026-01-09: Added Module 03b ALRA imputation (post-normalization)
# - 2026-01-09: ALRA auto-skips if SCTransform selected from benchmarking
# - 2026-01-09: Added CLTS clustering source selection (clts_clustering_source)
# - 2026-01-09: Added DE object source selection (de_object_sources)
# - 2026-01-09: Added cross-object comparison for DE (de_run_cross_object_comparison)
# - 2026-01-07: Added generic GROUP_ID/GROUP_LABEL filtering system
# - 2026-01-06: Moved functions to functions.R (load_sample_sheet, validate_params, etc.)
# - 2026-01-06: Added integration_selection_mode for balanced vs aggressive selection
# - 2026-01-06: Updated to use PREPROCESS_DIR and DOWNSTREAM_DIR environment vars
# - 2026-01-06: Compatible with new Output_dir_<dataset>/Single_cell_* structure
# - 2026-01-05: Made output path relative to script/project location
# - 2026-01-05: Added flexible input path configuration
# - 2026-01-05: Added comprehensive filtering control flags
#
# ==============================================================================


# ==============================================================================
# SECTION 0: PROJECT ROOT AND ENVIRONMENT DETECTION
# ==============================================================================
# The directories are automatically determined based on:
#   1. Environment variables set by submit_all_jobs.sh:
#      - PROJECT_ROOT, PREPROCESS_DIR, DOWNSTREAM_DIR, DATASET_NAME
#      - GROUP_ID, GROUP_LABEL (new generic filtering)
#      - VENTRICLE_FILTER (backward compatible)
#   2. Location of this params.R file (fallback)
#   3. Current working directory (ultimate fallback)
#
# This makes the pipeline portable - just copy the project folder and run.
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("LOADING PIPELINE CONFIGURATION\n")
cat("================================================================================\n\n")

# ------------------------------------------------------------------------------
# Set global R options for parallel processing (future package)
# ------------------------------------------------------------------------------
# SCTransform and other functions use the future package for parallelization.
# Large datasets can exceed the default 500 MB globals limit, causing failures.
# This must be set early since params.R is sourced by all modules.
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB

# ------------------------------------------------------------------------------
# 0A. Read Environment Variables
# ------------------------------------------------------------------------------

env_project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
env_preprocess_dir <- Sys.getenv("PREPROCESS_DIR", unset = "")
env_downstream_dir <- Sys.getenv("DOWNSTREAM_DIR", unset = "")
env_dataset_name <- Sys.getenv("DATASET_NAME", unset = "")
env_ventricle_filter <- Sys.getenv("VENTRICLE_FILTER", unset = "")
env_group_id <- Sys.getenv("GROUP_ID", unset = "")
env_group_label <- Sys.getenv("GROUP_LABEL", unset = "")

cat("Environment variables detected:\n")
cat("  PROJECT_ROOT:     ", if (env_project_root != "") env_project_root else "<not set>", "\n")
cat("  PREPROCESS_DIR:   ", if (env_preprocess_dir != "") env_preprocess_dir else "<not set>", "\n")
cat("  DOWNSTREAM_DIR:   ", if (env_downstream_dir != "") env_downstream_dir else "<not set>", "\n")
cat("  DATASET_NAME:     ", if (env_dataset_name != "") env_dataset_name else "<not set>", "\n")
cat("  GROUP_ID:         ", if (env_group_id != "") env_group_id else "<not set>", "\n")
cat("  GROUP_LABEL:      ", if (env_group_label != "") env_group_label else "<not set>", "\n")
cat("  VENTRICLE_FILTER: ", if (env_ventricle_filter != "") env_ventricle_filter else "<not set>", "\n")
cat("\n")

# ------------------------------------------------------------------------------
# 0B. Determine Project Root
# ------------------------------------------------------------------------------

project_root <- env_project_root

if (project_root == "" || !dir.exists(project_root)) {
  # Fallback 1: Derive from this script's location
  # params.R is at: {project}/Scripts/scrnaseq_pipeline/config/params.R
  # So project root is 4 levels up

  script_path <- tryCatch({
    if (sys.nframe() > 0) {
      for (i in seq_len(sys.nframe())) {
        ofile <- sys.frame(i)$ofile
        if (!is.null(ofile) && grepl("params\\.R$", ofile)) {
          normalizePath(ofile)
        }
      }
    }
    NULL
  }, error = function(e) NULL)

  if (!is.null(script_path) && file.exists(script_path)) {
    project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))
  } else {
    # Fallback 2: Try common relative paths from current directory
    possible_roots <- c(
      getwd(),
      file.path(getwd(), ".."),
      file.path(getwd(), "..", ".."),
      file.path(getwd(), "..", "..", "..")
    )

    for (root in possible_roots) {
      if (dir.exists(root)) {
        if (file.exists(file.path(root, "submit_all_jobs.sh")) ||
            file.exists(file.path(root, "samplesheet.csv"))) {
          project_root <- normalizePath(root)
          break
        }
      }
    }

    # Ultimate fallback: use current working directory
    if (project_root == "" || !dir.exists(project_root)) {
      project_root <- getwd()
      cat("WARNING: Could not detect project root. Using current directory.\n")
    }
  }
}

cat("Project root:", project_root, "\n")

# ------------------------------------------------------------------------------
# 0C. Determine Dataset Name
# ------------------------------------------------------------------------------

dataset_name <- env_dataset_name

if (dataset_name == "") {
  # Try to read from samplesheet
  samplesheet_path <- file.path(project_root, "samplesheet.csv")
  if (file.exists(samplesheet_path)) {
    tryCatch({
      ss <- read.csv(samplesheet_path, stringsAsFactors = FALSE, nrows = 2)
      if ("dataset_name" %in% colnames(ss) && nrow(ss) > 0) {
        dataset_name <- ss$dataset_name[1]
        cat("Dataset name from samplesheet:", dataset_name, "\n")
      }
    }, error = function(e) NULL)
  }
}

if (dataset_name == "") {
  # Default dataset name
  dataset_name <- "default"
  cat("WARNING: Dataset name not found, using 'default'\n")
}

# ------------------------------------------------------------------------------
# 0D. Source utility functions (needed for load_sample_sheet, validate_params, etc.)
# ------------------------------------------------------------------------------

scripts_dir <- file.path(project_root, "Scripts")
pipeline_scripts_dir <- file.path(scripts_dir, "scrnaseq_pipeline")

utils_file <- file.path(pipeline_scripts_dir, "utils", "functions.R")
if (file.exists(utils_file)) {
  source(utils_file)
} else {
  # Try alternative paths
  alt_paths <- c(
    file.path(getwd(), "utils", "functions.R"),
    file.path(getwd(), "..", "utils", "functions.R"),
    file.path(dirname(sys.frame(1)$ofile), "..", "utils", "functions.R")
  )
  for (alt_path in alt_paths) {
    if (file.exists(alt_path)) {
      source(alt_path)
      break
    }
  }
}


# ==============================================================================
# SECTION 1: INPUT FILES CONFIGURATION
# ==============================================================================
# Configure where your input files are located and how they are named.
# This is the MOST IMPORTANT section to configure for a new experiment.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1A. INPUT DIRECTORY
# ------------------------------------------------------------------------------
# Priority:
#   1. PREPROCESS_DIR environment variable + step folder
#   2. Derived from project_root + dataset_name
#   3. Legacy hardcoded path (for backward compatibility)
#
# UPDATED 2026-01-15: Now reads from Step 8 (DecontX) instead of Step 7 (scCDC)
# The DecontX output contains scCDC + DecontX corrected counts in the
# 'scCDC_corrected' layer (same layer name, updated contents).

if (env_preprocess_dir != "" && dir.exists(env_preprocess_dir)) {
  # Use environment variable from submit_all_jobs.sh
  input_dir <- file.path(env_preprocess_dir, "8_DecontX_correction")
  cat("Input directory from PREPROCESS_DIR:", input_dir, "\n")
} else if (dataset_name != "" && dataset_name != "default") {
  # Derive from dataset name
  input_dir <- file.path(project_root, paste0("Output_dir_", dataset_name),
                         "Single_cell_preprocessed", "8_DecontX_correction")
  cat("Input directory derived from dataset_name:", input_dir, "\n")
} else {
  # Fallback to legacy hardcoded path for backward compatibility
  input_dir <- "/scicore/home/doetsch/kaiser0001/Revision_NatureComm_Sex/Single_cell_data_expression_analysis/CHOIR_Output/scCDC_Counts/All_Separate"  # <<CFG:legacy_input_dir>>
  cat("WARNING: Using legacy hardcoded input directory\n")
  cat("Input directory:", input_dir, "\n")
}

# ------------------------------------------------------------------------------
# 1B. FILE NAMING PATTERN
# ------------------------------------------------------------------------------
# How are your files named? Use {sample_name} as a placeholder.
# The {sample_name} will be replaced with values from your sample sheet.
#
# For new pipeline structure (Steps 1-8):
#   Files are in: 8_DecontX_correction/<sample_name>/<sample_name>_decontX_corrected.rds
#
# For legacy structure:
#   Files are in: input_dir/<sample_name>_qClus_CHOIR_scCDC_corrected.rds
#
# UPDATED 2026-01-15: New pipeline now reads from DecontX output (Step 8)

# Detect which structure we're using
if (env_preprocess_dir != "" || (dataset_name != "" && dataset_name != "default")) {
  # New pipeline structure - DecontX output
  input_file_pattern <- "{sample_name}_decontX_corrected.rds"
  files_in_subdirectories <- TRUE
  cat("Using NEW pipeline file structure (files in sample subdirectories)\n")
} else {
  # Legacy structure
  input_file_pattern <- "{sample_name}_qClus_CHOIR_scCDC_corrected.rds"  # <<CFG:legacy_input_file_pattern>>
  files_in_subdirectories <- FALSE  # <<CFG:legacy_files_in_subdirectories>>
  cat("Using LEGACY file structure (all files in one directory)\n")
}

# ------------------------------------------------------------------------------
# 1C. COUNTS LAYER TO USE
# ------------------------------------------------------------------------------
# Which layer in the Seurat object contains the counts for analysis?
#
#   "auto"             -> Auto-detect (tries scCDC_corrected, then counts, then data)
#   "scCDC_corrected"  -> Use ambient RNA corrected counts (scCDC + DecontX)
#   "counts"           -> Use raw counts
#
# NOTE: The DecontX output stores scCDC + DecontX corrected counts in the
# 'scCDC_corrected' layer (same layer name as before, but now with additional
# DecontX ambient RNA removal applied).

counts_layer_to_use <- "auto"  # <<CFG:counts_layer_to_use>>


# ==============================================================================
# SECTION 2: OUTPUT CONFIGURATION
# ==============================================================================
# Output is placed based on DOWNSTREAM_DIR or derived from project structure.
# ==============================================================================

# ------------------------------------------------------------------------------
# 2A. OUTPUT BASE DIRECTORY
# ------------------------------------------------------------------------------

if (env_downstream_dir != "") {
  # Use environment variable from submit_all_jobs.sh
  output_base_dir <- env_downstream_dir
  cat("Output base from DOWNSTREAM_DIR:", output_base_dir, "\n")
} else if (dataset_name != "" && dataset_name != "default") {
  # Derive from dataset name
  output_base_dir <- file.path(project_root, paste0("Output_dir_", dataset_name),
                                "Single_cell_clustering")
  cat("Output base derived from dataset_name:", output_base_dir, "\n")
} else {
  # Fallback to project root
  output_base_dir <- project_root
  cat("WARNING: Using project_root as output base\n")
}

# ------------------------------------------------------------------------------
# 2B. GROUP-BASED OR VENTRICLE-SPECIFIC ANALYSIS
# ------------------------------------------------------------------------------
# The pipeline supports two filtering systems:
#   1. NEW: Generic group-based filtering (GROUP_ID + GROUP_LABEL)
#   2. LEGACY: Ventricle-specific filtering (VENTRICLE_FILTER)
#
# Group-based filtering takes priority when GROUP_LABEL is set.

group_id <- env_group_id
group_label <- env_group_label
ventricle_filter <- env_ventricle_filter

# Determine the active filter and output label
if (group_label != "") {
  # New generic system: use group_label for output directory naming
  analysis_label <- group_label
  cat(">>> GROUP-BASED ANALYSIS: Group", group_id, "(", group_label, ") <<<\n")
} else if (ventricle_filter != "") {
  # Backward compatibility: validate and use ventricle for output directory
  if (!ventricle_filter %in% c("LV", "4V", "ALL")) {
    stop("Invalid VENTRICLE_FILTER: '", ventricle_filter, "'. Must be 'LV', '4V', or 'ALL'")
  }
  analysis_label <- ventricle_filter
  cat(">>> VENTRICLE-SPECIFIC ANALYSIS:", ventricle_filter, "<<<\n")
} else {
  # No filter: analyze all samples
  analysis_label <- ""
  cat(">>> FULL DATASET ANALYSIS (all samples) <<<\n")
}

# Set output directory based on analysis label
if (analysis_label != "") {
  out_root <- file.path(output_base_dir, paste0("10_Downstream_Analysis_", analysis_label))
} else {
  out_root <- file.path(output_base_dir, "10_Downstream_Analysis")
}

cat("Output directory:", out_root, "\n\n")


# ==============================================================================
# SECTION 3: SAMPLE SHEET CONFIGURATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 3A. SAMPLE SHEET LOCATION
# ------------------------------------------------------------------------------
# Search for sample sheet in common locations

# Search order for sample sheet
sample_sheet_candidates <- c(
  file.path(pipeline_scripts_dir, "config", "sample_sheet.csv"),
  file.path(project_root, "samplesheet.csv"),
  file.path(project_root, "sample_sheet.csv"),
  file.path(project_root, "config", "sample_sheet.csv")
)

sample_sheet_path <- NULL
for (candidate in sample_sheet_candidates) {
  if (file.exists(candidate)) {
    sample_sheet_path <- candidate
    cat("Found sample sheet:", sample_sheet_path, "\n")
    break
  }
}

if (is.null(sample_sheet_path)) {
  sample_sheet_path <- sample_sheet_candidates[1]  # Default location
  cat("WARNING: Sample sheet not found. Expected at:", sample_sheet_path, "\n")
}

# ------------------------------------------------------------------------------
# 3B. LOAD SAMPLE SHEET
# ------------------------------------------------------------------------------
# Note: load_sample_sheet() function is now in functions.R
# If functions.R doesn't have group filtering support, we handle it here

sample_metadata <- tryCatch({
  # Try loading with the updated function that supports group filtering
  if (exists("load_sample_sheet")) {
    loaded_data <- load_sample_sheet(
      sample_sheet_path = sample_sheet_path,
      input_dir = input_dir,
      input_file_pattern = input_file_pattern,
      files_in_subdirectories = files_in_subdirectories,
      ventricle_filter = ventricle_filter
    )

    # Apply group-based filtering if GROUP_ID is set and not already filtered
    if (group_id != "" && "group_id" %in% colnames(loaded_data)) {
      loaded_data <- loaded_data[as.character(loaded_data$group_id) == group_id, ]
      cat("Applied group_id filter:", group_id, "->", nrow(loaded_data), "samples\n")
    }

    loaded_data
  } else {
    stop("load_sample_sheet function not found")
  }
}, error = function(e) {
  cat("WARNING: Could not load sample sheet:", e$message, "\n")
  cat("Using default Vandebroucke sample metadata...\n")

  # Fallback metadata - MODIFY THIS for your specific dataset
  all_samples <- data.frame(
    sample_name = c("M_22w_LV", "M_22w_4V", "F_7w_LV", "F_7w_4V"),
    sex = c("Male", "Male", "Female", "Female"),
    age = c("22w", "22w", "7w", "7w"),
    batch = c("Batch1", "Batch1", "Batch1", "Batch1"),
    ventricle = c("LV", "4V", "LV", "4V"),
    group_id = c("1", "2", "1", "2"),
    group_label = c("Lateral_Ventricle", "Fourth_Ventricle", "Lateral_Ventricle", "Fourth_Ventricle"),
    include = rep(TRUE, 4),
    stringsAsFactors = FALSE
  )

  # Apply group-based filtering (new generic system)
  if (group_id != "" && "group_id" %in% colnames(all_samples)) {
    all_samples <- all_samples[as.character(all_samples$group_id) == group_id, ]
    cat("Applied group_id filter:", group_id, "->", nrow(all_samples), "samples\n")
  } else if (ventricle_filter != "" && ventricle_filter != "ALL") {
    # Backward compatibility: filter by ventricle
    all_samples <- all_samples[all_samples$ventricle == ventricle_filter, ]
    cat("Applied ventricle filter:", ventricle_filter, "->", nrow(all_samples), "samples\n")
  }

  all_samples$input_file <- sapply(all_samples$sample_name, function(s) {
    filename <- gsub("\\{sample_name\\}", s, input_file_pattern)
    if (files_in_subdirectories) {
      file.path(input_dir, s, filename)
    } else {
      file.path(input_dir, filename)
    }
  })

  all_samples
})


# ==============================================================================
# SECTION 4: MODULE 01 - QC AND FILTERING PARAMETERS
# ==============================================================================

# ------------------------------------------------------------------------------
# 4A. MASTER FILTERING CONTROL
# ------------------------------------------------------------------------------
skip_all_filtering <- FALSE  # <<CFG:skip_all_filtering>>

# ------------------------------------------------------------------------------
# 4B. CELL QC FILTERING
# ------------------------------------------------------------------------------
apply_qc_filtering <- TRUE  # <<CFG:apply_qc_filtering>>
filter_by_min_features <- TRUE  # <<CFG:filter_by_min_features>>
filter_by_max_features <- TRUE  # <<CFG:filter_by_max_features>>
filter_by_percent_mt <- TRUE  # <<CFG:filter_by_percent_mt>>

# ------------------------------------------------------------------------------
# 4C. GENE FILTERING
# ------------------------------------------------------------------------------
filter_genes_by_min_cells <- TRUE  # <<CFG:filter_genes_by_min_cells>>

# ------------------------------------------------------------------------------
# 4D. DOUBLET FILTERING
# ------------------------------------------------------------------------------
filter_doublets <- TRUE  # <<CFG:filter_doublets>>
doublet_column <- "doublet_consensus"  # <<CFG:doublet_column>>
doublet_value_to_remove <- "Doublet"  # <<CFG:doublet_value_to_remove>>
use_doublet_vote_threshold <- FALSE  # <<CFG:use_doublet_vote_threshold>>
doublet_vote_column <- "doublet_votes"  # <<CFG:doublet_vote_column>>
doublet_vote_threshold <- 2  # <<CFG:doublet_vote_threshold>>

# ------------------------------------------------------------------------------
# 4E. HEMOGLOBIN FILTERING
# ------------------------------------------------------------------------------
filter_hemoglobin <- TRUE  # <<CFG:filter_hemoglobin>>
max_percent_hb <- 2  # <<CFG:max_percent_hb>>
hemoglobin_pattern <- "^HB[AB]-"  # <<CFG:hemoglobin_pattern>>

# ------------------------------------------------------------------------------
# 4F. QC THRESHOLDS
# ------------------------------------------------------------------------------
min_features <- 200  # <<CFG:min_features>>
max_features <- 10000  # <<CFG:max_features>>
max_percent_mt <- 25  # <<CFG:max_percent_mt>>
min_cells_per_gene <- 3  # <<CFG:min_cells_per_gene>>


# ==============================================================================
# SECTION 5: MAIN PARAMETERS LIST
# ==============================================================================
# Parameters are organized by pipeline module

params <- list(

  # ===========================================================================
  # PROJECT PATHS AND CONFIGURATION
  # ===========================================================================
  project_root = project_root,
  dataset_name = dataset_name,
  input_dir = input_dir,
  input_file_pattern = input_file_pattern,
  files_in_subdirectories = files_in_subdirectories,
  counts_layer_to_use = counts_layer_to_use,
  output_base_dir = output_base_dir,
  out_root = out_root,
  scripts_dir = scripts_dir,
  pipeline_scripts_dir = pipeline_scripts_dir,
  sample_sheet_path = sample_sheet_path,

  # Environment variables (for reference)
  env_preprocess_dir = env_preprocess_dir,
  env_downstream_dir = env_downstream_dir,

  # Legacy compatibility
  base_dir = output_base_dir,
  sccdc_input_dir = input_dir,

  # ===========================================================================
  # GROUP-BASED ANALYSIS (NEW GENERIC SYSTEM)
  # ===========================================================================
  group_id = group_id,
  group_label = group_label,
  analysis_label = analysis_label,

  # Ventricle analysis (backward compatible)
  ventricle_analysis = if (analysis_label != "") analysis_label else "All",
  ventricle_filter = ventricle_filter,

  # ===========================================================================
  # SAMPLE CONFIGURATION
  # ===========================================================================
  sample_metadata = sample_metadata,

  # ===========================================================================
  # COUNTS LAYER
  # ===========================================================================
  counts_layer = counts_layer_to_use,

  # ===========================================================================
  # MODULE 01: QC AND FILTERING PARAMETERS
  # ===========================================================================
  skip_all_filtering = skip_all_filtering,
  apply_qc_filtering = apply_qc_filtering,
  filter_by_min_features = filter_by_min_features,
  filter_by_max_features = filter_by_max_features,
  filter_by_percent_mt = filter_by_percent_mt,
  filter_genes_by_min_cells = filter_genes_by_min_cells,
  filter_doublets = filter_doublets,
  doublet_column = doublet_column,
  doublet_value_to_remove = doublet_value_to_remove,
  use_doublet_vote_threshold = use_doublet_vote_threshold,
  doublet_vote_column = doublet_vote_column,
  doublet_vote_threshold = doublet_vote_threshold,
  filter_hemoglobin = filter_hemoglobin,
  max_percent_hb = max_percent_hb,
  hemoglobin_pattern = hemoglobin_pattern,
  min_features = min_features,
  max_features = max_features,
  max_percent_mt = max_percent_mt,
  min_cells_per_gene = min_cells_per_gene,

  # ===========================================================================
  # MODULE 02b & 03b: IMPUTATION PARAMETERS
  # ===========================================================================
  # Two imputation methods are available:
  #
  # 1. afMF (Module 02b): Counts-based imputation BEFORE normalization
  #    - Uses Low-rank Full Matrix Factorization
  #    - Works on raw counts
  #    - Output: "imputed" assay
  #
  # 2. ALRA (Module 03b): Normalized-data imputation AFTER normalization
  #    - Uses Adaptively-thresholded Low Rank Approximation
  #    - Works on log-normalized data (LogNormalize or scran only)
  #    - Output: "ALRA" assay
  #    - AUTOMATICALLY SKIPPED if SCTransform wins benchmarking
  #
  # IMPUTATION METHOD SELECTION:
  #   "none"  - No imputation
  #   "afmf"  - Only afMF (Module 02b)
  #   "alra"  - Only ALRA (Module 03b, if compatible normalization)
  #   "both"  - Both methods (DEFAULT)
  # ===========================================================================

  imputation_method = "both",  # Options: "none", "afmf", "alra", "both"  # <<CFG:imputation_method>>

  # ---------------------------------------------------------------------------
  # afMF Parameters (Module 02b - counts-based)
  # ---------------------------------------------------------------------------
  # afMF: Low-rank Full Matrix Factorization for dropout imputation
  # Reference: https://github.com/GO3295/SCImputation
  #
  # afMF runs on RAW COUNTS before normalization.
  # The imputed object is created when imputation_method = "afmf" or "both".
  # use_afmf_for_normalization controls which counts Module 03 uses.
  # ---------------------------------------------------------------------------

  # Use afMF-imputed counts for normalization in Module 03?
  # FALSE = use original scCDC counts (DEFAULT, recommended for DE)
  # TRUE = use afMF-imputed counts
  use_afmf_for_normalization = FALSE,  # <<CFG:use_afmf_for_normalization>>

  # Python environment for afMF (dedicated conda env)
  afmf_python = file.path(Sys.getenv("HOME"), ".conda/envs/afMF_SCImputation_env/bin/python"),  # <<CFG:afmf_python>>

  # afMF algorithm parameters
  afmf_max_iter = 100,  # Maximum iterations for convergence  # <<CFG:afmf_max_iter>>
  afmf_tol = 0.00001,  # Convergence tolerance  # <<CFG:afmf_tol>>
  afmf_min_cells_expressing = 10,  # Min cells expressing gene for imputation  # <<CFG:afmf_min_cells_expressing>>

  # ---------------------------------------------------------------------------
  # ALRA Parameters (Module 03b - normalized-data-based)
  # ---------------------------------------------------------------------------
  # ALRA: Adaptively-thresholded Low Rank Approximation
  # Reference: https://github.com/KlugerLab/ALRA
  #
  # ALRA runs on LOG-NORMALIZED DATA after Module 03.
  # ONLY compatible with LogNormalize and scran (NOT SCTransform).
  # If SCTransform wins normalization benchmarking, ALRA is automatically skipped.
  #
  # The ALRA-imputed data is stored in "ALRA" assay.
  # use_alra_for_downstream controls which data Module 04+ uses.
  # ---------------------------------------------------------------------------

  # Use ALRA-imputed data for downstream analysis (Module 04+)?
  # FALSE = use original normalized data (DEFAULT, recommended for DE)
  # TRUE = use ALRA-imputed data
  use_alra_for_downstream = FALSE,  # <<CFG:use_alra_for_downstream>>

  # ALRA algorithm parameters
  alra_k = NULL,  # SVD rank (NULL = auto-detect via choose_k())  # <<CFG:alra_k>>
  alra_q = 10,  # Number of power iterations for randomized SVD  # <<CFG:alra_q>>
  alra_quantile_prob = 0.001,  # Quantile probability for thresholding  # <<CFG:alra_quantile_prob>>

  # ALRA normalization compatibility
  # ALRA only works with these normalization methods (NOT SCTransform):
  alra_compatible_methods = c("LogNormalize", "scran"),  # <<CFG:alra_compatible_methods>>

  # ===========================================================================
  # BACKWARD COMPATIBILITY - Deprecated parameters (DO NOT USE)
  # ===========================================================================
  # These are kept for backward compatibility but should not be used.
  # Use imputation_method instead.
  run_imputation = NULL,            # Deprecated: use imputation_method
  use_imputed_counts = NULL,        # Deprecated: use use_afmf_for_normalization

  # ===========================================================================
  # MODULE 03: NORMALIZATION PARAMETERS
  # ===========================================================================
  run_sctransform = TRUE,  # <<CFG:run_sctransform>>
  run_scran = TRUE,  # <<CFG:run_scran>>
  run_lognorm = TRUE,  # <<CFG:run_lognorm>>
  run_sckwarn = TRUE,  # <<CFG:run_sckwarn>>
  integration_normalization_method = "auto",  # <<CFG:integration_normalization_method>>
  sct_vars_to_regress = c("percent.mt"),  # <<CFG:sct_vars_to_regress>>
  sct_method = "glmGamPoi",  # <<CFG:sct_method>>
  scran_min_mean = 0.1,  # <<CFG:scran_min_mean>>
  run_normalization_benchmarking = TRUE,  # <<CFG:run_normalization_benchmarking>>
  norm_benchmark_max_cells = 5000,  # <<CFG:norm_benchmark_max_cells>>

  # ---------------------------------------------------------------------------
  # Normalization Benchmarking Integration Methods (NEW 2026-02-05)
  # ---------------------------------------------------------------------------
  # Which integration methods to use when benchmarking normalizations in Module 03.
  # Each normalization is tested with each integration method, then the best
  # normalization is selected by majority vote across integration methods.
  #
  # With only 1 method, majority vote cannot activate and the pipeline falls
  # back to a simple composite score. Use 2+ methods for proper voting.
  #
  # Available options (must also be installed):
  #   "harmony"  - Harmony (R package)
  #   "mnn"      - FastMNN via batchelor + SeuratWrappers (R packages)
  #   "scvi"     - scVI (Python, requires scvi-tools)
  #   "sccobra"  - scCobra (Python)
  #   "concord"  - Concord (Python)
  #
  # Default: c("harmony", "mnn") - enables majority vote with two R-based methods
  # ---------------------------------------------------------------------------
  norm_benchmark_integration_methods = c("harmony", "mnn", "rpca", "cca"),  # <<CFG:norm_benchmark_integration_methods>>

  # ===========================================================================
  # MODULE 04: INTEGRATION PARAMETERS
  # ===========================================================================
  run_batch_integration = TRUE,  # <<CFG:run_batch_integration>>
  batch_variable = "sample_name",  # Instead of "batch",  # <<CFG:batch_variable>>
  vars_to_regress = c("percent.mt"),  # <<CFG:vars_to_regress>>

  # ---------------------------------------------------------------------------
  # Integration Method Selection Mode
  # ---------------------------------------------------------------------------
  # This controls HOW the best integration method is selected from benchmarking.
  #
  # Options:
  #   "batch_removal"  - Minimizes batch_variance (most aggressive correction)
  #                      Criterion: which.min(batch_variance)
  #                      Use when: Batches are TECHNICAL replicates only
  #                      (same biology, different processing)
  #                      Example: Same cell line processed on different days
  #                      WARNING: May remove real biological signal!
  #
  #   "balanced"       - Maximizes composite score across ALL batch metrics
  #                      Criterion: which.max(mean(batch_var_norm, asw_norm, lisi_norm))
  #                      Use when: Batches have BIOLOGICAL meaning
  #                      (different conditions, sexes, timepoints, treatments)
  #                      Example: Male vs Female, Treatment vs Control
  #                      >>> RECOMMENDED FOR MOST BIOLOGICAL STUDIES <<<
  #
  #   "conservative"   - Prioritizes LISI score (moderate mixing)
  #                      Criterion: which.max(lisi_norm)
  #                      Use when: Preserving subtle biological differences is critical
  #                      Example: Rare cell types, subtle disease phenotypes
  #
  # For your Male vs Female choroid plexus comparison, use "balanced"!
  # ---------------------------------------------------------------------------
  integration_selection_mode = "balanced",  # <<CFG:integration_selection_mode>>
  celltype_column = "cluster_name_MapMyCells",  # <<CFG:celltype_column>>
  
  # ---------------------------------------------------------------------------
  # Biological vs Batch Weight for Composite Scoring
  # ---------------------------------------------------------------------------
  # Controls the relative importance of biological conservation vs batch
  # correction in the scIB-style composite score used for method selection.
  #
  # Default: 0.4 bio + 0.6 batch (scIB paper recommendation)
  #
  # Applied in both Module 03 (normalization) and Module 04 (integration).
  # ---------------------------------------------------------------------------
  bio_weight = 0.4,  # <<CFG:bio_weight>>
  batch_weight = 0.6,  # <<CFG:batch_weight>>

  # Integration method settings
  nfeatures_integration = 3000,  # <<CFG:nfeatures_integration>>
  dims_use = 30,  # <<CFG:dims_use>>
  integration_methods_r = c("harmony", "cca", "rpca"),  # <<CFG:integration_methods_r>>
  run_integration_benchmarking = TRUE,  # <<CFG:run_integration_benchmarking>>

  # ===========================================================================
  # MODULE 04: PYTHON INTEGRATION METHODS
  # ===========================================================================
  unified_python = "/scicore/home/doetsch/kaiser0001/miniforge3/envs/kaiser_test_py3.11/bin/python",  # <<CFG:unified_python>>
  run_python_integrations = TRUE,  # <<CFG:run_python_integrations>>
  integration_methods_python = c("scvi", "scanorama", "bbknn"),  # <<CFG:integration_methods_python>>

  # ===========================================================================
  # MODULE 04b: CONCORD PARAMETERS (Nature Biotechnology 2025)
  # ===========================================================================
  run_concord = TRUE,  # <<CFG:run_concord>>
  concord_n_top_features = 2000,  # <<CFG:concord_n_top_features>>
  concord_n_latent = 30,  # <<CFG:concord_n_latent>>
  concord_max_epochs = 200,  # <<CFG:concord_max_epochs>>
  concord_batch_size = 256,  # <<CFG:concord_batch_size>>
  concord_lr = 0.001,  # <<CFG:concord_lr>>
  concord_early_stopping_patience = 15,  # <<CFG:concord_early_stopping_patience>>
  concord_preload_dense = TRUE,  # <<CFG:concord_preload_dense>>
  concord_device = "auto",  # <<CFG:concord_device>>

  # ===========================================================================
  # MODULE 05: CHOIR CLUSTERING PARAMETERS
  # ===========================================================================
  run_choir_clustering = TRUE,  # <<CFG:run_choir_clustering>>
  choir_alpha = 0.05,  # <<CFG:choir_alpha>>
  choir_use_assay = "RNA",  # <<CFG:choir_use_assay>>

  # ===========================================================================
  # MODULE 05: scAURA CLUSTERING PARAMETERS
  # ===========================================================================
  # scAURA: Graph debiased contrastive learning for unsupervised cell type discovery
  # Reference: https://github.com/bozdaglab/scAURA
  #
  # scAURA runs alongside CHOIR in Module 05. Both methods can be enabled

  # independently. Use scice_clustering_source to select which feeds into
  # downstream analysis (Module 06 scICE subclustering).
  # ===========================================================================

  run_scaura_clustering = TRUE,  # <<CFG:run_scaura_clustering>>

  # Python environment for scAURA
  scaura_python = "/scicore/home/doetsch/kaiser0001/miniforge3/envs/kaiser_test_py3.11/bin/python",  # <<CFG:scaura_python>>

  # scAURA repository path
  scaura_repo_path = "/scicore/home/doetsch/kaiser0001/GITHUB_repositories/scAURA",  # <<CFG:scaura_repo_path>>

  # Number of clusters (K_CLUSTERS in scAURA)
  scaura_k_clusters = 7,  # <<CFG:scaura_k_clusters>>

  # Adaptive k-NN parameters
  scaura_kmax = 40,  # <<CFG:scaura_kmax>>

  # GCN architecture
  scaura_hidden_dim = 64,  # <<CFG:scaura_hidden_dim>>

  # Contrastive learning parameters
  scaura_tau = 0.7,  # <<CFG:scaura_tau>>
  scaura_tau_plus = 0.1,  # <<CFG:scaura_tau_plus>>
  scaura_pe = 0.3,  # <<CFG:scaura_pe>>
  scaura_pf = 0.3,  # <<CFG:scaura_pf>>

  # Training parameters
  scaura_epochs = 100,  # <<CFG:scaura_epochs>>
  scaura_lr = 0.001,  # <<CFG:scaura_lr>>

  # Self-supervised refinement (SSC)
  scaura_self_train = TRUE,  # <<CFG:scaura_self_train>>

  # Input: number of top HVGs to use
  scaura_n_top_genes = 2000,  # <<CFG:scaura_n_top_genes>>

  # Use GPU if available
  scaura_use_gpu = TRUE,  # <<CFG:scaura_use_gpu>>

  # ===========================================================================
  # CLUSTERING SOURCE SELECTION FOR DOWNSTREAM ANALYSIS
  # ===========================================================================
  # Controls which clustering method is used for scICE subclustering (Module 06)
  # and as the primary "seurat_clusters" identity.
  #
  # Options:
  #   "choir"  - Use CHOIR clusters (default if CHOIR succeeds)
  #   "scaura" - Use scAURA clusters
  #   "auto"   - Automatically select: CHOIR > scAURA > none
  #   "both"   - subcluster BOTH CHOIR and scAURA independently
  #              Creates suffixed columns: scice_subcluster_choir, scice_subcluster_scaura, etc.
  #              Also creates unsuffixed columns from primary source for backward compat
  # ===========================================================================
  scice_clustering_source = "both",  # <<CFG:scice_clustering_source>>

  # ===========================================================================
  # MODULE 05b: scCobra PARAMETERS
  # ===========================================================================
  run_sccobra = TRUE,  # <<CFG:run_sccobra>>

  # ===========================================================================
  # MODULE 06: scICE SUBCLUSTERING PARAMETERS
  # ===========================================================================
  run_scice_subclustering = TRUE,  # <<CFG:run_scice_subclustering>>
  scice_target_clusters = NULL,  # <<CFG:scice_target_clusters>>
  scice_k_min = 2,  # <<CFG:scice_k_min>>
  scice_k_max = 15,  # <<CFG:scice_k_max>>
  scice_ic_threshold = 1.005,  # <<CFG:scice_ic_threshold>>
  scice_min_cells = 100,  # <<CFG:scice_min_cells>>
  julia_bin = file.path(Sys.getenv("HOME"), "julia", "bin", "julia"),  # <<CFG:julia_bin>>
  scice_env = file.path(Sys.getenv("HOME"), "julia", "julia_envs", "scICE_env"),  # <<CFG:scice_env>>
  scice_pkg_dir = file.path(Sys.getenv("HOME"), "julia", "julia_envs", "scICE_env", "scICE"),  # <<CFG:scice_pkg_dir>>

  # ===========================================================================
  # MODULE 06: IDclust SUBCLUSTERING PARAMETERS
  # ===========================================================================
  # IDclust: Iterative Differential Clustering
  # Reference: Prompsy et al., NAR Genomics and Bioinformatics, 2024
  #
  # IDclust recursively splits clusters, validating each split with

  # differential expression testing. Only biologically meaningful splits
  # (with sufficient DEGs) are retained.
  #
  # Can run alongside or instead of scICE subclustering.
  # Install with: devtools::install_github("vallotlab/IDclust")
  # ===========================================================================

  run_idclust_subclustering = TRUE,  # <<CFG:run_idclust_subclustering>>
  idclust_target_clusters = NULL,  # NULL = all clusters meeting min_cells  # <<CFG:idclust_target_clusters>>
  idclust_logFC_th = log2(1.5),  # Log2 fold-change threshold for DE validation  # <<CFG:idclust_logFC_th>>
  idclust_qval_th = 0.01,  # Adjusted p-value threshold  # <<CFG:idclust_qval_th>>
  idclust_min_DEGs = 5,  # Minimum DEGs required to validate a split  # <<CFG:idclust_min_DEGs>>
  idclust_max_depth = 10,  # Maximum recursion depth  # <<CFG:idclust_max_depth>>
  idclust_min_frac_assigned = 0.1,  # Minimum fraction of cells assigned  # <<CFG:idclust_min_frac_assigned>>
  idclust_n_dims = 50,  # PCA dimensions for subclustering  # <<CFG:idclust_n_dims>>
  idclust_starting_resolution = 0.1,  # Initial Louvain resolution  # <<CFG:idclust_starting_resolution>>
  idclust_resolution = 0.8,  # Subsequent resolution  # <<CFG:idclust_resolution>>
  idclust_starting_k = 100,  # Initial k-nearest neighbors  # <<CFG:idclust_starting_k>>
  idclust_k = 100,  # Subsequent k neighbors  # <<CFG:idclust_k>>
  idclust_min_cells = 100,  # Minimum cells per cluster to subcluster  # <<CFG:idclust_min_cells>>
  idclust_plotting = TRUE,  # Generate IDclust internal plots  # <<CFG:idclust_plotting>>

  # ===========================================================================
  # UNIFIED SUBCLUSTERING CONTROL (Module 06)
  # ===========================================================================
  # Controls which subclustering methods run and their input source.
  #
  # subclustering_methods: character vector of methods to run
  #   c("scice", "idclust") - run both (default)
  #   c("scice")            - only scICE
  #   c("idclust")          - only IDclust
  #
  # subclustering_source: which primary clustering to subcluster
  #   "auto"   - auto-detect from Module 05 (default)
  #   "choir"  - use CHOIR clusters
  #   "scaura" - use scAURA clusters
  # ===========================================================================
  subclustering_methods = c("scice", "idclust"),  # <<CFG:subclustering_methods>>
  subclustering_source = "both",  # <<CFG:subclustering_source>>

  # ===========================================================================
  # MODULE 07: LEIDEN CLUSTERING PARAMETERS
  # ===========================================================================
  run_leiden_clustering = TRUE,  # <<CFG:run_leiden_clustering>>
  leiden_resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1),  # <<CFG:leiden_resolutions>>
  leiden_algorithm = 4,  # <<CFG:leiden_algorithm>>
  leiden_n_neighbors = 20,  # <<CFG:leiden_n_neighbors>>
  final_resolution = 0.5,  # <<CFG:final_resolution>>
  run_clustering_quality = TRUE,  # <<CFG:run_clustering_quality>>

  # ===========================================================================
  # MODULE 07b: CLTS RE-NORMALIZATION PARAMETERS
  # ===========================================================================
  # CLTS (Count based on Linearized Transcriptome Size) preserves biological
  # transcriptome size variation across cell types while correcting technical
  # sequencing depth differences between samples.
  #
  # Module 07b can run after Module 06 (scICE) or Module 07 (Leiden)
  #
  # Reference: Lu et al., Nature Communications 2025
  # DOI: 10.1038/s41467-025-56623-1
  # ===========================================================================

  run_clts_renormalization = TRUE,  # <<CFG:run_clts_renormalization>>

  # ---------------------------------------------------------------------------
  # CLTS CLUSTERING SOURCE - Which clustering to use for CLTS normalization
  # ---------------------------------------------------------------------------
  # Options:
  #   "scice"  - Apply CLTS using scICE subclusters (DEFAULT)
  #   "leiden" - Apply CLTS using Leiden clusters
  #   "choir"  - Apply CLTS using CHOIR clusters
  #   "both"   - Apply CLTS to BOTH scICE and Leiden (creates two _redeconv objects)
  #   "all"    - Apply CLTS to scICE, Leiden, AND CHOIR
  #
  # Output files generated:
  #   "scice"  -> scice_subclustered_object_redeconv.rds
  #   "leiden" -> leiden_clustered_object_redeconv.rds
  #   "choir"  -> choir_clustered_object_redeconv.rds
  #   "both"   -> scice + leiden _redeconv.rds files
  #   "all"    -> all three _redeconv.rds files
  # ---------------------------------------------------------------------------
  clts_clustering_source = "scice",  # <<CFG:clts_clustering_source>>

  # Which cluster column to use (auto-detected if NULL)
  clts_cluster_column = NULL,  # <<CFG:clts_cluster_column>>

  # Minimum cells per cluster to include in regression
  clts_min_cells_per_cluster = 50,  # <<CFG:clts_min_cells_per_cluster>>

  # Baseline sample selection: "auto" (highest correlation) or specific sample name
  clts_baseline_sample = "auto",  # <<CFG:clts_baseline_sample>>

  # Marker detection thresholds for benchmark comparison
  clts_marker_logfc_threshold = 0.5,  # <<CFG:clts_marker_logfc_threshold>>
  clts_marker_pval_threshold = 0.05,  # <<CFG:clts_marker_pval_threshold>>
  clts_marker_min_pct = 0.25,  # <<CFG:clts_marker_min_pct>>

  # Run marker benchmark comparing original vs CLTS normalization
  clts_run_benchmark = TRUE,  # <<CFG:clts_run_benchmark>>

  # ===========================================================================
  # MODULE 08: DIFFERENTIAL EXPRESSION PARAMETERS
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # DE OBJECT SOURCES - Which objects to run DE analysis on
  # ---------------------------------------------------------------------------
  # This is a CHARACTER VECTOR - you can specify multiple objects to run
  # DE analysis on in parallel. Results will be saved in separate subdirectories.
  #
  # Valid options (can combine multiple):
  #   "scice"          - Standard scICE subclustered object
  #   "scice_redeconv" - CLTS-normalized scICE object
  #   "leiden"         - Standard Leiden clustered object
  #   "leiden_redeconv"- CLTS-normalized Leiden object
  #   "choir"          - Standard CHOIR clustered object
  #   "choir_redeconv" - CLTS-normalized CHOIR object
  #
  # DEFAULT: Run DE on both standard scICE and CLTS-normalized scICE
  #
  # Examples:
  #   c("scice", "scice_redeconv")           - Compare standard vs CLTS on scICE
  #   c("leiden", "leiden_redeconv")         - Compare standard vs CLTS on Leiden
  #   c("scice", "scice_redeconv", "leiden") - Three parallel analyses
  #   c("scice_redeconv")                    - Only CLTS-normalized scICE
  # ---------------------------------------------------------------------------
  de_object_sources = c("scice", "scice_redeconv"),  # <<CFG:de_object_sources>>

  # Cross-object comparison (when multiple de_object_sources specified)
  de_run_cross_object_comparison = TRUE,  # <<CFG:de_run_cross_object_comparison>>

  # Comparison settings
  de_comparison_variable = "sex",  # <<CFG:de_comparison_variable>>
  de_group1 = "Male",  # <<CFG:de_group1>>
  de_group2 = "Female",  # <<CFG:de_group2>>
  de_comparison_scope = "whole_dataset",  # <<CFG:de_comparison_scope>>
  de_target_clusters = NULL,  # <<CFG:de_target_clusters>>
  de_covariates = c("percent.mt", "nFeature_RNA", "batch"),  # <<CFG:de_covariates>>

  # DE method flags
  run_mast = TRUE,  # <<CFG:run_mast>>
  run_dream = TRUE,  # <<CFG:run_dream>>
  run_pseudobulk_edger = TRUE,  # <<CFG:run_pseudobulk_edger>>
  run_pseudobulk_deseq2 = TRUE,  # <<CFG:run_pseudobulk_deseq2>>
  run_permutation = TRUE,  # <<CFG:run_permutation>>
  run_negative_control = TRUE,  # <<CFG:run_negative_control>>

  # DE thresholds
  de_logfc_threshold = 0.5,  # <<CFG:de_logfc_threshold>>
  de_pval_threshold = 0.05,  # <<CFG:de_pval_threshold>>
  de_min_pct = 0.1,  # <<CFG:de_min_pct>>

  # ===========================================================================
  # MODULE 09: VISUALIZATION AND GENES OF INTEREST
  # ===========================================================================
  sex_marker_genes = c("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty"),  # <<CFG:sex_marker_genes>>
  cp_marker_genes = c("Ttr", "Folr1", "Aqp1", "Cldn1", "Otx2", "Enpp2"),  # <<CFG:cp_marker_genes>>
  fibroblast_markers = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra", "Vim"),  # <<CFG:fibroblast_markers>>
  genes_of_interest = c("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty",
                        "Ttr", "Folr1", "Aqp1", "Cldn1", "Otx2",
                        "Col1a1", "Col3a1", "Dcn", "Pdgfra"),  # <<CFG:genes_of_interest>>

  # ===========================================================================
  # PLOT SETTINGS
  # ===========================================================================
  plot_formats = c("png", "pdf"),  # <<CFG:plot_formats>>
  plot_width_in = 7,  # <<CFG:plot_width_in>>
  plot_height_in = 5,  # <<CFG:plot_height_in>>
  plot_dpi = 300,  # <<CFG:plot_dpi>>
  tiff_compression = "lzw",  # <<CFG:tiff_compression>>

  # ===========================================================================
  # MISC SETTINGS
  # ===========================================================================
  max_cells_plot = 50000,  # <<CFG:max_cells_plot>>
  random_seed = 42  # <<CFG:random_seed>>
)


# ==============================================================================
# SECTION 6: DERIVE SAMPLES TO ANALYZE
# ==============================================================================
# Note: get_samples_to_analyze() function is now in functions.R

samples_to_analyze <- get_samples_to_analyze(params)
params$samples_to_analyze <- samples_to_analyze

params$analysis_metadata <- params$sample_metadata[
  params$sample_metadata$sample_name %in% samples_to_analyze,
]


# ==============================================================================
# SECTION 7: GENERATE INPUT FILE PATHS
# ==============================================================================
# Note: get_input_paths() function is now in functions.R

params$input_paths <- get_input_paths(params)


# ==============================================================================
# SECTION 8: PARAMETER VALIDATION
# ==============================================================================
# Note: validate_params() function is now in functions.R
# Additional validation for new parameters

# Validate CLTS and DE parameters
validate_clts_de_params <- function(params) {
  # Validate clts_clustering_source
  valid_clts_sources <- c("scice", "leiden", "choir", "both", "all")
  if (!params$clts_clustering_source %in% valid_clts_sources) {
    cat("  WARNING: Invalid clts_clustering_source '", params$clts_clustering_source,
        "'. Setting to 'scice'.\n", sep = "")
    params$clts_clustering_source <- "scice"
  }

  # Validate de_object_sources
  valid_de_sources <- c("scice", "scice_redeconv", "leiden", "leiden_redeconv",
                        "choir", "choir_redeconv")
  invalid_sources <- setdiff(params$de_object_sources, valid_de_sources)
  if (length(invalid_sources) > 0) {
    cat("  WARNING: Invalid de_object_sources: ", paste(invalid_sources, collapse = ", "),
        ". Removing.\n", sep = "")
    params$de_object_sources <- intersect(params$de_object_sources, valid_de_sources)
  }

  if (length(params$de_object_sources) == 0) {
    cat("  WARNING: No valid de_object_sources. Setting to 'scice'.\n")
    params$de_object_sources <- "scice"
  }

  # Warn about redeconv sources that may not be available
  clts_source <- params$clts_clustering_source
  de_sources <- params$de_object_sources

  if ("scice_redeconv" %in% de_sources && !clts_source %in% c("scice", "both", "all")) {
    cat("  NOTE: 'scice_redeconv' requested for DE but clts_clustering_source='",
        clts_source, "'.\n", sep = "")
    cat("        scice_subclustered_object_redeconv.rds may not be generated.\n")
  }
  if ("leiden_redeconv" %in% de_sources && !clts_source %in% c("leiden", "both", "all")) {
    cat("  NOTE: 'leiden_redeconv' requested for DE but clts_clustering_source='",
        clts_source, "'.\n", sep = "")
    cat("        leiden_clustered_object_redeconv.rds may not be generated.\n")
  }
  if ("choir_redeconv" %in% de_sources && !clts_source %in% c("choir", "all")) {
    cat("  NOTE: 'choir_redeconv' requested for DE but clts_clustering_source='",
        clts_source, "'.\n", sep = "")
    cat("        choir_clustered_object_redeconv.rds may not be generated.\n")
  }

  return(params)
}

# Validate imputation parameters
validate_imputation_params <- function(params) {
  # Validate imputation_method
  valid_methods <- c("none", "afmf", "alra", "both")
  if (is.null(params$imputation_method) || !params$imputation_method %in% valid_methods) {
    cat("  WARNING: Invalid imputation_method. Setting to 'both'.\n")
    params$imputation_method <- "both"
  }

  # Handle backward compatibility for deprecated parameters
  if (!is.null(params$run_imputation) && is.null(params$imputation_method)) {
    if (isTRUE(params$run_imputation)) {
      params$imputation_method <- "afmf"
      cat("  NOTE: Migrated deprecated run_imputation=TRUE to imputation_method='afmf'\n")
    }
  }

  # Check afMF Python path if afMF is enabled
  if (params$imputation_method %in% c("afmf", "both")) {
    if (!file.exists(params$afmf_python)) {
      cat("  WARNING: afMF Python not found at:", params$afmf_python, "\n")
      cat("           Run setup_afMF_env.sh to create the environment.\n")
      cat("           afMF imputation will be skipped if environment is not available.\n")
    }
  }

  # Validate afMF parameters
  if (is.null(params$afmf_max_iter) || params$afmf_max_iter < 1) {
    params$afmf_max_iter <- 100
  }
  if (is.null(params$afmf_tol) || params$afmf_tol <= 0) {
    params$afmf_tol <- 1e-5
  }
  if (is.null(params$afmf_min_cells_expressing) || params$afmf_min_cells_expressing < 1) {
    params$afmf_min_cells_expressing <- 10
  }

  # Validate ALRA parameters
  if (!is.null(params$alra_k) && (!is.numeric(params$alra_k) || params$alra_k < 1)) {
    cat("  WARNING: Invalid alra_k. Setting to NULL (auto-detect).\n")
    params$alra_k <- NULL
  }
  if (is.null(params$alra_q) || params$alra_q < 1) {
    params$alra_q <- 10
  }
  if (is.null(params$alra_quantile_prob) || params$alra_quantile_prob <= 0 || params$alra_quantile_prob >= 1) {
    params$alra_quantile_prob <- 0.001
  }

  # Ensure ALRA compatible methods is set
  if (is.null(params$alra_compatible_methods)) {
    params$alra_compatible_methods <- c("LogNormalize", "scran")
  }

  return(params)
}

params <- validate_params(params)
params <- validate_clts_de_params(params)
params <- validate_imputation_params(params)


# ==============================================================================
# SECTION 9: CONFIGURATION SUMMARY
# ==============================================================================

cat("================================================================================\n")
cat("PIPELINE CONFIGURATION SUMMARY\n")
cat("================================================================================\n\n")

cat("Project root:", params$project_root, "\n")
cat("Dataset name:", params$dataset_name, "\n")

# Show group or ventricle analysis scope
if (params$group_label != "") {
  cat("Analysis scope: Group", params$group_id, "(", params$group_label, ")\n")
} else if (params$ventricle_filter != "") {
  cat("Analysis scope: Ventricle", params$ventricle_filter, "\n")
} else {
  cat("Analysis scope: All samples\n")
}

cat("Sample sheet:", params$sample_sheet_path, "\n")

cat("\nSamples to analyze (", length(params$samples_to_analyze), "):\n", sep = "")
for (s in params$samples_to_analyze) {
  meta <- params$sample_metadata[params$sample_metadata$sample_name == s, ]
  vent <- if ("ventricle" %in% colnames(meta)) meta$ventricle else "-"
  grp <- if ("group_id" %in% colnames(meta) && !is.na(meta$group_id)) paste0("G", meta$group_id) else "-"
  cat("  ", s, " (", meta$sex, ", ", meta$batch, ", ", vent, ", ", grp, ")\n", sep = "")
}

cat("\nInput:", params$input_dir, "\n")
cat("Output:", params$out_root, "\n")

cat("\n--- Imputation Settings ---\n")
cat("  imputation_method:", params$imputation_method, "\n")
if (params$imputation_method %in% c("afmf", "both")) {
  cat("  afMF (Module 02b): ENABLED\n")
  cat("    use_afmf_for_normalization:", params$use_afmf_for_normalization, "\n")
} else {
  cat("  afMF (Module 02b): DISABLED\n")
}
if (params$imputation_method %in% c("alra", "both")) {
  cat("  ALRA (Module 03b): ENABLED (if scran or LogNormalize selected)\n")
  cat("    use_alra_for_downstream:", params$use_alra_for_downstream, "\n")
  cat("    NOTE: ALRA auto-skips if SCTransform wins benchmarking\n")
} else {
  cat("  ALRA (Module 03b): DISABLED\n")
}

cat("\n--- Module 03: Normalization ---\n")
cat("  run_sctransform:", params$run_sctransform, "\n")
cat("  run_scran:", params$run_scran, "\n")
cat("  run_lognorm:", params$run_lognorm, "\n")
cat("  run_sckwarn:", params$run_sckwarn, "\n")
cat("  norm_benchmark_max_cells:",
    if (is.null(params$norm_benchmark_max_cells)) "ALL" else params$norm_benchmark_max_cells, "\n")
cat("  norm_benchmark_integration_methods:",
    paste(params$norm_benchmark_integration_methods, collapse = ", "), "\n")

cat("\n--- Module 04: Integration ---\n")
cat("  Integration selection mode:", params$integration_selection_mode, "\n")
cat("  Batch variable:", params$batch_variable, "\n")
cat("  run_concord:", params$run_concord, "\n")

cat("\n--- Module 05-07: Clustering ---\n")
cat("  run_choir_clustering:", params$run_choir_clustering, "\n")
cat("  run_sccobra:", params$run_sccobra, "\n")
cat("  run_scice_subclustering:", params$run_scice_subclustering, "\n")
cat("  run_idclust_subclustering:", params$run_idclust_subclustering, "\n")
cat("  subclustering_methods:", paste(params$subclustering_methods, collapse = ", "), "\n")
cat("  subclustering_source:", params$subclustering_source, "\n")
cat("  run_leiden_clustering:", params$run_leiden_clustering, "\n")

cat("\n--- Module 07b: CLTS Re-normalization ---\n")
cat("  run_clts_renormalization:", params$run_clts_renormalization, "\n")
cat("  clts_clustering_source:", params$clts_clustering_source, "\n")
cat("  clts_min_cells_per_cluster:", params$clts_min_cells_per_cluster, "\n")
cat("  clts_baseline_sample:", params$clts_baseline_sample, "\n")
cat("  clts_run_benchmark:", params$clts_run_benchmark, "\n")

cat("\n--- Module 08: Differential Expression ---\n")
cat("  de_object_sources:", paste(params$de_object_sources, collapse = ", "), "\n")
cat("  Comparison:", params$de_group1, "vs", params$de_group2, "\n")
cat("  run_mast:", params$run_mast, "\n")
cat("  run_dream:", params$run_dream, "\n")
cat("  run_pseudobulk_edger:", params$run_pseudobulk_edger, "\n")
cat("  run_pseudobulk_deseq2:", params$run_pseudobulk_deseq2, "\n")
cat("  run_permutation:", params$run_permutation, "\n")
cat("  de_run_cross_object_comparison:", params$de_run_cross_object_comparison, "\n")

cat("\n================================================================================\n")
