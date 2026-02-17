# ==============================================================================
# MODULE 00: ENVIRONMENT SETUP (UNIFIED MULTI-SAMPLE PIPELINE)
# ==============================================================================
#
# This module sets up the computing environment for the downstream analysis
# pipeline that processes scCDC-corrected files from preprocessing.
#
# FEATURES:
# 1. Pins Python environment before any reticulate usage
# 2. Loads all required R packages
# 3. Checks Julia/scICE availability and instantiates packages
# 4. Checks afMF Python environment for imputation (Module 02b)
# 5. Checks ALRA R package for imputation (Module 03b)
# 6. Checks scKWARN R package for normalization (Module 03)
# 7. Creates output directory structure
# 8. Validates parameters and input files
# 9. Supports multi-sample analysis with configurable sample selection
#
# INPUT: Configuration from params.R
# OUTPUT: Pipeline environment saved for subsequent modules
#
# UPDATES:
# - 2026-01-09: Added ALRA package availability check for Module 03b
# - 2026-01-09: Added imputation_method parameter support
# - 2026-01-09: Added 03b_ALRA_Imputation output directory
# - 2026-01-09: Added afMF Python environment check for imputation (Module 02b)
# - 2026-01-09: Added 02b_Imputation output directory
# - 2026-01-09: Added 07b_CLTS_Normalization output directory
# - 2026-01-09: Added CLTS and DE object source parameters to summary
# - 2026-01-09: Added variancePartition check for DREAM DE
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("CHOROID PLEXUS scRNA-seq ANALYSIS PIPELINE\n")
cat("MODULE 00: ENVIRONMENT SETUP\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ==============================================================================
# Get script directory and source configuration
# ==============================================================================
script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  # If running interactively, use current directory
  "."
})

# Check for PIPELINE_DIR environment variable
pipeline_dir <- Sys.getenv("PIPELINE_DIR", unset = "")
if (pipeline_dir == "") {
  # Try to determine from script location
  if (file.exists(file.path(script_dir, "..", "config", "params.R"))) {
    pipeline_dir <- normalizePath(file.path(script_dir, ".."))
  } else if (file.exists("config/params.R")) {
    pipeline_dir <- normalizePath(".")
  } else {
    stop("Cannot determine PIPELINE_DIR. Set environment variable or run from pipeline directory.")
  }
  Sys.setenv(PIPELINE_DIR = pipeline_dir)
}

cat("PIPELINE_DIR:", pipeline_dir, "\n")

# Source parameters
params_file <- file.path(pipeline_dir, "config", "params.R")
if (!file.exists(params_file)) {
  stop("Configuration file not found: ", params_file)
}
source(params_file)

cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Samples to analyze:", length(params$samples_to_analyze), "\n")
cat("Batch variable:", params$batch_variable, "\n\n")

# ==============================================================================
# PIN PYTHON ENVIRONMENT BEFORE ANY RETICULATE USAGE
# ==============================================================================
cat("--- Pinning Python Environment ---\n")

Sys.setenv(RETICULATE_PYTHON = params$unified_python)
cat("RETICULATE_PYTHON set to:", params$unified_python, "\n\n")

# ==============================================================================
# LOAD RETICULATE AND VERIFY PYTHON
# ==============================================================================
suppressPackageStartupMessages(library(reticulate))

reticulate::use_python(params$unified_python, required = TRUE)

cat("reticulate::py_config():\n")
print(reticulate::py_config())
# FORCE Python initialization NOW - before any other package can override
invisible(reticulate::py_eval("1"))  # This locks in the Python version

# ==============================================================================
# VERIFY PYTHON PACKAGES
# ==============================================================================
cat("\n--- Verifying Python packages ---\n")

python_packages_ok <- TRUE

# Check torch
tryCatch({
  torch <- reticulate::import("torch")
  cuda_available <- torch$cuda$is_available()
  cat("  [OK] torch (CUDA available:", cuda_available, ")\n")
  if (cuda_available) {
    cat("       GPU:", torch$cuda$get_device_name(0L), "\n")
  }
}, error = function(e) {
  cat("  [FAILED] torch:", e$message, "\n")
  python_packages_ok <<- FALSE
})

# Check scvi
tryCatch({
  scvi <- reticulate::import("scvi")
  cat("  [OK] scvi version:", scvi$`__version__`, "\n")
}, error = function(e) {
  cat("  [FAILED] scvi:", e$message, "\n")
  python_packages_ok <<- FALSE
})

# Check scanorama
tryCatch({
  scanorama <- reticulate::import("scanorama")
  cat("  [OK] scanorama\n")
}, error = function(e) {
  cat("  [FAILED] scanorama:", e$message, "\n")
  python_packages_ok <<- FALSE
})

# Check bbknn
tryCatch({
  bbknn <- reticulate::import("bbknn")
  cat("  [OK] bbknn\n")
}, error = function(e) {
  cat("  [FAILED] bbknn:", e$message, "\n")
  python_packages_ok <<- FALSE
})

# Check scCobra - use consistent path from GITHUB_repositories
tryCatch({
  # First try direct import
  sccobra <- reticulate::import("scCobra")
  cat("  [OK] scCobra\n")
}, error = function(e) {
  # Try with path insertion - use GITHUB_repositories path consistently
  sccobra_path <- "/scicore/home/doetsch/kaiser0001/GITHUB_repositories/scCobra"
  if (dir.exists(sccobra_path)) {
    sys <- reticulate::import("sys")
    sys$path$insert(0L, sccobra_path)
    tryCatch({
      sccobra <- reticulate::import("scCobra")
      cat("  [OK] scCobra (from", sccobra_path, ")\n")
    }, error = function(e2) {
      cat("  [WARNING] scCobra not available from", sccobra_path, ":", e2$message, "\n")
    })
  } else {
    cat("  [WARNING] scCobra path not found:", sccobra_path, "\n")
  }
})

# Check CONCORD (Nature Biotechnology 2025)
tryCatch({
  concord <- reticulate::import("concord")
  cat("  [OK] concord\n")
}, error = function(e) {
  cat("  [WARNING] concord not available:", e$message, "\n")
  cat("       Install with: pip install concord-sc\n")
})

if (python_packages_ok) {
  cat("\n>>> ALL CORE PYTHON IMPORTS SUCCESSFUL <<<\n")
} else {
  cat("\n>>> WARNING: Some Python packages failed to import <<<\n")
}

# ==============================================================================
# VERIFY afMF PYTHON ENVIRONMENT (FOR MODULE 02b IMPUTATION)
# ==============================================================================
cat("\n--- Verifying afMF Python Environment (Module 02b) ---\n")

has_afmf <- FALSE
afmf_python <- NULL

# Get afMF Python path from params
if (!is.null(params$afmf_python) && params$afmf_python != "") {
  afmf_python <- params$afmf_python
} else {
  # Default path
  afmf_python <- file.path(Sys.getenv("HOME"), ".conda/envs/afMF_SCImputation_env/bin/python")
}

cat("  afMF Python path:", afmf_python, "\n")

if (file.exists(afmf_python)) {
  cat("  [OK] afMF Python binary found\n")

  # Test if afMF package is importable
  # We need to test in the afMF environment, not the current one
  afmf_test <- tryCatch({
    result <- system2(afmf_python,
                      args = c("-c", "\"from afMF.runafMF import afMF; print('AFMF_OK')\""),
                      stdout = TRUE, stderr = TRUE,
                      timeout = 30)
    any(grepl("AFMF_OK", result))
  }, error = function(e) {
    FALSE
  })

  if (afmf_test) {
    has_afmf <- TRUE
    cat("  [OK] afMF package importable\n")
  } else {
    cat("  [WARNING] afMF package not importable\n")
    cat("       To install, run:\n")
    cat("       conda activate afMF_SCImputation_env\n")
    cat("       cd ~/GITHUB_repositories/SCImputation/afMF\n")
    cat("       pip install .\n")
  }
} else {
  cat("  [NOT FOUND] afMF Python environment not found\n")
  cat("       Expected at:", afmf_python, "\n")
  cat("       To create the environment, run:\n")
  cat("       conda create -n afMF_SCImputation_env python=3.10 -y\n")
  cat("       conda activate afMF_SCImputation_env\n")
  cat("       pip install numpy scipy pandas scikit-learn\n")
  cat("       cd ~/GITHUB_repositories/SCImputation/afMF\n")
  cat("       pip install .\n")
}

# Check if afMF imputation is enabled but afMF not available
if (params$imputation_method %in% c("afmf", "both") && !has_afmf) {
  cat("\n  WARNING: imputation_method includes afMF but afMF is not available\n")
  cat("           Module 02b will be skipped unless afMF is properly installed\n")
}

if (has_afmf) {
  cat("\n>>> afMF ENVIRONMENT READY <<<\n")
} else {
  cat("\n>>> afMF NOT CONFIGURED (afMF imputation will be skipped) <<<\n")
}

# ==============================================================================
# VERIFY ALRA R PACKAGE (FOR MODULE 03b IMPUTATION)
# ==============================================================================
cat("\n--- Verifying ALRA R Package (Module 03b) ---\n")

has_alra <- FALSE

# Check if ALRA package is available
if (requireNamespace("ALRA", quietly = TRUE)) {
  has_alra <- TRUE
  cat("  [OK] ALRA package installed\n")
  
  # Try to load and check version
  tryCatch({
    library(ALRA, quietly = TRUE)
    cat("  [OK] ALRA package loaded successfully\n")
    
    # Check for required functions
    if (exists("alra") && exists("choose_k")) {
      cat("  [OK] ALRA core functions (alra, choose_k) available\n")
    } else {
      cat("  [WARNING] ALRA functions not found - package may be incomplete\n")
      has_alra <- FALSE
    }
  }, error = function(e) {
    cat("  [WARNING] Failed to load ALRA:", e$message, "\n")
    has_alra <- FALSE
  })
} else {
  cat("  [NOT FOUND] ALRA package not installed\n")
  cat("       To install ALRA, run:\n")
  cat("       # Option 1: From CRAN (if available)\n")
  cat("       install.packages('ALRA')\n")
  cat("       \n")
  cat("       # Option 2: From GitHub\n")
  cat("       devtools::install_github('KlugerLab/ALRA')\n")
}

# Check if ALRA imputation is enabled but ALRA not available
if (params$imputation_method %in% c("alra", "both") && !has_alra) {
  cat("\n  WARNING: imputation_method includes ALRA but ALRA package is not available\n")
  cat("           Module 03b will be skipped unless ALRA is properly installed\n")
}

# Additional note about SCTransform incompatibility
if (params$imputation_method %in% c("alra", "both")) {
  cat("\n  NOTE: ALRA only works with LogNormalize and scran normalization\n")
  cat("        If SCTransform wins benchmarking, Module 03b will be automatically skipped\n")
}

if (has_alra) {
  cat("\n>>> ALRA PACKAGE READY <<<\n")
} else {
  cat("\n>>> ALRA NOT CONFIGURED (ALRA imputation will be skipped) <<<\n")
}

# ==============================================================================
# VERIFY scKWARN R PACKAGE (FOR MODULE 03 NORMALIZATION)
# ==============================================================================
cat("\n--- Verifying scKWARN R Package (Module 03) ---\n")

has_sckwarn <- FALSE

# Check if scKWARN package is available
if (requireNamespace("scKWARN", quietly = TRUE)) {
  has_sckwarn <- TRUE
  cat("  [OK] scKWARN package installed\n")

  # Try to load and check for LocASN function

  tryCatch({
    library(scKWARN, quietly = TRUE)
    cat("  [OK] scKWARN package loaded successfully\n")

    # Check for required function
    if (exists("LocASN")) {
      cat("  [OK] LocASN function available\n")
    } else {
      cat("  [WARNING] LocASN function not found - package may be incomplete\n")
      has_sckwarn <- FALSE
    }
  }, error = function(e) {
    cat("  [WARNING] Failed to load scKWARN:", e$message, "\n")
    has_sckwarn <- FALSE
  })
} else {
  cat("  [NOT FOUND] scKWARN package not installed\n")
  cat("       To install scKWARN, run:\n")
  cat("       devtools::install_github('cyhsuTN/scKWARN')\n")
}

# Check if scKWARN is enabled but not available
if (isTRUE(params$run_sckwarn) && !has_sckwarn) {
  cat("\n  WARNING: run_sckwarn is TRUE but scKWARN package is not available\n")
  cat("           scKWARN normalization will be skipped in Module 03\n")
}

if (has_sckwarn) {
  cat("\n>>> scKWARN PACKAGE READY <<<\n")
} else {
  cat("\n>>> scKWARN NOT CONFIGURED (scKWARN normalization will be skipped) <<<\n")
}

# ==============================================================================
# VERIFY JULIA AND scICE ENVIRONMENT
# ==============================================================================
cat("\n--- Verifying Julia/scICE Environment ---\n")

# Get Julia paths from params or use defaults
julia_bin <- if (!is.null(params$julia_bin)) {
  params$julia_bin
} else {
  file.path(Sys.getenv("HOME"), "julia", "bin", "julia")
}

# scICE package directory (the scLENS package directory containing Project.toml)
scice_pkg_dir <- if (!is.null(params$scice_pkg_dir)) {
  params$scice_pkg_dir
} else {
  file.path(Sys.getenv("HOME"), "julia", "julia_envs", "scICE_env", "scICE")
}

# Legacy: scice_env points to parent (for backward compatibility)
scice_env <- dirname(scice_pkg_dir)

scice_source <- file.path(scice_pkg_dir, "src", "scICE.jl")
sclens_project <- file.path(scice_pkg_dir, "Project.toml")

# Check Julia binary
has_julia <- FALSE
julia_version <- NULL

if (file.exists(julia_bin)) {
  julia_version <- tryCatch({
    result <- system2(julia_bin, "--version", stdout = TRUE, stderr = TRUE)
    if (length(result) > 0) result[1] else NULL
  }, error = function(e) NULL)

  if (!is.null(julia_version)) {
    has_julia <- TRUE
    cat("  [OK] Julia found:", julia_version, "\n")
    cat("       Binary:", julia_bin, "\n")
  } else {
    cat("  [FAILED] Julia binary exists but cannot execute\n")
    cat("       Path:", julia_bin, "\n")
  }
} else {
  cat("  [NOT FOUND] Julia binary not found at:", julia_bin, "\n")
  cat("       To install Julia, run:\n")
  cat("       wget https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.5-linux-x86_64.tar.gz -P ~/julia\n")
  cat("       cd ~/julia && tar -xzf julia-1.11.5-linux-x86_64.tar.gz\n")
}

# Check scICE package directory (contains Project.toml)
has_scice_env <- FALSE

if (file.exists(sclens_project)) {
  cat("  [OK] scLENS Project.toml found:", sclens_project, "\n")
  has_scice_env <- TRUE
} else {
  cat("  [NOT FOUND] scLENS Project.toml not found at:", sclens_project, "\n")
  cat("       Expected scICE/scLENS package at:", scice_pkg_dir, "\n")
}

# Check scICE source file
has_scice_source <- FALSE

if (file.exists(scice_source)) {
  cat("  [OK] scICE.jl source found:", scice_source, "\n")
  has_scice_source <- TRUE
} else {
  cat("  [NOT FOUND] scICE.jl not found at:", scice_source, "\n")
  if (has_scice_env) {
    cat("       To install, run:\n")
    cat("       cd", scice_env, "\n")
    cat("       wget https://github.com/Mathbiomed/scICE/archive/refs/heads/main.zip -O scICE.zip\n")
    cat("       unzip scICE.zip && mv scICE-main scICE && rm scICE.zip\n")
  }
}

# ==============================================================================
# INSTANTIATE JULIA PACKAGES (Critical for SLURM jobs) - SIMPLIFIED VERSION
# ==============================================================================
# This checks if packages are installed by looking for Manifest.toml
# Only runs Pkg.instantiate() if Manifest is missing

has_scice_packages <- FALSE

if (has_julia && has_scice_env) {
  cat("\n--- Checking Julia packages ---\n")

  # Pin depot so login + SLURM jobs see the same Julia packages/artifacts
  julia_depot <- Sys.getenv("JULIA_DEPOT_PATH", unset = file.path(Sys.getenv("HOME"), ".julia"))
  Sys.setenv(JULIA_DEPOT_PATH = julia_depot)
  Sys.setenv(JULIA_PROJECT = scice_pkg_dir)

  cat("  JULIA_DEPOT_PATH:", julia_depot, "\n")
  cat("  JULIA_PROJECT:", scice_pkg_dir, "\n")

  # Check if Manifest.toml exists (indicates packages were resolved previously)
  manifest_file <- file.path(scice_pkg_dir, "Manifest.toml")

  if (file.exists(manifest_file)) {
    manifest_lines <- length(readLines(manifest_file))
    cat("  Manifest.toml:", manifest_lines, "lines\n")

    if (manifest_lines > 100) {
      # Packages are installed - do a quick startup test
      cat("  Testing Julia startup...")

      startup_test <- tryCatch({
        system2(julia_bin,
                args = c(paste0("--project=", scice_pkg_dir),
                         "-e", "println(\"STARTUP_OK\")"),
                stdout = TRUE, stderr = TRUE,
                timeout = 60)
      }, error = function(e) NULL)

      if (!is.null(startup_test) && any(grepl("STARTUP_OK", startup_test))) {
        has_scice_packages <- TRUE
        cat(" [OK]\n")
      } else {
        cat(" [WARNING] Startup test failed\n")
        cat("  Will attempt to continue anyway\n")
        has_scice_packages <- TRUE  # Assume OK since Manifest exists
      }
    } else {
      cat("  [WARNING] Manifest.toml seems incomplete (", manifest_lines, " lines)\n")
    }
  } else {
    cat("  [WARNING] Manifest.toml not found\n")
    cat("  Packages have never been resolved.\n")
    cat("  To fix, run on login node:\n")
    cat("    ", julia_bin, " --project=", scice_pkg_dir, "\n", sep="")
    cat("    julia> using Pkg; Pkg.instantiate(); Pkg.precompile()\n")
  }
} else {
  cat("  Skipping Julia check (Julia or scICE env not found)\n")
}

# Overall scICE readiness
has_scice <- has_julia && has_scice_env && has_scice_source && has_scice_packages


if (has_scice) {
  cat("\n>>> scICE ENVIRONMENT READY <<<\n")
} else {
  cat("\n>>> scICE NOT FULLY CONFIGURED <<<\n")
  if (isTRUE(params$run_scice_subclustering)) {
    cat("    WARNING: run_scice_subclustering is TRUE but scICE is not ready\n")
    cat("    Module 06 will be skipped unless scICE is properly installed\n")
  }
}

# ==============================================================================
# LOAD R PACKAGES
# ==============================================================================
cat("\n--- Loading R packages ---\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(Matrix)
  library(future)
  library(MAST)
  library(edgeR)
  library(limma)
  library(scran)
  library(scuttle)
  library(SingleCellExperiment)
  library(reshape2)
  library(viridis)
  library(RColorBrewer)
  library(pheatmap)
  library(ggrepel)
  library(scales)
  library(cluster)
})

# Check optional packages
has_SeuratIntegrate <- requireNamespace("SeuratIntegrate", quietly = TRUE)
has_CHOIR <- requireNamespace("CHOIR", quietly = TRUE)
has_glmGamPoi <- requireNamespace("glmGamPoi", quietly = TRUE)
has_DESeq2 <- requireNamespace("DESeq2", quietly = TRUE)
has_harmony <- requireNamespace("harmony", quietly = TRUE)
has_lisi <- requireNamespace("lisi", quietly = TRUE)
has_kBET <- requireNamespace("kBET", quietly = TRUE)
has_batchelor <- requireNamespace("batchelor", quietly = TRUE)
has_SeuratWrappers <- requireNamespace("SeuratWrappers", quietly = TRUE)
has_MAST <- requireNamespace("MAST", quietly = TRUE)
has_clustree <- requireNamespace("clustree", quietly = TRUE)
has_bluster <- requireNamespace("bluster", quietly = TRUE)
has_variancePartition <- requireNamespace("variancePartition", quietly = TRUE)

# Note: has_sckwarn already checked above with detailed verification

cat("\nR Package availability:\n")
cat("  SeuratIntegrate:", has_SeuratIntegrate, "\n")
cat("  CHOIR:", has_CHOIR, "\n")
cat("  glmGamPoi:", has_glmGamPoi, "\n")
cat("  DESeq2:", has_DESeq2, "\n")
cat("  harmony:", has_harmony, "\n")
cat("  lisi:", has_lisi, "\n")
cat("  kBET:", has_kBET, "\n")
cat("  batchelor (for MNN):", has_batchelor, "\n")
cat("  SeuratWrappers (for MNN):", has_SeuratWrappers, "\n")
cat("  MAST:", has_MAST, "\n")
cat("  clustree:", has_clustree, "\n")
cat("  bluster:", has_bluster, "\n")
cat("  variancePartition (for DREAM):", has_variancePartition, "\n")
cat("  ALRA:", has_alra, "\n")
cat("  scKWARN:", has_sckwarn, "\n")

# Load SeuratIntegrate if available
if (has_SeuratIntegrate) {
  suppressPackageStartupMessages(library(SeuratIntegrate))
  cat("\nSeuratIntegrate loaded.\n")
}

# Load CHOIR if available
if (has_CHOIR) {
  suppressPackageStartupMessages(library(CHOIR))
  cat("CHOIR loaded.\n")
}

# Load DESeq2 if available
if (has_DESeq2) {
  suppressPackageStartupMessages(library(DESeq2))
  cat("DESeq2 loaded.\n")
}

# Load harmony if available
if (has_harmony) {
  suppressPackageStartupMessages(library(harmony))
  cat("harmony loaded.\n")
}

# Load SeuratWrappers if available (for FastMNN)
if (has_SeuratWrappers) {
  suppressPackageStartupMessages(library(SeuratWrappers))
  cat("SeuratWrappers loaded.\n")
}

# Load variancePartition if available (for DREAM)
if (has_variancePartition) {
  suppressPackageStartupMessages(library(variancePartition))
  cat("variancePartition loaded (DREAM DE).\n")
}

# Load ALRA if available
if (has_alra) {
  suppressPackageStartupMessages(library(ALRA))
  cat("ALRA loaded.\n")
}

# ==============================================================================
# CREATE OUTPUT DIRECTORIES
# ==============================================================================
cat("\n--- Creating output directories ---\n")

out_base <- params$out_root
cat("Output root:", out_base, "\n")

# Main output directories
output_dirs <- list(
  root = out_base,
  plots = file.path(out_base, "plots"),
  tables = file.path(out_base, "tables"),
  objects = file.path(out_base, "objects"),
  qc = file.path(out_base, "01_QC"),
  normalization = file.path(out_base, "02_Normalization"),
  imputation_afmf = file.path(out_base, "02b_Imputation_afMF"),
  imputation_alra = file.path(out_base, "03b_Imputation_ALRA"),
  integration = file.path(out_base, "03_Integration"),
  benchmarking = file.path(out_base, "04_Benchmarking"),
  sccobra = file.path(out_base, "05_scCobra"),
  choir_clustering = file.path(out_base, "05_CHOIR_Clustering"),
  scice_subclustering = file.path(out_base, "05_CHOIR_Clustering", "scICE_subclustering"),
  leiden_clustering = file.path(out_base, "07_Leiden_Clustering"),
  clts_normalization = file.path(out_base, "07b_CLTS_Normalization"),
  de = file.path(out_base, "08_Differential_Expression"),
  gene_viz = file.path(out_base, "09_Gene_Visualization"),
  reports = file.path(out_base, "reports")
)

# Create all directories
for (dir_name in names(output_dirs)) {
  dir.create(output_dirs[[dir_name]], showWarnings = FALSE, recursive = TRUE)
  cat("  ", dir_name, ":", output_dirs[[dir_name]], "\n")
}

# Subdirectories
subdirs <- list(
  # Plot subdirectories
  plots_qc = file.path(output_dirs$plots, "qc"),
  plots_normalization = file.path(output_dirs$plots, "normalization"),
  plots_integration = file.path(output_dirs$plots, "integration"),
  plots_clustering = file.path(output_dirs$plots, "clustering"),
  plots_leiden = file.path(output_dirs$plots, "leiden_clustering"),
  plots_scice = file.path(output_dirs$plots, "scice_subclustering"),
  plots_clts = file.path(output_dirs$plots, "clts_normalization"),
  plots_de = file.path(output_dirs$plots, "de"),
  plots_gene_expression = file.path(output_dirs$plots, "gene_expression"),
  plots_sccobra = file.path(output_dirs$plots, "sccobra"),
  plots_imputation_afmf = file.path(output_dirs$plots, "imputation_afmf"),
  plots_imputation_alra = file.path(output_dirs$plots, "imputation_alra"),

  # Normalization subdirectories
  norm_sct = file.path(output_dirs$normalization, "SCTransform"),
  norm_scran = file.path(output_dirs$normalization, "scran"),
  norm_lognorm = file.path(output_dirs$normalization, "LogNormalize"),
  norm_benchmark = file.path(output_dirs$normalization, "benchmarking"),

  # afMF Imputation subdirectories
  afmf_tables = file.path(output_dirs$imputation_afmf, "tables"),
  afmf_plots = file.path(output_dirs$imputation_afmf, "plots"),
  afmf_benchmark = file.path(output_dirs$imputation_afmf, "benchmark"),

  # ALRA Imputation subdirectories
  alra_tables = file.path(output_dirs$imputation_alra, "tables"),
  alra_plots = file.path(output_dirs$imputation_alra, "plots"),
  alra_benchmark = file.path(output_dirs$imputation_alra, "benchmark"),

  # Clustering subdirectories
  resolution_testing = file.path(output_dirs$leiden_clustering, "resolution_testing"),
  quality_metrics = file.path(output_dirs$leiden_clustering, "quality_metrics"),

  # CLTS subdirectories (per clustering type)
  clts_scice = file.path(output_dirs$clts_normalization, "scice"),
  clts_leiden = file.path(output_dirs$clts_normalization, "leiden"),
  clts_choir = file.path(output_dirs$clts_normalization, "choir"),

  # DE subdirectories (legacy - now dynamic per object source)
  de_mast = file.path(output_dirs$de, "MAST"),
  de_edger = file.path(output_dirs$de, "pseudobulk_edgeR"),
  de_deseq2 = file.path(output_dirs$de, "pseudobulk_DESeq2"),

  # Tables subdirectories
  de_results = file.path(output_dirs$tables, "de_results")
)

# Create all subdirs
for (dir_name in names(subdirs)) {
  dir.create(subdirs[[dir_name]], showWarnings = FALSE, recursive = TRUE)
}

cat("\nSubdirectories created.\n")

# ==============================================================================
# LOAD UTILITY FUNCTIONS
# ==============================================================================
utils_file <- file.path(pipeline_dir, "utils", "functions.R")
if (file.exists(utils_file)) {
  source(utils_file)
  cat("Utility functions loaded from:", utils_file, "\n")
} else {
  cat("WARNING: Utility functions file not found:", utils_file, "\n")
}

# ==============================================================================
# VALIDATE INPUT FILES
# ==============================================================================
cat("\n--- Validating Input Files ---\n")

input_validation <- data.frame(
  sample = character(),
  path = character(),
  exists = logical(),
  size_mb = numeric(),
  stringsAsFactors = FALSE
)

for (sample in params$samples_to_analyze) {
  path <- params$input_paths[[sample]]
  exists <- file.exists(path)
  size_mb <- if (exists) round(file.info(path)$size / (1024^2), 2) else NA

  input_validation <- rbind(input_validation, data.frame(
    sample = sample,
    path = path,
    exists = exists,
    size_mb = size_mb,
    stringsAsFactors = FALSE
  ))

  status <- if (exists) paste0("[OK] ", size_mb, " MB") else "[MISSING]"
  cat("  ", sample, ":", status, "\n")
}

# Check for missing files
missing_samples <- input_validation$sample[!input_validation$exists]
if (length(missing_samples) > 0) {
  cat("\n  WARNING: Missing input files for samples:", paste(missing_samples, collapse = ", "), "\n")
  cat("  Pipeline will fail when trying to load these samples.\n")
  cat("  Either run preprocessing first or update samples_to_exclude in params.R\n")
}

# Save validation results
write.csv(input_validation, file.path(output_dirs$tables, "input_file_validation.csv"), row.names = FALSE)

# ==============================================================================
# SAVE ENVIRONMENT FOR SUBSEQUENT MODULES
# ==============================================================================
env_file <- file.path(output_dirs$objects, "pipeline_environment.RData")

save(
  params,
  output_dirs,
  subdirs,
  out_base,
  pipeline_dir,
  input_validation,
  # R package availability
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
  has_sckwarn,
  # Python status
  python_packages_ok,
  # afMF status (for Module 02b imputation)
  has_afmf,
  afmf_python,
  # Julia/scICE status
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
  file = env_file
)

cat("\nEnvironment saved to:", env_file, "\n")

# ==============================================================================
# WRITE SAMPLE METADATA TO FILE
# ==============================================================================
sample_meta_file <- file.path(output_dirs$tables, "sample_metadata.csv")
write.csv(params$analysis_metadata, sample_meta_file, row.names = FALSE)
cat("Sample metadata saved to:", sample_meta_file, "\n")

# ==============================================================================
# PRINT SUMMARY
# ==============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("ENVIRONMENT SETUP SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Analysis configuration:\n")
cat("  Samples to analyze:", length(params$samples_to_analyze), "\n")
for (s in params$samples_to_analyze) {
  meta <- params$sample_metadata[params$sample_metadata$sample_name == s, ]
  cat("    -", s, "(", meta$sex, ",", meta$batch, ")\n")
}
cat("  Batch variable:", params$batch_variable, "\n")
cat("  Output directory:", out_base, "\n")

cat("\nPython environment:", params$unified_python, "\n")
cat("Python packages OK:", python_packages_ok, "\n")

cat("\n--- Imputation Settings ---\n")
cat("  imputation_method:", params$imputation_method, "\n")

cat("\nafMF environment (Module 02b - counts-based):\n")
cat("  afMF Python path:", afmf_python, "\n")
cat("  afMF available:", has_afmf, "\n")
if (params$imputation_method %in% c("afmf", "both")) {
  cat("  afMF imputation: ENABLED\n")
  cat("  use_afmf_for_normalization:", isTRUE(params$use_afmf_for_normalization), "\n")
} else {
  cat("  afMF imputation: DISABLED\n")
}

cat("\nALRA package (Module 03b - normalized-data-based):\n")
cat("  ALRA available:", has_alra, "\n")
if (params$imputation_method %in% c("alra", "both")) {
  cat("  ALRA imputation: ENABLED (if scran or LogNormalize selected)\n")
  cat("  use_alra_for_downstream:", isTRUE(params$use_alra_for_downstream), "\n")
  cat("  NOTE: ALRA auto-skips if SCTransform wins benchmarking\n")
} else {
  cat("  ALRA imputation: DISABLED\n")
}

cat("\nJulia/scICE environment:\n")
cat("  Julia installed:", has_julia, "\n")
if (has_julia) cat("  Julia version:", julia_version, "\n")
cat("  scICE pkg dir:", scice_pkg_dir, "\n")
cat("  scICE env ready:", has_scice_env, "\n")
cat("  scICE source found:", has_scice_source, "\n")
cat("  scICE packages OK:", has_scice_packages, "\n")
cat("  scICE READY:", has_scice, "\n")

cat("\nKey analysis flags:\n")
cat("  run_sctransform:", params$run_sctransform, "\n")
cat("  run_scran:", params$run_scran, "\n")
cat("  run_lognorm:", params$run_lognorm, "\n")
cat("  run_sckwarn:", isTRUE(params$run_sckwarn), "(available:", has_sckwarn, ")\n")
cat("  run_batch_integration:", params$run_batch_integration, "\n")
cat("  run_integration_benchmarking:", params$run_integration_benchmarking, "\n")
cat("  run_sccobra:", params$run_sccobra, "\n")
cat("  run_concord:", params$run_concord, "\n")
cat("  run_choir_clustering:", params$run_choir_clustering, "\n")
cat("  run_scice_subclustering:", params$run_scice_subclustering, "\n")
cat("  run_leiden_clustering:", params$run_leiden_clustering, "\n")

cat("\nCLTS Re-normalization (Module 07b):\n")
cat("  run_clts_renormalization:", isTRUE(params$run_clts_renormalization), "\n")
if (isTRUE(params$run_clts_renormalization)) {
  cat("  clts_clustering_source:", params$clts_clustering_source, "\n")
  cat("  clts_min_cells_per_cluster:", params$clts_min_cells_per_cluster, "\n")
  cat("  clts_baseline_sample:", params$clts_baseline_sample, "\n")
  cat("  clts_run_benchmark:", params$clts_run_benchmark, "\n")
}

cat("\nDifferential Expression (Module 08):\n")
de_sources <- if (!is.null(params$de_object_sources)) params$de_object_sources else "scice"
cat("  de_object_sources:", paste(de_sources, collapse = ", "), "\n")
cat("  run_mast:", params$run_mast, "\n")
cat("  run_dream:", isTRUE(params$run_dream), "\n")
cat("  run_pseudobulk_edger:", params$run_pseudobulk_edger, "\n")
cat("  run_pseudobulk_deseq2:", params$run_pseudobulk_deseq2, "\n")
cat("  run_permutation:", isTRUE(params$run_permutation), "\n")
cat("  de_run_cross_object_comparison:", isTRUE(params$de_run_cross_object_comparison), "\n")

cat("\nFastMNN availability:\n")
cat("  batchelor:", has_batchelor, "\n")
cat("  SeuratWrappers:", has_SeuratWrappers, "\n")
cat("  FastMNN ready:", has_batchelor && has_SeuratWrappers, "\n")

cat("\nDREAM availability:\n")
cat("  variancePartition:", has_variancePartition, "\n")

cat("\n>>> MODULE 00 COMPLETE <<<\n")
