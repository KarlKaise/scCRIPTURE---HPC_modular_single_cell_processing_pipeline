# ==============================================================================
# install_r_packages.R
# ==============================================================================
# Post-install script for R packages NOT available via conda-forge or bioconda.
#
# Run AFTER creating the conda environment:
#   conda activate kaiser_test_py3.11
#   Rscript install_r_packages.R
#
# These packages are installed from CRAN or GitHub into the conda R library.
# ==============================================================================

message("=============================================================")
message("  R Post-Install: Packages not available in conda channels   ")
message("=============================================================")
message(paste("R version:", R.version.string))
message(paste("Library:", .libPaths()[1]))
message("")

# --- Helper: install from CRAN if not already installed -----------------------
install_cran <- function(pkg, ...) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("[CRAN] Installing:", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", ...)
  } else {
    message(paste("[CRAN] Already installed:", pkg, packageVersion(pkg)))
  }
}

# --- Helper: install from GitHub if not already installed ---------------------
install_gh <- function(repo, pkg = NULL, ...) {
  if (is.null(pkg)) pkg <- basename(repo)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("[GitHub] Installing:", repo))
    devtools::install_github(repo, upgrade = "never", ...)
  } else {
    message(paste("[GitHub] Already installed:", pkg, packageVersion(pkg)))
  }
}

# --- Helper: install from Bioconductor if not already installed ---------------
install_bioc <- function(pkg, ...) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("[Bioc] Installing:", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE, ...)
  } else {
    message(paste("[Bioc] Already installed:", pkg, packageVersion(pkg)))
  }
}

# ==============================================================================
# 1. Verify build tools
# ==============================================================================
message("\n--- Checking build tools ---")
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
library(devtools)
library(BiocManager)

# ==============================================================================
# 2. CRAN packages not reliably in conda-forge
# ==============================================================================
message("\n--- CRAN packages ---")

install_cran("tictoc")            # Timing benchmarks
install_cran("scCustomize")       # Seurat QC/visualization wrappers
install_cran("countsplit")        # Count splitting for cluster validation
install_cran("qs")                # Fast serialization (used by some pipelines)

# ==============================================================================
# 3. Bioconductor packages (verify / fallback)
# ==============================================================================
message("\n--- Bioconductor packages (verify/fallback) ---")

# These should already be installed via bioconda, but verify:
install_bioc("DESeq2")              # Differential expression (bulk / pseudo-bulk)
install_bioc("DropletUtils")
install_bioc("SingleCellExperiment")
install_bioc("scDblFinder")
install_bioc("decontX")
install_bioc("MAST")
install_bioc("zellkonverter")
install_bioc("scRNAseq")
install_bioc("scran")
install_bioc("scater")
install_bioc("glmGamPoi")
install_bioc("ComplexHeatmap")
install_bioc("limma")
install_bioc("edgeR")

# ==============================================================================
# 4. GitHub-only packages
# ==============================================================================
message("\n--- GitHub packages ---")

# DropletQC - Nuclear fraction-based QC
install_gh("powellgenomicslab/DropletQC")

# DoubletFinder - Doublet detection
install_gh("chris-mcginnis-ucsf/DoubletFinder")

# scCDC - Contamination detection/correction
install_gh("lydiaMyr/scCDC")

# CHOIR - Clustering Hierarchy Optimization for Integrated Resolution
install_gh("corceslab/CHOIR")

# ALRA - Adaptively-thresholded Low Rank Approximation (imputation)
install_gh("KlugerLab/ALRA")

# IDclust - Iterative Differential Clustering
install_gh("vallotlab/IDclust")

# dreamlet - Differential expression for multi-sample scRNA-seq
install_gh("DiseaseNeurogenomics/dreamlet")

# ==============================================================================
# 5. Verify all installations
# ==============================================================================
message("\n===========================================================")
message("  Verification: Loading all required packages               ")
message("===========================================================\n")

required_packages <- c(
  # Seurat ecosystem
  "Seurat", "SeuratObject",
  # Core R
  "ggplot2", "dplyr", "tidyr", "Matrix", "patchwork", "reshape2",
  "future", "cluster", "reticulate", "ragg", "ggrepel", "ggnewscale",
  # Bioconductor
  "DESeq2",
  "DropletUtils", "SingleCellExperiment", "scDblFinder", "decontX",
  "MAST", "zellkonverter", "scRNAseq", "scran", "scater",
  "glmGamPoi", "ComplexHeatmap", "limma", "edgeR",
  # Specialized
  "harmony", "scCustomize", "tictoc", "countsplit",
  # GitHub
  "DropletQC", "DoubletFinder", "scCDC", "CHOIR", "ALRA", "IDclust", "dreamlet"
)

results <- data.frame(
  Package = required_packages,
  Installed = sapply(required_packages, function(p) {
    requireNamespace(p, quietly = TRUE)
  }),
  Version = sapply(required_packages, function(p) {
    tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
  }),
  stringsAsFactors = FALSE
)

# Print results
message(sprintf("%-25s %-10s %s", "Package", "Status", "Version"))
message(paste(rep("-", 50), collapse = ""))
for (i in seq_len(nrow(results))) {
  status <- if (results$Installed[i]) "OK" else "MISSING"
  ver <- if (is.na(results$Version[i])) "-" else results$Version[i]
  message(sprintf("%-25s %-10s %s", results$Package[i], status, ver))
}

n_ok <- sum(results$Installed)
n_total <- nrow(results)
message(paste0("\n", n_ok, "/", n_total, " packages installed successfully."))

if (any(!results$Installed)) {
  message("\nMISSING packages:")
  message(paste("  -", results$Package[!results$Installed], collapse = "\n"))
  message("\nTry installing missing packages manually.")
}

message("\n=== R post-install complete ===")
