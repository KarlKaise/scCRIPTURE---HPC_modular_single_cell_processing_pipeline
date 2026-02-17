#!/usr/bin/env Rscript
# ==============================================================================
# MODULE 06: SUBCLUSTERING (scICE + IDclust)
# ==============================================================================
#
# This module performs subclustering using one or both of:
# - scICE: Julia-based scLENS for information-theoretic subclustering
# - IDclust: Iterative differential clustering with DE-validated splits
#
# Both methods operate per parent cluster (CHOIR/scAURA), subsetting cells
# and finding biologically meaningful subclusters within each.
#
# CONTROL VIA params.R:
#   subclustering_methods = c("scice", "idclust")  # which to run
#   subclustering_source  = "auto"                  # which primary clustering
#   run_scice_subclustering  = TRUE                 # individual toggle
#   run_idclust_subclustering = TRUE                # individual toggle
#
# SUBCLUSTERING SOURCE OPTIONS:
#   "auto"   - auto-detect from Module 05 primary method
#   "choir"  - subcluster CHOIR clusters only
#   "scaura" - subcluster scAURA clusters only
#   "both"   - subcluster BOTH CHOIR and scAURA clusters independently
#              Creates separate metadata columns per source:
#                scice_subcluster_choir, scice_subcluster_scaura
#                idclust_subcluster_choir, idclust_subcluster_scaura
#
# INPUT: Clustered Seurat object from Module 05 (CHOIR/scAURA)
# OUTPUT: Subclustered Seurat object with scice_subcluster and/or
#         idclust_subcluster metadata columns (with source suffix if "both")
#
# CRITICAL ENVIRONMENT NOTES (scICE):
# - Julia must run with sanitized LD_LIBRARY_PATH (no CUDA module, no conda libs)
# - LD_PRELOAD must be unset for Julia (conflicts with Julia's libstdc++)
# - UMAP.jl and Python umap-learn have namespace collision - handled by isolation
#
# UPDATES:
# - 2026-02-05: Added "both" subclustering_source for CHOIR+scAURA subclustering
# - 2026-02-05: Added IDclust subclustering alongside scICE
# - 2026-02-05: Flexible clustering source (auto/choir/scaura/both)
# - 2026-02-05: Unified subclustering_methods control parameter
# - 2026-01-15: CRITICAL FIX - Updated UMAP.jl API for v0.1.x
# - 2026-01-14: CRITICAL FIX - Proper Julia environment isolation
# - 2026-01-06: Added verbose output, enhanced Julia checks
# - 2026-01-04: Changed output directory to separate 06_scICE_subclustering
# - 2026-01-03: Added JoinLayers before count extraction (Seurat v5)
#
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MODULE 06: SUBCLUSTERING (scICE + IDclust)\n")
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
})

out_base <- params$out_root
load(file.path(out_base, "objects", "pipeline_environment.RData"))

# ==============================================================================
# Determine which subclustering methods to run
# ==============================================================================
cat("--- Subclustering method selection ---\n")

subclustering_methods <- if (!is.null(params$subclustering_methods)) {
  params$subclustering_methods
} else {
  methods <- c()
  if (isTRUE(params$run_scice_subclustering)) methods <- c(methods, "scice")
  if (isTRUE(params$run_idclust_subclustering)) methods <- c(methods, "idclust")
  methods
}

run_scice <- "scice" %in% subclustering_methods && isTRUE(params$run_scice_subclustering)
run_idclust <- "idclust" %in% subclustering_methods && isTRUE(params$run_idclust_subclustering)

has_IDclust <- FALSE
if (run_idclust) {
  has_IDclust <- requireNamespace("IDclust", quietly = TRUE)
  if (!has_IDclust) {
    cat("WARNING: IDclust requested but package not installed.\n")
    cat("  Install with: devtools::install_github('vallotlab/IDclust')\n")
    run_idclust <- FALSE
  }
}

subclustering_source <- if (!is.null(params$subclustering_source)) {
  params$subclustering_source
} else {
  "auto"
}

cat("  scICE:  ", run_scice, "\n")
cat("  IDclust:", run_idclust, if (run_idclust) paste0("(v", packageVersion("IDclust"), ")") else "", "\n")
cat("  Subclustering source:", subclustering_source, "\n\n")

# ==============================================================================
# Skip if no methods enabled
# ==============================================================================
if (!run_scice && !run_idclust) {
  cat("No subclustering methods enabled. Skipping Module 06.\n")
  scice_success <- FALSE
  scice_results <- NULL
  idclust_success <- FALSE
  idclust_results <- NULL
  save(scice_success, scice_results, idclust_success, idclust_results,
       file = file.path(output_dirs$objects, "06_scice_data.RData"))
  cat("\n>>> MODULE 06 SKIPPED <<<\n")
  quit(save = "no", status = 0)
}

# ==============================================================================
# Load clustered object (flexible source)
# ==============================================================================
cat("--- Loading clustered object ---\n")

# For "both" mode, we need the unified object that has both cluster columns
# For single mode, we can use source-specific objects

actual_source <- subclustering_source

if (actual_source == "auto") {
  clustering_data_file <- file.path(output_dirs$objects, "05_clustering_data.RData")
  if (file.exists(clustering_data_file)) {
    temp_env <- new.env()
    load(clustering_data_file, envir = temp_env)
    if (exists("clustering_method_used", envir = temp_env)) {
      actual_source <- get("clustering_method_used", envir = temp_env)
      cat("Auto-detected primary clustering method:", actual_source, "\n")
    }
  }
  if (actual_source == "auto" || actual_source == "none") {
    actual_source <- "choir"
    cat("Falling back to 'choir' as source\n")
  }
}

# For "both" mode, prefer the unified clustered object
if (subclustering_source == "both") {
  object_candidates <- c(
    file.path(output_dirs$objects, "clustered_object.rds"),
    file.path(output_dirs$objects, "choir_clustered_object.rds"),
    file.path(output_dirs$objects, "scaura_clustered_object.rds")
  )
} else {
  object_candidates <- c(
    file.path(output_dirs$objects, "clustered_object.rds"),
    file.path(output_dirs$objects, paste0(actual_source, "_clustered_object.rds")),
    file.path(output_dirs$objects, "choir_clustered_object.rds"),
    file.path(output_dirs$objects, "scaura_clustered_object.rds")
  )
}

choir_obj <- NULL
loaded_from <- NULL
for (candidate in object_candidates) {
  if (file.exists(candidate)) {
    choir_obj <- readRDS(candidate)
    loaded_from <- candidate
    cat("Loaded:", basename(candidate), "\n")
    break
  }
}

if (is.null(choir_obj)) {
  cat("ERROR: No clustered object found. Tried:\n")
  for (cand in object_candidates) cat("  ", cand, "\n")
  scice_success <- FALSE
  scice_results <- NULL
  idclust_success <- FALSE
  idclust_results <- NULL
  save(scice_success, scice_results, idclust_success, idclust_results,
       file = file.path(output_dirs$objects, "06_scice_data.RData"))
  quit(save = "no", status = 1)
}

cat("Object:", ncol(choir_obj), "cells,", nrow(choir_obj), "genes\n")
cat("Source:", loaded_from, "\n")
cat("Subclustering source mode:", subclustering_source, "\n\n")

# ==============================================================================
# Ensure layers are joined for count extraction (Seurat v5 compatibility)
# ==============================================================================
cat("--- Preparing object for count extraction ---\n")

tryCatch({
  if ("RNA" %in% names(choir_obj@assays)) {
    rna_layers <- Layers(choir_obj[["RNA"]])
    counts_layers <- grep("^counts", rna_layers, value = TRUE)

    if (length(counts_layers) > 1) {
      cat("Found", length(counts_layers), "counts layers - joining for export...\n")
      choir_obj[["RNA"]] <- JoinLayers(choir_obj[["RNA"]])
      cat("Successfully joined RNA layers\n")
    } else {
      cat("RNA layers already unified (", length(counts_layers), " counts layer)\n", sep = "")
    }
  }
}, error = function(e) {
  cat("Note: Layer preparation -", e$message, "\n")
})

# ==============================================================================
# Determine sources to process
# ==============================================================================
cat("\n--- Determining sources to process ---\n")

# Helper: find cluster column for a given source
get_cluster_col_for_source <- function(source_name, meta_cols) {
  if (source_name == "choir") {
    choir_pattern_cols <- grep("^CHOIR_clusters", meta_cols, value = TRUE)
    candidates <- c(
      choir_pattern_cols,
      "choir_clusters",
      "seurat_clusters"
    )
  } else if (source_name == "scaura") {
    candidates <- c("scaura_clusters", "scaura_ssc_clusters", "scaura_kmeans_clusters")
  } else {
    choir_pattern_cols <- grep("^CHOIR_clusters", meta_cols, value = TRUE)
    candidates <- c(
      choir_pattern_cols,
      "seurat_clusters",
      "choir_clusters",
      "scaura_clusters",
      "leiden_clusters"
    )
  }
  for (col in candidates) {
    if (col %in% meta_cols) return(col)
  }
  return(NULL)
}

meta_cols <- colnames(choir_obj@meta.data)

if (subclustering_source == "both") {
  # Check which sources actually have cluster columns in the object
  sources_to_process <- c()

  choir_col <- get_cluster_col_for_source("choir", meta_cols)
  if (!is.null(choir_col)) {
    sources_to_process <- c(sources_to_process, "choir")
    cat("  CHOIR clusters found: column '", choir_col, "'\n", sep = "")
  } else {
    cat("  WARNING: No CHOIR cluster column found in object - skipping CHOIR\n")
  }

  scaura_col <- get_cluster_col_for_source("scaura", meta_cols)
  if (!is.null(scaura_col)) {
    sources_to_process <- c(sources_to_process, "scaura")
    cat("  scAURA clusters found: column '", scaura_col, "'\n", sep = "")
  } else {
    cat("  WARNING: No scAURA cluster column found in object - skipping scAURA\n")
  }

  if (length(sources_to_process) == 0) {
    cat("ERROR: 'both' requested but no cluster columns found for CHOIR or scAURA\n")
    scice_success <- FALSE
    scice_results <- NULL
    idclust_success <- FALSE
    idclust_results <- NULL
    save(scice_success, scice_results, idclust_success, idclust_results,
         file = file.path(output_dirs$objects, "06_scice_data.RData"))
    quit(save = "no", status = 1)
  }

  is_multi_source <- length(sources_to_process) > 1
  cat("  Processing", length(sources_to_process), "source(s):", paste(sources_to_process, collapse = ", "), "\n")

} else {
  sources_to_process <- c(actual_source)
  is_multi_source <- FALSE
  cat("  Single source:", actual_source, "\n")
}

cat("\n")

# ==============================================================================
# Initialize global result accumulators
# ==============================================================================
all_scice_results <- list()
all_idclust_results <- list()
scice_success <- FALSE
idclust_success <- FALSE
scice_cols_added <- c()
idclust_cols_added <- c()


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                    scICE ONE-TIME SETUP (before source loop)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

scice_julia_ok <- FALSE
julia_depot <- NULL
n_cores <- NULL
scice_k_min <- NULL
scice_k_max <- NULL
scice_ic_threshold <- NULL
min_cells_subcluster <- NULL
julia_template <- NULL

if (run_scice) {

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("scICE ONE-TIME SETUP\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ==============================================================================
# CRITICAL: Julia Environment Setup Functions
# ==============================================================================

sanitize_ld_library_path_for_julia <- function(ld_path) {
  if (is.na(ld_path) || ld_path == "" || is.null(ld_path)) {
    return("")
  }
  parts <- strsplit(ld_path, ":", fixed = TRUE)[[1]]
  parts <- parts[!grepl("/scicore/soft/easybuild/apps/CUDA/", parts)]
  parts <- parts[!grepl("/miniforge3/", parts)]
  parts <- parts[!grepl("/conda/", parts)]
  parts <- parts[nzchar(parts)]
  parts <- parts[!duplicated(parts)]
  paste(parts, collapse = ":")
}

get_julia_env <- function(n_cores, julia_depot, julia_project) {
  ld_orig <- Sys.getenv("LD_LIBRARY_PATH", unset = "")
  ld_clean <- sanitize_ld_library_path_for_julia(ld_orig)
  env_vars <- c(
    paste0("JULIA_DEPOT_PATH=", julia_depot),
    paste0("JULIA_PROJECT=", julia_project),
    paste0("JULIA_NUM_THREADS=", n_cores),
    "JULIA_CUDA_MEMORY_POOL=none",
    "JULIA_PKG_OFFLINE=true",
    paste0("LD_LIBRARY_PATH=", ld_clean),
    "LD_PRELOAD="
  )
  return(env_vars)
}

run_julia_isolated <- function(julia_bin, args, julia_depot, julia_project,
                                n_cores = 4, timeout = 7200) {
  env_vars <- get_julia_env(n_cores, julia_depot, julia_project)
  cat("  Running Julia with isolated environment:\n")
  cat("    JULIA_PROJECT:", julia_project, "\n")
  cat("    JULIA_NUM_THREADS:", n_cores, "\n")
  cat("    LD_LIBRARY_PATH: [sanitized - no CUDA module, no conda]\n")
  cat("    LD_PRELOAD: [unset]\n\n")
  tryCatch({
    system2(
      julia_bin,
      args = args,
      stdout = TRUE,
      stderr = TRUE,
      timeout = timeout,
      env = env_vars
    )
  }, error = function(e) {
    cat("ERROR running Julia:", e$message, "\n")
    NULL
  })
}

# ==============================================================================
# Check Julia availability
# ==============================================================================
cat("--- Checking Julia availability ---\n")

julia_bin <- if (!is.null(params$julia_bin)) {
  params$julia_bin
} else {
  file.path(Sys.getenv("HOME"), "julia", "bin", "julia")
}

scice_pkg_dir <- if (!is.null(params$scice_pkg_dir)) {
  params$scice_pkg_dir
} else {
  file.path(Sys.getenv("HOME"), "julia", "julia_envs", "scICE_env", "scICE")
}

scice_source <- file.path(scice_pkg_dir, "src", "scICE.jl")
sclens_project <- file.path(scice_pkg_dir, "Project.toml")

scice_julia_ok <- TRUE

if (!file.exists(julia_bin)) {
  cat("ERROR: Julia not found at:", julia_bin, "\n")
  scice_julia_ok <- FALSE
}

if (scice_julia_ok && !file.exists(sclens_project)) {
  cat("ERROR: scLENS Project.toml not found at:", sclens_project, "\n")
  scice_julia_ok <- FALSE
}

if (scice_julia_ok && !file.exists(scice_source)) {
  cat("ERROR: scICE.jl not found at:", scice_source, "\n")
  scice_julia_ok <- FALSE
}

if (!scice_julia_ok) {
  cat("WARNING: Julia/scICE not available. scICE will be skipped.\n\n")
} else {

cat("Julia binary:", julia_bin, "\n")
cat("scICE/scLENS package:", scice_pkg_dir, "\n")

julia_version <- run_julia_isolated(
  julia_bin,
  args = c("--version"),
  julia_depot = Sys.getenv("JULIA_DEPOT_PATH",
                           unset = file.path(Sys.getenv("HOME"), ".julia")),
  julia_project = scice_pkg_dir,
  n_cores = 1,
  timeout = 30
)
if (!is.null(julia_version)) {
  cat("Julia version:", julia_version[1], "\n\n")
} else {
  cat("Julia version: [could not determine]\n\n")
}

# ==============================================================================
# Julia depot and project configuration
# ==============================================================================
cat("--- Configuring Julia environment ---\n")

julia_depot <- Sys.getenv("JULIA_DEPOT_PATH",
                          unset = file.path(Sys.getenv("HOME"), ".julia"))

cat("Julia project:", scice_pkg_dir, "\n")
cat("Julia depot  :", julia_depot, "\n")

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
cat("Julia threads:", n_cores, "\n\n")

# ==============================================================================
# Check package availability
# ==============================================================================
manifest_file <- file.path(scice_pkg_dir, "Manifest.toml")
if (file.exists(manifest_file)) {
  manifest_lines <- length(readLines(manifest_file))
  cat("Manifest.toml found:", manifest_lines, "lines\n")
  if (manifest_lines > 100) {
    cat("[OK] Packages appear to be resolved - skipping instantiation\n")
    cat("     (Run scice_readycheck.sbatch if you encounter issues)\n\n")
  } else {
    cat("[WARNING] Manifest.toml seems incomplete\n")
    cat("          Consider running: sbatch scice_readycheck.sbatch\n\n")
  }
} else {
  cat("WARNING: No Manifest.toml found!\n")
  cat("Run on a GPU node first:\n")
  cat("  sbatch Scripts/Slurm_scripts/scice_readycheck.sbatch\n\n")
}

cat("Quick Julia startup test...")
startup_test <- run_julia_isolated(
  julia_bin,
  args = c("--startup-file=no", paste0("--project=", scice_pkg_dir),
           "-e", "println(\"OK\")"),
  julia_depot = julia_depot,
  julia_project = scice_pkg_dir,
  n_cores = 1,
  timeout = 60
)
if (!is.null(startup_test) && any(grepl("OK", startup_test))) {
  cat(" [OK]\n\n")
} else {
  cat(" [WARNING] Julia startup test failed\n")
  cat("Output:", paste(startup_test, collapse = "\n"), "\n\n")
}

# ==============================================================================
# Get scICE parameters
# ==============================================================================
scice_k_min <- ifelse(is.null(params$scice_k_min), 2, params$scice_k_min)
scice_k_max <- ifelse(is.null(params$scice_k_max), 15, params$scice_k_max)
scice_ic_threshold <- ifelse(is.null(params$scice_ic_threshold), 1.005, params$scice_ic_threshold)
min_cells_subcluster <- ifelse(is.null(params$scice_min_cells), 100, params$scice_min_cells)

cat("scICE Parameters:\n")
cat("  K range:", scice_k_min, "-", scice_k_max, "\n")
cat("  IC threshold:", scice_ic_threshold, "\n")
cat("  Min cells:", min_cells_subcluster, "\n\n")

# ==============================================================================
# Julia script template (defined once, used per source/cluster)
# ==============================================================================
julia_template <- '
println("="^70)
println("scICE Analysis - Cluster {{CLUSTER_ID}}")
println("="^70)

ENV["NUM_CORES"] = "{{N_CORES}}"

scice_pkg_dir = "{{SCICE_PKG_DIR}}"
input_file = "{{INPUT_FILE}}"
output_dir = "{{OUTPUT_DIR}}"

mkpath(output_dir)

println("[1/10] Loading packages...")
println("  Activating: ", scice_pkg_dir)
using Pkg
Pkg.activate(scice_pkg_dir)

println("  Loading core packages...")
using CUDA, CSV, DataFrames

println("  Loading scLENS...")
using scLENS

println("  Loading scICE functions...")
include(joinpath(scice_pkg_dir, "src", "scICE.jl"))

println("  Loading CairoMakie...")
using CairoMakie
CairoMakie.activate!(type="png")
println("  Packages loaded successfully")

println("[2/10] Device selection...")
cur_dev = CUDA.functional() ? "gpu" : "cpu"
println("  Using: ", cur_dev)
if cur_dev == "gpu"
    CUDA.versioninfo()
end

println("[3/10] Loading data...")
println("  File: ", input_file)
ndf = scLENS.read_file(input_file)
println("  Loaded: ", size(ndf, 1), " cells x ", size(ndf, 2)-1, " genes")
println("  First column: ", names(ndf)[1])

cell_names_original = copy(ndf.cell)
println("  Saved ", length(cell_names_original), " cell names")

println("[4/10] Preprocessing...")
pre_df = scLENS.preprocess(ndf, mito_percent=20., min_genes_per_cell=100)
println("  After QC: ", size(pre_df, 1), " cells x ", size(pre_df, 2), " features")

println("[5/10] Creating scLENS embedding...")
println("  Device: ", cur_dev)
println("  This may take 10-30 minutes...")
t_start = time()
sclens_embedding = scLENS.sclens(pre_df, device_=cur_dev)
t_elapsed = round((time() - t_start) / 60, digits=1)
println("  Embedding completed in ", t_elapsed, " minutes")
println("  Keys: ", keys(sclens_embedding))

n_cells_after = size(pre_df, 1)
if n_cells_after == length(cell_names_original)
    sclens_embedding[:cell_id] = cell_names_original
    println("  Stored original cell names")
else
    println("  WARNING: Cell count changed (", length(cell_names_original), " -> ", n_cells_after, ")")
    sclens_embedding[:cell_id] = ["cell_$i" for i in 1:n_cells_after]
end

CSV.write(joinpath(output_dir, "pca.csv"), sclens_embedding[:pca_n1])
println("  PCA saved")

println("[6/10] UMAP transformation...")
scLENS.apply_umap!(sclens_embedding)
CSV.write(joinpath(output_dir, "umap.csv"), DataFrame(sclens_embedding[:umap], :auto))
println("  UMAP saved")

println("[7/10] UMAP visualization...")
try
    panel_0 = scLENS.plot_embedding(sclens_embedding)
    save(joinpath(output_dir, "umap_dist.png"), panel_0)
    println("  Plot saved")
catch e
    println("  Plot failed: ", e)
end

println("[8/10] scICE clustering...")
println("  Range: [{{K_MIN}}, {{K_MAX}}]")
println("  This may take 20-60 minutes...")
t_start = time()
clustering_success = false

try
    clustering!(sclens_embedding, [{{K_MIN}}, {{K_MAX}}])
    global clustering_success = true
    t_elapsed = round((time() - t_start) / 60, digits=1)
    println("  Clustering completed in ", t_elapsed, " minutes")
catch e
    println("  Range [{{K_MIN}},{{K_MAX}}] failed: ", e)
    println("  Trying default range [1,20]...")
    try
        clustering!(sclens_embedding)
        global clustering_success = true
        println("  Default range succeeded")
    catch e2
        println("  All clustering attempts failed: ", e2)
    end
end

if haskey(sclens_embedding, :n_cluster)
    println("  Candidate clusters: ", sclens_embedding[:n_cluster])
else
    println("  WARNING: No :n_cluster key - no valid clusters found")
end

println("[9/10] IC plot...")
if clustering_success && haskey(sclens_embedding, :n_cluster) && length(sclens_embedding[:n_cluster]) > 0
    try
        panel_1 = plot_ic(sclens_embedding)
        save(joinpath(output_dir, "ic_plot.png"), panel_1)
        println("  IC plot saved")
    catch e
        println("  IC plot failed: ", e)
    end
else
    println("  Skipped - no valid clustering results")
end

println("[10/10] Extracting consistent labels...")
label_out = nothing

if clustering_success && haskey(sclens_embedding, :n_cluster) && length(sclens_embedding[:n_cluster]) > 0
    if haskey(sclens_embedding, :pca) && !("cell" in names(sclens_embedding[:pca]))
        if haskey(sclens_embedding, :cell_id)
            pca_df = sclens_embedding[:pca]
            if size(pca_df, 1) == length(sclens_embedding[:cell_id])
                insertcols!(pca_df, 1, :cell => sclens_embedding[:cell_id])
                sclens_embedding[:pca] = pca_df
                println("  Added :cell column to :pca DataFrame")
            end
        end
    end

    for th in [{{IC_THRESHOLD}}, 1.01, 1.02, 1.05, 1.1]
        println("  Trying IC threshold: ", th)
        try
            global label_out = get_rlabel!(sclens_embedding, th)
            if !isnothing(label_out) && size(label_out, 2) > 1
                n_solutions = size(label_out, 2) - 1
                CSV.write(joinpath(output_dir, "consistent_labels.csv"), label_out)
                println("  SUCCESS: Found ", n_solutions, " solutions at threshold ", th)
                break
            else
                println("    No valid labels at this threshold")
            end
        catch e
            println("    Failed: ", e)
        end
    end
else
    println("  Skipped - no valid clustering results")
end

if isnothing(label_out)
    println()
    println("  Using k-means fallback...")
    using Clustering

    umap_coords = sclens_embedding[:umap]
    label_dict = Dict{String, Vector{Int}}()

    for k in [3, 5, 7]
        try
            km = kmeans(umap_coords\', k; maxiter=100)
            label_dict["l_$k"] = km.assignments
            println("    k=", k, " complete")
        catch e
            println("    k=", k, " failed")
        end
    end

    if !isempty(label_dict)
        cell_ids = haskey(sclens_embedding, :cell_id) ? sclens_embedding[:cell_id] : ["cell_$i" for i in 1:size(umap_coords, 1)]
        global label_out = DataFrame(:cell => cell_ids)
        for (k, v) in label_dict
            label_out[!, k] = v
        end
        CSV.write(joinpath(output_dir, "kmeans_labels.csv"), label_out)
        println("  K-means labels saved")
    end
end

if !isnothing(label_out)
    println()
    println("Creating cluster visualizations...")
    cluster_cols = filter(n -> startswith(string(n), "l_"), names(label_out))
    for col in cluster_cols[1:min(2, length(cluster_cols))]
        try
            panel = scLENS.plot_embedding(sclens_embedding, label_out[!, col])
            save(joinpath(output_dir, "umap_$(col).png"), panel)
            println("  Saved: umap_$(col).png")
        catch e
            println("  Failed to visualize ", col)
        end
    end
end

println()
println("="^70)
println("CLUSTER {{CLUSTER_ID}} COMPLETE")
println("="^70)
println()
println("Output directory: ", output_dir)
println("Files created:")
for f in readdir(output_dir)
    fpath = joinpath(output_dir, f)
    size_mb = round(filesize(fpath) / 1024^2, digits=2)
    println("  ", f, " (", size_mb, " MB)")
end
'

}  # end scice_julia_ok check (inner)
}  # end run_scice (one-time setup)


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                    MAIN SOURCE LOOP
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

for (current_source in sources_to_process) {

cat("\n", paste(rep("#", 80), collapse = ""), "\n")
cat("PROCESSING SOURCE:", toupper(current_source), "\n")
cat(paste(rep("#", 80), collapse = ""), "\n\n")

# ==============================================================================
# Determine column suffix and output directories for this source
# ==============================================================================
col_suffix <- if (is_multi_source) paste0("_", current_source) else ""
source_label <- toupper(current_source)

# Output directories
source_scice_base <- if (is_multi_source) {
  file.path(out_base, "06_scICE_subclustering", current_source)
} else {
  file.path(out_base, "06_scICE_subclustering")
}

source_idclust_base <- if (is_multi_source) {
  file.path(out_base, "06_IDclust_subclustering", current_source)
} else {
  file.path(out_base, "06_IDclust_subclustering")
}

if (run_scice) {
  dir.create(source_scice_base, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(source_scice_base, "input"), showWarnings = FALSE)
  dir.create(file.path(source_scice_base, "output"), showWarnings = FALSE)
  dir.create(file.path(source_scice_base, "scripts"), showWarnings = FALSE)
}
if (run_idclust) {
  dir.create(source_idclust_base, showWarnings = FALSE, recursive = TRUE)
}

cat("Source:", current_source, "\n")
cat("Column suffix:", if (col_suffix == "") "(none)" else col_suffix, "\n")
if (run_scice) cat("scICE output:", source_scice_base, "\n")
if (run_idclust) cat("IDclust output:", source_idclust_base, "\n")
cat("\n")

# ==============================================================================
# Identify clusters for this source
# ==============================================================================
cat("--- Identifying clusters for", source_label, "---\n")

cluster_col <- get_cluster_col_for_source(current_source, colnames(choir_obj@meta.data))

if (is.null(cluster_col)) {
  cat("WARNING: No cluster column found for source '", current_source, "'\n", sep = "")
  cat("  Available columns:", paste(head(colnames(choir_obj@meta.data), 20), collapse = ", "), "\n")
  cat("  Skipping this source.\n\n")
  next
}

cat("Using cluster column:", cluster_col, "\n")

clusters <- as.character(choir_obj@meta.data[[cluster_col]])
cluster_sizes <- table(clusters)
cat("Cluster distribution:\n")
print(cluster_sizes)
cat("\n")

# ==============================================================================
# Initialize per-source result variables
# ==============================================================================
source_scice_results <- list()
source_scice_success <- FALSE
source_idclust_results <- list()
source_idclust_success <- FALSE


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#               scICE SUBCLUSTERING (per source)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if (run_scice && scice_julia_ok) {

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("scICE SUBCLUSTERING [", source_label, "]\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Determine target clusters for this source
scice_target_clusters <- params$scice_target_clusters

if (is.null(scice_target_clusters)) {
  scice_target_clusters <- names(cluster_sizes)[cluster_sizes >= min_cells_subcluster]
  cat("Auto-selected clusters with >=", min_cells_subcluster, "cells:",
      paste(scice_target_clusters, collapse = ", "), "\n")
}

if (length(scice_target_clusters) == 0) {
  cat("WARNING: No clusters meet minimum cell requirement for scICE\n")
} else {

# ==============================================================================
# Export RAW COUNT MATRICES per cluster
# ==============================================================================
cat("\n--- Exporting raw count matrices ---\n")
cat("NOTE: scLENS requires CSV with rows=cells, cols=genes, first col='cell'\n\n")

counts_matrix <- GetAssayData(choir_obj, assay = "RNA", layer = "counts")

if (is.null(counts_matrix) || ncol(counts_matrix) == 0) {
  cat("ERROR: Could not extract counts matrix\n")
} else {

cat("Full counts matrix:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "cells\n")
cat("Matrix class:", class(counts_matrix)[1], "\n")
cat("Matrix sum:", format(sum(counts_matrix), big.mark = ","), "\n\n")

for (clust in scice_target_clusters) {
  cat("Exporting cluster", clust, "...")
  cells_in_cluster <- colnames(choir_obj)[clusters == clust]
  n_cells <- length(cells_in_cluster)
  cluster_counts <- as.matrix(counts_matrix[, cells_in_cluster])
  cluster_counts_t <- t(cluster_counts)
  count_df <- data.frame(
    cell = rownames(cluster_counts_t),
    cluster_counts_t,
    check.names = FALSE
  )
  output_file <- file.path(source_scice_base, "input", paste0("cluster_", clust, "_counts.csv"))
  write.csv(count_df, output_file, row.names = FALSE)
  file_size <- file.size(output_file) / (1024^2)
  cat(" saved", nrow(count_df), "cells x", ncol(count_df)-1, "genes",
      "(", round(file_size, 1), "MB)\n", sep = " ")
}

# ==============================================================================
# Run scICE on each cluster
# ==============================================================================
cat("\n--- Running scICE subclustering [", source_label, "] ---\n\n", sep = "")

for (clust in scice_target_clusters) {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("scICE [", source_label, "]: Processing cluster:", clust, "\n", sep = "")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  input_file <- file.path(source_scice_base, "input", paste0("cluster_", clust, "_counts.csv"))
  output_subdir <- file.path(source_scice_base, "output", paste0("cluster_", clust))
  dir.create(output_subdir, showWarnings = FALSE, recursive = TRUE)

  julia_script <- file.path(source_scice_base, "scripts", paste0("run_scice_cluster_", clust, ".jl"))

  julia_code <- julia_template
  julia_code <- gsub("\\{\\{CLUSTER_ID\\}\\}", clust, julia_code)
  julia_code <- gsub("\\{\\{N_CORES\\}\\}", n_cores, julia_code)
  julia_code <- gsub("\\{\\{SCICE_PKG_DIR\\}\\}", scice_pkg_dir, julia_code)
  julia_code <- gsub("\\{\\{INPUT_FILE\\}\\}", input_file, julia_code)
  julia_code <- gsub("\\{\\{OUTPUT_DIR\\}\\}", output_subdir, julia_code)
  julia_code <- gsub("\\{\\{K_MIN\\}\\}", scice_k_min, julia_code)
  julia_code <- gsub("\\{\\{K_MAX\\}\\}", scice_k_max, julia_code)
  julia_code <- gsub("\\{\\{IC_THRESHOLD\\}\\}", format(scice_ic_threshold, nsmall = 4), julia_code)

  writeLines(julia_code, julia_script)
  cat("Created Julia script:", julia_script, "\n")
  cat("Running scICE (expect 30-60 min per cluster)...\n\n")

  julia_result <- run_julia_isolated(
    julia_bin,
    args = c("--startup-file=no", paste0("--project=", scice_pkg_dir), julia_script),
    julia_depot = julia_depot,
    julia_project = scice_pkg_dir,
    n_cores = n_cores,
    timeout = 7200
  )

  if (!is.null(julia_result)) {
    cat(paste(julia_result, collapse = "\n"), "\n")
  }

  labels_file <- file.path(output_subdir, "consistent_labels.csv")
  kmeans_file <- file.path(output_subdir, "kmeans_labels.csv")

  if (file.exists(labels_file)) {
    labels_df <- read.csv(labels_file)
    source_scice_results[[clust]] <- list(
      cluster = clust, labels = labels_df, success = TRUE,
      method = "scICE", n_cells = nrow(labels_df)
    )
    cat("\n>>> scICE SUCCESS: cluster", clust, "(", nrow(labels_df), "cells)\n")
  } else if (file.exists(kmeans_file)) {
    labels_df <- read.csv(kmeans_file)
    source_scice_results[[clust]] <- list(
      cluster = clust, labels = labels_df, success = TRUE,
      method = "kmeans", n_cells = nrow(labels_df)
    )
    cat("\n>>> scICE SUCCESS: K-means fallback for cluster", clust, "(", nrow(labels_df), "cells)\n")
  } else {
    source_scice_results[[clust]] <- list(cluster = clust, success = FALSE)
    cat("\n>>> scICE WARNING: No labels for cluster", clust, "\n")
  }

  cat("\n")
}

# ==============================================================================
# Integrate scICE results for this source
# ==============================================================================
cat("\n--- Integrating scICE results [", source_label, "] ---\n", sep = "")

scice_col_name <- paste0("scice_subcluster", col_suffix)
choir_obj[[scice_col_name]] <- as.character(choir_obj@meta.data[[cluster_col]])
successful_clusters <- names(source_scice_results)[sapply(source_scice_results, function(x) x$success)]

cat("Column:", scice_col_name, "\n")
cat("scICE successful:", length(successful_clusters), "/", length(scice_target_clusters), "\n")

for (clust in successful_clusters) {
  result <- source_scice_results[[clust]]
  labels_df <- result$labels

  label_cols <- grep("^l_", colnames(labels_df), value = TRUE)
  if (length(label_cols) > 0) {
    best_col <- label_cols[which.max(sapply(label_cols, function(c) sum(!is.na(labels_df[[c]]))))]

    cell_col <- intersect(c("cell", "cell_id"), colnames(labels_df))
    if (length(cell_col) == 0) cell_col <- colnames(labels_df)[1]
    cell_ids <- labels_df[[cell_col[1]]]

    subcluster_labels <- labels_df[[best_col]]
    new_labels <- paste0(clust, ".", subcluster_labels)

    cells_in_cluster <- colnames(choir_obj)[clusters == clust]
    matched_idx <- match(cell_ids, cells_in_cluster)
    valid <- !is.na(matched_idx)

    if (sum(valid) > 0) {
      update_cells <- cells_in_cluster[matched_idx[valid]]
      choir_obj[[scice_col_name]][match(update_cells, colnames(choir_obj)), 1] <- new_labels[valid]
      cat("  Cluster", clust, ":", sum(valid), "cells assigned (", result$method, ")\n")
    }
  }
}

source_scice_success <- length(successful_clusters) > 0
if (source_scice_success) {
  scice_success <- TRUE
  scice_cols_added <- c(scice_cols_added, scice_col_name)
}
all_scice_results[[current_source]] <- source_scice_results

}  # end counts_matrix check
}  # end scice_target_clusters check
}  # end run_scice && scice_julia_ok


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#               IDclust SUBCLUSTERING (per source)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if (run_idclust) {

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("IDclust SUBCLUSTERING [", source_label, "]\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

library(IDclust)

# ==============================================================================
# IDclust parameters
# ==============================================================================
idclust_logFC_th <- if (!is.null(params$idclust_logFC_th)) params$idclust_logFC_th else log2(1.5)
idclust_qval_th <- if (!is.null(params$idclust_qval_th)) params$idclust_qval_th else 0.01
idclust_min_DEGs <- if (!is.null(params$idclust_min_DEGs)) params$idclust_min_DEGs else 5
idclust_max_depth <- if (!is.null(params$idclust_max_depth)) params$idclust_max_depth else 10
idclust_min_frac <- if (!is.null(params$idclust_min_frac_assigned)) params$idclust_min_frac_assigned else 0.1
idclust_n_dims <- if (!is.null(params$idclust_n_dims)) params$idclust_n_dims else 50
idclust_starting_res <- if (!is.null(params$idclust_starting_resolution)) params$idclust_starting_resolution else 0.1
idclust_resolution <- if (!is.null(params$idclust_resolution)) params$idclust_resolution else 0.8
idclust_starting_k <- if (!is.null(params$idclust_starting_k)) params$idclust_starting_k else 100
idclust_k <- if (!is.null(params$idclust_k)) params$idclust_k else 100
idclust_min_cells <- if (!is.null(params$idclust_min_cells)) params$idclust_min_cells else 100
idclust_plotting <- if (!is.null(params$idclust_plotting)) params$idclust_plotting else TRUE

cat("IDclust Parameters:\n")
cat("  logFC threshold:      ", idclust_logFC_th, "(log2 scale)\n")
cat("  q-value threshold:    ", idclust_qval_th, "\n")
cat("  Min DEGs per split:   ", idclust_min_DEGs, "\n")
cat("  Max recursion depth:  ", idclust_max_depth, "\n")
cat("  Min fraction assigned:", idclust_min_frac, "\n")
cat("  PCA dimensions:       ", idclust_n_dims, "\n")
cat("  Starting resolution:  ", idclust_starting_res, "\n")
cat("  Subsequent resolution:", idclust_resolution, "\n")
cat("  Starting K neighbors: ", idclust_starting_k, "\n")
cat("  Subsequent K:         ", idclust_k, "\n")
cat("  Min cells per cluster:", idclust_min_cells, "\n")
cat("  Internal plotting:    ", idclust_plotting, "\n")

# ==============================================================================
# Determine target clusters for IDclust
# ==============================================================================
idclust_target_clusters <- if (!is.null(params$idclust_target_clusters)) {
  params$idclust_target_clusters
} else {
  names(cluster_sizes)[cluster_sizes >= idclust_min_cells]
}

cat("\nIDclust target clusters (", length(idclust_target_clusters), "):",
    paste(idclust_target_clusters, collapse = ", "), "\n\n")

if (length(idclust_target_clusters) == 0) {
  cat("WARNING: No clusters meet minimum cell requirement for IDclust\n")
} else {

# ==============================================================================
# Run IDclust per cluster
# ==============================================================================

for (clust in idclust_target_clusters) {
  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat("IDclust [", source_label, "]: Processing cluster ", clust, "\n", sep = "")
  cat(paste(rep("-", 70), collapse = ""), "\n")

  cells <- colnames(choir_obj)[clusters == clust]
  n_cells <- length(cells)

  if (n_cells < idclust_min_cells) {
    cat("  Skipping: only", n_cells, "cells (min:", idclust_min_cells, ")\n\n")
    source_idclust_results[[clust]] <- list(cluster = clust, success = FALSE, reason = "too_few_cells")
    next
  }

  cat("  Cells:", n_cells, "\n")

  cluster_out_dir <- file.path(source_idclust_base, paste0("cluster_", clust))
  dir.create(cluster_out_dir, showWarnings = FALSE, recursive = TRUE)

  tryCatch({
    subset_obj <- subset(choir_obj, cells = cells)

    tryCatch({
      if ("RNA" %in% names(subset_obj@assays)) {
        sub_rna_layers <- Layers(subset_obj[["RNA"]])
        sub_counts_layers <- grep("^counts", sub_rna_layers, value = TRUE)
        if (length(sub_counts_layers) > 1) {
          subset_obj[["RNA"]] <- JoinLayers(subset_obj[["RNA"]])
          cat("  Joined RNA layers for subset\n")
        }
      }
    }, error = function(e) cat("  Note: Layer join -", e$message, "\n"))

    # Convert Assay5 -> Assay for IDclust compatibility (requires @counts slot)
    tryCatch({
      if (inherits(subset_obj[["RNA"]], "Assay5")) {
        subset_obj[["RNA"]] <- as(subset_obj[["RNA"]], "Assay")
        cat("  Converted Assay5 -> Assay for IDclust compatibility\n")
      }
    }, error = function(e) cat("  WARNING: Assay conversion failed -", e$message, "\n"))

    DefaultAssay(subset_obj) <- "RNA"

    cat("  Running IDclust::iterative_differential_clustering...\n")
    cat("  (This may take several minutes per cluster)\n")

    t_start <- Sys.time()

    subset_obj <- IDclust::iterative_differential_clustering(
      subset_obj,
      output_dir = cluster_out_dir,
      plotting = isTRUE(idclust_plotting),
      saving = TRUE,
      n_dims = idclust_n_dims,
      dim_red = "pca",
      vizualization_dim_red = "umap",
      logFC.th = idclust_logFC_th,
      qval.th = idclust_qval_th,
      min_frac_cell_assigned = idclust_min_frac,
      limit = idclust_max_depth,
      starting.resolution = idclust_starting_res,
      resolution = idclust_resolution,
      starting.k = idclust_starting_k,
      k = idclust_k
    )

    t_elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 1)
    cat("  IDclust completed in", as.numeric(t_elapsed), "minutes\n")

    if ("IDcluster" %in% colnames(subset_obj@meta.data)) {
      labels <- as.character(subset_obj$IDcluster)
      n_subclusters <- length(unique(labels))
      cat("  SUCCESS: Found", n_subclusters, "hierarchical subclusters\n")
      cat("  Subcluster distribution:\n")
      print(table(labels))

      source_idclust_results[[clust]] <- list(
        cluster = clust,
        success = TRUE,
        labels = data.frame(
          cell = colnames(subset_obj),
          IDcluster = labels,
          stringsAsFactors = FALSE
        ),
        n_subclusters = n_subclusters,
        n_cells = n_cells,
        elapsed_min = as.numeric(t_elapsed)
      )

      idc_summary_file <- file.path(cluster_out_dir, "IDC_summary.qs")
      if (file.exists(idc_summary_file)) {
        tryCatch({
          idc_summary <- qs::qread(idc_summary_file)
          source_idclust_results[[clust]]$idc_summary <- idc_summary
          cat("  IDC_summary loaded (", nrow(idc_summary), " rows)\n", sep = "")
        }, error = function(e) {
          cat("  Note: Could not read IDC_summary:", e$message, "\n")
        })
      }

      tryCatch({
        idc_summary_for_plot <- source_idclust_results[[clust]]$idc_summary
        p_network <- IDclust::plot_cluster_network(
          subset_obj,
          IDC_summary = idc_summary_for_plot
        )
        network_file <- file.path(cluster_out_dir, paste0("cluster_", clust, "_network.png"))
        ggplot2::ggsave(network_file, p_network, width = 10, height = 8, dpi = 150)
        cat("  Network plot saved:", basename(network_file), "\n")
      }, error = function(e) {
        cat("  Note: Network plot failed:", e$message, "\n")
      })

    } else {
      cat("  WARNING: No IDcluster column in result\n")
      cat("  Available columns:", paste(head(colnames(subset_obj@meta.data), 10), collapse = ", "), "\n")
      source_idclust_results[[clust]] <- list(cluster = clust, success = FALSE, reason = "no_IDcluster_column")
    }

  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    source_idclust_results[[clust]] <<- list(cluster = clust, success = FALSE, reason = e$message)
  })

  cat("\n")
}

# ==============================================================================
# Integrate IDclust results for this source
# ==============================================================================
cat("\n--- Integrating IDclust results [", source_label, "] ---\n", sep = "")

idclust_col_name <- paste0("idclust_subcluster", col_suffix)
choir_obj[[idclust_col_name]] <- as.character(choir_obj@meta.data[[cluster_col]])
idclust_successful <- names(source_idclust_results)[sapply(source_idclust_results, function(x) isTRUE(x$success))]

cat("Column:", idclust_col_name, "\n")
cat("IDclust successful:", length(idclust_successful), "/", length(idclust_target_clusters), "\n")

for (clust in idclust_successful) {
  result <- source_idclust_results[[clust]]
  labels_df <- result$labels

  matched <- match(labels_df$cell, colnames(choir_obj))
  valid <- !is.na(matched)

  if (sum(valid) > 0) {
    new_labels <- paste0(clust, "_IDC_", labels_df$IDcluster[valid])
    choir_obj[[idclust_col_name]][matched[valid], 1] <- new_labels
    cat("  Cluster", clust, ":", sum(valid), "cells ->", result$n_subclusters, "subclusters\n")
  }
}

source_idclust_success <- length(idclust_successful) > 0
if (source_idclust_success) {
  idclust_success <- TRUE
  idclust_cols_added <- c(idclust_cols_added, idclust_col_name)
}
all_idclust_results[[current_source]] <- source_idclust_results

if (source_idclust_success) {
  idclust_tab <- table(choir_obj[[idclust_col_name]])
  idclust_summary_df <- data.frame(
    subcluster = names(idclust_tab),
    n_cells = as.vector(idclust_tab)
  )
  idclust_summary_df <- idclust_summary_df[order(idclust_summary_df$subcluster), ]
  write.csv(idclust_summary_df, file.path(source_idclust_base, "idclust_subcluster_summary.csv"), row.names = FALSE)
  cat("  Saved IDclust subcluster summary\n")
}

}  # end idclust_target_clusters check
}  # end run_idclust

}  # ===== END SOURCE LOOP =====


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#          BACKWARD COMPATIBILITY: Create unsuffixed columns for "both"
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# When "both" mode is used, also create unsuffixed columns using the primary
# method's results for backward compatibility with downstream modules
if (is_multi_source) {
  cat("\n--- Creating backward-compatible unsuffixed columns ---\n")

  # Determine primary source for unsuffixed column (prefer choir)
  primary_for_compat <- if ("choir" %in% sources_to_process) "choir" else sources_to_process[1]
  cat("Primary source for unsuffixed columns:", primary_for_compat, "\n")

  primary_suffix <- paste0("_", primary_for_compat)

  if (scice_success) {
    suffixed_scice_col <- paste0("scice_subcluster", primary_suffix)
    if (suffixed_scice_col %in% colnames(choir_obj@meta.data)) {
      choir_obj$scice_subcluster <- choir_obj@meta.data[[suffixed_scice_col]]
      scice_cols_added <- c(scice_cols_added, "scice_subcluster")
      cat("  scice_subcluster <- ", suffixed_scice_col, "\n", sep = "")
    }
  }

  if (idclust_success) {
    suffixed_idclust_col <- paste0("idclust_subcluster", primary_suffix)
    if (suffixed_idclust_col %in% colnames(choir_obj@meta.data)) {
      choir_obj$idclust_subcluster <- choir_obj@meta.data[[suffixed_idclust_col]]
      idclust_cols_added <- c(idclust_cols_added, "idclust_subcluster")
      cat("  idclust_subcluster <- ", suffixed_idclust_col, "\n", sep = "")
    }
  }
}

# Flatten results for backward compatibility
scice_results <- do.call(c, unname(all_scice_results))
idclust_results <- do.call(c, unname(all_idclust_results))


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                       VISUALIZATION AND SAVE
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# ==============================================================================
# Comparison visualization
# ==============================================================================
cat("\n--- Creating visualizations ---\n")

# Base output dir for plots
plot_output_dir <- file.path(out_base, "06_scICE_subclustering")
dir.create(plot_output_dir, showWarnings = FALSE, recursive = TRUE)

if (!"umap" %in% names(choir_obj@reductions)) {
  cat("Running UMAP for visualization...\n")
  reduction_to_use <- "pca"
  for (red in c("integrated.mnn", "harmony", "integrated.cca", "pca")) {
    if (red %in% names(choir_obj@reductions)) {
      reduction_to_use <- red
      break
    }
  }
  choir_obj <- RunUMAP(choir_obj, reduction = reduction_to_use, dims = 1:30, verbose = FALSE)
}

tryCatch({
  plots_list <- list()

  # Original cluster plots (one per source processed)
  for (src in sources_to_process) {
    src_col <- get_cluster_col_for_source(src, colnames(choir_obj@meta.data))
    if (!is.null(src_col) && src_col %in% colnames(choir_obj@meta.data)) {
      p <- DimPlot(choir_obj, reduction = "umap", group.by = src_col,
                    label = TRUE, label.size = 3) +
        ggtitle(paste0("Original (", toupper(src), ")")) + NoLegend()
      plots_list[[paste0("original_", src)]] <- p
    }
  }

  # scICE subcluster plots
  for (col_name in scice_cols_added) {
    if (col_name %in% colnames(choir_obj@meta.data)) {
      p <- DimPlot(choir_obj, reduction = "umap", group.by = col_name,
                    label = TRUE, label.size = 2) +
        ggtitle(paste0("scICE: ", col_name)) +
        theme(legend.text = element_text(size = 6))
      plots_list[[col_name]] <- p
    }
  }

  # IDclust subcluster plots
  for (col_name in idclust_cols_added) {
    if (col_name %in% colnames(choir_obj@meta.data)) {
      p <- DimPlot(choir_obj, reduction = "umap", group.by = col_name,
                    label = TRUE, label.size = 2) +
        ggtitle(paste0("IDclust: ", col_name)) +
        theme(legend.text = element_text(size = 6))
      plots_list[[col_name]] <- p
    }
  }

  if (length(plots_list) > 0) {
    n_plots <- length(plots_list)
    ncol_plots <- min(n_plots, 3)
    p_combined <- wrap_plots(plots_list, ncol = ncol_plots)

    plot_file <- file.path(plot_output_dir, "subclustering_comparison_plot.png")
    ggsave(plot_file, p_combined, width = 7 * ncol_plots,
           height = 7 * ceiling(n_plots / ncol_plots), dpi = 150)
    cat("Saved comparison plot:", plot_file, "\n")
  }

  # Detailed per-column plots
  for (col_name in scice_cols_added) {
    if (col_name %in% colnames(choir_obj@meta.data)) {
      p_full <- DimPlot(choir_obj, reduction = "umap", group.by = col_name,
                         label = TRUE, label.size = 3) +
        ggtitle(paste0(col_name, " (n=", length(unique(choir_obj@meta.data[[col_name]])), ")"))
      save_plot_multi(p_full, paste0(col_name, "_UMAP"),
                      output_dir = plot_output_dir, width = 12, height = 8)
    }
  }

  for (col_name in idclust_cols_added) {
    if (col_name %in% colnames(choir_obj@meta.data)) {
      idclust_plot_dir <- file.path(out_base, "06_IDclust_subclustering")
      dir.create(idclust_plot_dir, showWarnings = FALSE, recursive = TRUE)

      p_full <- DimPlot(choir_obj, reduction = "umap", group.by = col_name,
                         label = TRUE, label.size = 3) +
        ggtitle(paste0(col_name, " (n=", length(unique(choir_obj@meta.data[[col_name]])), ")"))
      save_plot_multi(p_full, paste0(col_name, "_UMAP"),
                      output_dir = idclust_plot_dir, width = 12, height = 8)

      if ("sex" %in% colnames(choir_obj@meta.data)) {
        p_sex <- DimPlot(choir_obj, reduction = "umap", group.by = col_name,
                          split.by = "sex", label = TRUE, label.size = 2) +
          ggtitle(paste0(col_name, " by Sex"))
        save_plot_multi(p_sex, paste0(col_name, "_by_sex"),
                        output_dir = idclust_plot_dir, width = 18, height = 8)
      }
    }
  }

}, error = function(e) {
  cat("WARNING: Could not create plots:", e$message, "\n")
})

# ==============================================================================
# Save results
# ==============================================================================
cat("\n--- Saving results ---\n")

scice_rds <- file.path(output_dirs$objects, "scice_subclustered_object.rds")
saveRDS(choir_obj, scice_rds)
cat("Saved:", scice_rds, "\n")

subclustered_rds <- file.path(output_dirs$objects, "subclustered_object.rds")
saveRDS(choir_obj, subclustered_rds)
cat("Saved:", subclustered_rds, "\n")

save(scice_success, scice_results, idclust_success, idclust_results,
     all_scice_results, all_idclust_results,
     scice_cols_added, idclust_cols_added,
     is_multi_source, sources_to_process,
     choir_obj,
     file = file.path(output_dirs$objects, "06_scice_data.RData"))
cat("Saved: 06_scice_data.RData\n")

# Save per-source scICE summaries
for (src in names(all_scice_results)) {
  src_results <- all_scice_results[[src]]
  src_successful <- names(src_results)[sapply(src_results, function(x) isTRUE(x$success))]
  if (length(src_successful) > 0) {
    scice_col <- paste0("scice_subcluster", if (is_multi_source) paste0("_", src) else "")
    if (scice_col %in% colnames(choir_obj@meta.data)) {
      subcluster_summary <- table(choir_obj@meta.data[[scice_col]])
      summary_df <- data.frame(
        subcluster = names(subcluster_summary),
        n_cells = as.vector(subcluster_summary)
      )
      summary_df <- summary_df[order(summary_df$subcluster), ]
      src_dir <- if (is_multi_source) {
        file.path(out_base, "06_scICE_subclustering", src)
      } else {
        file.path(out_base, "06_scICE_subclustering")
      }
      write.csv(summary_df, file.path(src_dir, "scice_subcluster_summary.csv"), row.names = FALSE)
    }
  }
}

write_readme(plot_output_dir, "scICE Subclustering",
             paste0("Julia-based scLENS/scICE subclustering\n",
                    "Subclustering source: ", subclustering_source, "\n",
                    "Sources processed: ", paste(sources_to_process, collapse = ", "), "\n",
                    if (run_scice && !is.null(scice_k_min)) paste0(
                      "K range: ", scice_k_min, "-", scice_k_max, "\n",
                      "IC threshold: ", scice_ic_threshold, "\n",
                      "Min cells: ", min_cells_subcluster, "\n"
                    ) else "scICE was not run\n"),
             list("input/" = "Count matrices per cluster",
                  "output/" = "scICE results per cluster",
                  "scripts/" = "Julia scripts",
                  "scice_subcluster_summary.csv" = "Subcluster cell counts"))

idclust_readme_dir <- file.path(out_base, "06_IDclust_subclustering")
if (run_idclust && dir.exists(idclust_readme_dir)) {
  write_readme(idclust_readme_dir, "IDclust Subclustering",
               paste0("IDclust: Iterative Differential Clustering\n",
                      "Reference: Prompsy et al., NAR Genomics and Bioinformatics, 2024\n\n",
                      "Subclustering source: ", subclustering_source, "\n",
                      "Sources processed: ", paste(sources_to_process, collapse = ", "), "\n\n",
                      "logFC.th: ", if (exists("idclust_logFC_th")) idclust_logFC_th else "N/A",
                      "  qval.th: ", if (exists("idclust_qval_th")) idclust_qval_th else "N/A", "\n",
                      "max_depth: ", if (exists("idclust_max_depth")) idclust_max_depth else "N/A",
                      "  min_frac: ", if (exists("idclust_min_frac")) idclust_min_frac else "N/A", "\n",
                      "n_dims: ", if (exists("idclust_n_dims")) idclust_n_dims else "N/A",
                      "  min_cells: ", if (exists("idclust_min_cells")) idclust_min_cells else "N/A", "\n"),
               list("cluster_*/" = "Per-cluster IDclust outputs",
                    "idclust_subcluster_summary.csv" = "Subcluster cell counts"))
}

# ==============================================================================
# Final summary
# ==============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat(">>> MODULE 06 COMPLETE <<<\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Subclustering source:", subclustering_source, "\n")
cat("Sources processed:", paste(sources_to_process, collapse = ", "), "\n")
cat("Multi-source mode:", is_multi_source, "\n")
cat("Object loaded from:", basename(loaded_from), "\n\n")

cat("--- scICE ---\n")
cat("  Enabled:", run_scice, "| Julia OK:", scice_julia_ok, "| Success:", scice_success, "\n")
if (scice_success) {
  cat("  Columns added:", paste(scice_cols_added, collapse = ", "), "\n")
  for (col_name in scice_cols_added) {
    if (col_name %in% colnames(choir_obj@meta.data)) {
      cat("    ", col_name, ": ", length(unique(choir_obj@meta.data[[col_name]])), " unique values\n", sep = "")
    }
  }
}

cat("\n--- IDclust ---\n")
cat("  Enabled:", run_idclust, "| Success:", idclust_success, "\n")
if (idclust_success) {
  cat("  Columns added:", paste(idclust_cols_added, collapse = ", "), "\n")
  for (col_name in idclust_cols_added) {
    if (col_name %in% colnames(choir_obj@meta.data)) {
      cat("    ", col_name, ": ", length(unique(choir_obj@meta.data[[col_name]])), " unique values\n", sep = "")
    }
  }
  # Per-source timings
  for (src in names(all_idclust_results)) {
    src_results <- all_idclust_results[[src]]
    src_successful <- names(src_results)[sapply(src_results, function(x) isTRUE(x$success))]
    if (length(src_successful) > 0) {
      timings <- sapply(src_results[src_successful], function(x) x$elapsed_min)
      cat("  ", toupper(src), " per-cluster time (min): ", paste(round(timings, 1), collapse = ", "), "\n", sep = "")
    }
  }
}

cat("\nMetadata columns added:\n")
for (col_name in c(scice_cols_added, idclust_cols_added)) {
  cat("  ", col_name, "\n")
}

cat("\nOutput directories:\n")
cat("  scICE:  ", file.path(out_base, "06_scICE_subclustering"), "\n")
if (run_idclust) cat("  IDclust:", file.path(out_base, "06_IDclust_subclustering"), "\n")