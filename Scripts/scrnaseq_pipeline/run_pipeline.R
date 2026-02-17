#!/usr/bin/env Rscript
# ==============================================================================
# MAIN PIPELINE ORCHESTRATION SCRIPT (MULTI-SAMPLE VERSION)
# ==============================================================================
#
# This script runs the complete scRNA-seq downstream analysis pipeline.
# It processes scCDC-corrected files from preprocessing Step 7.
#
# Usage examples:
#   Rscript run_pipeline.R
#   Rscript run_pipeline.R --modules 0,1,2,3
#   Rscript run_pipeline.R --modules 0,1,2,2b,3,3b  # Include imputation modules
#   Rscript run_pipeline.R --start 0 --stop 7
#   Rscript run_pipeline.R --modules 10             # Run only HTML report
#
# Module mapping:
#   0   -> 00_setup_environment.R
#   1   -> 01_load_data.R
#   2   -> 02_qc_validation.R
#   2b  -> 02b_imputation.R           (optional, afMF counts-based imputation)
#   3   -> 03_normalization.R
#   3b  -> 03b_alra_imputation.R      (optional, ALRA normalized-data imputation)
#   4   -> 04_integration.R
#   5   -> 05_choir_clustering.R
#   6   -> 06_scice_subclustering.R
#   7   -> 07_leiden_clustering.R
#   7b  -> 07b_clts_renormalization.R (optional, CLTS re-normalization)
#   8   -> 08_differential_expression.R
#   9   -> 09_gene_visualization.R
#   10  -> 10_final_summary.R
#   11  -> 11_html_report.R
#
# Note: Modules 2b, 3b, and 7b are optional intermediate modules.
#       - 2b runs automatically after 2 if imputation_method="afmf" or "both"
#       - 3b runs automatically after 3 if imputation_method="alra" or "both"
#         AND selected_normalization_method is NOT SCTransform
#       - 7b runs automatically after 7 if run_clts_renormalization=TRUE
#       - Can also be explicitly requested via --modules 2b, 3b, or 7b
#
# ==============================================================================

cat("================================================================================\n")
cat("CHOROID PLEXUS scRNA-seq DOWNSTREAM ANALYSIS PIPELINE\n")
cat("================================================================================\n\n")

# ==============================================================================
# Determine PIPELINE_DIR FIRST (before anything else)
# ==============================================================================

get_script_dir <- function() {
  # Method 1: commandArgs (works with Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("--file=", "", file_arg[1]))
    return(dirname(script_path))
  }

  # Method 2: Environment variable
  env_dir <- Sys.getenv("PIPELINE_DIR", unset = "")
  if (env_dir != "" && dir.exists(env_dir)) {
    return(env_dir)
  }

  # Method 3: check current directory
  if (file.exists("run_pipeline.R") && dir.exists("modules") && dir.exists("config")) {
    return(normalizePath("."))
  }

  stop("Cannot determine script directory. Run from pipeline directory or use Rscript.")
}

PIPELINE_DIR <- get_script_dir()

# Ensure modules can find the pipeline root if they use Sys.getenv("PIPELINE_DIR")
Sys.setenv(PIPELINE_DIR = PIPELINE_DIR)

cat("PIPELINE_DIR:", PIPELINE_DIR, "\n\n")

# ==============================================================================
# Parse command line arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

modules_to_run <- NULL
start_module <- 0
stop_module <- 11  # HTML is module 11

i <- 1
while (i <= length(args)) {
  if (args[i] == "--modules" && i < length(args)) {
    # Accept: "1,2,3" or "11" or "2b,3,3b,7b" (mixed integer and string)
    raw_mods <- strsplit(args[i + 1], ",")[[1]]
    modules_to_run <- trimws(raw_mods)
    i <- i + 2
  } else if (args[i] == "--start" && i < length(args)) {
    start_module <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--stop" && i < length(args)) {
    stop_module <- as.numeric(args[i + 1])
    i <- i + 2
  } else {
    i <- i + 1
  }
}

cat("Configuration:\n")
cat(
  "  Modules:",
  if (is.null(modules_to_run) || length(modules_to_run) == 0) {
    paste0(start_module, " to ", stop_module, " (plus optional 2b, 3b, 7b if enabled)")
  } else {
    paste(modules_to_run, collapse = ", ")
  },
  "\n\n"
)

# ==============================================================================
# Setup paths
# ==============================================================================

modules_dir <- file.path(PIPELINE_DIR, "modules")
config_dir  <- file.path(PIPELINE_DIR, "config")

cat("Directories:\n")
cat("  Pipeline:", PIPELINE_DIR, "\n")
cat("  Modules:", modules_dir, "\n")
cat("  Config:",  config_dir, "\n\n")

if (!dir.exists(modules_dir)) stop("Modules directory not found: ", modules_dir)

# ==============================================================================
# Define modules
# ==============================================================================
# Module IDs use numeric sort order:
#   - Main modules: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
#   - Sub-modules: 2.5 (displayed as 2b), 3.5 (displayed as 3b), 7.5 (displayed as 7b)
#
# Sub-modules (2b, 3b, 7b) are:
#   - Automatically included when running a range that spans them
#   - Controlled by params (imputation_method, run_clts_renormalization)
#   - Can be explicitly requested via --modules
# ==============================================================================

module_list <- list(
  list(num = 0,    id = "0",   name = "Environment Setup",        file = "00_setup_environment.R",       required = TRUE,  auto = TRUE),
  list(num = 1,    id = "1",   name = "Load Data",                file = "01_load_data.R",               required = TRUE,  auto = TRUE),
  list(num = 2,    id = "2",   name = "QC Validation",            file = "02_qc_validation.R",           required = TRUE,  auto = TRUE),
  list(num = 2.5,  id = "2b",  name = "Imputation (afMF)",        file = "02b_imputation.R",             required = FALSE, auto = FALSE),
  list(num = 3,    id = "3",   name = "Normalization",            file = "03_normalization.R",           required = TRUE,  auto = TRUE),
  list(num = 3.5,  id = "3b",  name = "Imputation (ALRA)",        file = "03b_alra_imputation.R",        required = FALSE, auto = FALSE),
  list(num = 4,    id = "4",   name = "Integration",              file = "04_integration.R",             required = FALSE, auto = TRUE),
  list(num = 5,    id = "5",   name = "CHOIR Clustering",         file = "05_choir_clustering.R",        required = FALSE, auto = TRUE),
  list(num = 6,    id = "6",   name = "scICE Subclustering",      file = "06_scice_subclustering.R",     required = FALSE, auto = TRUE),
  list(num = 7,    id = "7",   name = "Leiden Clustering",        file = "07_leiden_clustering.R",       required = TRUE,  auto = TRUE),
  list(num = 7.5,  id = "7b",  name = "CLTS Re-normalization",    file = "07b_clts_renormalization.R",   required = FALSE, auto = FALSE),
  list(num = 8,    id = "8",   name = "Differential Expression",  file = "08_differential_expression.R", required = FALSE, auto = TRUE),
  list(num = 9,    id = "9",   name = "Gene Visualization",       file = "09_gene_visualization.R",      required = FALSE, auto = TRUE),
  list(num = 10,   id = "10",  name = "Final Summary",            file = "10_final_summary.R",           required = TRUE,  auto = TRUE),
  list(num = 11,   id = "11",  name = "HTML Report",              file = "11_html_report.R",             required = FALSE, auto = TRUE)
)

# ==============================================================================
# Load params to check which optional modules are enabled
# ==============================================================================

params_file <- file.path(config_dir, "params.R")
params <- NULL

if (file.exists(params_file)) {
  # Source params in a temporary environment to get settings
  temp_env <- new.env()
  tryCatch({
    source(params_file, local = temp_env)
    if (exists("params", envir = temp_env)) {
      params <- get("params", envir = temp_env)
    }
  }, error = function(e) {
    cat("Note: Could not pre-load params.R for module selection:", e$message, "\n")
  })
}

# Determine if optional modules should run
# afMF imputation (Module 2b): runs if imputation_method = "afmf" or "both"
run_afmf <- if (!is.null(params) && !is.null(params$imputation_method)) {
  params$imputation_method %in% c("afmf", "both")
} else {
  FALSE
}

# ALRA imputation (Module 3b): runs if imputation_method = "alra" or "both"
# Note: Module 3b will also check at runtime if SCTransform was selected and skip if so
run_alra <- if (!is.null(params) && !is.null(params$imputation_method)) {
  params$imputation_method %in% c("alra", "both")
} else {
  FALSE
}

# CLTS re-normalization (Module 7b)
run_clts <- if (!is.null(params) && !is.null(params$run_clts_renormalization)) {
  isTRUE(params$run_clts_renormalization)
} else {
  FALSE
}

cat("Optional module settings (from params.R):\n")
cat("  imputation_method:", if (!is.null(params$imputation_method)) params$imputation_method else "not set", "\n")
cat("  Module 2b (afMF):", run_afmf, "\n")
cat("  Module 3b (ALRA):", run_alra, "(will skip if SCTransform selected)\n")
cat("  Module 7b (CLTS):", run_clts, "\n\n")

# ==============================================================================
# Determine which modules to run
# ==============================================================================

# Helper function to parse module ID (handles "2b", "3b", "7b", etc.)
parse_module_id <- function(id_str) {
  id_str <- trimws(tolower(id_str))

  # Check for special IDs
  if (id_str == "2b") return(2.5)
  if (id_str == "3b") return(3.5)
  if (id_str == "7b") return(7.5)

  # Try numeric conversion
  num <- suppressWarnings(as.numeric(id_str))
  if (!is.na(num)) return(num)

  return(NA)
}

# Helper to get display ID
get_display_id <- function(num) {
  if (num == 2.5) return("2b")
  if (num == 3.5) return("3b")
  if (num == 7.5) return("7b")
  return(sprintf("%d", as.integer(num)))
}

if (!is.null(modules_to_run) && length(modules_to_run) > 0) {
  # Explicit module list provided
  modules_selected <- sapply(modules_to_run, parse_module_id)
  modules_selected <- modules_selected[!is.na(modules_selected)]
  modules_selected <- sort(unique(modules_selected))
} else {
  # Range-based selection
  modules_selected <- c()

  for (mod in module_list) {
    # Check if module is in range
    in_range <- mod$num >= start_module & mod$num <= stop_module

    if (in_range) {
      # For auto modules, always include if in range
      if (isTRUE(mod$auto)) {
        modules_selected <- c(modules_selected, mod$num)
      }
      # For sub-modules, include if enabled in params
      else if (mod$id == "2b" && run_afmf) {
        modules_selected <- c(modules_selected, mod$num)
      }
      else if (mod$id == "3b" && run_alra) {
        modules_selected <- c(modules_selected, mod$num)
      }
      else if (mod$id == "7b" && run_clts) {
        modules_selected <- c(modules_selected, mod$num)
      }
    }
  }

  modules_selected <- sort(unique(modules_selected))
}

cat("Modules to run:\n")
for (mod in module_list) {
  if (mod$num %in% modules_selected) {
    status <- if (file.exists(file.path(modules_dir, mod$file))) "[FOUND]" else "[MISSING]"
    cat(sprintf("  %3s  %-25s %s\n", mod$id, mod$name, status))
  }
}
cat("\n")

# ==============================================================================
# Run modules
# ==============================================================================

pipeline_start <- Sys.time()
module_times <- list()
failed_modules <- c()

for (mod in module_list) {
  if (mod$num %in% modules_selected) {
    module_path <- file.path(modules_dir, mod$file)

    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat(sprintf("RUNNING MODULE %s: %s\n", mod$id, mod$name))
    cat(paste(rep("=", 80), collapse = ""), "\n\n")

    if (!file.exists(module_path)) {
      cat("WARNING: Module file not found:", module_path, "\n")
      if (isTRUE(mod$required)) {
        stop("Required module missing.")
      } else {
        cat("Skipping optional module.\n\n")
        next
      }
    }

    module_start <- Sys.time()

    result <- tryCatch({
      # Capture and print the call stack at the moment an error occurs.
      # This is much more informative than traceback() after tryCatch.
      withCallingHandlers(
        source(module_path, local = FALSE),
        error = function(e) {
          cat("\n*** ERROR DURING EXECUTION OF MODULE", mod$id, "(", mod$name, ") ***\n")
          cat("Error message: ", conditionMessage(e), "\n\n", sep = "")

          cat("--- CALL STACK (most recent last) ---\n")
          calls <- sys.calls()
          # Print a compact stack; still very readable in slurm logs
          for (k in seq_along(calls)) {
            cat(sprintf("  [%02d] %s\n", k, deparse(calls[[k]])))
          }
          cat("\n--- SESSION INFO (abbrev) ---\n")
          cat("R: ", getRversion(), "\n", sep = "")
          cat("Module file: ", module_path, "\n\n", sep = "")

          # Re-throw so outer tryCatch marks the module as FAILED
          stop(e)
        }
      )

      "SUCCESS"
    }, error = function(e) {
      # Outer handler: single place that converts errors to FAILED status
      paste("FAILED:", conditionMessage(e))
    })

    module_end <- Sys.time()
    module_times[[mod$id]] <- as.numeric(difftime(module_end, module_start, units = "mins"))

    if (grepl("^FAILED", result)) {
      failed_modules <- c(failed_modules, mod$id)
      if (isTRUE(mod$required)) {
        cat("\nRequired module failed. Stopping pipeline.\n")
        break
      }
    }

    cat(sprintf("\nModule %s completed in %.2f minutes.\n\n",
                mod$id, module_times[[mod$id]]))
  }
}


# ==============================================================================
# Summary
# ==============================================================================

pipeline_duration <- difftime(Sys.time(), pipeline_start, units = "mins")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("PIPELINE SUMMARY\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Total runtime:", round(as.numeric(pipeline_duration), 2), "minutes\n\n")

cat("Module timings:\n")
for (mod in module_list) {
  if (mod$id %in% names(module_times)) {
    status <- if (mod$id %in% failed_modules) "FAILED" else "OK"
    cat(sprintf("  %3s  %-25s %.2f min [%s]\n",
                mod$id, mod$name, module_times[[mod$id]], status))
  }
}

if (length(failed_modules) > 0) {
  cat("\nFailed modules:", paste(failed_modules, collapse = ", "), "\n")
  cat("\n*** PIPELINE COMPLETED WITH ERRORS ***\n")
  quit(status = 1)
} else {
  cat("\n*** PIPELINE COMPLETED SUCCESSFULLY ***\n")
  quit(status = 0)
}
