#!/usr/bin/env Rscript
# ==============================================================================
# Transfer MapMyCells annotations to preprocessed Seurat objects
# ==============================================================================
# This script transfers cell-type labels from MapMyCells hierarchical mapping
# output to:
#   1. Per-sample DecontX-corrected Seurat objects (8_DecontX_correction/)
#   2. Integrated CHOIR Seurat object (9_CHOIR_integration/)
#
# MapMyCells cell_id format: BARCODE_SAMPLEID_all
# e.g., AAACCCAGTGACTCTA_c15_12_all
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

cat("==============================================================================\n")
cat("  MapMyCells Label Transfer to Seurat Objects (R)\n")
cat("==============================================================================\n\n")

# --- Configuration -----------------------------------------------------------
base_dir <- "/scicore/home/doetsch/kaiser0001/Single_cell_paper/Datasets/Human_Covid_LV_ChP_PMID_34153974"

mapmycells_csv <- file.path(base_dir,
  "MapMyCells_Human_Covid/All_combined_CHOIR_integrated_10xWholeHumanBrain(CCN202210140)_HierarchicalMapping_UTC_1770459613213.csv")

decontx_dir <- file.path(base_dir,
  "Output_dir_Human_Covid_LV_ChP/Single_cell_preprocessed/8_DecontX_correction")

choir_dir <- file.path(base_dir,
  "Output_dir_Human_Covid_LV_ChP/Single_cell_preprocessed/9_CHOIR_integration/All_combined")

sample_ids <- c("c15_12", "c15_17", "c16_18", "CP_90", "CP_91")

# Columns to transfer (will be suffixed with _MapMyCells)
annotation_cols <- c(
  "supercluster_label", "supercluster_name", "supercluster_bootstrapping_probability",
  "cluster_label", "cluster_name", "cluster_bootstrapping_probability",
  "subcluster_label", "subcluster_name", "subcluster_alias",
  "subcluster_bootstrapping_probability"
)

# --- Load MapMyCells annotations ----------------------------------------------
cat("[1] Loading MapMyCells annotations...\n")
mmc <- read.csv(mapmycells_csv, stringsAsFactors = FALSE)
cat(sprintf("    Total annotated cells: %d\n", nrow(mmc)))
cat(sprintf("    Columns: %s\n", paste(colnames(mmc), collapse = ", ")))

# Parse cell_id into barcode and sample_id
# Format: BARCODE_SAMPLEID_all
# Sample IDs can contain underscores (e.g., c15_12), so we strip the trailing "_all"
# and then match against known sample IDs
mmc$cell_id_clean <- sub("_all$", "", mmc$cell_id)

# Extract sample ID by matching known sample IDs from the end of cell_id_clean
mmc$parsed_sample <- NA_character_
mmc$parsed_barcode <- NA_character_

for (sid in sample_ids) {
  pattern <- paste0("_", sid, "$")
  mask <- grepl(pattern, mmc$cell_id_clean)
  mmc$parsed_sample[mask] <- sid
  mmc$parsed_barcode[mask] <- sub(pattern, "", mmc$cell_id_clean[mask])
}

cat("\n    Parsed sample distribution:\n")
sample_counts <- table(mmc$parsed_sample, useNA = "ifany")
for (s in names(sample_counts)) {
  cat(sprintf("      %-10s : %d cells\n", ifelse(is.na(s), "UNMATCHED", s), sample_counts[s]))
}

n_unmatched_parse <- sum(is.na(mmc$parsed_sample))
if (n_unmatched_parse > 0) {
  cat(sprintf("\n    WARNING: %d cells could not be parsed to any sample ID!\n", n_unmatched_parse))
  cat("    Example unmatched cell_ids:\n")
  unmatched_examples <- head(mmc$cell_id[is.na(mmc$parsed_sample)], 5)
  for (ex in unmatched_examples) cat(sprintf("      %s\n", ex))
}

# --- Helper: print barcode inspection -----------------------------------------
print_barcode_inspection <- function(obj, label, n = 5) {
  barcodes <- colnames(obj)
  cat(sprintf("\n    First %d barcodes in %s:\n", min(n, length(barcodes)), label))
  for (i in seq_len(min(n, length(barcodes)))) {
    bc <- barcodes[i]
    # Check which MapMyCells columns exist
    has_labels <- any(grepl("_MapMyCells$", colnames(obj@meta.data)))
    if (has_labels) {
      sc_name <- obj@meta.data[bc, "supercluster_name_MapMyCells"]
      cl_name <- obj@meta.data[bc, "cluster_name_MapMyCells"]
      sub_name <- obj@meta.data[bc, "subcluster_name_MapMyCells"]
      cat(sprintf("      [%d] %s => supercluster=%s | cluster=%s | subcluster=%s\n",
                  i, bc,
                  ifelse(is.na(sc_name), "NA", sc_name),
                  ifelse(is.na(cl_name), "NA", cl_name),
                  ifelse(is.na(sub_name), "NA", sub_name)))
    } else {
      cat(sprintf("      [%d] %s => NO MapMyCells columns found!\n", i, bc))
    }
  }
}

# ==============================================================================
# PART A: Per-sample DecontX objects
# ==============================================================================
cat("\n==============================================================================\n")
cat("  PART A: Transferring labels to per-sample DecontX objects\n")
cat("==============================================================================\n")

overall_stats <- data.frame(
  sample = character(),
  total_cells_obj = integer(),
  total_cells_mmc = integer(),
  matched = integer(),
  unmatched_obj = integer(),
  unmatched_mmc = integer(),
  match_pct = numeric(),
  stringsAsFactors = FALSE
)

for (sid in sample_ids) {
  cat(sprintf("\n--- Processing sample: %s ---\n", sid))

  rds_path <- file.path(decontx_dir, sid, paste0(sid, "_decontX_corrected.rds"))

  if (!file.exists(rds_path)) {
    cat(sprintf("    SKIP: File not found: %s\n", rds_path))
    next
  }

  # Load Seurat object
  cat(sprintf("    Loading: %s\n", basename(rds_path)))
  obj <- readRDS(rds_path)
  obj_barcodes <- colnames(obj)
  cat(sprintf("    Cells in object: %d\n", length(obj_barcodes)))

  # Get MapMyCells annotations for this sample
  mmc_sample <- mmc[mmc$parsed_sample == sid & !is.na(mmc$parsed_sample), ]
  cat(sprintf("    Cells in MapMyCells for %s: %d\n", sid, nrow(mmc_sample)))

  if (nrow(mmc_sample) == 0) {
    cat("    WARNING: No MapMyCells annotations found for this sample!\n")
    next
  }

  # Inspect barcode formats
  cat(sprintf("    Object barcode examples: %s\n",
              paste(head(obj_barcodes, 3), collapse = ", ")))
  cat(sprintf("    MapMyCells barcode examples: %s\n",
              paste(head(mmc_sample$parsed_barcode, 3), collapse = ", ")))

  # Try direct match first
  direct_match <- sum(mmc_sample$parsed_barcode %in% obj_barcodes)
  cat(sprintf("    Direct barcode match: %d / %d\n", direct_match, nrow(mmc_sample)))

  # If direct match is poor, try stripping -1 suffix from object barcodes
  if (direct_match < nrow(mmc_sample) * 0.5) {
    obj_barcodes_stripped <- sub("-1$", "", obj_barcodes)
    stripped_match <- sum(mmc_sample$parsed_barcode %in% obj_barcodes_stripped)
    cat(sprintf("    After stripping '-1' from object barcodes: %d / %d\n",
                stripped_match, nrow(mmc_sample)))

    if (stripped_match > direct_match) {
      cat("    Using stripped barcodes for matching.\n")
      # Create lookup: stripped -> original
      barcode_lookup <- setNames(obj_barcodes, obj_barcodes_stripped)
      match_type <- "stripped"
    } else {
      barcode_lookup <- setNames(obj_barcodes, obj_barcodes)
      match_type <- "direct"
    }
  } else {
    barcode_lookup <- setNames(obj_barcodes, obj_barcodes)
    match_type <- "direct"
  }

  # Build annotation dataframe indexed by object barcodes
  lookup_keys <- names(barcode_lookup)
  mmc_sample_indexed <- mmc_sample
  rownames(mmc_sample_indexed) <- mmc_sample_indexed$parsed_barcode

  matched_keys <- intersect(lookup_keys, mmc_sample_indexed$parsed_barcode)
  matched_obj_barcodes <- barcode_lookup[matched_keys]

  n_matched <- length(matched_keys)
  n_total_obj <- length(obj_barcodes)
  n_total_mmc <- nrow(mmc_sample)
  match_pct <- round(n_matched / n_total_obj * 100, 2)

  cat(sprintf("\n    MATCHING SUMMARY (%s):\n", match_type))
  cat(sprintf("      Cells in object:         %d\n", n_total_obj))
  cat(sprintf("      Cells in MapMyCells:      %d\n", n_total_mmc))
  cat(sprintf("      Successfully matched:     %d (%.2f%% of object)\n", n_matched, match_pct))
  cat(sprintf("      Unmatched in object:      %d\n", n_total_obj - n_matched))
  cat(sprintf("      Unmatched in MapMyCells:  %d\n", n_total_mmc - n_matched))

  # Add annotation columns
  for (col in annotation_cols) {
    new_col <- paste0(col, "_MapMyCells")
    obj@meta.data[[new_col]] <- NA
    if (n_matched > 0) {
      obj@meta.data[matched_obj_barcodes, new_col] <- mmc_sample_indexed[matched_keys, col]
    }
  }

  # Show annotation distribution
  cat("\n    Supercluster distribution (matched cells):\n")
  sc_dist <- table(obj@meta.data$supercluster_name_MapMyCells, useNA = "ifany")
  for (nm in names(sc_dist)) {
    cat(sprintf("      %-30s : %d\n", ifelse(is.na(nm), "NA (unmatched)", nm), sc_dist[nm]))
  }

  cat("\n    Cluster distribution (top 10, matched cells):\n")
  cl_dist <- sort(table(obj@meta.data$cluster_name_MapMyCells, useNA = "no"), decreasing = TRUE)
  for (i in seq_len(min(10, length(cl_dist)))) {
    cat(sprintf("      %-30s : %d\n", names(cl_dist)[i], cl_dist[i]))
  }

  # Save
  cat(sprintf("\n    Saving to: %s\n", rds_path))
  saveRDS(obj, rds_path)

  # Reload and verify
  cat("    Reloading to verify...\n")
  obj_verify <- readRDS(rds_path)

  mmc_cols_found <- grep("_MapMyCells$", colnames(obj_verify@meta.data), value = TRUE)
  cat(sprintf("    MapMyCells columns found after reload: %d\n", length(mmc_cols_found)))
  cat(sprintf("      %s\n", paste(mmc_cols_found, collapse = ", ")))

  n_annotated <- sum(!is.na(obj_verify@meta.data$supercluster_name_MapMyCells))
  cat(sprintf("    Cells with annotations after reload: %d / %d\n",
              n_annotated, ncol(obj_verify)))

  print_barcode_inspection(obj_verify, paste0(sid, " (reloaded)"))

  overall_stats <- rbind(overall_stats, data.frame(
    sample = sid,
    total_cells_obj = n_total_obj,
    total_cells_mmc = n_total_mmc,
    matched = n_matched,
    unmatched_obj = n_total_obj - n_matched,
    unmatched_mmc = n_total_mmc - n_matched,
    match_pct = match_pct,
    stringsAsFactors = FALSE
  ))

  rm(obj, obj_verify, mmc_sample, mmc_sample_indexed)
  gc(verbose = FALSE)
}

# ==============================================================================
# PART B: Integrated CHOIR object (.rds)
# ==============================================================================
cat("\n==============================================================================\n")
cat("  PART B: Transferring labels to integrated CHOIR object (.rds)\n")
cat("==============================================================================\n")

choir_rds <- file.path(choir_dir, "All_combined_CHOIR_integrated.rds")

if (file.exists(choir_rds)) {
  cat(sprintf("\n    Loading: %s\n", basename(choir_rds)))
  obj_int <- readRDS(choir_rds)
  int_barcodes <- colnames(obj_int)
  cat(sprintf("    Cells in integrated object: %d\n", length(int_barcodes)))
  cat(sprintf("    Integrated barcode examples: %s\n",
              paste(head(int_barcodes, 5), collapse = ", ")))
  cat(sprintf("    MapMyCells cell_id examples (clean): %s\n",
              paste(head(mmc$cell_id_clean, 5), collapse = ", ")))

  # Try matching strategies
  # Strategy 1: Direct match cell_id_clean vs integrated barcodes
  direct_match <- sum(mmc$cell_id_clean %in% int_barcodes)
  cat(sprintf("\n    Strategy 1 - Direct match (cell_id_clean): %d / %d\n",
              direct_match, nrow(mmc)))

  # Strategy 2: cell_id (with _all) vs integrated barcodes
  full_match <- sum(mmc$cell_id %in% int_barcodes)
  cat(sprintf("    Strategy 2 - Full cell_id match: %d / %d\n",
              full_match, nrow(mmc)))

  # Strategy 3: Strip -1 from integrated barcodes
  int_barcodes_stripped <- sub("-1$", "", int_barcodes)
  stripped_match <- sum(mmc$cell_id_clean %in% int_barcodes_stripped)
  cat(sprintf("    Strategy 3 - Stripped '-1' from integrated: %d / %d\n",
              stripped_match, nrow(mmc)))

  # Strategy 4: Try matching with sample suffix patterns
  # Integrated barcodes might be BARCODE-1_SAMPLEID or BARCODE_SAMPLEID
  # Build a flexible lookup
  # Extract the base barcode + sample from integrated barcodes
  int_parsed <- data.frame(
    original = int_barcodes,
    stripped = int_barcodes_stripped,
    stringsAsFactors = FALSE
  )

  # Choose best strategy
  strategies <- c(direct_match, full_match, stripped_match)
  best_strat <- which.max(strategies)
  cat(sprintf("\n    Best matching strategy: %d (matched %d cells)\n",
              best_strat, max(strategies)))

  if (best_strat == 1) {
    mmc_key <- mmc$cell_id_clean
    int_key <- int_barcodes
  } else if (best_strat == 2) {
    mmc_key <- mmc$cell_id
    int_key <- int_barcodes
  } else {
    mmc_key <- mmc$cell_id_clean
    int_key <- int_barcodes_stripped
  }

  # Create lookup from mmc_key
  rownames(mmc) <- mmc_key
  int_lookup <- setNames(int_barcodes, int_key)

  matched_keys <- intersect(int_key, mmc_key)
  matched_int_barcodes <- int_lookup[matched_keys]

  n_matched_int <- length(matched_keys)
  match_pct_int <- round(n_matched_int / length(int_barcodes) * 100, 2)

  cat(sprintf("\n    INTEGRATED MATCHING SUMMARY:\n"))
  cat(sprintf("      Cells in integrated object: %d\n", length(int_barcodes)))
  cat(sprintf("      Cells in MapMyCells:        %d\n", nrow(mmc)))
  cat(sprintf("      Successfully matched:       %d (%.2f%% of integrated)\n",
              n_matched_int, match_pct_int))
  cat(sprintf("      Unmatched in integrated:    %d\n",
              length(int_barcodes) - n_matched_int))
  cat(sprintf("      Unmatched in MapMyCells:    %d\n",
              nrow(mmc) - n_matched_int))

  # Add annotations
  for (col in annotation_cols) {
    new_col <- paste0(col, "_MapMyCells")
    obj_int@meta.data[[new_col]] <- NA
    if (n_matched_int > 0) {
      obj_int@meta.data[matched_int_barcodes, new_col] <- mmc[matched_keys, col]
    }
  }

  # Show distributions
  cat("\n    Supercluster distribution (integrated, matched cells):\n")
  sc_dist_int <- table(obj_int@meta.data$supercluster_name_MapMyCells, useNA = "ifany")
  for (nm in names(sc_dist_int)) {
    cat(sprintf("      %-30s : %d\n", ifelse(is.na(nm), "NA (unmatched)", nm), sc_dist_int[nm]))
  }

  cat("\n    Cluster distribution (integrated, top 15):\n")
  cl_dist_int <- sort(table(obj_int@meta.data$cluster_name_MapMyCells, useNA = "no"), decreasing = TRUE)
  for (i in seq_len(min(15, length(cl_dist_int)))) {
    cat(sprintf("      %-30s : %d\n", names(cl_dist_int)[i], cl_dist_int[i]))
  }

  # Save
  cat(sprintf("\n    Saving integrated object to: %s\n", choir_rds))
  saveRDS(obj_int, choir_rds)

  # Reload and verify
  cat("    Reloading integrated object to verify...\n")
  obj_int_verify <- readRDS(choir_rds)

  mmc_cols_int <- grep("_MapMyCells$", colnames(obj_int_verify@meta.data), value = TRUE)
  cat(sprintf("    MapMyCells columns found after reload: %d\n", length(mmc_cols_int)))

  n_annotated_int <- sum(!is.na(obj_int_verify@meta.data$supercluster_name_MapMyCells))
  cat(sprintf("    Cells with annotations after reload: %d / %d\n",
              n_annotated_int, ncol(obj_int_verify)))

  print_barcode_inspection(obj_int_verify, "Integrated CHOIR (reloaded)")

  rm(obj_int, obj_int_verify)
  gc(verbose = FALSE)
} else {
  cat(sprintf("    SKIP: Integrated RDS not found: %s\n", choir_rds))
}

# ==============================================================================
# OVERALL SUMMARY
# ==============================================================================
cat("\n==============================================================================\n")
cat("  OVERALL SUMMARY\n")
cat("==============================================================================\n\n")

cat("  Per-sample transfer results:\n")
cat(sprintf("  %-10s | %8s | %8s | %8s | %8s | %8s | %6s\n",
            "Sample", "Obj", "MMC", "Matched", "Unm.Obj", "Unm.MMC", "Pct"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (i in seq_len(nrow(overall_stats))) {
  r <- overall_stats[i, ]
  cat(sprintf("  %-10s | %8d | %8d | %8d | %8d | %8d | %5.1f%%\n",
              r$sample, r$total_cells_obj, r$total_cells_mmc,
              r$matched, r$unmatched_obj, r$unmatched_mmc, r$match_pct))
}
cat(paste(rep("-", 75), collapse = ""), "\n")
cat(sprintf("  %-10s | %8d | %8d | %8d | %8d | %8d | %5.1f%%\n",
            "TOTAL",
            sum(overall_stats$total_cells_obj),
            sum(overall_stats$total_cells_mmc),
            sum(overall_stats$matched),
            sum(overall_stats$unmatched_obj),
            sum(overall_stats$unmatched_mmc),
            round(sum(overall_stats$matched) / sum(overall_stats$total_cells_obj) * 100, 1)))

cat("\n[DONE] R label transfer complete.\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
