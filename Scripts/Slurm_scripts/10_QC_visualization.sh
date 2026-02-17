#!/bin/bash
#SBATCH --job-name=QC_viz
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --output=logs/qc_viz/qc_viz_%j.out
#SBATCH --error=logs/qc_viz/qc_viz_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ==============================================================================
# STEP 10: QC Visualization
# ==============================================================================
#
# UNIVERSAL MODULAR VERSION - Reads samples from samplesheet.csv
#
# Generates comprehensive QC visualization plots from all pipeline steps.
#
# Usage:
#   sbatch 10_QC_visualization.sh
#
# ==============================================================================

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

if [[ -n "${PROJECT_ROOT:-}" ]]; then
    BASE_DIR="${PROJECT_ROOT}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    BASE_DIR="${SLURM_SUBMIT_DIR}"
else
    BASE_DIR="/scicore/home/doetsch/kaiser0001/Revision_NatureComm_Sex/Vandebroucke_fibroblast_paper"
fi

SAMPLESHEET="${SAMPLESHEET:-${BASE_DIR}/samplesheet.csv}"

get_dataset_name() {
    tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f12
}

if [[ -z "${PREPROCESS_DIR:-}" ]]; then
    DATASET_NAME=$(get_dataset_name)
    if [[ -n "$DATASET_NAME" ]]; then
        PREPROCESS_DIR="${BASE_DIR}/Output_dir_${DATASET_NAME}/Single_cell_preprocessed"
    else
        PREPROCESS_DIR="${BASE_DIR}"
    fi
fi

# Directories - UPDATED for new numbering
DROPLETQC_DIR="${PREPROCESS_DIR}/2_DropletQC_output"
QCLUS_DIR="${PREPROCESS_DIR}/3_qClus_empty_droplets"
VAEDA_DIR="${PREPROCESS_DIR}/4_Vaeda_doublet_detection"
SEURAT_PRE_DIR="${PREPROCESS_DIR}/5_Seurat_conversion"
DOUBLET_DIR="${PREPROCESS_DIR}/6_Doublet_consensus"
SCCDC_DIR="${PREPROCESS_DIR}/7_scCDC_correction"
DECONTX_DIR="${PREPROCESS_DIR}/8_DecontX_correction"
CHOIR_DIR="${PREPROCESS_DIR}/9_CHOIR_integration"

OUTPUT_DIR="${PREPROCESS_DIR}/10_QC_visualization"
PLOT_DIR="${OUTPUT_DIR}/plots"
SAMPLE_PLOT_DIR="${PLOT_DIR}/per_sample"
CHOIR_PLOT_DIR="${PLOT_DIR}/CHOIR"
SAMPLE_DATA_DIR="${OUTPUT_DIR}/per_sample"
SCRIPT_DIR="${BASE_DIR}/Scripts/R_scripts"
LOG_DIR="${BASE_DIR}/logs/qc_viz"
README_FILE="${OUTPUT_DIR}/README.txt"

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

error_exit() {
    log_msg "ERROR: $1"
    exit 1
}

get_col_index() {
    local file=$1
    local col_name=$2
    head -1 "$file" | tr ',' '\n' | grep -n "^${col_name}$" | cut -d: -f1
}

get_unique_samples() {
    local col=$(get_col_index "$SAMPLESHEET" "sample_id")
    tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | cut -d',' -f"$col" | sort -u
}

# ==============================================================================
# INITIALIZE
# ==============================================================================

log_msg "============================================================================"
log_msg "QC Visualization Pipeline - Step 10"
log_msg "============================================================================"
log_msg "Job ID:        ${SLURM_JOB_ID:-local}"
log_msg "Date:          $(date)"
log_msg "Host:          $(hostname)"
log_msg "Base Dir:      ${BASE_DIR}"
log_msg "============================================================================"
echo ""

if [[ ! -f "$SAMPLESHEET" ]]; then
    error_exit "Samplesheet not found: $SAMPLESHEET"
fi

log_msg "Samplesheet: $SAMPLESHEET"

SAMPLES=($(get_unique_samples))
log_msg "Samples to process: ${SAMPLES[*]}"
echo ""

mkdir -p "${OUTPUT_DIR}" "${PLOT_DIR}" "${SAMPLE_PLOT_DIR}" "${CHOIR_PLOT_DIR}" "${SAMPLE_DATA_DIR}"
mkdir -p "${SCRIPT_DIR}" "${LOG_DIR}"

# ==============================================================================
# Initialize README
# ==============================================================================
cat > "${README_FILE}" << EOF
================================================================================
STEP 10: QC Visualization
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')
Pipeline: Universal Modular scRNA-seq Pipeline

SAMPLESHEET: ${SAMPLESHEET}
SAMPLES: ${SAMPLES[*]}

INPUT SOURCES:
  - 2_DropletQC_output/<SAMPLE>/*_dropletqc_summary.csv
  - 3_qClus_empty_droplets/<SAMPLE>/*_qclus_dropletqc_status.csv
  - 5_Seurat_conversion/<SAMPLE>/*_qClus_dropletqc_vaeda.rds
  - 6_Doublet_consensus/<SAMPLE>/*_doublets_removed.rds
  - 7_scCDC_correction/<SAMPLE>/*_scCDC_corrected.rds
  - 8_DecontX_correction/<SAMPLE>/*_decontX_corrected.rds
  - 9_CHOIR_integration/**/**.rds

OUTPUT:
  - 10_QC_visualization/
    - combined_filtering_summary.csv
    - filtering_flow_summary.csv
    - plots/
================================================================================
EOF

# ==============================================================================
# Environment Setup
# ==============================================================================

log_msg "Setting up environment..."

source /scicore/home/doetsch/kaiser0001/miniforge3/etc/profile.d/conda.sh
conda activate /scicore/home/doetsch/kaiser0001/miniforge3/envs/R_4_5

export LANG=C.UTF-8
export LC_ALL=C.UTF-8

log_msg "R: $(which R)"
R --version | head -n 1
echo ""

# ==============================================================================
# Write R Script
# ==============================================================================

R_SCRIPT="${SCRIPT_DIR}/run_QC_visualization.R"
SAMPLES_CSV=$(IFS=','; echo "${SAMPLES[*]}")

cat > "${R_SCRIPT}" << 'RSCRIPT_EOF'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(viridis)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: run_QC_visualization.R <base_dir> <output_dir> <samples_csv>\n")

BASE_DIR   <- args[[1]]
OUTPUT_DIR <- args[[2]]
SAMPLES    <- strsplit(args[[3]], ",")[[1]]

cat("\nSamples:", paste(SAMPLES, collapse=", "), "\n\n")

PLOT_DIR        <- file.path(OUTPUT_DIR, "plots")
SAMPLE_PLOT_DIR <- file.path(PLOT_DIR, "per_sample")
CHOIR_PLOT_DIR  <- file.path(PLOT_DIR, "CHOIR")
SAMPLE_DATA_DIR <- file.path(OUTPUT_DIR, "per_sample")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SAMPLE_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHOIR_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SAMPLE_DATA_DIR, recursive = TRUE, showWarnings = FALSE)

DROPLETQC_DIR  <- file.path(BASE_DIR, "2_DropletQC_output")
QCLUS_DIR      <- file.path(BASE_DIR, "3_qClus_empty_droplets")
VAEDA_DIR      <- file.path(BASE_DIR, "4_Vaeda_doublet_detection")
SEURAT_PRE_DIR <- file.path(BASE_DIR, "5_Seurat_conversion")
DOUBLET_DIR    <- file.path(BASE_DIR, "6_Doublet_consensus")
SCCDC_DIR      <- file.path(BASE_DIR, "7_scCDC_correction")
DECONTX_DIR    <- file.path(BASE_DIR, "8_DecontX_correction")
CHOIR_DIR      <- file.path(BASE_DIR, "9_CHOIR_integration")

cat("================================================================================\n")
cat("QC Visualization Pipeline\n")
cat("================================================================================\n\n")

first_existing <- function(paths) {
  paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  ok <- file.exists(paths)
  if (any(ok)) return(paths[which(ok)[1]])
  NA_character_
}

read_csv_safe <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}

ensure_umap <- function(seu, dims = 1:30, seed = 1) {
  rednames <- Reductions(seu)
  if ("umap" %in% rednames || "UMAP" %in% rednames) return(seu)
  set.seed(seed)
  DefaultAssay(seu) <- "RNA"
  if (!("pca" %in% rednames)) {
    seu <- NormalizeData(seu, verbose = FALSE)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    seu <- ScaleData(seu, verbose = FALSE)
    seu <- RunPCA(seu, npcs = max(dims), verbose = FALSE)
  }
  seu <- FindNeighbors(seu, dims = dims, verbose = FALSE)
  if (!("seurat_clusters" %in% colnames(seu[[]]))) {
    seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
  }
  seu <- RunUMAP(seu, dims = dims, reduction = "pca", verbose = FALSE)
  seu
}

pick_umap_name <- function(seu) {
  if ("umap" %in% Reductions(seu)) return("umap")
  if ("UMAP" %in% Reductions(seu)) return("UMAP")
  NA_character_
}

ensure_meta_cols <- function(seu, cols) {
  md <- seu[[]]
  missing <- setdiff(cols, colnames(md))
  if (length(missing) > 0) {
    for (m in missing) md[[m]] <- NA
    seu <- AddMetaData(seu, metadata = md)
  }
  seu
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

plot_group_if_present <- function(seu, red, group_by, title, out_png, width=12, height=10, pt=0.2) {
  md <- colnames(seu[[]])
  if (!(group_by %in% md)) return(FALSE)
  p <- DimPlot(seu, reduction = red, group.by = group_by, pt.size = pt, raster = FALSE) +
    ggtitle(title) + theme(legend.position = "right")
  ggsave(out_png, p, width = width, height = height, dpi = 150)
  TRUE
}

all_stats <- tibble()

for (sample in SAMPLES) {
  cat(paste0("Processing: ", sample, "\n"))

  dropletqc_summary <- first_existing(c(
    file.path(DROPLETQC_DIR, sample, paste0(sample, "_dropletqc_summary.csv")),
    file.path(DROPLETQC_DIR, sample, paste0(sample, "_dropletqc_summary.csv.gz"))
  ))
  dqc <- read_csv_safe(dropletqc_summary)

  qclus_status <- first_existing(c(
    file.path(QCLUS_DIR, sample, paste0(sample, "_qclus_dropletqc_status.csv")),
    file.path(QCLUS_DIR, sample, paste0(sample, "_qclus_dropletqc_status.csv.gz"))
  ))
  qcs <- read_csv_safe(qclus_status)

  pre_rds <- first_existing(c(
    file.path(SEURAT_PRE_DIR, sample, paste0(sample, "_qClus_dropletqc_vaeda.rds")),
    file.path(SEURAT_PRE_DIR, sample, paste0(sample, "_qclus_dropletqc_vaeda.rds"))
  ))
  post_rds <- first_existing(c(
    file.path(DOUBLET_DIR, sample, paste0(sample, "_doublets_removed.rds"))
  ))

  decontx_summary <- first_existing(c(
    file.path(DECONTX_DIR, sample, paste0(sample, "_decontX_summary.csv")),
    file.path(DECONTX_DIR, sample, paste0(sample, "_decontX_summary.csv.gz"))
  ))
  decontx <- read_csv_safe(decontx_summary)

  dbl_sum <- first_existing(c(
    file.path(DOUBLET_DIR, sample, paste0(sample, "_doublet_summary.csv")),
    file.path(DOUBLET_DIR, sample, paste0(sample, "_doublet_summary.csv.gz"))
  ))
  dbl <- read_csv_safe(dbl_sum)

  st <- list(sample = sample)

  if (!is.null(dqc)) {
    colmap <- list(
      total   = c("total_barcodes", "Total_BC", "total", "Total"),
      empty   = c("empty_droplets", "Empty", "empty", "Empty_BC"),
      damaged = c("damaged_cells", "Damaged", "damaged"),
      intact  = c("intact_cells", "Intact", "intact")
    )
    get1 <- function(df, keys) {
      k <- keys[keys %in% colnames(df)]
      if (length(k) == 0) return(NA_real_)
      as.numeric(df[[k[1]]][1])
    }
    st$cellranger_total  <- get1(dqc, colmap$total)
    st$dropletqc_empty   <- get1(dqc, colmap$empty)
    st$dropletqc_damaged <- get1(dqc, colmap$damaged)
    st$dropletqc_intact  <- get1(dqc, colmap$intact)
  } else if (!is.null(qcs)) {
    st$cellranger_total  <- nrow(qcs)
    st$dropletqc_empty   <- sum(qcs$dropletqc_empty == TRUE, na.rm = TRUE)
    st$dropletqc_damaged <- sum(qcs$dropletqc_damaged == TRUE, na.rm = TRUE)
    st$dropletqc_intact  <- sum(qcs$dropletqc_intact == TRUE, na.rm = TRUE)
  } else {
    st$cellranger_total  <- NA_real_
    st$dropletqc_empty   <- NA_real_
    st$dropletqc_damaged <- NA_real_
    st$dropletqc_intact  <- NA_real_
  }

  if (!is.null(qcs)) {
    in_input <- if ("in_qclus_input" %in% colnames(qcs)) qcs$in_qclus_input else rep(NA, nrow(qcs))
    in_kept  <- if ("in_qclus_kept"  %in% colnames(qcs)) qcs$in_qclus_kept  else rep(NA, nrow(qcs))
    st$qclus_input    <- sum(in_input == TRUE, na.rm = TRUE)
    st$qclus_kept     <- sum(in_kept  == TRUE, na.rm = TRUE)
    st$qclus_filtered <- sum(in_input == TRUE & in_kept == FALSE, na.rm = TRUE)
    st$pass_combined  <- if ("pass_combined" %in% colnames(qcs)) sum(qcs$pass_combined == TRUE, na.rm = TRUE) else NA_real_
  } else {
    st$qclus_input    <- NA_real_
    st$qclus_kept     <- NA_real_
    st$qclus_filtered <- NA_real_
    st$pass_combined  <- NA_real_
  }

  vaeda_summary <- first_existing(c(
    file.path(VAEDA_DIR, sample, paste0(sample, "_vaeda_summary.csv")),
    file.path(VAEDA_DIR, sample, paste0(sample, "_vaeda_summary.csv.gz"))
  ))
  vaeda <- read_csv_safe(vaeda_summary)
  st$vaeda_doublets <- if (!is.null(vaeda) && "n_doublets" %in% colnames(vaeda)) as.numeric(vaeda$n_doublets[1]) else NA_real_
  st$vaeda_singlets <- if (!is.null(vaeda) && "n_singlets" %in% colnames(vaeda)) as.numeric(vaeda$n_singlets[1]) else NA_real_

  if (!is.null(dbl)) {
    pick <- function(df, keys) {
      k <- keys[keys %in% colnames(df)]
      if (length(k) == 0) return(NA_real_)
      as.numeric(df[[k[1]]][1])
    }
    st$scDblFinder_doublets   <- pick(dbl, c("scDblFinder_doublets","scDblFinder_n_doublets","n_scDblFinder_doublets"))
    st$DoubletFinder_doublets <- pick(dbl, c("DoubletFinder_doublets","DoubletFinder_n_doublets","n_DoubletFinder_doublets"))
    st$consensus_doublets     <- pick(dbl, c("consensus_doublets","n_consensus_doublets"))
    st$final_cells_reported   <- pick(dbl, c("final_cells","n_final_cells"))
    st$input_cells_reported   <- pick(dbl, c("input_cells","n_input_cells"))
  } else {
    st$scDblFinder_doublets   <- NA_real_
    st$DoubletFinder_doublets <- NA_real_
    st$consensus_doublets     <- NA_real_
    st$final_cells_reported   <- NA_real_
    st$input_cells_reported   <- NA_real_
  }

  if (!is.null(decontx)) {
    get_metric <- function(df, metric_name) {
      idx <- which(df$metric == metric_name)
      if (length(idx) > 0) return(as.numeric(df$value[idx[1]]))
      NA_real_
    }
    st$decontx_mean_contamination <- get_metric(decontx, "mean_contamination")
    st$decontx_median_contamination <- get_metric(decontx, "median_contamination")
    st$decontx_umi_reduction_pct <- get_metric(decontx, "umi_reduction_pct")
  } else {
    st$decontx_mean_contamination <- NA_real_
    st$decontx_median_contamination <- NA_real_
    st$decontx_umi_reduction_pct <- NA_real_
  }

  st$pre_doublet_cells  <- if (!is.na(pre_rds)  && file.exists(pre_rds))  ncol(readRDS(pre_rds))  else NA_real_
  st$post_singlet_cells <- if (!is.na(post_rds) && file.exists(post_rds)) ncol(readRDS(post_rds)) else NA_real_
  st$doublets_removed_inferred <- if (!is.na(st$pre_doublet_cells) && !is.na(st$post_singlet_cells)) {
    st$pre_doublet_cells - st$post_singlet_cells
  } else {
    NA_real_
  }

  all_stats <- bind_rows(all_stats, as_tibble(st))
}

stats_file <- file.path(OUTPUT_DIR, "combined_filtering_summary.csv")
write.csv(all_stats, stats_file, row.names = FALSE)
cat(paste0("\nSaved: ", stats_file, "\n\n"))

if (nrow(all_stats) > 0) {
  tot <- summarise(all_stats,
                   cellranger_total    = sum(cellranger_total, na.rm = TRUE),
                   dropletqc_empty     = sum(dropletqc_empty,  na.rm = TRUE),
                   dropletqc_damaged   = sum(dropletqc_damaged,na.rm = TRUE),
                   qclus_filtered      = sum(qclus_filtered,   na.rm = TRUE),
                   pre_doublet_cells   = sum(pre_doublet_cells,na.rm = TRUE),
                   post_singlet_cells  = sum(post_singlet_cells,na.rm = TRUE),
                   doublets_removed    = sum(doublets_removed_inferred, na.rm = TRUE))

  flow_summary <- tibble(
    Step = c("Cell Ranger Input","DropletQC Empty","DropletQC Damaged","QClus Filtered",
             "Pre-doublet Cells","Doublets Removed","Final Singlets"),
    Cells = c(tot$cellranger_total,tot$dropletqc_empty,tot$dropletqc_damaged,
              tot$qclus_filtered,tot$pre_doublet_cells,tot$doublets_removed,tot$post_singlet_cells),
    Type = c("Input","Removed","Removed","Removed","Checkpoint","Removed","Final")
  )
  write.csv(flow_summary, file.path(OUTPUT_DIR, "filtering_flow_summary.csv"), row.names = FALSE)

  plot_data <- all_stats %>%
    select(sample, cellranger_total, dropletqc_empty, dropletqc_damaged,
           qclus_filtered, pre_doublet_cells, doublets_removed_inferred, post_singlet_cells) %>%
    pivot_longer(cols = -sample, names_to = "Metric", values_to = "Cells") %>%
    mutate(Metric = factor(Metric, levels = c("cellranger_total","dropletqc_empty","dropletqc_damaged",
                                             "qclus_filtered","pre_doublet_cells","doublets_removed_inferred","post_singlet_cells")))

  metric_labels <- c(cellranger_total="CR Input",dropletqc_empty="Empty",dropletqc_damaged="Damaged",
                     qclus_filtered="QClus",pre_doublet_cells="Pre-doublet",
                     doublets_removed_inferred="Doublets",post_singlet_cells="Final")

  p1 <- ggplot(plot_data, aes(x = Metric, y = Cells, fill = sample)) +
    geom_col(position = "dodge") +
    scale_x_discrete(labels = metric_labels) +
    labs(title = "Cells across pipeline checkpoints", x = NULL, y = "Count", fill = "Sample") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(PLOT_DIR, "barplot_cells_key_checkpoints.png"), p1, width = 13, height = 8, dpi = 150)

  retention <- all_stats %>%
    mutate(
      Overall_retention    = post_singlet_cells / cellranger_total * 100,
      DropletQC_retention  = (cellranger_total - dropletqc_empty - dropletqc_damaged) / cellranger_total * 100,
      QClus_retention      = qclus_kept / qclus_input * 100,
      Doublet_retention    = post_singlet_cells / pre_doublet_cells * 100
    ) %>%
    select(sample, DropletQC_retention, QClus_retention, Doublet_retention, Overall_retention) %>%
    pivot_longer(cols = -sample, names_to = "Step", values_to = "Retention") %>%
    mutate(Step = factor(Step,
                         levels = c("DropletQC_retention","QClus_retention","Doublet_retention","Overall_retention"),
                         labels = c("DropletQC","QClus","Doublet","Overall")))

  p2 <- ggplot(retention, aes(x = Step, y = sample, fill = Retention)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.finite(Retention), sprintf("%.1f%%", Retention), "NA")), size = 3) +
    scale_fill_viridis(option = "plasma", limits = c(0, 100), oob = squish) +
    labs(title = "Retention rates", x = NULL, y = NULL, fill = "Retention (%)") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(PLOT_DIR, "heatmap_retention_rates.png"), p2, width = 10, height = 7, dpi = 150)

  if (any(!is.na(all_stats$decontx_mean_contamination))) {
    p_decontx <- ggplot(all_stats, aes(x = sample, y = decontx_mean_contamination * 100, fill = sample)) +
      geom_col() +
      geom_text(aes(label = sprintf("%.2f%%", decontx_mean_contamination * 100)), vjust = -0.5, size = 3) +
      labs(title = "DecontX Mean Contamination", x = NULL, y = "Contamination (%)", fill = "Sample") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(PLOT_DIR, "barplot_decontx_contamination.png"), p_decontx, width = 10, height = 6, dpi = 150)
  }

  combined_plot <- (p1 / p2) + plot_annotation(title = "QC Summary", subtitle = paste0("Generated: ", Sys.time()))
  ggsave(file.path(PLOT_DIR, "combined_QC_summary.png"), combined_plot, width = 14, height = 14, dpi = 150)
}

cat("\nGenerating per-sample UMAP panels...\n")

for (sample in SAMPLES) {
  cat(paste0("  -> ", sample, "\n"))
  pre_rds <- first_existing(c(
    file.path(SEURAT_PRE_DIR, sample, paste0(sample, "_qClus_dropletqc_vaeda.rds")),
    file.path(SEURAT_PRE_DIR, sample, paste0(sample, "_qclus_dropletqc_vaeda.rds"))
  ))
  if (is.na(pre_rds) || !file.exists(pre_rds)) {
    cat("     WARNING: missing RDS, skipping.\n")
    next
  }
  seu <- readRDS(pre_rds)
  expected <- c("dropletqc_nuclear_fraction","dropletqc_status","filter_status",
                "fraction_unspliced","pct_counts_MT","scDblFinder_class","DoubletFinder_class",
                "doublet_votes","doublet_consensus")
  seu <- ensure_meta_cols(seu, expected)
  seu <- ensure_umap(seu, dims = 1:30, seed = 1)
  umap_name <- pick_umap_name(seu)
  if (is.na(umap_name)) {
    cat("     WARNING: no UMAP, skipping.\n")
    rm(seu); gc()
    next
  }
  md <- seu[[]]
  for (v in c("dropletqc_status","filter_status","scDblFinder_class","DoubletFinder_class","doublet_votes","doublet_consensus")) {
    if (v %in% colnames(md)) md[[v]] <- as.factor(md[[v]])
  }
  seu <- AddMetaData(seu, md)

  pA <- DimPlot(seu, reduction = umap_name, group.by = "dropletqc_status", pt.size = 0.25, raster = FALSE) +
    ggtitle(paste0(sample, " | DropletQC")) + theme(legend.position = "right")
  pB <- DimPlot(seu, reduction = umap_name, group.by = "filter_status", pt.size = 0.25, raster = FALSE) +
    ggtitle(paste0(sample, " | QClus")) + theme(legend.position = "right")
  pC <- DimPlot(seu, reduction = umap_name, group.by = "scDblFinder_class", pt.size = 0.25, raster = FALSE) +
    ggtitle(paste0(sample, " | scDblFinder")) + theme(legend.position = "right")
  pD <- DimPlot(seu, reduction = umap_name, group.by = "DoubletFinder_class", pt.size = 0.25, raster = FALSE) +
    ggtitle(paste0(sample, " | DoubletFinder")) + theme(legend.position = "right")
  pE <- DimPlot(seu, reduction = umap_name, group.by = "doublet_votes", pt.size = 0.25, raster = FALSE) +
    ggtitle(paste0(sample, " | Votes")) + theme(legend.position = "right")
  pF <- DimPlot(seu, reduction = umap_name, group.by = "doublet_consensus", pt.size = 0.25, raster = FALSE) +
    ggtitle(paste0(sample, " | Consensus")) + theme(legend.position = "right")

  panel_cat <- (pA | pB) / (pC | pD) / (pE | pF) +
    plot_annotation(title = paste0(sample, " - QC categorical"))
  ggsave(file.path(SAMPLE_PLOT_DIR, paste0(sample, "_UMAP_QC_categorical.png")),
         panel_cat, width = 14, height = 18, dpi = 150)

  pG <- FeaturePlot(seu, reduction = umap_name, features = "dropletqc_nuclear_fraction", pt.size = 0.25, raster = FALSE) +
    ggtitle(paste0(sample, " | nuclear_fraction"))
  pH <- FeaturePlot(seu, reduction = umap_name, features = "fraction_unspliced", pt.size = 0.25, raster = FALSE) +
    ggtitle(paste0(sample, " | unspliced"))
  pI <- FeaturePlot(seu, reduction = umap_name, features = "pct_counts_MT", pt.size = 0.25, raster = FALSE) +
    ggtitle(paste0(sample, " | MT%"))
  panel_feat <- (pG | pH | pI) + plot_annotation(title = paste0(sample, " - QC continuous"))
  ggsave(file.path(SAMPLE_PLOT_DIR, paste0(sample, "_UMAP_QC_continuous.png")),
         panel_feat, width = 18, height = 6.5, dpi = 150)

  write.csv(tibble(column = colnames(seu[[]])),
            file.path(SAMPLE_DATA_DIR, paste0(sample, "_seurat_metadata_columns.csv")), row.names = FALSE)
  rm(seu); gc()
}

cat("\nGenerating CHOIR plots...\n")

choir_integrated_files <- unique(c(
  list.files(CHOIR_DIR, pattern = "_CHOIR_integrated\\.rds$", recursive = TRUE, full.names = TRUE),
  list.files(CHOIR_DIR, pattern = "_CHOIR_integrated\\.rd$",  recursive = TRUE, full.names = TRUE)
))
choir_individual_files <- list.files(file.path(CHOIR_DIR, "Individual"),
                                     pattern = "_CHOIR\\.rds$", recursive = TRUE, full.names = TRUE)
all_choir_files <- unique(c(choir_integrated_files, choir_individual_files))
all_choir_files <- all_choir_files[file.exists(all_choir_files)]

if (length(all_choir_files) == 0) {
  cat("  No CHOIR files found.\n")
} else {
  cat("  Found ", length(all_choir_files), " CHOIR object(s)\n")
}

for (f in all_choir_files) {
  rel <- sub(paste0("^", normalizePath(CHOIR_DIR, winslash="/"), "/?"), "", normalizePath(f, winslash="/"))
  tag <- safe_filename(gsub("\\.rds$|\\.rd$", "", rel, ignore.case = TRUE))
  out_dir <- file.path(CHOIR_PLOT_DIR, tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat(paste0("  -> CHOIR: ", rel, "\n"))
  seu <- tryCatch({ readRDS(f) }, error = function(e) { cat("     ERROR\n"); NULL })
  if (is.null(seu)) next
  seu <- ensure_umap(seu, dims = 1:30, seed = 1)
  red <- pick_umap_name(seu)
  if (is.na(red)) { rm(seu); gc(); next }
  plot_group_if_present(seu, red, "seurat_clusters", paste0("CHOIR | ", rel),
    file.path(out_dir, paste0(tag, "_UMAP_by_clusters.png")))
  plot_group_if_present(seu, red, "Sample_ID", paste0("CHOIR | ", rel),
    file.path(out_dir, paste0(tag, "_UMAP_by_Sample_ID.png")))
  plot_group_if_present(seu, red, "Sex", paste0("CHOIR | ", rel),
    file.path(out_dir, paste0(tag, "_UMAP_by_Sex.png")))
  plot_group_if_present(seu, red, "Ventricle", paste0("CHOIR | ", rel),
    file.path(out_dir, paste0(tag, "_UMAP_by_Ventricle.png")))
  write.csv(tibble(column = colnames(seu[[]])), file.path(out_dir, paste0(tag, "_metadata_columns.csv")), row.names = FALSE)
  rm(seu); gc()
}

for (s in SAMPLES) {
  write.csv(all_stats %>% filter(sample == !!s),
            file.path(SAMPLE_DATA_DIR, paste0(s, "_filtering_summary.csv")), row.names = FALSE)
}

cat("\n================================================================================\n")
cat("QC Visualization Pipeline COMPLETE\n")
cat("================================================================================\n")
RSCRIPT_EOF

chmod +x "${R_SCRIPT}"
log_msg "R script written to: ${R_SCRIPT}"

# ==============================================================================
# Run R script
# ==============================================================================

log_msg "============================================================================"
log_msg "Running QC Visualization Pipeline"
log_msg "============================================================================"
echo ""

Rscript "${R_SCRIPT}" "${PREPROCESS_DIR}" "${OUTPUT_DIR}" "${SAMPLES_CSV}"

R_EXIT=$?

cat >> "${README_FILE}" << EOF

PROCESSING COMPLETE: $(date '+%Y-%m-%d %H:%M:%S')
================================================================================
EOF

echo ""
log_msg "============================================================================"
log_msg "QC Visualization Complete"
log_msg "============================================================================"
log_msg "Output: ${OUTPUT_DIR}"
echo ""
find "${OUTPUT_DIR}" -type f \( -name "*.csv" -o -name "*.png" \) | head -20
echo ""
log_msg "============================================================================"

exit $R_EXIT
