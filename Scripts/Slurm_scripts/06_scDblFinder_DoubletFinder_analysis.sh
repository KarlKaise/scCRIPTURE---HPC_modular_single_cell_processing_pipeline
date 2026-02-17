#!/bin/bash
#SBATCH --job-name=doublet
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --output=logs/doublet/doublet_%a_%A.out
#SBATCH --error=logs/doublet/doublet_%a_%A.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=karol.kaiser@unibas.ch

# ============================================================================
# STEP 6: scDblFinder + DoubletFinder Consensus Doublet Analysis
# ============================================================================
#
# UNIVERSAL MODULAR VERSION - Reads samples from samplesheet.csv
#
# Input:  5_Seurat_conversion/<SAMPLE>/<SAMPLE>_qClus_dropletqc_vaeda.rds
# Output: 6_Doublet_consensus/<SAMPLE>/<SAMPLE>_doublets_removed.rds
#
# Consensus strategy: doublet if >= 2 of 3 (VAEDA + scDblFinder + DoubletFinder)
#
# Usage:
#   sbatch --array=1-N 06_scDblFinder_DoubletFinder_analysis.sh
#
# ============================================================================

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Base directory (from environment or default)
if [[ -n "${PROJECT_ROOT:-}" ]]; then
    BASE_DIR="${PROJECT_ROOT}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    BASE_DIR="${SLURM_SUBMIT_DIR}"
else
    BASE_DIR="/scicore/home/doetsch/kaiser0001/Revision_NatureComm_Sex/Vandebroucke_fibroblast_paper"
fi

# Samplesheet
SAMPLESHEET="${SAMPLESHEET:-${BASE_DIR}/samplesheet.csv}"

# Get dataset name from samplesheet (column 12)
get_dataset_name() {
    tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f12
}

# Set PREPROCESS_DIR from environment or derive from samplesheet
if [[ -z "${PREPROCESS_DIR:-}" ]]; then
    DATASET_NAME=$(get_dataset_name)
    if [[ -n "$DATASET_NAME" ]]; then
        PREPROCESS_DIR="${BASE_DIR}/Output_dir_${DATASET_NAME}/Single_cell_preprocessed"
    else
        PREPROCESS_DIR="${BASE_DIR}"
    fi
fi

# Directories
INPUT_ROOT="${PREPROCESS_DIR}/5_Seurat_conversion"
OUTPUT_ROOT="${PREPROCESS_DIR}/6_Doublet_consensus"
SCRIPT_DIR="${BASE_DIR}/Scripts/R_scripts"
LOG_DIR="${BASE_DIR}/logs/doublet"
README_FILE="${OUTPUT_ROOT}/README.txt"

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

get_sample_count() {
    get_unique_samples | wc -l
}

get_sample_by_index() {
    local idx=$1
    get_unique_samples | sed -n "${idx}p"
}

get_sample_field() {
    local sample=$1
    local field=$2
    local col=$(get_col_index "$SAMPLESHEET" "$field")
    if [[ -n "$col" ]]; then
        grep "^${sample}," "$SAMPLESHEET" | head -1 | cut -d',' -f"$col"
    fi
}

# ==============================================================================
# INITIALIZE
# ==============================================================================

log_msg "============================================================================"
log_msg "Doublet Consensus Analysis - Modular Pipeline"
log_msg "============================================================================"
log_msg "Job ID:        ${SLURM_JOB_ID:-local}"
log_msg "Array Task:    ${SLURM_ARRAY_TASK_ID:-1}"
log_msg "Date:          $(date)"
log_msg "Host:          $(hostname)"
log_msg "Base Dir:      ${BASE_DIR}"
log_msg "============================================================================"
echo ""

# Validate samplesheet
if [[ ! -f "$SAMPLESHEET" ]]; then
    error_exit "Samplesheet not found: $SAMPLESHEET"
fi

log_msg "Samplesheet: $SAMPLESHEET"

# Get sample
N_SAMPLES=$(get_sample_count)
TASK_IDX=${SLURM_ARRAY_TASK_ID:-1}

if [[ $TASK_IDX -gt $N_SAMPLES ]]; then
    error_exit "Array task $TASK_IDX exceeds sample count $N_SAMPLES"
fi

SAMPLE=$(get_sample_by_index $TASK_IDX)

if [[ -z "$SAMPLE" ]]; then
    error_exit "Could not determine sample for task index $TASK_IDX"
fi

log_msg "Processing: $SAMPLE (task $TASK_IDX of $N_SAMPLES)"

# Get sample metadata
SAMPLE_SEX=$(get_sample_field "$SAMPLE" "sex")
SAMPLE_BATCH=$(get_sample_field "$SAMPLE" "batch")
SAMPLE_VENTRICLE=$(get_sample_field "$SAMPLE" "ventricle")

log_msg "Sample metadata: Sex=${SAMPLE_SEX:-N/A}, Batch=${SAMPLE_BATCH:-N/A}, Ventricle=${SAMPLE_VENTRICLE:-N/A}"
echo ""

# Paths
INPUT_RDS="${INPUT_ROOT}/${SAMPLE}/${SAMPLE}_qClus_dropletqc_vaeda.rds"
SAMPLE_OUTPUT_DIR="${OUTPUT_ROOT}/${SAMPLE}"
OUTPUT_RDS="${SAMPLE_OUTPUT_DIR}/${SAMPLE}_doublets_removed.rds"

START_TIME=$(date '+%H:%M:%S')

# Create directories
mkdir -p "${SAMPLE_OUTPUT_DIR}"
mkdir -p "${SCRIPT_DIR}"
mkdir -p "${LOG_DIR}"

# ==============================================================================
# INITIALIZE README (first task only)
# ==============================================================================

if [[ ${SLURM_ARRAY_TASK_ID:-1} -eq 1 ]]; then
    cat > "${README_FILE}" << EOF
================================================================================
STEP 6: scDblFinder + DoubletFinder Consensus Doublet Analysis
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')
Pipeline: Universal Modular scRNA-seq Pipeline

DESCRIPTION:
  Runs three independent doublet detection methods and applies consensus voting:
    1. VAEDA (from Step 4)
    2. scDblFinder (Bioconductor)
    3. DoubletFinder (Seurat-based)

  Consensus: Cell is doublet if >= 2 of 3 methods agree.

SAMPLESHEET: ${SAMPLESHEET}
INPUT: 5_Seurat_conversion/<SAMPLE>/<SAMPLE>_qClus_dropletqc_vaeda.rds

SAMPLES PROCESSED:
EOF
fi

# ==============================================================================
# VERIFY INPUT
# ==============================================================================

if [[ ! -f "${INPUT_RDS}" ]]; then
    error_exit "Input file not found: ${INPUT_RDS}"
fi

INPUT_SIZE=$(du -h "${INPUT_RDS}" | cut -f1)
log_msg "Input: ${INPUT_RDS} (${INPUT_SIZE})"
log_msg "Output: ${OUTPUT_RDS}"
echo ""

# ==============================================================================
# ACTIVATE R ENVIRONMENT
# ==============================================================================

log_msg "Activating R environment..."

source /scicore/home/doetsch/kaiser0001/miniforge3/etc/profile.d/conda.sh
conda activate /scicore/home/doetsch/kaiser0001/miniforge3/envs/R_4_5

log_msg "R: $(which R)"
R --version | head -n 1
echo ""

# ==============================================================================
# CREATE R SCRIPT
# ==============================================================================

R_SCRIPT="${SCRIPT_DIR}/run_doublet_consensus_${SAMPLE}.R"

cat > "${R_SCRIPT}" << 'REOF'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(DoubletFinder)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: run_doublet_consensus.R <input_rds> <output_rds> <output_dir>\n")
input_rds <- normalizePath(args[[1]], mustWork = TRUE)
output_rds <- args[[2]]
output_dir <- args[[3]]

dir.create(dirname(output_rds), recursive=TRUE, showWarnings=FALSE)

# ============================
# PARAMETERS
# ============================
MIN_UMI <- 200
MIN_GENE_CELLS <- 5
MIN_GENES_PER_CELL <- 10

NPCS <- 30
PN <- 0.25
EXPECTED_DOUBLET_RATE <- 0.075
CLUSTER_RES <- 0.5

CONSENSUS_K <- 2  # doublet if >=2 of 3 (VAEDA + scDblFinder + DoubletFinder)

# ============================
# HELPERS
# ============================
pkgver <- function(p) tryCatch(as.character(packageVersion(p)), error=function(e) "NA")

nnz <- function(m) {
  if (is.null(m)) return(NA_integer_)
  if (inherits(m, "dgCMatrix")) return(length(m@x))
  if (inherits(m, "Matrix")) return(Matrix::nnzero(m))
  suppressWarnings(sum(m != 0))
}

is_empty_mat <- function(m) {
  is.null(m) || any(dim(m) == 0) || (inherits(m, "dgCMatrix") && length(m@x) == 0)
}

get_counts <- function(seu) {
  m <- tryCatch(GetAssayData(seu, assay="RNA", layer="counts"), error=function(e) NULL)
  if (!is_empty_mat(m)) return(m)
  m <- tryCatch(GetAssayData(seu, assay="RNA", slot="counts"), error=function(e) NULL)
  if (!is_empty_mat(m)) return(m)
  stop("Counts matrix is missing/empty in RNA assay.")
}

df_get_fun <- function(candidates) {
  ns <- asNamespace("DoubletFinder")
  for (nm in candidates) {
    f <- get0(nm, envir = ns, inherits = FALSE)
    if (is.function(f)) return(list(name=nm, fun=f))
  }
  stop("None of these DoubletFinder functions exist: ", paste(candidates, collapse=", "))
}

coerce_vec_n <- function(x, n) {
  if (is.null(x)) return(rep(NA, n))
  if (is.data.frame(x)) x <- x[[1]]
  if (is.matrix(x)) x <- x[,1,drop=TRUE]
  if (is.list(x) && !is.atomic(x)) x <- unlist(x, use.names=FALSE)
  if (length(x)==1 && n>1) x <- rep(x, n)
  if (length(x) < n) x <- c(x, rep(NA, n-length(x)))
  if (length(x) > n) x <- x[seq_len(n)]
  x
}

as_is_doublet_generic <- function(x, n) {
  x <- coerce_vec_n(x, n)
  if (is.factor(x)) x <- as.character(x)
  x <- tolower(trimws(as.character(x)))
  x %in% c("doublet","doublets","dbl","true","1","yes")
}

as_is_doublet_df <- function(x, n) {
  x <- coerce_vec_n(x, n)
  if (is.factor(x)) x <- as.character(x)
  x <- tolower(trimws(as.character(x)))
  x %in% c("doublet")
}

df_paramSweep <- df_get_fun(c("paramSweep", "paramSweep_v3"))
df_summarize  <- df_get_fun(c("summarizeSweep"))
df_findpk     <- df_get_fun(c("find.pK", "findpK", "find_pK"))
df_doublet    <- df_get_fun(c("doubletFinder", "doubletFinder_v3"))

cat("\n[DoubletFinder API]\n")
cat("  paramSweep:    ", df_paramSweep$name, "\n", sep="")
cat("  summarizeSweep:", df_summarize$name, "\n", sep="")
cat("  find.pK:       ", df_findpk$name, "\n", sep="")
cat("  doubletFinder: ", df_doublet$name, "\n", sep="")
cat("  DoubletFinder version: ", pkgver("DoubletFinder"), "\n", sep="")

# ============================
# MAIN ANALYSIS
# ============================

cat("\n================================================================================\n")
cat("INPUT:  ", input_rds, "\n", sep="")
cat("OUTPUT: ", output_rds, "\n", sep="")
cat("================================================================================\n")

seu <- readRDS(input_rds)
DefaultAssay(seu) <- "RNA"

counts <- get_counts(seu)

# Store input cell count
n0_cells <- ncol(seu)
n0_genes <- nrow(seu)

# ---------- QC filters ----------
umi <- Matrix::colSums(counts)
genes_per_cell <- Matrix::colSums(counts > 0)
gene_cells <- Matrix::rowSums(counts > 0)

keep_cells <- (umi >= MIN_UMI) & (genes_per_cell >= MIN_GENES_PER_CELL)
keep_genes <- (gene_cells >= MIN_GENE_CELLS)

seu <- subset(seu, cells = colnames(seu)[keep_cells])
seu <- subset(seu, features = rownames(seu)[keep_genes])
n1_cells <- ncol(seu); n1_genes <- nrow(seu)

cat(sprintf("[QC] Start: cells=%d genes=%d | After QC: cells=%d genes=%d\n", n0_cells, n0_genes, n1_cells, n1_genes))

# ---------- Preprocessing ----------
seu <- NormalizeData(seu, verbose=FALSE)
seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=2000, verbose=FALSE)
seu <- ScaleData(seu, verbose=FALSE)
seu <- RunPCA(seu, npcs=NPCS, verbose=FALSE)
seu <- RunUMAP(seu, dims=1:NPCS, verbose=FALSE)
seu <- FindNeighbors(seu, dims=1:NPCS, verbose=FALSE)
seu <- FindClusters(seu, resolution=CLUSTER_RES, verbose=FALSE)

# ---------- scDblFinder ----------
sce <- as.SingleCellExperiment(seu)
samples_arg <- if ("orig.ident" %in% colnames(SummarizedExperiment::colData(sce))) "orig.ident" else NULL
sce <- scDblFinder::scDblFinder(sce, samples=samples_arg, verbose=FALSE)
cd <- as.data.frame(SummarizedExperiment::colData(sce))
seu$scDblFinder_score <- cd$scDblFinder.score
seu$scDblFinder_class <- cd$scDblFinder.class

# ---------- DoubletFinder: pK selection ----------
sweep.res.list <- df_paramSweep$fun(seu, PCs = 1:NPCS, sct = FALSE)
sweep.stats <- df_summarize$fun(sweep.res.list, GT = FALSE)
bcmvn <- df_findpk$fun(sweep.stats)

best_idx <- which.max(bcmvn$BCmetric)
best_pK <- as.numeric(as.character(bcmvn$pK[best_idx]))
cat(sprintf("[DoubletFinder] best pK=%s\n", as.character(best_pK)))

# ---------- Homotypic adjust ----------
annotations <- seu@meta.data$seurat_clusters
homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
nExp_poi <- round(EXPECTED_DOUBLET_RATE * nrow(seu@meta.data))
nExp_adj <- round(nExp_poi * (1 - homotypic.prop))

cat(sprintf("[DoubletFinder] n_cells=%d expected=%d homotypic=%.3f adjusted=%d\n",
            ncol(seu), nExp_poi, homotypic.prop, nExp_adj))

# Run DoubletFinder
seu <- df_doublet$fun(
  seu, PCs = 1:NPCS, pN = PN, pK = best_pK, nExp = nExp_adj,
  reuse.pANN = NULL, sct = FALSE
)

# Capture newest DF.classifications column
md_cols <- colnames(seu[[]])
df_class_cols <- grep("^DF\\.classifications", md_cols, value=TRUE)
if (length(df_class_cols) == 0) stop("DoubletFinder classification column not found.")
newest_class <- df_class_cols[length(df_class_cols)]
seu$DoubletFinder_class <- seu[[newest_class]][,1]

# ---------- Consensus (VAEDA + scDblFinder + DoubletFinder) ----------
n <- ncol(seu)
vaeda_calls <- if ("vaeda_prediction" %in% colnames(seu[[]])) seu$vaeda_prediction else NULL

vaeda_is <- as_is_doublet_generic(vaeda_calls, n)
scdbl_is <- as_is_doublet_generic(seu$scDblFinder_class, n)
df_is    <- as_is_doublet_df(seu$DoubletFinder_class, n)

votes <- as.integer(vaeda_is) + as.integer(scdbl_is) + as.integer(df_is)

seu$doublet_vaeda <- vaeda_is
seu$doublet_scDblFinder <- scdbl_is
seu$doublet_DoubletFinder <- df_is
seu$doublet_votes <- votes
seu$doublet_consensus <- votes >= CONSENSUS_K

n_cons  <- sum(seu$doublet_consensus, na.rm=TRUE)
n_vaeda <- sum(seu$doublet_vaeda, na.rm=TRUE)
n_scdbl <- sum(seu$doublet_scDblFinder, na.rm=TRUE)
n_df    <- sum(seu$doublet_DoubletFinder, na.rm=TRUE)

cat(sprintf("\n[Doublets] VAEDA=%d scDblFinder=%d DoubletFinder=%d | Consensus(>=%d)=%d\n",
            n_vaeda, n_scdbl, n_df, CONSENSUS_K, n_cons))

# ---------- Remove doublets ----------
seu_filt <- subset(seu, cells = colnames(seu)[!seu$doublet_consensus])
saveRDS(seu_filt, output_rds)

# ---------- Save summary ----------
summary_df <- data.frame(
  sample = basename(dirname(output_rds)),
  input_cells = n0_cells,
  after_qc = n1_cells,
  vaeda_doublets = n_vaeda,
  scDblFinder_doublets = n_scdbl,
  DoubletFinder_doublets = n_df,
  consensus_doublets = n_cons,
  final_cells = ncol(seu_filt)
)
summary_file <- file.path(output_dir, paste0(basename(dirname(output_rds)), "_doublet_summary.csv"))
write.csv(summary_df, summary_file, row.names = FALSE)

cat(sprintf("\n[Summary] Input: %d cells | After QC: %d | Final: %d cells\n",
            n0_cells, n1_cells, ncol(seu_filt)))

cat("\nDONE\n")
REOF

chmod +x "${R_SCRIPT}"

# ==============================================================================
# RUN DOUBLET ANALYSIS
# ==============================================================================

log_msg "============================================================================"
log_msg "Running Doublet Consensus Analysis"
log_msg "============================================================================"
echo ""

Rscript "${R_SCRIPT}" "${INPUT_RDS}" "${OUTPUT_RDS}" "${SAMPLE_OUTPUT_DIR}"

R_EXIT=$?
END_TIME=$(date '+%H:%M:%S')

# ==============================================================================
# UPDATE README
# ==============================================================================

update_readme() {
    local sample=$1
    local status=$2

    (
        flock -x 200
        if [[ -f "${SAMPLE_OUTPUT_DIR}/${sample}_doublet_summary.csv" ]]; then
            INPUT_CELLS=$(awk -F',' 'NR==2 {print $2}' "${SAMPLE_OUTPUT_DIR}/${sample}_doublet_summary.csv")
            FINAL_CELLS=$(awk -F',' 'NR==2 {print $8}' "${SAMPLE_OUTPUT_DIR}/${sample}_doublet_summary.csv")
            CONSENSUS=$(awk -F',' 'NR==2 {print $7}' "${SAMPLE_OUTPUT_DIR}/${sample}_doublet_summary.csv")
            echo "  ${sample}: ${status} | Input: ${INPUT_CELLS}, Final: ${FINAL_CELLS}, Doublets: ${CONSENSUS}" >> "${README_FILE}"
        else
            echo "  ${sample}: ${status}" >> "${README_FILE}"
        fi
    ) 200>"${README_FILE}.lock"
}

if [[ $R_EXIT -eq 0 ]]; then
    update_readme "${SAMPLE}" "SUCCESS"
    log_msg "Status: SUCCESS"
else
    update_readme "${SAMPLE}" "FAILED"
    error_exit "Doublet analysis failed"
fi

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

echo ""
log_msg "============================================================================"
log_msg "Doublet Analysis Complete - ${SAMPLE}"
log_msg "============================================================================"
log_msg "Output: ${OUTPUT_RDS}"
echo ""
ls -lh "${SAMPLE_OUTPUT_DIR}/"
echo ""
log_msg "============================================================================"
