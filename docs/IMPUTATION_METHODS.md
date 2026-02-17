# Imputation Methods in the scRNA-seq Pipeline

This document describes the two imputation strategies available in the pipeline — **afMF** (Module 02b) and **ALRA** (Module 03b) — along with the complementary **CLTS re-normalisation** (Module 07b). It covers algorithm design, placement in the pipeline, configuration, compatibility constraints, and practical recommendations.

---

## Table of Contents

- [Background: The Dropout Problem](#background-the-dropout-problem)
- [Method Overview](#method-overview)
- [Module 02b: afMF Imputation (Counts-Based)](#module-02b-afmf-imputation-counts-based)
  - [Algorithm](#afmf-algorithm)
  - [Pipeline Placement](#afmf-pipeline-placement)
  - [Input / Output](#afmf-input--output)
  - [Parameters](#afmf-parameters)
  - [Environment Requirements](#afmf-environment-requirements)
- [Module 03b: ALRA Imputation (Normalised-Data)](#module-03b-alra-imputation-normalised-data)
  - [Algorithm](#alra-algorithm)
  - [Pipeline Placement](#alra-pipeline-placement)
  - [Input / Output](#alra-input--output)
  - [Parameters](#alra-parameters)
  - [Normalisation Compatibility](#alra-normalisation-compatibility)
- [Module 07b: CLTS Re-normalisation](#module-07b-clts-re-normalisation)
- [Unified Configuration](#unified-configuration)
  - [The `imputation_method` Parameter](#the-imputation_method-parameter)
  - [Execution Matrix](#execution-matrix)
  - [Full params.R Reference](#full-paramsr-reference)
- [Recommended Settings](#recommended-settings)
  - [For Differential Expression Analysis](#for-differential-expression-analysis)
  - [For Visualisation and Exploration](#for-visualisation-and-exploration)
  - [For Trajectory Inference](#for-trajectory-inference)
- [Method Comparison](#method-comparison)
- [Downstream Compatibility](#downstream-compatibility)
- [Troubleshooting](#troubleshooting)
- [References](#references)

---

## Background: The Dropout Problem

Single-cell RNA-seq data contain a high proportion of zero values — often 80–95% of the expression matrix. These zeros arise from two distinct sources: **biological zeros** (genes genuinely not expressed in a given cell) and **technical zeros** or "dropouts" (genes that are expressed but fail to be captured due to inefficient mRNA capture, reverse transcription, or amplification). Dropouts introduce sparsity that can bias downstream analyses including clustering, differential expression, and trajectory inference.

Imputation methods attempt to recover the expression signal lost to technical dropout while preserving true biological zeros. The two methods implemented in this pipeline — afMF and ALRA — are both matrix-theory-based approaches that exploit the low-rank structure of expression matrices. They were chosen because benchmarking studies have demonstrated that matrix-factorisation methods provide stable, consistent improvements across diverse downstream tasks without the overfitting and data distortion risks associated with more complex deep-learning or smoothing-based alternatives (Huang et al. 2025, Clin Transl Med).

Critically, both methods are provided as **optional complements** to the standard analysis. By default, imputed data are stored in separate assays (`"imputed"` for afMF, `"ALRA"` for ALRA) and are **not** used for downstream differential expression, which operates on raw or normalised counts.

---

## Method Overview

| Property | afMF (Module 02b) | ALRA (Module 03b) |
|----------|-------------------|-------------------|
| Full name | Adaptive low-rank Full Matrix Factorisation | Adaptively-thresholded Low-Rank Approximation |
| Category | Matrix-theory (iterative gradient descent) | Matrix-theory (randomised SVD + thresholding) |
| Operates on | Raw counts (before normalisation) | Log-normalised data (after normalisation) |
| Pipeline position | Between Module 02 (QC) and Module 03 (Normalisation) | Between Module 03 (Normalisation) and Module 04 (Integration) |
| Output assay | `"imputed"` | `"ALRA"` |
| SCTransform compatible | Yes | **No** (automatically skipped) |
| Language | Python | R |
| Key reference | Huang et al. 2025 | Linderman et al. 2022 |

---

## Module 02b: afMF Imputation (Counts-Based)

### afMF Algorithm

afMF (adaptive low-rank Full Matrix Factorisation) is an improved matrix-theory-based algorithm that builds upon ALRA. While ALRA employs randomised SVD for a single-pass low-rank approximation, afMF uses an **iterative gradient descent process** to optimise two low-rank factor matrices (W and H, where the expression matrix X ≈ W × H). This iterative refinement allows afMF to capture more signal structure than a single SVD pass.

The algorithm works as follows:

1. **Gene filtering**: Remove genes expressed in fewer than `min_cells_expressing` cells to reduce noise dimensions.
2. **Rank estimation**: Determine the intrinsic rank of the data from the singular value spectrum.
3. **Iterative factorisation**: Decompose the counts matrix into two low-rank matrices via gradient descent, iterating until convergence (controlled by `max_iter` and `tol`).
4. **Reconstruction**: Multiply the factor matrices to produce the imputed count matrix.
5. **Zero preservation**: Values corresponding to genes not expected to be expressed in a cell type are handled conservatively to avoid false positives.

afMF operates on **raw counts** (before any normalisation), which is conceptually appropriate because the dropout process occurs at the count level. Benchmarking has shown that afMF ranks among the top methods for cell type annotation, trajectory inference (DPT), and AUCell/SCENIC regulatory analysis, with improvements in differential expression when using MAST or Wilcoxon rank-sum tests (Huang et al. 2025).

### afMF Pipeline Placement

```
Module 02 (QC Validation)
    │
    ▼
Module 02b (afMF Imputation)  ◄── operates on merged_object.rds (raw counts)
    │
    ▼
Module 03 (Normalisation)     ◄── uses ORIGINAL counts by default
```

afMF runs **after** QC filtering but **before** normalisation. By default, the downstream normalisation in Module 03 uses the **original** (non-imputed) counts. Set `use_afmf_for_normalization = TRUE` to normalise the imputed counts instead (not recommended for DE analysis).

### afMF Input / Output

**Input**: `objects/merged_object.rds` — the merged, QC-filtered Seurat object from Module 02 containing raw counts in the `RNA` assay.

**Output**: `objects/merged_object_imputed.rds` — the same Seurat object with an additional `"imputed"` assay containing the afMF-imputed count matrix.

**Additional outputs** (in `02b_Imputation_afMF/`):
- `tables/afmf_imputation_summary.csv` — summary statistics (genes imputed, sparsity before/after, convergence metrics)
- `plots/afmf_sparsity_comparison.png` — per-sample sparsity before and after imputation
- `plots/afmf_gene_recovery.png` — distribution of recovered expression values

### afMF Parameters

```r
# In params.R
afmf_python = file.path(Sys.getenv("HOME"),
                        ".conda/envs/afMF_SCImputation_env/bin/python"),
afmf_max_iter = 100,              # Maximum optimisation iterations
afmf_tol = 1e-5,                  # Convergence tolerance (stop when ΔL < tol)
afmf_min_cells_expressing = 10,   # Minimum cells expressing a gene to retain it
use_afmf_for_normalization = FALSE # Use imputed counts in Module 03?
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `afmf_python` | `~/.conda/envs/afMF_SCImputation_env/bin/python` | Path to the Python interpreter in the afMF conda environment |
| `afmf_max_iter` | `100` | Maximum number of gradient descent iterations. Increase for very large or complex datasets |
| `afmf_tol` | `1e-5` | Convergence tolerance. The algorithm stops when the change in loss falls below this threshold |
| `afmf_min_cells_expressing` | `10` | Genes expressed in fewer cells than this are excluded from imputation |
| `use_afmf_for_normalization` | `FALSE` | If `TRUE`, Module 03 normalises the imputed counts instead of the original counts. **Not recommended for DE** |

### afMF Environment Requirements

afMF requires a dedicated Python environment with specific dependencies:

```bash
# Create the afMF environment
conda create -n afMF_SCImputation_env python=3.11
conda activate afMF_SCImputation_env
pip install numpy scipy scikit-learn anndata h5py
```

Module 00 (Environment Setup) automatically verifies the availability of the afMF Python environment and sets the internal flag `has_afmf`. If the environment is not found, Module 02b is skipped with a warning.

---

## Module 03b: ALRA Imputation (Normalised-Data)

### ALRA Algorithm

ALRA (Adaptively-thresholded Low-Rank Approximation) exploits the fact that the true expression matrix is both **non-negative** and **low-rank** — it can be well-approximated by a matrix of much lower rank than its full dimensions. The key innovation is a thresholding step that preserves biological zeros after the low-rank reconstruction.

The algorithm proceeds in three steps:

1. **Rank-k approximation via randomised SVD**: Compute the rank-k approximation of the log-normalised expression matrix using randomised singular value decomposition. The rank k is automatically determined by identifying where the singular value spectrum transitions from signal to noise — specifically, k is chosen as the first index where the gap between successive singular values exceeds μ + 6σ (where μ and σ are the mean and standard deviation of the tail gaps).

2. **Adaptive thresholding to preserve biological zeros**: After SVD reconstruction, every zero in the original matrix is replaced by a non-zero value. However, elements corresponding to biological zeros are approximately symmetrically distributed around zero (because the true value is zero plus noise). For each gene, ALRA thresholds at the magnitude of the most negative value (or a specified quantile), setting all values below this threshold back to zero. This preserves biological zeros while retaining imputed technical zeros.

3. **Rescaling**: The imputed matrix is rescaled to restore the original library-size normalisation properties.

ALRA is conservative by design — it preferentially preserves biological zeros at the cost of potentially leaving some technical zeros unimputed. This is preferable for downstream analyses where false positive expression is more problematic than missing data.

### ALRA Pipeline Placement

```
Module 03 (Normalisation)
    │
    ▼
Module 03b (ALRA Imputation)  ◄── operates on log-normalised data
    │
    ▼
Module 04 (Integration)       ◄── uses ORIGINAL normalised data by default
```

ALRA runs **after** normalisation because it requires log-normalised expression values — the symmetry argument for biological zeros depends on the log transformation. By default, Module 04 proceeds with the **original** normalised data. Set `use_alra_for_downstream = TRUE` to use ALRA-imputed data for integration (not recommended for DE).

### ALRA Input / Output

**Input**: The normalised Seurat object from Module 03 (e.g., `merged_normalized/merged_LogNormalize_unintegrated.rds`).

**Output**: `merged_normalized/merged_LogNormalize_alra.rds` (or equivalent for scran) — the Seurat object with an additional `"ALRA"` assay containing the imputed log-normalised expression values.

**Additional outputs** (in `03b_Imputation_ALRA/`):
- `tables/alra_imputation_summary.csv` — summary (chosen rank k, genes with recovered expression, sparsity metrics)
- `plots/alra_singular_values.png` — singular value spectrum with the chosen k marked
- `plots/alra_sparsity_comparison.png` — sparsity before and after imputation

### ALRA Parameters

```r
# In params.R
alra_k = NULL,                    # SVD rank (NULL = auto-detect from singular value gaps)
alra_q = 10,                      # Number of power iterations for randomised SVD
alra_quantile_prob = 0.001,       # Quantile for thresholding (controls zero-preservation stringency)
alra_compatible_methods = c("LogNormalize", "scran"),  # Normalisation methods compatible with ALRA
use_alra_for_downstream = FALSE   # Use ALRA-imputed data in Module 04+?
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alra_k` | `NULL` (auto) | Rank of the SVD approximation. `NULL` uses the automatic detection based on singular value gaps. Manual specification useful if the auto-detection picks a suboptimal rank |
| `alra_q` | `10` | Power iterations for the randomised SVD, improving accuracy of the low-rank approximation. Higher values are more accurate but slower |
| `alra_quantile_prob` | `0.001` | Quantile of the negative value distribution used for thresholding. Smaller values preserve more zeros (more conservative); larger values impute more aggressively |
| `alra_compatible_methods` | `c("LogNormalize", "scran")` | Which normalisation methods produce data suitable for ALRA. **Do not add SCTransform** |
| `use_alra_for_downstream` | `FALSE` | If `TRUE`, downstream modules use the ALRA-imputed assay. **Not recommended for DE** |

### ALRA Normalisation Compatibility

ALRA's thresholding logic requires **log-normalised** expression values where zero-valued elements in the original matrix map to symmetrically-distributed values around zero after SVD reconstruction. This assumption holds for:

| Normalisation Method | Compatible | Reason |
|---------------------|-----------|--------|
| **LogNormalise** | ✅ Yes | Produces log-normalised values; zeros map to zero in log-space |
| **scran** | ✅ Yes | Produces log-normalised values via deconvolution size factors |
| **SCTransform** | ❌ No | Produces Pearson residuals (centred around zero for all genes), not log-normalised counts. The symmetry assumption for biological zeros does not hold |

**Automatic skip logic**: If Module 03 selects SCTransform as the best normalisation method (via benchmarking), Module 03b is automatically skipped with the following warning:

```
⚠ ALRA imputation skipped: SCTransform produces Pearson residuals,
  not log-normalised data. ALRA requires log-normalised input.
  Set imputation_method = "afmf" to use counts-based imputation only.
```

This applies even when `imputation_method = "both"` — afMF (Module 02b) still runs on raw counts, but ALRA is skipped.

---

## Module 07b: CLTS Re-normalisation

CLTS (Cell-type Level Total Sum) re-normalisation is conceptually distinct from imputation but is grouped with the optional modules because it addresses a related bias: **cell-type composition effects** on normalisation.

Standard library-size normalisation assumes that total mRNA content is similar across cells. When cell types with very different transcriptome sizes (e.g., large secretory cells vs small immune cells) are present in the same dataset, this assumption is violated. Highly expressed genes in large cells can suppress the apparent expression of genes in small cells after normalisation, distorting differential expression results.

CLTS addresses this by re-normalising within each cell-type cluster independently, then rescaling so that between-cluster comparisons remain valid.

**Configuration**:

```r
run_clts_renormalization = FALSE,   # Set TRUE to enable
clts_clustering_source = "scice",   # Which clustering to use: "scice", "leiden", "choir"
clts_min_cells_per_cluster = 50,    # Minimum cluster size for re-normalisation
clts_baseline_sample = NULL,        # Optional reference sample
clts_run_benchmark = TRUE           # Compare pre- vs post-CLTS DE results
```

**Pipeline placement**: Module 07b runs after clustering (Modules 05–07) and before differential expression (Module 08). Its output is an additional Seurat object (`<source>_clts_renormalized_object.rds`) that Module 08 can use alongside or instead of the standard object.

---

## Unified Configuration

### The `imputation_method` Parameter

A single parameter in `params.R` controls both imputation modules:

```r
imputation_method = "none"   # Options: "none", "afmf", "alra", "both"
```

| Value | Module 02b (afMF) | Module 03b (ALRA) |
|-------|-------------------|-------------------|
| `"none"` | Skip | Skip |
| `"afmf"` | **Run** | Skip |
| `"alra"` | Skip | **Run** (if normalisation is compatible) |
| `"both"` | **Run** | **Run** (if normalisation is compatible) |

### Execution Matrix

The actual execution depends on both `imputation_method` and the normalisation method selected by Module 03:

| `imputation_method` | Normalisation = LogNorm | Normalisation = scran | Normalisation = SCTransform |
|---------------------|------------------------|----------------------|-----------------------------|
| `"none"` | — | — | — |
| `"afmf"` | afMF only | afMF only | afMF only |
| `"alra"` | ALRA only | ALRA only | **Skipped** (warning) |
| `"both"` | afMF + ALRA | afMF + ALRA | **afMF only** (ALRA skipped with warning) |

### Full params.R Reference

```r
# ==============================================================================
# IMPUTATION SETTINGS (UNIFIED)
# ==============================================================================

# Main control — determines which imputation modules run
imputation_method = "both",           # "none", "afmf", "alra", "both"

# ==============================================================================
# MODULE 02b: afMF IMPUTATION (COUNTS-BASED)
# ==============================================================================
use_afmf_for_normalization = FALSE,   # Use imputed counts for Module 03?
afmf_python = file.path(Sys.getenv("HOME"),
                        ".conda/envs/afMF_SCImputation_env/bin/python"),
afmf_max_iter = 100,                  # Maximum optimisation iterations
afmf_tol = 1e-5,                      # Convergence tolerance
afmf_min_cells_expressing = 10,       # Gene filter: min cells expressing

# ==============================================================================
# MODULE 03b: ALRA IMPUTATION (NORMALISED-DATA)
# ==============================================================================
use_alra_for_downstream = FALSE,      # Use ALRA data for Module 04+?
alra_k = NULL,                        # SVD rank (NULL = auto-detect)
alra_q = 10,                          # Power iterations for rSVD
alra_quantile_prob = 0.001,           # Thresholding quantile
alra_compatible_methods = c("LogNormalize", "scran"),

# ==============================================================================
# MODULE 07b: CLTS RE-NORMALISATION
# ==============================================================================
run_clts_renormalization = FALSE,     # Set TRUE to enable
clts_clustering_source = "scice",     # Clustering source: "scice", "leiden", "choir"
clts_min_cells_per_cluster = 50,      # Min cluster size
clts_baseline_sample = NULL,          # Reference sample (optional)
clts_run_benchmark = TRUE             # Pre- vs post-CLTS benchmark
```

---

## Recommended Settings

### For Differential Expression Analysis

DE analysis is sensitive to imputation artefacts. Both afMF and ALRA can introduce correlated expression patterns that inflate test statistics. The recommended approach is to **run both methods** for completeness but **not use imputed data** for the DE tests themselves:

```r
imputation_method = "both"            # Generate both imputed objects
use_afmf_for_normalization = FALSE    # Normalise from original counts
use_alra_for_downstream = FALSE       # Integrate from original normalised data
run_clts_renormalization = TRUE        # CLTS genuinely improves DE accuracy
```

Imputed assays remain available in the Seurat object for inspection and visualisation, while DE operates on the standard pipeline output. CLTS re-normalisation (Module 07b) is recommended because it addresses a genuine normalisation bias rather than imputing missing values.

Note that imputation is incompatible with pseudobulk DE analysis via limma-trend — pseudobulk aggregation already smooths dropout effects, and prior imputation introduces redundant or conflicting corrections (Huang et al. 2025).

### For Visualisation and Exploration

Imputed data produce cleaner UMAP embeddings, sharper marker gene dot plots, and more interpretable heatmaps:

```r
imputation_method = "both"
# After pipeline completion, switch to imputed assays for plotting:
# DefaultAssay(seurat_obj) <- "imputed"   # afMF counts
# DefaultAssay(seurat_obj) <- "ALRA"      # ALRA normalised
```

### For Trajectory Inference

Both afMF and ALRA improve trajectory inference via diffusion pseudotime (DPT). However, imputation is **not compatible** with Slingshot-based trajectory methods — prior imputation either shows no improvement or generates false positive trajectory branches (Huang et al. 2025).

```r
# For DPT-based trajectories: imputation helps
imputation_method = "afmf"            # afMF ranked highest for DPT

# For Slingshot: do NOT impute
imputation_method = "none"
```

---

## Method Comparison

| Aspect | afMF | ALRA |
|--------|------|------|
| **Approach** | Iterative gradient descent on full matrix factorisation (X ≈ W × H) | Single-pass randomised SVD + adaptive thresholding |
| **Relationship** | Extended version of ALRA; replaces SVD with iterative optimisation | Original method; afMF was built upon this |
| **Operating data** | Raw counts | Log-normalised expression |
| **Zero preservation** | Via factorisation constraints | Via symmetry-based thresholding per gene |
| **Conservatism** | Moderate — recovers more zeros than ALRA | Conservative — prioritises preserving biological zeros |
| **Speed** | Slower (iterative; Python GPU-capable) | Fast (single SVD; R) |
| **Memory** | Moderate–high | Low–moderate |
| **SCTransform** | Compatible | Incompatible |
| **DE improvement** | Improves MAST / Wilcoxon; incompatible with pseudobulk | Similar profile to afMF |
| **Best for** | Cell type annotation, AUCell/SCENIC, DPT trajectories | General visualisation, marker gene recovery, quick exploration |

---

## Downstream Compatibility

The following table summarises which downstream analyses benefit from or are harmed by imputation, based on the benchmarking results from Huang et al. (2025):

| Downstream Analysis | afMF | ALRA | Recommendation |
|---------------------|------|------|----------------|
| Cell clustering (Leiden / Louvain) | Slight improvement | Slight improvement | Either; improvement modest |
| Automatic cell type annotation | **Good improvement** | Good improvement | afMF preferred |
| DE analysis (MAST / Wilcoxon) | **Improvement** | Improvement | Use with caution; keep `use_*_for_downstream = FALSE` |
| DE analysis (pseudobulk / limma) | **Incompatible** | **Incompatible** | Do not impute |
| GSEA | Improvement | Improvement | Either |
| AUCell / SCENIC | **Good improvement** | Good improvement | afMF preferred |
| Trajectory (DPT) | **Good improvement** | Improvement | afMF preferred |
| Trajectory (Slingshot) | **Incompatible** | **Incompatible** | Do not impute |
| Cell-cell communication | **Incompatible** | **Incompatible** | Do not impute |
| Visualisation (UMAP, heatmaps) | Good improvement | Good improvement | Either; cleaner plots |

---

## Troubleshooting

**Module 02b is skipped even though `imputation_method = "afmf"`**

Check that the afMF Python environment exists and is accessible:
```bash
ls ~/.conda/envs/afMF_SCImputation_env/bin/python
# If missing, install:
conda create -n afMF_SCImputation_env python=3.11
conda activate afMF_SCImputation_env
pip install numpy scipy scikit-learn anndata h5py
```
Also verify that Module 00 set the `has_afmf` flag in `objects/pipeline_environment.RData`.

**Module 03b is skipped with "SCTransform incompatible" warning**

This is expected behaviour. ALRA cannot operate on Pearson residuals. Either accept SCTransform without ALRA, or force a different normalisation method in `params.R`:
```r
normalization_methods = c("LogNormalize")  # Skip SCTransform benchmarking
```

**ALRA auto-detected rank k seems too low / too high**

Override with a manual value:
```r
alra_k = 40  # Set manually instead of NULL
```
Check the singular value spectrum plot in `03b_Imputation_ALRA/plots/` to guide your choice — k should correspond to the "elbow" where signal transitions to noise.

**afMF does not converge within `max_iter` iterations**

Increase the iteration limit or relax the tolerance:
```r
afmf_max_iter = 200
afmf_tol = 1e-4
```
Also ensure that very lowly-expressed genes are filtered (increase `afmf_min_cells_expressing`).

**Memory errors during imputation**

Both methods create dense matrices. For large datasets (>50,000 cells), ensure adequate memory allocation in SLURM:
```bash
#SBATCH --mem=128G
```
For afMF, the Python process may benefit from GPU acceleration if available.

---

## References

**afMF**:
Huang Y, et al. "An in-depth benchmark framework for evaluating single cell RNA-seq dropout imputation methods and the development of an improved algorithm afMF." *Clinical and Translational Medicine* 15(3):e70283 (2025). DOI: [10.1002/ctm2.70283](https://doi.org/10.1002/ctm2.70283)

**ALRA**:
Linderman GC, Zhao J, Rber E, Fong R, Kluger Y. "Zero-preserving imputation of single-cell RNA-seq data." *Nature Communications* 13:192 (2022). DOI: [10.1038/s41467-021-27729-z](https://doi.org/10.1038/s41467-021-27729-z)

**ALRA (Seurat integration)**:
Integrated into Seurat v3+ via `SeuratWrappers::RunALRA()`. GitHub: [KlugerLab/ALRA](https://github.com/KlugerLab/ALRA)

**General imputation benchmarking**:
Hou W, Ji Z, Ji H, Hicks SC. "A systematic evaluation of single-cell RNA-sequencing imputation methods." *Genome Biology* 21:218 (2020). DOI: [10.1186/s13059-020-02132-x](https://doi.org/10.1186/s13059-020-02132-x)
