# scCRIPTURE — Modular HPC Single-Cell scRNA-seq Bioinformatics Pipeline

> **End-to-end, HPC-native scRNA-seq pipeline for 10x Genomics Cell Ranger data**: FASTQ → Cell Ranger → multi-layered QC → doublet consensus → contamination correction → normalisation benchmarking (SCTransform / scran / LogNormalise / scKWARN) → multi-method integration benchmarking (Harmony, CCA, RPCA, FastMNN, scVI, Scanorama, BBKNN, scCobra, CONCORD) → CHOIR & scAURA clustering → scICE & IDclust subclustering → differential expression (DESeq2, MAST, Dream, EdgeR) → interactive reporting. Fully modular, samplesheet-driven, zero hardcoded paths.

<p align="center">
  <img src="scheme.png" alt="scCRIPTURE pipeline schematic" width="100%"/>
</p>

---

## Table of Contents

- [Overview](#overview)
- [Key Design Principles](#key-design-principles)
- [Installation](#installation)
- [Pipeline Architecture](#pipeline-architecture)
  - [Phase 1: Preprocessing (Steps 1–10)](#phase-1-preprocessing-steps-110)
  - [Phase 2: Downstream Analysis (Steps 11–21)](#phase-2-downstream-analysis-steps-1121)
  - [Optional Modules](#optional-modules)
  - [MapMyCells Cell Type Annotation Transfer](#mapmycells-cell-type-annotation-transfer)
- [Configuration](#configuration)
  - [Samplesheet System](#samplesheet-system)
  - [Pipeline Configurator GUI](#pipeline-configurator-gui)
  - [params.R Configuration](#paramsr-configuration)
- [Usage](#usage)
- [Directory Structure](#directory-structure)
  - [Script Organisation](#script-organisation)
  - [Output Folder Structure](#output-folder-structure)
- [Environment Variables](#environment-variables)
- [Module Reference](#module-reference)
- [Citation](#citation)
- [License](#license)

---

## Overview

**scCRIPTURE** processes single-cell RNA-seq data through a two-phase workflow. **Phase 1** (preprocessing, Steps 1–10) handles per-sample processing as SLURM array jobs — alignment, empty droplet detection, doublet removal, and contamination correction. **Phase 2** (downstream analysis, Steps 11–21) performs group-specific multi-sample integration, clustering, differential expression, and visualisation using a modular R framework.

The pipeline is species- and tissue-agnostic. All configuration is driven by a central `samplesheet.csv` and a `params.R` file, with no hardcoded sample names or paths. A browser-based configurator GUI is provided for guided parameter setup.

## Key Design Principles

- **Modular**: Each step is self-contained; modules can be run independently or as a chain with automatic dependency tracking.
- **Universal**: All sample metadata is read from `samplesheet.csv` — no hardcoded identifiers anywhere.
- **Array Jobs**: Preprocessing steps 1–8 run as SLURM array jobs (one task per sample) for efficient parallelisation.
- **Group-Based Design**: Downstream analysis supports arbitrary experimental groupings via `group_id` / `group_label` columns, with full backward compatibility for legacy `--ventricle` syntax.
- **Optional Modules**: Imputation (afMF, ALRA) and CLTS re-normalisation can be toggled on/off via `params.R`.
- **Environment Variables**: Configuration propagated through `PROJECT_ROOT`, `PREPROCESS_DIR`, `DOWNSTREAM_DIR`, `DATASET_NAME`, `SAMPLESHEET`, `GROUP_ID`, `GROUP_LABEL`.

---

## Installation

The pipeline requires Python 3.11, R ≥ 4.5, and Julia (for scICE). A four-step installation process sets up the complete environment:

```bash
# Step 1: Create the conda environment (Python + R + Julia + conda packages)
conda env create -f kaiser_test_py3.11.yml

# Step 2: Install Python pip packages (using uv for speed)
pip install uv
uv pip install -r requirements.txt

# Step 3: Install GitHub-only R packages
Rscript install_r_packages.R

# Step 4: Install additional tools (scAURA + scICE)
bash install_tools.sh
```

Steps 2 and 4 can be combined by running `install_tools.sh`, which handles `uv` installation and `requirements.txt` automatically. Some packages (scanpy, leidenalg, igraph, celltypist, gseapy) install more reliably via conda — see `kaiser_test_py3.11.yml` for the recommended conda-based installation.

**GPU requirements**: PyTorch, DGL, and RAPIDS packages require NVIDIA GPU with CUDA 12.x drivers. Remove the GPU sections from `requirements.txt` if no GPU is available.

**Conda environments used by the pipeline**:

| Environment | Used by |
|---|---|
| `R_4_5` | R-based steps (2, 5–8, 10, 11–21) |
| `qclus` | Step 3 (QClus) |
| `vaeda_env` | Step 4 (VAEDA) |
| `afMF_SCImputation_env` | Module 02b (afMF imputation) |
| `kaiser_test_py3.11` | Python scripts, MapMyCells preparation |

---

## Pipeline Architecture

### Phase 1: Preprocessing (Steps 1–10)

Preprocessing runs as SLURM array jobs (one task per sample). Each step reads from the preceding step's output and writes to a numbered subdirectory under `Single_cell_preprocessed/`.

| Step | Module | Description | Key Tools |
|------|--------|-------------|-----------|
| 1 | Cell Ranger Count | Align reads, call cells, generate count matrix | Cell Ranger |
| 2 | DropletQC | Identify empty droplets and damaged cells via nuclear fraction analysis | DropletQC (R) |
| 3 | QClus | Additional empty droplet filtering using the QClus algorithm | QClus (Python) |
| 4 | VAEDA | Variational autoencoder-based doublet detection (labels only, no removal) | VAEDA (Python) |
| 5 | Format Conversion | Convert AnnData (h5ad) → Seurat (RDS), preserving all QC metadata | anndata2ri / SeuratDisk |
| 6 | Consensus Doublets | Run scDblFinder + DoubletFinder, then majority-vote with VAEDA (≥ 2/3 → remove) | scDblFinder, DoubletFinder (R) |
| 7 | scCDC Correction | Detect and correct ambient RNA contamination; produces the primary input for downstream analysis | scCDC (R) |
| 8 | DecontX Correction | Additional ambient RNA decontamination | DecontX / celda (R) |
| 9 | CHOIR Integration | Preliminary sample integration and clustering across all samples and per-group | CHOIR (R) |
| 10 | QC Visualisation | Generate comprehensive QC plots and filtering-flow summaries | R / ggplot2 |

**Data flow**: Raw FASTQ → Cell Ranger → DropletQC → QClus → VAEDA → Seurat conversion → Consensus doublet removal → scCDC correction → DecontX → CHOIR integration → QC plots. The primary input for downstream analysis is the scCDC-corrected RDS from Step 7.

### Phase 2: Downstream Analysis (Steps 11–21)

Downstream analysis runs as a single SLURM job per experimental group. The R pipeline (`run_pipeline.R`) executes modules sequentially, controlled by `params.R`.

| Module | Description | Key Methods |
|--------|-------------|-------------|
| 00 | Environment Setup | Load libraries, validate paths, detect optional modules |
| 01 | Load Data | Read scCDC-corrected RDS files, merge samples by group |
| 02 | QC Validation | Apply additional QC filters (nFeature, nCount, mito%, ribo%) |
| 02b | **[Optional]** afMF Imputation | Low-rank matrix factorisation on raw counts (before normalisation) |
| 03 | Normalisation | Benchmark SCTransform vs scran vs LogNormalise vs scKWARN; select best method |
| 03b | **[Optional]** ALRA Imputation | Impute dropouts on log-normalised data (incompatible with SCTransform) |
| 04 | Integration & Benchmarking | Multi-method integration benchmarking (Harmony, CCA, RPCA, FastMNN, scVI, Scanorama, BBKNN, scCobra, CONCORD) |
| 05 | CHOIR & scAURA Clustering | Hierarchical clustering (CHOIR) and/or autoencoder-based clustering (scAURA) |
| 06 | Subclustering | Refine clusters via scICE (Julia) and/or IDclust (R) |
| 07 | Leiden Clustering | Standard Leiden community detection at multiple resolutions |
| 07b | **[Optional]** CLTS Re-normalisation | Remove cell-type composition bias for improved DE accuracy |
| 08 | Differential Expression | DE testing via DESeq2, MAST, Dream, and EdgeR |
| 09 | Gene Visualisation | Marker gene expression plots (dot plots, feature plots, violin plots) |
| 10 | Final Summary | Compile summary statistics and metadata tables |
| 11 | HTML Report | Generate interactive HTML report with all results |

### Optional Modules

Three optional modules can be enabled via `params.R`:

**afMF Imputation (Module 02b)** — operates on raw counts *before* normalisation. Controlled by `imputation_method = "afmf"` or `"both"`. Requires a dedicated Python environment (`afMF_SCImputation_env`).

**ALRA Imputation (Module 03b)** — operates on log-normalised data *after* normalisation. Controlled by `imputation_method = "alra"` or `"both"`. Compatible with LogNormalise, scran, and scKWARN only; automatically skipped if SCTransform is selected.

**CLTS Re-normalisation (Module 07b)** — removes cell-type-specific composition bias. Controlled by `run_clts_renormalization = TRUE`. Improves differential expression accuracy between conditions.

| Method | Module | Stage | Input | Output Assay | SCTransform Compatible |
|--------|--------|-------|-------|-------------|----------------------|
| afMF | 02b | Before norm | Raw counts | `imputed` | Yes |
| ALRA | 03b | After norm | Log-normalised | `ALRA` | No (auto-skipped) |

### MapMyCells Cell Type Annotation Transfer

The pipeline includes utility scripts for transferring cell type labels from the [Allen Institute MapMyCells](https://knowledge.brain-map.org/mapmycells/process) service. MapMyCells maps user-submitted single-cell transcriptomic data onto Allen Institute reference taxonomies (Whole Mouse Brain, SEA-AD Human MTG, Whole Human Brain), providing hierarchical cell type annotations at the class, subclass, supertype, and cluster levels.

**Workflow**:

1. **Prepare data for upload** (`Python_scripts/MapMyCells_preparation_script.ipynb`): Extract the expression matrix from the Seurat/AnnData object, format it as a compressed h5ad file (cells × genes, ≤ 2 GB), and validate gene identifiers against the target taxonomy.

2. **Run MapMyCells**: Upload the prepared h5ad file to the [MapMyCells web interface](https://knowledge.brain-map.org/mapmycells/process) or use the [MapMyCells Python API](https://github.com/AllenInstitute/cell_type_mapper) for programmatic access. Select the appropriate reference taxonomy and mapping algorithm (hierarchical correlation, flat correlation, or deep generative mapping).

3. **Transfer labels back** to the Seurat object:
   - **Python** (`Python_scripts/transfer_mapmycells_labels.py`): Reads the MapMyCells CSV output, merges hierarchical annotations (class, subclass, supertype, cluster) and confidence scores into the AnnData object, and exports the updated metadata.
   - **R** (`R_scripts/transfer_mapmycells_labels.R`): Reads the MapMyCells CSV output, joins annotations onto the Seurat object's metadata, and saves the updated RDS. Annotations are added as metadata columns (`mapmycells_class`, `mapmycells_subclass`, `mapmycells_supertype`, `mapmycells_cluster`) alongside their corresponding confidence scores.

4. **SLURM submission** (`Slurm_scripts/MapMyCells_transfer_for_preprocessing.sh`): Wrapper script for batch execution of the label transfer on the HPC cluster.

**Available reference taxonomies**: Whole Mouse Brain (10x, ~4M cells, 5,322 clusters), SEA-AD Human MTG (aged human cohort, 139 supertypes), Whole Human Brain (3M+ cells, 3,313 subclusters), and additional species-specific BICAN atlases.

---

## Configuration

### Samplesheet System

The pipeline uses two CSV samplesheets:

**`samplesheet.csv`** (primary metadata — used by all steps):

| Column | Name | Required | Description |
|--------|------|----------|-------------|
| 1 | `sample_id` | Yes | Unique sample identifier |
| 2 | `sex` | Yes | M or F |
| 3 | `age` | No | Sample age (e.g., P60, 8w) |
| 4 | `batch` | Yes | Batch identifier for integration |
| 5 | `ventricle` | No | Anatomical region (backward compat) |
| 6 | `condition` | No | Treatment condition |
| 7 | `SRA_accession` | No | SRA accession number |
| 8–10 | Paths | No | FTP, local FASTQ, genome reference |
| 11 | `include` | No | TRUE/FALSE to include sample |
| 12 | `dataset_name` | Recommended | Dataset identifier for output naming |
| 13 | `group_id` | Recommended | Numeric group ID for filtering |
| 14 | `group_label` | Recommended | Human-readable group name |

**`fastq_samplesheet.csv`** (FASTQ file locations — used by Step 1 only): Maps original FASTQ filenames (e.g., SRA accessions) to Cell Ranger naming convention. Supports multi-lane samples via multiple rows.

### Pipeline Configurator GUI

A self-contained, browser-based HTML application (`params_scrna_pipeline_configurator.html`) provides a guided wizard for configuring all pipeline parameters. It replaces manual editing of `params.R` with a step-by-step interface across 15 steps in two phases:

**Phase A (Sample Sheet, steps A1–A5)**: Load, assign column roles, inspect experimental design, validate, and export the sample sheet.

**Phase B (params.R Configuration, steps B0–B9)**: Configure project metadata, QC thresholds, normalisation methods, imputation strategy, batch integration, clustering methods, differential expression settings, and visualisation options. A live R-code preview shows the generated configuration.

Key features include tag-based template injection (`<<CFG:key>>` tags), per-parameter approve/modify workflow, input validation with contextual warnings (e.g., SCTransform ↔ ALRA incompatibility), session save/restore (JSON), and progress tracking. Runs entirely client-side from a single HTML file.

### params.R Configuration

The `params.R` file controls all downstream analysis behaviour. Key parameter groups:

- **Project metadata**: paths, dataset name, species, tissue
- **QC thresholds**: nFeature, nCount, mitochondrial %, ribosomal %
- **Normalisation**: method selection and benchmarking (SCTransform / scran / LogNormalise / scKWARN)
- **Imputation**: `imputation_method` ∈ {`"none"`, `"afmf"`, `"alra"`, `"both"`}
- **Integration**: method benchmarking weights (`bio_weight`, `batch_weight`), methods (Harmony, CCA, RPCA, FastMNN, scVI, Scanorama, BBKNN, scCobra, CONCORD)
- **Clustering**: CHOIR, scAURA, Leiden, scICE, IDclust, scCobra parameters
- **Differential expression**: DESeq2, MAST, Dream, EdgeR; thresholds, covariates
- **Visualisation**: output paths, plot dimensions, gene lists

---

## Usage

```bash
# Validate setup before running
bash validate_pipeline.sh

# List available experimental groups
./submit_all_jobs.sh --list-groups

# Run full preprocessing (Steps 1–10, array jobs)
./submit_all_jobs.sh --preprocess

# Run specific preprocessing steps
./submit_all_jobs.sh --only 1          # Just Cell Ranger
./submit_all_jobs.sh --only 2,3        # DropletQC + QClus
./submit_all_jobs.sh --from 4 --to 7   # Steps 4–7

# Downstream analysis by group
./submit_all_jobs.sh --downstream --1              # Group 1
./submit_all_jobs.sh --downstream --2              # Group 2
./submit_all_jobs.sh --downstream --group 1        # Explicit group ID

# Downstream with optional modules
./submit_all_jobs.sh --downstream --1 --modules 2b,3       # afMF + normalisation
./submit_all_jobs.sh --downstream --1 --modules 3,3b,4     # Norm + ALRA + integration
./submit_all_jobs.sh --downstream --1 --modules 7b,8       # CLTS + DE

# Legacy ventricle syntax (still supported)
./submit_all_jobs.sh --downstream --ventricle LV

# Dry run (show what would be submitted)
./submit_all_jobs.sh --preprocess --dry-run
```

---

## Directory Structure

### Script Organisation

```
Scripts/
├── Slurm_scripts/                          # SLURM job submission scripts
│   ├── 01_cellranger_count.sh
│   ├── 02_DropletQC_nuclear_fraction.sh
│   ├── 03_qClus_empty_droplets.sh
│   ├── 04_Vaeda_doublet_analysis.sh
│   ├── 05_Vaeda_anndata_to_seurat.sh
│   ├── 06_scDblFinder_DoubletFinder_analysis.sh
│   ├── 07_scCDC_correction.sh
│   ├── 08_DecontX_correction.sh
│   ├── 09_CHOIR_integration.sh
│   ├── 10_QC_visualization.sh
│   ├── 11_downstream_analysis.sh
│   └── MapMyCells_transfer_for_preprocessing.sh
│
├── Python_scripts/                         # MapMyCells utilities
│   ├── MapMyCells_preparation_script.ipynb #   Data preparation notebook
│   └── transfer_mapmycells_labels.py       #   Label transfer (Python)
│
├── R_scripts/                              # MapMyCells utilities
│   └── transfer_mapmycells_labels.R        #   Label transfer (R/Seurat)
│
├── Analytical_block/                       # Post-pipeline analytical modules
│   ├── ANALYTICAL_BLOCK_ARCHITECTURE.md
│   ├── Python_scripts/
│   │   └── AB01_General_overview/          #   Overview analysis (Python)
│   ├── R_scripts/
│   │   └── AB01_General_overview/          #   Overview analysis (R)
│   │       ├── ab_01_cluster_overview.R
│   │       └── ab_01_config.R
│   └── Slurm_scripts/
│       └── Ab01_general_overview.sh
│
└── scrnaseq_pipeline/                      # Core downstream R pipeline
    ├── run_pipeline.R                      #   Main pipeline runner
    ├── config/
    │   ├── params.R                        #   Pipeline parameters
    │   ├── params_scrna_pipeline_configurator.html  # GUI configurator
    │   └── sample_sheet.csv                #   Downstream sample inclusion
    ├── modules/
    │   ├── 00_setup_environment.R
    │   ├── 01_load_data.R
    │   ├── 02_qc_validation.R
    │   ├── 02b_imputation.R                #   [Optional] afMF
    │   ├── 03_normalization.R
    │   ├── 03b_alra_imputation.R           #   [Optional] ALRA
    │   ├── 04_integration.R
    │   ├── 05_choir_clustering.R
    │   ├── 06_scice_subclustering.R
    │   ├── 07_leiden_clustering.R
    │   ├── 07b_clts_renormalization.R      #   [Optional] CLTS
    │   ├── 08_differential_expression.R
    │   ├── 09_gene_visualization.R
    │   ├── 10_final_summary.R
    │   └── 11_html_report.R
    ├── docs/
    │   └── IMPUTATION_METHODS.md
    └── utils/
        ├── functions.R                     #   Shared helper functions
        └── python/
            └── run_afmf_imputation.py      #   afMF Python wrapper
```

### Output Folder Structure

```
project_root/
├── samplesheet.csv                         # Primary sample metadata
├── fastq_samplesheet.csv                   # FASTQ file mappings
├── submit_all_jobs.sh                      # Master job submission
├── validate_pipeline.sh                    # Pre-run validation
├── Complete_FASTQ/                         # Raw FASTQ files
├── config/                                 # Samplesheet parsers
│   ├── parse_samplesheet.R
│   └── parse_samplesheet.py
├── logs/                                   # SLURM output logs
│   ├── cellranger/
│   ├── dropletqc/
│   ├── qclus/
│   ├── vaeda/
│   ├── h5ad2rds/
│   ├── doublets/
│   ├── scCDC/
│   ├── decontx/
│   ├── CHOIR/
│   ├── qc_viz/
│   └── downstream/
│
└── Output_dir_<dataset_name>/
    │
    ├── Single_cell_preprocessed/           # Steps 1–10 outputs
    │   ├── 1_CellRanger_output/
    │   │   └── {SAMPLE}/outs/
    │   │       ├── filtered_feature_bc_matrix/
    │   │       ├── raw_feature_bc_matrix/
    │   │       ├── possorted_genome_bam.bam
    │   │       └── web_summary.html
    │   ├── 2_DropletQC_output/{SAMPLE}/
    │   ├── 3_qClus_empty_droplets/{SAMPLE}/
    │   ├── 4_Vaeda_doublet_detection/{SAMPLE}/
    │   ├── 5_Seurat_conversion/{SAMPLE}/
    │   ├── 6_Doublet_consensus/{SAMPLE}/
    │   ├── 7_scCDC_correction/{SAMPLE}/    # ← Primary downstream input
    │   │   └── {SAMPLE}_scCDC_corrected.rds
    │   ├── 8_DecontX_correction/{SAMPLE}/
    │   ├── 9_CHOIR_integration/
    │   │   ├── Individual/{SAMPLE}/
    │   │   ├── {group_label}_only/
    │   │   └── All_combined/
    │   └── 10_QC_visualization/
    │
    └── Single_cell_clustering/             # Steps 11+ outputs
        ├── 11_Downstream_Analysis_<group_label_A>/
        │   ├── objects/
        │   │   ├── pipeline_environment.RData
        │   │   ├── merged_object.rds
        │   │   ├── merged_object_imputed.rds       # [if afMF enabled]
        │   │   ├── merged_normalized/
        │   │   │   ├── merged_LogNormalize_unintegrated.rds
        │   │   │   └── merged_LogNormalize_alra.rds # [if ALRA enabled]
        │   │   ├── scice_subclustered_object.rds
        │   │   └── final_clustered_object.rds
        │   ├── plots/
        │   │   ├── qc/
        │   │   ├── normalization/
        │   │   ├── integration/
        │   │   ├── clustering/
        │   │   ├── clts_normalization/             # [if CLTS enabled]
        │   │   ├── de/
        │   │   └── gene_expression/
        │   ├── tables/
        │   ├── 02b_Imputation_afMF/                # [if afMF enabled]
        │   ├── 03b_Imputation_ALRA/                # [if ALRA enabled]
        │   ├── 05_CHOIR_Clustering/
        │   ├── 07b_CLTS_Normalization/             # [if CLTS enabled]
        │   ├── 08_Differential_Expression/
        │   └── reports/
        │       └── Analysis_Report_MultiSample.html
        │
        └── 11_Downstream_Analysis_<group_label_B>/
            └── [same structure as above]
```

---

## Environment Variables

| Variable | Description |
|----------|-------------|
| `PROJECT_ROOT` | Base directory for the project |
| `SAMPLESHEET` | Path to `samplesheet.csv` |
| `DATASET_NAME` | Dataset identifier (from samplesheet column 12) |
| `PREPROCESS_DIR` | `Output_dir_<dataset_name>/Single_cell_preprocessed` |
| `DOWNSTREAM_DIR` | `Output_dir_<dataset_name>/Single_cell_clustering` |
| `PIPELINE_DIR` | `Scripts/scrnaseq_pipeline` |
| `GROUP_ID` | Numeric group ID for downstream filtering |
| `GROUP_LABEL` | Human-readable group label |
| `VENTRICLE_FILTER` | LV, 4V, or ALL (backward compatible) |

---

## Module Reference

| Module | Name | Optional | Controlled by |
|--------|------|----------|---------------|
| 00 | Environment Setup | No | Always runs |
| 01 | Load Data | No | Always runs |
| 02 | QC Validation | No | Always runs |
| 02b | Imputation (afMF) | **Yes** | `imputation_method ∈ {"afmf", "both"}` |
| 03 | Normalisation | No | Always runs |
| 03b | Imputation (ALRA) | **Yes** | `imputation_method ∈ {"alra", "both"}` AND norm ∈ {LogNorm, scran} |
| 04 | Integration | No | Always runs |
| 05 | CHOIR & scAURA Clustering | No | `run_choir_clustering` / scAURA settings |
| 06 | Subclustering | No | `subclustering_methods` / `run_scice_subclustering` / `run_idclust_subclustering` |
| 07 | Leiden Clustering | No | `run_leiden_clustering` |
| 07b | CLTS Re-normalisation | **Yes** | `run_clts_renormalization = TRUE` |
| 08 | Differential Expression | No | Always runs |
| 09 | Gene Visualisation | No | Always runs |
| 10 | Final Summary | No | Always runs |
| 11 | HTML Report | No | Always runs |

---

## Citation

If you use this pipeline, please cite the individual tools used in each step:

| Tool | Reference | DOI |
|------|-----------|-----|
| **Cell Ranger** | 10x Genomics | [`10.1038/ncomms14049`](https://doi.org/10.1038/ncomms14049) |
| **DropletQC** | Muskovic & Powell (2021) | [`10.1186/s13059-021-02547-0`](https://doi.org/10.1186/s13059-021-02547-0) |
| **QClus** | Schmauch et al. (2024) | [`10.1093/nar/gkae1145`](https://doi.org/10.1093/nar/gkae1145) |
| **VAEDA** | Doublet detection via variational autoencoder | [`10.1093/bioinformatics/btac720`](https://doi.org/10.1093/bioinformatics/btac720) |
| **scDblFinder** | Germain et al. (2022) | [`10.12688/f1000research.73600.2`](https://doi.org/10.12688/f1000research.73600.2) |
| **DoubletFinder** | McGinnis et al. (2019) | [`10.1016/j.cels.2019.03.003`](https://doi.org/10.1016/j.cels.2019.03.003) |
| **scCDC** | Contamination correction | [`10.1186/s13059-024-03284-w`](https://doi.org/10.1186/s13059-024-03284-w) |
| **DecontX** | Yang et al. (2020) | [`10.1186/s13059-020-1950-6`](https://doi.org/10.1186/s13059-020-1950-6) |
| **scKWARN** | Hsu et al. (2024) | [`10.1093/bioinformatics/btae008`](https://doi.org/10.1093/bioinformatics/btae008) |
| **CHOIR** | Hierarchical clustering | [`10.1038/s41588-025-02148-8`](https://doi.org/10.1038/s41588-025-02148-8) |
| **scICE** | Consistent subclustering | [`10.1038/s41467-025-60702-8`](https://doi.org/10.1038/s41467-025-60702-8) |
| **IDclust** | Prompsy et al. (2024) | [`10.1093/nargab/lqae174`](https://doi.org/10.1093/nargab/lqae174) |
| **Seurat** | Hao et al. (2024) | [`10.1038/s41587-023-01767-y`](https://doi.org/10.1038/s41587-023-01767-y) |
| **scVI-tools** | Gayoso et al. (2022) | [`10.1038/s41587-021-01206-w`](https://doi.org/10.1038/s41587-021-01206-w) |
| **CONCORD** | Zhu et al. (2026) | [`10.1038/s41587-025-02950-z`](https://doi.org/10.1038/s41587-025-02950-z) |
| **scCobra** | Zhao et al. (2025) | [`10.1038/s42003-025-07692-x`](https://doi.org/10.1038/s42003-025-07692-x) |
| **scAURA** | Rifat et al. (2026, preprint) | [`10.64898/2026.01.25.701579`](https://doi.org/10.64898/2026.01.25.701579) |
| **MapMyCells** | Yao et al. (2023) — [Allen Institute](https://knowledge.brain-map.org/mapmycells/process) | [`10.1038/s41586-023-06812-z`](https://doi.org/10.1038/s41586-023-06812-z) |
| **MAST** | Finak et al. (2015) | [`10.1186/s13059-015-0844-5`](https://doi.org/10.1186/s13059-015-0844-5) |
| **dreamlet** | Hoffman et al. (2023, preprint) | [`10.1101/2023.03.17.533005`](https://doi.org/10.1101/2023.03.17.533005) |

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
