# ==============================================================================
# UNIVERSAL MODULAR scRNA-seq PIPELINE - DOCUMENTATION
# ==============================================================================
#
# Project: Universal Modular scRNA-seq Pipeline (adaptable to any scRNA-seq project)
#
# ==============================================================================

================================================================================
TABLE OF CONTENTS
================================================================================

1. OVERVIEW
2. FOLDER STRUCTURE
3. SAMPLESHEET SYSTEM
4. PIPELINE STEPS (1-21)
5. SLURM SCRIPTS - MODULAR DESIGN
6. HELPER FUNCTIONS
7. SUBMIT_ALL_JOBS.SH - MASTER CONTROLLER
8. DOWNSTREAM ANALYSIS (STEPS 11-21)
9. ADAPTING FOR NEW DATASETS
10. TROUBLESHOOTING
11. QUICK REFERENCE

================================================================================
1. OVERVIEW
================================================================================

This pipeline processes single-cell RNA-seq data from raw FASTQ files through
Cell Ranger, quality control, doublet detection, batch correction, clustering,
and downstream analysis.

KEY DESIGN PRINCIPLES:
----------------------
- MODULAR: Each step is self-contained and reads configuration from samplesheets
- UNIVERSAL: No hardcoded sample names - all derived from samplesheet.csv
- ARRAY JOBS: Steps 1-8 run as SLURM array jobs (one task per sample)
- ENVIRONMENT VARIABLES: Configuration passed via PROJECT_ROOT, PREPROCESS_DIR,
  DOWNSTREAM_DIR, DATASET_NAME, SAMPLESHEET, GROUP_ID, GROUP_LABEL
- DATASET-SPECIFIC OUTPUT: All outputs go to Output_dir_<dataset_name>/
- GENERIC GROUP FILTERING: Downstream analysis can filter by any samplesheet-
  defined group using numeric shorthand (--1, --2) or group IDs
- BACKWARD COMPATIBLE: Legacy --ventricle syntax still fully supported
- OPTIONAL MODULES: Modules 02b (afMF), 03b (ALRA), and 07b (CLTS) can be
  enabled via params.R settings
- SUBCLUSTERING METHODS: Module 06 supports scICE and IDclust; enable/disable
  each method and choose the cluster source via params.R

DATASET CONFIGURATION:
----------------------
- Species: Configured in samplesheet.csv
- Tissue: Configured in samplesheet.csv
- Samples: Auto-detected from samplesheet.csv
- Reference: Specified per sample in samplesheet.csv
- Dataset Name: Specified in samplesheet.csv (column 12)
- Groups: Defined by group_id and group_label columns in samplesheet.csv

================================================================================
2. FOLDER STRUCTURE
================================================================================

project_root/
│
├── samplesheet.csv                 # PRIMARY: Universal sample metadata
├── fastq_samplesheet.csv           # FASTQ file locations for Cell Ranger
├── submit_all_jobs.sh              # Master job submission script
├── validate_pipeline.sh            # Pre-run validation script
├── README.txt                      # This file
│
├── Complete_FASTQ/                 # Raw FASTQ files
│   ├── {SRA_ID}_1.fastq           # Index read (I1)
│   ├── {SRA_ID}_2.fastq           # Read 1 (R1) - cell barcode + UMI
│   ├── {SRA_ID}_3.fastq           # Read 2 (R2) - cDNA
│   └── ...
│
├── config/                         # Samplesheet parsers
│   ├── parse_samplesheet.R         # R helper functions
│   └── parse_samplesheet.py        # Python helper functions (with group support)
│
├── Scripts/
│   ├── Slurm_scripts/              # SLURM job scripts (Steps 1-11)
│   │   ├── 01_cellranger_count.sh
│   │   ├── 02_DropletQC_nuclear_fraction.sh
│   │   ├── 03_qClus_empty_droplets.sh
│   │   ├── 04_Vaeda_doublet_analysis.sh
│   │   ├── 05_Vaeda_anndata_to_seurat.sh
│   │   ├── 06_scDblFinder_DoubletFinder_analysis.sh
│   │   ├── 07_scCDC_correction.sh
│   │   ├── 08_DecontX_correction.sh
│   │   ├── 09_CHOIR_integration.sh
│   │   ├── 10_QC_visualization.sh
│   │   └── 11_downstream_analysis.sh
│   │
│   ├── R_scripts/                  # Generated R scripts (per-sample)
│   ├── Python_scripts/             # Generated Python scripts (per-sample)
│   │
│   └── scrnaseq_pipeline/          # Downstream analysis modules
│       ├── run_pipeline.R          # Main pipeline runner
│       ├── Define_input_paths.txt  # Path configuration guide
│       ├── PIPELINE_PROCESSING_LOGIC_SUMMARY.txt  # Processing documentation
│       ├── config/
│       │   ├── params.R            # Pipeline parameters (reads env vars)
│       │   └── sample_sheet.csv    # Downstream sample inclusion
│       ├── modules/
│       │   ├── 00_setup_environment.R
│       │   ├── 01_load_data.R
│       │   ├── 02_qc_validation.R
│       │   ├── 02b_imputation.R        # [OPTIONAL] afMF counts-based imputation
│       │   ├── 03_normalization.R
│       │   ├── 03b_alra_imputation.R   # [OPTIONAL] ALRA normalized-data imputation
│       │   ├── 04_integration.R
│       │   ├── 05_choir_clustering.R
│       │   ├── 06_scice_subclustering.R
│       │   ├── 07_leiden_clustering.R
│       │   ├── 07b_clts_renormalization.R  # [OPTIONAL] CLTS re-normalization
│       │   ├── 08_differential_expression.R
│       │   ├── 09_gene_visualization.R
│       │   ├── 10_final_summary.R
│       │   └── 11_html_report.R
│       ├── docs/
│       │   └── IMPUTATION_METHODS.md   # Detailed imputation documentation
│       └── utils/
│           ├── functions.R
│           └── python/
│               └── run_afmf_imputation.py
│
├── logs/                           # SLURM output logs
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
└── Output_dir_<dataset_name>/      # DATASET-SPECIFIC OUTPUT ROOT
    │
    ├── Single_cell_preprocessed/   # Steps 1-10 outputs (PREPROCESS_DIR)
    │   ├── 1_CellRanger_output/
    │   │   └── {SAMPLE}/outs/
    │   ├── 2_DropletQC_output/
    │   │   └── {SAMPLE}/
    │   ├── 3_qClus_empty_droplets/
    │   │   └── {SAMPLE}/
    │   ├── 4_Vaeda_doublet_detection/
    │   │   └── {SAMPLE}/
    │   ├── 5_Seurat_conversion/
    │   │   └── {SAMPLE}/
    │   ├── 6_Doublet_consensus/
    │   │   └── {SAMPLE}/
    │   ├── 7_scCDC_correction/     # ← INPUT FOR DOWNSTREAM R PIPELINE
    │   │   └── {SAMPLE}/
    │   │       └── {SAMPLE}_scCDC_corrected.rds
    │   ├── 8_DecontX_correction/
    │   │   └── {SAMPLE}/
    │   ├── 9_CHOIR_integration/
    │   │   ├── Individual/
    │   │   │   └── {SAMPLE}/
    │   │   ├── {group_label}_only/
    │   │   └── All_combined/
    │   └── 10_QC_visualization/
    │
    └── Single_cell_clustering/     # Step 11+ outputs (DOWNSTREAM_DIR)
        ├── 11_Downstream_Analysis_<group_label_A>/   # Group 1 (--1)
        │   ├── objects/
        │   │   ├── pipeline_environment.RData
        │   │   ├── 01_loaded_data.RData
        │   │   ├── 02_qc_data.RData
        │   │   ├── merged_object.rds
        │   │   ├── merged_object_imputed.rds   # [if imputation_method includes afmf]
        │   │   ├── merged_normalized/
        │   │   │   ├── merged_LogNormalize_unintegrated.rds
        │   │   │   ├── merged_LogNormalize_alra.rds  # [if imputation_method includes alra]
        │   │   │   └── ...
        │   │   ├── 03_normalization_data.RData
        │   │   ├── individual_normalized/
        │   │   ├── scice_subclustered_object.rds
        │   │   ├── 06_scice_data.RData              # scICE + IDclust results
        │   │   ├── scice_clts_renormalized_object.rds  # [if run_clts]
        │   │   └── final_clustered_object.rds
        │   ├── plots/
        │   │   ├── qc/
        │   │   ├── normalization/
        │   │   ├── integration/
        │   │   ├── clustering/
        │   │   ├── clts_normalization/     # [if run_clts]
        │   │   ├── de/
        │   │   └── gene_expression/
        │   ├── tables/
        │   ├── 02b_Imputation_afMF/        # [if imputation_method includes afmf]
        │   │   ├── tables/
        │   │   └── plots/
        │   ├── 03b_Imputation_ALRA/        # [if imputation_method includes alra]
        │   │   ├── tables/
        │   │   └── plots/
        │   ├── 05_CHOIR_Clustering/
        │   ├── 07b_CLTS_Normalization/     # [if run_clts]
        │   ├── 08_Differential_Expression/
        │   ├── reports/
        │   │   └── Analysis_Report_MultiSample.html
        │   └── README.txt
        │
        └── 11_Downstream_Analysis_<group_label_B>/   # Group 2 (--2)
            └── [same structure as above]

================================================================================
3. SAMPLESHEET SYSTEM
================================================================================

The pipeline uses TWO samplesheets that work together:

------------------------------------------------------------------------------
3.1 samplesheet.csv - PRIMARY METADATA
------------------------------------------------------------------------------

Location: Project root
Used by: ALL pipeline steps
Purpose: Define samples and their metadata

REQUIRED COLUMNS:
  sample_id      : Unique sample identifier (e.g., sample_A, sample_B)
  sex            : M or F
  batch          : Batch identifier for integration

RECOMMENDED COLUMNS:
  age            : Sample age (e.g., P60, 8w)
  ventricle      : Anatomical region - for backward compatibility
  fastq_path     : Path to FASTQ directory
  cellranger_ref : Path to Cell Ranger reference
  species        : Species identifier (e.g., mouse, human)
  tissue         : Tissue type
  condition      : Treatment condition
  include        : TRUE/FALSE to include sample in analysis
  dataset_name   : Dataset identifier for Output_dir naming (COLUMN 12)
  group_id       : Numeric group identifier for filtering (COLUMN 13)
  group_label    : Human-readable group name (COLUMN 14)
  notes          : Free text notes

GROUP FILTERING COLUMNS:
-----------------------------------
  group_id       : Numeric identifier (1, 2, 3, ...) for shorthand filtering
  group_label    : Descriptive label used for output directory naming

These columns enable the generic filtering system:
  - Use --1 instead of --ventricle LV
  - Use --2 instead of --ventricle 4V
  - Output directories named by group_label (e.g., 11_Downstream_Analysis_<group_label>)

EXAMPLE (with group columns):
```
sample_id,sex,age,batch,ventricle,condition,SRA_accession,ftp_path,local_fastq_dir,genome,include,dataset_name,group_id,group_label
sample_A,M,P60,batch_A,region_A,control,SRR001,/path/to/ftp,/path/to/fastqs,path/to/reference,TRUE,my_dataset,1,Group_A
sample_B,M,P60,batch_A,region_B,control,SRR002,/path/to/ftp,/path/to/fastqs,path/to/reference,TRUE,my_dataset,2,Group_B
sample_C,F,P60,batch_A,region_A,control,SRR003,/path/to/ftp,/path/to/fastqs,path/to/reference,TRUE,my_dataset,1,Group_A
sample_D,F,P60,batch_A,region_B,control,SRR004,/path/to/ftp,/path/to/fastqs,path/to/reference,TRUE,my_dataset,2,Group_B
```

DATASET_NAME COLUMN (Column 12):
- Used by submit_all_jobs.sh to create Output_dir_<dataset_name>/
- All preprocessing outputs go to: Output_dir_<dataset_name>/Single_cell_preprocessed/
- All downstream outputs go to: Output_dir_<dataset_name>/Single_cell_clustering/
- Allows multiple datasets in same project directory with separate outputs

GROUP_ID COLUMN (Column 13):
- Numeric identifier for each analysis group
- Enables shorthand syntax: --1, --2, --3, etc.
- Samples with same group_id are analyzed together
- Must be positive integers

GROUP_LABEL COLUMN (Column 14):
- Human-readable name for each group
- Used for output directory naming: 11_Downstream_Analysis_<group_label>
- Should not contain spaces (use underscores)
- Examples: Lateral_Ventricle, Fourth_Ventricle, Treatment_Group, Control_Group

------------------------------------------------------------------------------
3.2 fastq_samplesheet.csv - FASTQ FILE LOCATIONS
------------------------------------------------------------------------------

Location: Project root
Used by: Step 1 (Cell Ranger) only
Purpose: Map FASTQ files to Cell Ranger naming convention

WHY NEEDED:
- Cell Ranger requires specific filename format: {Sample}_S{N}_L00{Lane}_R{Read}_001.fastq
- SRA downloads use different naming: SRR*_1.fastq, SRR*_2.fastq
- This samplesheet maps original → Cell Ranger names
- Multiple rows per sample = multiple lanes (auto-combined)

COLUMNS:
  sample_id      : Must match samplesheet.csv
  fastq_dir      : Directory containing FASTQ files
  fastq_prefix   : Original file prefix (e.g., SRR accession)
  lane           : Lane ID (L001, L002, etc.)
  read1_file     : Original R1 filename
  read2_file     : Original R2 filename
  index_file     : Optional I1 filename
  read1_renamed  : Cell Ranger-compatible R1 name (auto-generated if empty)
  read2_renamed  : Cell Ranger-compatible R2 name (auto-generated if empty)

EXAMPLE (with comments - supported):
```
# FASTQ SAMPLESHEET - Maps SRA files to Cell Ranger naming
sample_id,fastq_dir,fastq_prefix,lane,read1_file,read2_file,index_file,read1_renamed,read2_renamed
sample_A,/path/to/fastqs,SRR001,L001,SRR001_2.fastq,SRR001_3.fastq,,sample_A_S1_L001_R1_001.fastq,sample_A_S1_L001_R2_001.fastq
sample_A,/path/to/fastqs,SRR002,L002,SRR002_2.fastq,SRR002_3.fastq,,sample_A_S1_L002_R1_001.fastq,sample_A_S1_L002_R2_001.fastq
```

MULTI-LANE HANDLING:
- Multiple rows with same sample_id = multiple sequencing lanes
- Cell Ranger automatically combines all lanes for a sample
- Script creates symlinks with proper naming in scratch directory

================================================================================
4. PIPELINE STEPS (1-21)
================================================================================

PREPROCESSING (Steps 1-10) - Run as array jobs, outputs to Single_cell_preprocessed/
DOWNSTREAM (Steps 11-21) - Run per group, outputs to Single_cell_clustering/

------------------------------------------------------------------------------
STEP 1: Cell Ranger Count
------------------------------------------------------------------------------
Script:   01_cellranger_count.sh
Input:    Raw FASTQ files (from fastq_samplesheet.csv)
Output:   ${PREPROCESS_DIR}/1_CellRanger_output/{SAMPLE}/outs/
Purpose:  Align reads, call cells, generate count matrix
Time:     ~6-12 hours per sample
Memory:   64GB recommended

Key outputs:
  - filtered_feature_bc_matrix/  (cell-containing barcodes)
  - raw_feature_bc_matrix/       (all barcodes)
  - possorted_genome_bam.bam     (aligned reads)
  - web_summary.html             (QC report)

------------------------------------------------------------------------------
STEP 2: DropletQC - Nuclear Fraction
------------------------------------------------------------------------------
Script:   02_DropletQC_nuclear_fraction.sh
Input:    Cell Ranger BAM + filtered matrix
Output:   ${PREPROCESS_DIR}/2_DropletQC_output/{SAMPLE}/
Purpose:  Identify empty droplets and damaged cells using nuclear fraction
Time:     ~1-2 hours per sample
Memory:   32GB

Key outputs:
  - {SAMPLE}_nuclear_fraction.csv
  - {SAMPLE}_dropletqc_results.csv
  - {SAMPLE}_empty_barcodes.txt
  - {SAMPLE}_damaged_barcodes.txt

Metrics added:
  - dropletqc_nuclear_fraction
  - dropletqc_empty_droplet (boolean)
  - dropletqc_damaged_cell (boolean)
  - dropletqc_status (intact/empty/damaged)

------------------------------------------------------------------------------
STEP 3: QClus - Empty Droplet Detection
------------------------------------------------------------------------------
Script:   03_qClus_empty_droplets.sh
Input:    Cell Ranger + DropletQC results
Output:   ${PREPROCESS_DIR}/3_qClus_empty_droplets/{SAMPLE}/
Purpose:  Additional filtering using QClus algorithm
Time:     ~1-2 hours per sample
Memory:   32GB
Env:      qclus conda environment

Key outputs:
  - {SAMPLE}_qclus_dropletqc_filtered.h5ad (filtered AnnData)
  - {SAMPLE}_qclus_dropletqc_status.csv

Filtering logic:
  - pass_qclus: QClus-based filtering
  - pass_dropletqc: DropletQC-based filtering
  - pass_combined: Both filters pass (used for downstream)

------------------------------------------------------------------------------
STEP 4: VAEDA - Doublet Detection
------------------------------------------------------------------------------
Script:   04_Vaeda_doublet_analysis.sh
Input:    QClus-filtered h5ad
Output:   ${PREPROCESS_DIR}/4_Vaeda_doublet_detection/{SAMPLE}/
Purpose:  Variational autoencoder-based doublet detection
Time:     ~2-4 hours per sample
Memory:   64GB
Env:      vaeda_env conda environment

Key outputs:
  - {SAMPLE}_qClus_dropletqc_vaeda.h5ad
  - {SAMPLE}_vaeda_summary.csv

Metrics added:
  - vaeda_calls (singlet/doublet)
  - vaeda_scores (probability)
  - vaeda_prediction (alias for Step 6 compatibility)

NOTE: Does NOT remove doublets - only adds metadata for consensus voting.

------------------------------------------------------------------------------
STEP 5: AnnData to Seurat Conversion
------------------------------------------------------------------------------
Script:   05_Vaeda_anndata_to_seurat.sh
Input:    VAEDA-annotated h5ad
Output:   ${PREPROCESS_DIR}/5_Seurat_conversion/{SAMPLE}/
Purpose:  Convert to Seurat RDS for R-based analysis
Time:     ~30 minutes per sample
Memory:   32GB

Key outputs:
  - {SAMPLE}_qClus_dropletqc_vaeda.rds

Preserves:
  - All DropletQC metadata
  - All VAEDA metadata
  - Raw counts in counts layer

------------------------------------------------------------------------------
STEP 6: Consensus Doublet Detection
------------------------------------------------------------------------------
Script:   06_scDblFinder_DoubletFinder_analysis.sh
Input:    Seurat RDS from Step 5
Output:   ${PREPROCESS_DIR}/6_Doublet_consensus/{SAMPLE}/
Purpose:  Run additional doublet detectors, consensus voting
Time:     ~1-2 hours per sample
Memory:   32GB

Methods used:
  1. VAEDA (from Step 4)
  2. scDblFinder
  3. DoubletFinder

Consensus logic:
  - doublet_votes: Count of methods calling doublet (0-3)
  - doublet_consensus: TRUE if votes >= 2
  - Cells with consensus=TRUE are REMOVED

Key outputs:
  - {SAMPLE}_doublets_removed.rds (singlets only)
  - {SAMPLE}_doublet_summary.csv

------------------------------------------------------------------------------
STEP 7: scCDC - Contamination Correction
------------------------------------------------------------------------------
Script:   07_scCDC_correction.sh
Input:    Doublet-removed RDS
Output:   ${PREPROCESS_DIR}/7_scCDC_correction/{SAMPLE}/
Purpose:  Detect and correct ambient RNA contamination
Time:     ~1-2 hours per sample
Memory:   32GB

Key outputs:
  - {SAMPLE}_scCDC_corrected.rds    ← INPUT FOR DOWNSTREAM R PIPELINE
  - {SAMPLE}_GCGs.csv (globally contaminating genes)
  - {SAMPLE}_contamination_ratio.csv

Layers added:
  - RNA$scCDC_corrected (contamination-corrected counts)

------------------------------------------------------------------------------
STEP 8: DecontX - Ambient RNA Correction
------------------------------------------------------------------------------
Script:   08_DecontX_correction.sh
Input:    scCDC-corrected RDS files
Output:   ${PREPROCESS_DIR}/8_DecontX_correction/{SAMPLE}/
Purpose:  Additional ambient RNA decontamination using DecontX
Time:     ~1-2 hours per sample
Memory:   32GB

Key outputs:
  - {SAMPLE}_decontx_corrected.rds

------------------------------------------------------------------------------
STEP 9: CHOIR Integration
------------------------------------------------------------------------------
Script:   09_CHOIR_integration.sh
Input:    Corrected RDS files from preceding steps
Output:   ${PREPROCESS_DIR}/9_CHOIR_integration/
Purpose:  Integrate samples using CHOIR for preliminary clustering
Time:     ~2-4 hours
Memory:   64GB
Type:     SINGLE JOB (not array)

Key outputs:
  - Individual/{SAMPLE}/ (per-sample results)
  - {group_label}_only/ (group-specific integrations)
  - All_combined/ (all samples combined)

------------------------------------------------------------------------------
STEP 10: QC Visualization
------------------------------------------------------------------------------
Script:   10_QC_visualization.sh
Input:    All pipeline outputs
Output:   ${PREPROCESS_DIR}/10_QC_visualization/
Purpose:  Generate comprehensive QC plots and summaries
Time:     ~1-2 hours
Memory:   64GB
Type:     SINGLE JOB (not array)

Key outputs:
  - combined_filtering_summary.csv
  - filtering_flow_summary.csv
  - plots/barplot_cells_key_checkpoints.png
  - plots/heatmap_retention_rates.png
  - plots/per_sample/{SAMPLE}_UMAP_QC_*.png

------------------------------------------------------------------------------
STEPS 11-21: Downstream Analysis
------------------------------------------------------------------------------
Script:   11_downstream_analysis.sh
Input:    scCDC-corrected .rds files from ${PREPROCESS_DIR}/7_scCDC_correction/
Output:   ${DOWNSTREAM_DIR}/11_Downstream_Analysis_<group_label>/
Purpose:  Full downstream analysis pipeline
Time:     ~6-24 hours
Memory:   64GB + GPU
Type:     SINGLE JOB, group-specific

GROUP-BASED FILTERING:
-----------------------------------
  --1              → Analyze Group 1 (e.g., Group_A)
  --2              → Analyze Group 2 (e.g., Group_B)
  --group <N>      → Explicit group ID specification
  --group-label <L> → Direct label specification

VENTRICLE FILTERING (BACKWARD COMPATIBLE):
------------------------------------------
  --ventricle LV   → Analyze LV samples only (maps to group internally)
  --ventricle 4V   → Analyze 4V samples only (maps to group internally)
  --ventricle ALL  → Analyze all samples

Module mapping (Steps 11-21 → Modules 0-11, plus optional 2b, 3b, 7b):

  Step 11  → Module 0:  Setup environment
  Step 12  → Module 1:  Load data
  Step 13  → Module 2:  QC validation
           → Module 2b: Imputation (afMF) [OPTIONAL - counts-based]
  Step 14  → Module 3:  Normalization
           → Module 3b: Imputation (ALRA) [OPTIONAL - normalized-data]
  Step 15  → Module 4:  Integration
  Step 16  → Module 5:  CHOIR clustering
  Step 17  → Module 6:  Subclustering (scICE + IDclust)
  Step 18  → Module 7:  Leiden clustering
           → Module 7b: CLTS Re-normalization [OPTIONAL]
  Step 19  → Module 8:  Differential expression
  Step 20  → Module 9:  Gene visualization
  Step 21  → Modules 10+11: Final summary + HTML report

------------------------------------------------------------------------------
MODULE 06: SUBCLUSTERING (scICE + IDclust)
------------------------------------------------------------------------------

Module 06 refines the primary clustering by running one or both subclustering
methods:

1) scICE (Julia-based): Consistent subclustering using scICE.
2) IDclust (R-based): Iterative Differential Clustering (recursive splitting
   with DEG-based validation).

Both methods can be run on:
  - The selected primary clustering (auto / choir / scaura), OR
  - BOTH clustering sources independently (if subclustering_source="both").

Key outputs:
  - Seurat metadata columns:
      scice_subcluster            (if scICE runs; single-source mode)
      idclust_subcluster          (if IDclust runs; single-source mode)

    If subclustering_source="both", additional suffixed columns are created:
      scice_subcluster_choir,  scice_subcluster_scaura
      idclust_subcluster_choir, idclust_subcluster_scaura

    For backward compatibility, unsuffixed columns may be copied from the
    primary source (CHOIR preferred when available).

  - Objects:
      objects/scice_subclustered_object.rds  (updated object after Module 06)
      objects/06_scice_data.RData            (contains scICE + IDclust results)

  - Tables:
      tables/scice_subclustering_summary.csv
      tables/idclust_subclustering_summary.csv

  - Plots:
      plots/clustering/  (combined visualizations; method-specific outputs saved
                         under Module 06 subfolders)

OPTIONAL MODULES:
------------------------

The pipeline supports THREE optional imputation/normalization modules controlled
by params.R settings. See docs/IMPUTATION_METHODS.md for detailed documentation.

Module 2b - Imputation (afMF) - COUNTS-BASED:
  - Controlled by: imputation_method = "afmf" or "both"
  - Requires: afMF Python environment (~/.conda/envs/afMF_SCImputation_env/)
  - Purpose: Recover dropout events using low-rank matrix factorization
  - Operates on: RAW COUNTS (before normalization)
  - Input: merged_object.rds (from Module 02)
  - Output: merged_object_imputed.rds with "imputed" assay
  - Note: use_afmf_for_normalization=FALSE by default (recommended for DE)

Module 3b - Imputation (ALRA) - NORMALIZED-DATA:
  - Controlled by: imputation_method = "alra" or "both"
  - Requires: ALRA R package (CRAN or GitHub)
  - Purpose: Impute dropouts in log-normalized expression data
  - Operates on: LOG-NORMALIZED DATA (after normalization)
  - Input: Normalized object from Module 03
  - Output: *_alra.rds with "ALRA" assay
  - IMPORTANT: Only runs if scran or LogNormalize selected
  - AUTOMATICALLY SKIPPED if SCTransform wins benchmarking
  - Note: use_alra_for_downstream=FALSE by default (recommended for DE)

Module 7b - CLTS Re-normalization:
  - Controlled by: run_clts_renormalization = TRUE
  - Purpose: Remove cell-type-specific composition bias
  - Input: Clustered object (scice/leiden/choir based on clts_clustering_source)
  - Output: <source>_clts_renormalized_object.rds
  - Note: Improves DE analysis accuracy between conditions

IMPUTATION METHOD SELECTION:
-----------------------------------
The unified imputation_method parameter controls both afMF and ALRA:

  imputation_method = "none"  → No imputation (default behavior)
  imputation_method = "afmf"  → Only afMF (Module 2b)
  imputation_method = "alra"  → Only ALRA (Module 3b, if compatible)
  imputation_method = "both"  → Both methods (2b always, 3b if compatible)

ALRA COMPATIBILITY:
  - LogNormalize: YES (log-normalized data)
  - scran: YES (log-normalized data)
  - SCTransform: NO (uses Pearson residuals, not log-normalized)

If imputation_method="alra" or "both" but SCTransform wins normalization
benchmarking, Module 3b is automatically skipped with a warning.

To run optional modules explicitly:
  ./submit_all_jobs.sh --downstream --1 --modules 2b,3      # afMF + normalization
  ./submit_all_jobs.sh --downstream --1 --modules 3b,4      # ALRA + integration
  ./submit_all_jobs.sh --downstream --1 --modules 7b,8      # CLTS + DE

================================================================================
5. SLURM SCRIPTS - MODULAR DESIGN
================================================================================

All preprocessing scripts (01-10) follow identical modular pattern:

------------------------------------------------------------------------------
5.1 Configuration Section
------------------------------------------------------------------------------

```bash
# Base directory from environment or default
if [[ -n "${PROJECT_ROOT:-}" ]]; then
    BASE_DIR="${PROJECT_ROOT}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    BASE_DIR="${SLURM_SUBMIT_DIR}"
else
    BASE_DIR="/path/to/project"
fi

# Samplesheet
SAMPLESHEET="${SAMPLESHEET:-${BASE_DIR}/samplesheet.csv}"

# Dataset-specific output directory (from environment or samplesheet column 12)
if [[ -n "${PREPROCESS_DIR:-}" ]]; then
    OUTPUT_BASE="${PREPROCESS_DIR}"
else
    DATASET_NAME=$(tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f12)
    OUTPUT_BASE="${BASE_DIR}/Output_dir_${DATASET_NAME}/Single_cell_preprocessed"
fi
```

------------------------------------------------------------------------------
5.2 Helper Functions (identical in all scripts)
------------------------------------------------------------------------------

```bash
get_col_index()        # Get column number by name from CSV header
get_unique_samples()   # Extract unique sample IDs from samplesheet
get_sample_count()     # Count total samples
get_sample_by_index()  # Get nth sample (for array task mapping)
get_sample_field()     # Get metadata field value for specific sample
```

GROUP HELPER FUNCTIONS:
```bash
get_samples_by_group()     # Get samples filtered by group_id
get_group_label_by_id()    # Map group_id → group_label
get_group_id_by_ventricle()# Map ventricle → group_id (backward compat)
list_available_groups()    # Display all groups from samplesheet
```

------------------------------------------------------------------------------
5.3 Sample Selection (array job mapping)
------------------------------------------------------------------------------

```bash
N_SAMPLES=$(get_sample_count)
TASK_IDX=${SLURM_ARRAY_TASK_ID:-1}
SAMPLE=$(get_sample_by_index $TASK_IDX)

# Validate
if [[ $TASK_IDX -gt $N_SAMPLES ]]; then
    error_exit "Array task $TASK_IDX exceeds sample count $N_SAMPLES"
fi
```

------------------------------------------------------------------------------
5.4 Metadata Extraction
------------------------------------------------------------------------------

```bash
SAMPLE_SEX=$(get_sample_field "$SAMPLE" "sex")
SAMPLE_BATCH=$(get_sample_field "$SAMPLE" "batch")
SAMPLE_VENTRICLE=$(get_sample_field "$SAMPLE" "ventricle")
SAMPLE_GROUP_ID=$(get_sample_field "$SAMPLE" "group_id")
SAMPLE_GROUP_LABEL=$(get_sample_field "$SAMPLE" "group_label")
```

------------------------------------------------------------------------------
5.5 SBATCH Headers (relative paths, no array specification)
------------------------------------------------------------------------------

```bash
#SBATCH --job-name=STEP_NAME
#SBATCH --output=logs/step/step_%a_%A.out   # %a=array task, %A=job ID
#SBATCH --error=logs/step/step_%a_%A.err
# NOTE: --array=1-N is NOT in SBATCH header - set by submit_all_jobs.sh
```

================================================================================
6. HELPER FUNCTIONS - DETAILED
================================================================================

These functions are embedded in each script for portability:

------------------------------------------------------------------------------
get_col_index(file, col_name)
------------------------------------------------------------------------------
Purpose: Find column number for a given column name
Usage:   col=$(get_col_index "$SAMPLESHEET" "sample_id")
Returns: Integer column number (1-based)

```bash
get_col_index() {
    local file=$1
    local col_name=$2
    head -1 "$file" | tr ',' '\n' | grep -n "^${col_name}$" | cut -d: -f1
}
```

------------------------------------------------------------------------------
get_unique_samples()
------------------------------------------------------------------------------
Purpose: Get list of all unique sample IDs
Returns: Newline-separated sample IDs, sorted

```bash
get_unique_samples() {
    local col=$(get_col_index "$SAMPLESHEET" "sample_id")
    tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | \
        cut -d',' -f"$col" | sort -u
}
```

------------------------------------------------------------------------------
get_sample_count()
------------------------------------------------------------------------------
Purpose: Count total number of samples
Returns: Integer

```bash
get_sample_count() {
    get_unique_samples | wc -l
}
```

------------------------------------------------------------------------------
get_sample_by_index(idx)
------------------------------------------------------------------------------
Purpose: Map array task index to sample name
Usage:   SAMPLE=$(get_sample_by_index $SLURM_ARRAY_TASK_ID)

```bash
get_sample_by_index() {
    local idx=$1
    get_unique_samples | sed -n "${idx}p"
}
```

------------------------------------------------------------------------------
get_sample_field(sample, field)
------------------------------------------------------------------------------
Purpose: Get metadata value for a sample
Usage:   SEX=$(get_sample_field "sample_A" "sex")  # Returns "M"

```bash
get_sample_field() {
    local sample=$1
    local field=$2
    local col=$(get_col_index "$SAMPLESHEET" "$field")
    if [[ -n "$col" ]]; then
        grep "^${sample}," "$SAMPLESHEET" | head -1 | cut -d',' -f"$col"
    fi
}
```

------------------------------------------------------------------------------
get_samples_by_group(group_id)
------------------------------------------------------------------------------
Purpose: Get samples filtered by group_id
Usage:   SAMPLES=$(get_samples_by_group "1")

```bash
get_samples_by_group() {
    local group_id=$1
    local col_id=$(get_col_index "$SAMPLESHEET" "sample_id")
    local col_gid=$(get_col_index "$SAMPLESHEET" "group_id")
    tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | \
        awk -F',' -v col="$col_gid" -v gid="$group_id" -v id="$col_id" \
        '$col == gid {print $id}'
}
```

------------------------------------------------------------------------------
get_group_label_by_id(group_id)
------------------------------------------------------------------------------
Purpose: Get group_label for a given group_id
Usage:   LABEL=$(get_group_label_by_id "1")  # Returns "Group_A"

```bash
get_group_label_by_id() {
    local group_id=$1
    local col_gid=$(get_col_index "$SAMPLESHEET" "group_id")
    local col_glabel=$(get_col_index "$SAMPLESHEET" "group_label")
    tail -n +2 "$SAMPLESHEET" | grep -v '^#' | grep -v '^$' | \
        awk -F',' -v col="$col_gid" -v gid="$group_id" -v lc="$col_glabel" \
        '$col == gid {print $lc; exit}'
}
```

================================================================================
7. SUBMIT_ALL_JOBS.SH - MASTER CONTROLLER
================================================================================

The master script handles:
- Samplesheet validation
- Dynamic array job sizing
- Job dependencies
- Environment variable passing (including GROUP_ID, GROUP_LABEL)
- Dataset-specific output directory creation
- Group listing and selection
- Backward compatibility with --ventricle syntax
- Dry-run mode
- Optional module support (2b, 3b, 7b)

------------------------------------------------------------------------------
7.1 Key Features
------------------------------------------------------------------------------

VALIDATION:
```bash
./submit_all_jobs.sh --validate
```
Checks: samplesheets exist, required columns present, FASTQ files exist

DRY RUN:
```bash
./submit_all_jobs.sh --preprocess --dry-run
```
Shows what would be submitted without actually submitting

PREPROCESSING (Steps 1-10):
```bash
./submit_all_jobs.sh --preprocess
```

SPECIFIC STEPS:
```bash
./submit_all_jobs.sh --only 1      # Just Cell Ranger
./submit_all_jobs.sh --only 2,3    # DropletQC and QClus
./submit_all_jobs.sh --from 4 --to 7   # Steps 4-7
```

LIST AVAILABLE GROUPS:
```bash
./submit_all_jobs.sh --list-groups
```
Output:
```
Available groups (use --<group_id> or --group <group_id>):
  --1   →  Group_A              (2 samples: sample_A sample_C)
  --2   →  Group_B              (2 samples: sample_B sample_D)

Backward compatible (--ventricle):
  --ventricle LV  →  maps to group 1
  --ventricle 4V  →  maps to group 2
```

DOWNSTREAM WITH GROUP SHORTHAND:
```bash
./submit_all_jobs.sh --downstream --1          # Group 1 (Group_A)
./submit_all_jobs.sh --downstream --2          # Group 2 (Group_B)
./submit_all_jobs.sh --downstream --group 1    # Explicit group ID
./submit_all_jobs.sh --downstream --1 --modules 3,4,7  # With specific modules
```

DOWNSTREAM WITH OPTIONAL MODULES:
```bash
./submit_all_jobs.sh --downstream --1 --modules 2b,3       # afMF + normalization
./submit_all_jobs.sh --downstream --1 --modules 3,3b,4     # Norm + ALRA + integration
./submit_all_jobs.sh --downstream --1 --modules 7b,8       # CLTS + DE
./submit_all_jobs.sh --downstream --1 --modules 0,1,2,2b,3,3b,4,5,6,7,7b,8,9,10,11  # Full with optionals
```

DOWNSTREAM WITH VENTRICLE (BACKWARD COMPATIBLE):
```bash
./submit_all_jobs.sh --downstream --ventricle LV   # Still works
./submit_all_jobs.sh --downstream --ventricle 4V   # Still works
```

------------------------------------------------------------------------------
7.2 Environment Variable Exports
------------------------------------------------------------------------------

submit_all_jobs.sh exports these environment variables to all jobs:

```bash
# Read dataset name from samplesheet column 12
DATASET_NAME=$(tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f12)

# Construct dataset-specific output directories
export PROJECT_ROOT="${BASE_DIR}"
export DATASET_NAME="${DATASET_NAME}"
export PREPROCESS_DIR="${BASE_DIR}/Output_dir_${DATASET_NAME}/Single_cell_preprocessed"
export DOWNSTREAM_DIR="${BASE_DIR}/Output_dir_${DATASET_NAME}/Single_cell_clustering"
export SAMPLESHEET="${SAMPLESHEET}"
export PIPELINE_DIR="${BASE_DIR}/Scripts/scrnaseq_pipeline"

# Group-based filtering
export GROUP_ID="${GROUP_ID}"           # Numeric group ID (1, 2, ...)
export GROUP_LABEL="${GROUP_LABEL}"     # Human-readable label

# Backward compatible
export VENTRICLE_FILTER="${VENTRICLE}"  # LV, 4V, or ALL

# Create output directories
mkdir -p "${PREPROCESS_DIR}"
mkdir -p "${DOWNSTREAM_DIR}"
```

------------------------------------------------------------------------------
7.3 Array Job Submission
------------------------------------------------------------------------------

```bash
# Auto-detect sample count
N_SAMPLES=$(tail -n +2 "$SAMPLESHEET" | grep -v '^$' | wc -l)

# Submit array job with all environment variables
sbatch --array=1-${N_SAMPLES} \
       --export=ALL,SAMPLESHEET=${SAMPLESHEET},PROJECT_ROOT=${BASE_DIR},\
PREPROCESS_DIR=${PREPROCESS_DIR},DOWNSTREAM_DIR=${DOWNSTREAM_DIR},\
DATASET_NAME=${DATASET_NAME} \
       Scripts/Slurm_scripts/01_cellranger_count.sh
```

------------------------------------------------------------------------------
7.4 Job Dependencies
------------------------------------------------------------------------------

```bash
# Step 2 depends on Step 1
JOB_STEP1=$(sbatch --array=1-${N} ... 01_cellranger_count.sh | awk '{print $4}')
JOB_STEP2=$(sbatch --array=1-${N} --dependency=afterok:${JOB_STEP1} ... 02_*.sh)

# Step 10 (single job) depends on all Step 8 array tasks
JOB_STEP10=$(sbatch --dependency=afterok:${JOB_STEP8} ... 10_QC_visualization.sh)

# Downstream depends on preprocessing completion
JOB_DOWNSTREAM=$(sbatch --dependency=afterok:${JOB_STEP10} \
    --export=ALL,GROUP_ID=1,GROUP_LABEL=Group_A,... \
    11_downstream_analysis.sh)
```

------------------------------------------------------------------------------
7.5 Group ID Resolution
------------------------------------------------------------------------------

```bash
# Parse shorthand syntax (--1, --2, etc.)
if [[ "$arg" =~ ^--([0-9]+)$ ]]; then
    GROUP_ID="${BASH_REMATCH[1]}"
    GROUP_LABEL=$(get_group_label_by_id "$GROUP_ID")
fi

# Map --ventricle to group for backward compatibility
if [[ -n "$VENTRICLE" && -z "$GROUP_ID" ]]; then
    GROUP_ID=$(get_group_id_by_ventricle "$VENTRICLE")
    GROUP_LABEL=$(get_group_label_by_id "$GROUP_ID")
    echo "Note: Mapped --ventricle $VENTRICLE to --group $GROUP_ID ($GROUP_LABEL)"
fi
```

------------------------------------------------------------------------------
7.6 Optional Module Display
------------------------------------------------------------------------------

When running downstream analysis, optional modules are displayed in the
execution plan:

```
Modules to run:
   0 - Environment Setup
   1 - Load Data
   2 - QC Validation
  (2b - Imputation (afMF)) [if enabled in params.R]
   3 - Normalization
  (3b - Imputation (ALRA)) [if enabled in params.R]
   4 - Integration
   5 - CHOIR Clustering
   6 - Subclustering (scICE + IDclust)
   7 - Leiden Clustering
  (7b - CLTS Re-normalization) [if enabled in params.R]
   8 - Differential Expression
   9 - Gene Visualization
  10 - Final Summary
  11 - HTML Report

Note: Modules 2b, 3b, and 7b run automatically if enabled in params.R
```

================================================================================
8. DOWNSTREAM ANALYSIS (STEPS 11-21)
================================================================================

------------------------------------------------------------------------------
8.1 Group-Based Analysis
------------------------------------------------------------------------------

The downstream pipeline filters samples by group_id and uses group_label for
output directory naming:

```bash
# Group 1
./submit_all_jobs.sh --downstream --1
# OR: sbatch 11_downstream_analysis.sh --group-label Group_A
# Processes: sample_A, sample_C
# Input: ${PREPROCESS_DIR}/7_scCDC_correction/{sample_A,sample_C}/
# Output: ${DOWNSTREAM_DIR}/11_Downstream_Analysis_Group_A/

# Group 2
./submit_all_jobs.sh --downstream --2
# OR: sbatch 11_downstream_analysis.sh --group-label Group_B
# Processes: sample_B, sample_D
# Input: ${PREPROCESS_DIR}/7_scCDC_correction/{sample_B,sample_D}/
# Output: ${DOWNSTREAM_DIR}/11_Downstream_Analysis_Group_B/
```

------------------------------------------------------------------------------
8.2 Ventricle-Specific Analysis (BACKWARD COMPATIBLE)
------------------------------------------------------------------------------

Legacy --ventricle syntax is automatically mapped to groups:

```bash
# Maps to corresponding group internally
sbatch 11_downstream_analysis.sh --ventricle LV
# Output: ${DOWNSTREAM_DIR}/11_Downstream_Analysis_LV/

sbatch 11_downstream_analysis.sh --ventricle 4V
# Output: ${DOWNSTREAM_DIR}/11_Downstream_Analysis_4V/

# All samples (no group filter)
sbatch 11_downstream_analysis.sh --ventricle ALL
# Output: ${DOWNSTREAM_DIR}/11_Downstream_Analysis_ALL/
```

------------------------------------------------------------------------------
8.3 params.R Configuration
------------------------------------------------------------------------------

The downstream R pipeline reads these environment variables:

```R
# Set by bash script, read in params.R
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
PREPROCESS_DIR <- Sys.getenv("PREPROCESS_DIR")
DOWNSTREAM_DIR <- Sys.getenv("DOWNSTREAM_DIR")
DATASET_NAME <- Sys.getenv("DATASET_NAME")

# Group-based filtering
GROUP_ID <- Sys.getenv("GROUP_ID")           # "1", "2", etc.
GROUP_LABEL <- Sys.getenv("GROUP_LABEL")     # "Group_A", etc.

# Backward compatible
VENTRICLE_FILTER <- Sys.getenv("VENTRICLE_FILTER")

# Input directory (scCDC-corrected files)
input_dir <- file.path(PREPROCESS_DIR, "7_scCDC_correction")

# Determine output label (GROUP_LABEL takes priority)
if (GROUP_LABEL != "") {
    analysis_label <- GROUP_LABEL
} else if (VENTRICLE_FILTER != "") {
    analysis_label <- VENTRICLE_FILTER
} else {
    analysis_label <- "All"
}

# Output directory
out_root <- file.path(DOWNSTREAM_DIR,
                      paste0("11_Downstream_Analysis_", analysis_label))

# Filter samples by group_id (priority) or ventricle (fallback)
if (GROUP_ID != "" && "group_id" %in% colnames(samples)) {
    samples <- samples[samples$group_id == GROUP_ID, ]
} else if (VENTRICLE_FILTER != "" && VENTRICLE_FILTER != "ALL") {
    samples <- samples[samples$ventricle == VENTRICLE_FILTER, ]
}
```

------------------------------------------------------------------------------
8.4 Imputation Parameters
------------------------------------------------------------------------------

params.R contains unified imputation settings:

```R
# ==============================================================================
# IMPUTATION SETTINGS (UNIFIED)
# ==============================================================================

# Main control - determines which imputation modules run
imputation_method = "both",       # Options: "none", "afmf", "alra", "both"

# ==============================================================================
# MODULE 02b: afMF IMPUTATION (COUNTS-BASED)
# ==============================================================================
use_afmf_for_normalization = FALSE,  # Use imputed counts in Module 03?
afmf_python = file.path(Sys.getenv("HOME"),
                        ".conda/envs/afMF_SCImputation_env/bin/python"),
afmf_max_iter = 100,              # Maximum iterations
afmf_tol = 1e-5,                  # Convergence tolerance
afmf_min_cells_expressing = 10,   # Gene filter threshold

# ==============================================================================
# MODULE 03b: ALRA IMPUTATION (NORMALIZED-DATA)
# ==============================================================================
use_alra_for_downstream = FALSE,  # Use ALRA data in Module 04+?
alra_k = NULL,                    # SVD rank (NULL = auto-detect)
alra_q = 10,                      # Power iterations for randomized SVD
alra_quantile_prob = 0.001,       # Quantile for thresholding
alra_compatible_methods = c("LogNormalize", "scran"),  # Compatible normalizations

# ==============================================================================
# MODULE 07b: CLTS RE-NORMALIZATION
# ==============================================================================
run_clts_renormalization = FALSE, # Set TRUE to enable Module 07b
clts_clustering_source = "scice", # Which clustering to use: "scice", "leiden", "choir"
clts_min_cells_per_cluster = 50,  # Minimum cluster size
clts_baseline_sample = NULL,      # Reference sample (optional)
clts_run_benchmark = TRUE,        # Run benchmark comparison
```

IMPUTATION EXECUTION MATRIX:

| imputation_method | Module 02b (afMF) | Module 03b (ALRA) |
|-------------------|-------------------|-------------------|
| "none"            | SKIP              | SKIP              |
| "afmf"            | RUN               | SKIP              |
| "alra"            | SKIP              | RUN if compatible |
| "both"            | RUN               | RUN if compatible |

* "compatible" means selected_normalization_method ∈ {LogNormalize, scran}

IMPORTANT NOTES:

1. afMF Imputation (Module 02b):
   - Requires afMF Python environment to be installed
   - use_afmf_for_normalization=FALSE recommended for DE analysis
   - Imputed counts useful for visualization and trajectory analysis
   - Module 00 verifies afMF availability and sets has_afmf flag

2. ALRA Imputation (Module 03b):
   - Requires ALRA R package (install.packages("ALRA"))
   - AUTOMATICALLY SKIPPED if SCTransform wins normalization benchmarking
   - use_alra_for_downstream=FALSE recommended for DE analysis
   - Creates "ALRA" assay with imputed log-normalized data
   - Module 00 verifies ALRA availability and sets has_alra flag

3. CLTS Re-normalization (Module 07b):
   - Requires clustering to be completed first (Module 05, 06, or 07)
   - Improves DE analysis by removing composition bias
   - Multiple DE object sources can be analyzed (de_object_sources param)

------------------------------------------------------------------------------
8.5 Path Resolution Priority
------------------------------------------------------------------------------

params.R resolves paths in this order:

1. Environment variables (PREPROCESS_DIR, DOWNSTREAM_DIR) - highest priority
2. Dataset name derivation from DATASET_NAME or samplesheet column 12
3. Legacy fallback paths - lowest priority (for backward compatibility)

```R
# Priority 1: Environment variables
if (env_preprocess_dir != "" && dir.exists(env_preprocess_dir)) {
    input_dir <- file.path(env_preprocess_dir, "7_scCDC_correction")
    files_in_subdirectories <- TRUE
}

# Priority 2: Dataset name derivation
else if (dataset_name != "" && dataset_name != "default") {
    input_dir <- file.path(project_root,
                           paste0("Output_dir_", dataset_name),
                           "Single_cell_preprocessed",
                           "7_scCDC_correction")
}

# Priority 3: Legacy fallback
else {
    input_dir <- "/legacy/hardcoded/path"
}
```

------------------------------------------------------------------------------
8.6 Module Selection
------------------------------------------------------------------------------

Run specific modules:
```bash
# Only clustering modules
./submit_all_jobs.sh --downstream --1 --modules 5,6,7

# Translate step numbers to modules
./submit_all_jobs.sh --downstream --2 --steps 16,17,18

# Direct sbatch with modules
sbatch 11_downstream_analysis.sh --group-label Group_A --modules 5,6,7

# Include optional modules explicitly
./submit_all_jobs.sh --downstream --1 --modules 2b,3        # afMF + norm
./submit_all_jobs.sh --downstream --1 --modules 3,3b,4      # Norm + ALRA + integration
./submit_all_jobs.sh --downstream --1 --modules 7b,8        # CLTS + DE
./submit_all_jobs.sh --downstream --1 --modules 2,2b,3,3b,7,7b,8  # Selected with optionals
```

------------------------------------------------------------------------------
8.7 Module Execution Order
------------------------------------------------------------------------------

Standard pipeline:
  0 → 1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 → 9 → 10 → 11

With optional modules:
  0 → 1 → 2 → [2b] → 3 → [3b] → 4 → 5 → 6 → 7 → [7b] → 8 → 9 → 10 → 11

Module dependencies:
  Module 02b requires: Module 02 (merged_object.rds), afMF environment
  Module 03b requires: Module 03 (normalized object), ALRA package,
                       normalization method ∈ {LogNormalize, scran}
  Module 07b requires: Module 06 or 07 (clustered object)

================================================================================
9. ADAPTING FOR NEW DATASETS
================================================================================

To use this pipeline for a different dataset:

------------------------------------------------------------------------------
Step 1: Create samplesheet.csv with group columns
------------------------------------------------------------------------------

```csv
sample_id,sex,age,batch,ventricle,condition,SRA_accession,ftp_path,local_fastq_dir,genome,include,dataset_name,group_id,group_label
sample_A,M,P60,batch_A,region_A,treatment,SRR123,,/path/to/fastqs,path/to/reference,TRUE,my_dataset,1,Treatment_Group
sample_B,F,P60,batch_A,region_A,control,SRR124,,/path/to/fastqs,path/to/reference,TRUE,my_dataset,2,Control_Group
sample_C,M,P60,batch_A,region_A,treatment,SRR125,,/path/to/fastqs,path/to/reference,TRUE,my_dataset,1,Treatment_Group
sample_D,F,P60,batch_A,region_A,control,SRR126,,/path/to/fastqs,path/to/reference,TRUE,my_dataset,2,Control_Group
```

Required columns: sample_id, sex, batch
Recommended: group_id, group_label, dataset_name (column 12)
Optional: ventricle (for backward compatibility)

The dataset_name value determines output directory:
  - Output_dir_my_dataset/Single_cell_preprocessed/
  - Output_dir_my_dataset/Single_cell_clustering/

The group_id and group_label enable:
  - ./submit_all_jobs.sh --downstream --1   (Treatment_Group)
  - ./submit_all_jobs.sh --downstream --2   (Control_Group)
  - Output: 11_Downstream_Analysis_Treatment_Group/
  - Output: 11_Downstream_Analysis_Control_Group/

------------------------------------------------------------------------------
Step 2: Create fastq_samplesheet.csv
------------------------------------------------------------------------------

For each sample, list FASTQ files:
```csv
sample_id,fastq_dir,fastq_prefix,lane,read1_file,read2_file,index_file,read1_renamed,read2_renamed
sample_A,/path/to/fastqs,SRR123,L001,SRR123_R1.fastq.gz,SRR123_R2.fastq.gz,,sample_A_S1_L001_R1_001.fastq.gz,sample_A_S1_L001_R2_001.fastq.gz
```

------------------------------------------------------------------------------
Step 3: Update Cell Ranger Reference
------------------------------------------------------------------------------

Edit samplesheet.csv genome column or update default in scripts.

For human: /path/to/refdata-gex-GRCh38-2024-A
For mouse: /path/to/refdata-gex-GRCm39-2024-A

------------------------------------------------------------------------------
Step 4: Define Groups
------------------------------------------------------------------------------

Decide how to group your samples for downstream analysis:

EXAMPLES:
- By treatment:    group_id=1 → Treatment, group_id=2 → Control
- By timepoint:    group_id=1 → Day0, group_id=2 → Day7, group_id=3 → Day14
- By region:       group_id=1 → Cortex, group_id=2 → Hippocampus
- By genotype:     group_id=1 → Wildtype, group_id=2 → Knockout

Add group_id and group_label columns to samplesheet.csv accordingly.

------------------------------------------------------------------------------
Step 5: Configure Optional Modules
------------------------------------------------------------------------------

Edit params.R to enable optional modules if desired:

```R
# ==============================================================================
# MODULE 06: SUBCLUSTERING (scICE + IDclust)
# ==============================================================================
# Choose which subclustering methods to run (can run both)
subclustering_methods = c("scice", "idclust"),   # Options: "scice", "idclust"
subclustering_source = "auto",                  # "auto", "choir", "scaura", "both"

# scICE toggle (Julia/scICE required)
run_scice_subclustering = TRUE,

# IDclust toggle + thresholds (R/IDclust required)
run_idclust_subclustering = TRUE,
idclust_target_clusters = NULL,      # NULL = all clusters meeting min_cells
idclust_logFC_th = log2(1.5),
idclust_qval_th = 0.01,
idclust_min_DEGs = 5,
idclust_max_depth = 10,
idclust_min_frac_assigned = 0.1,
idclust_n_dims = 50,
idclust_starting_resolution = 0.1,
idclust_resolution = 0.8,
idclust_starting_k = 100,
idclust_k = 100,
idclust_min_cells = 100,
idclust_plotting = TRUE,

# Enable imputation (unified parameter)
imputation_method = "both",       # "none", "afmf", "alra", or "both"

# afMF settings (Module 02b)
use_afmf_for_normalization = FALSE,   # Keep FALSE for DE analysis

# ALRA settings (Module 03b)
use_alra_for_downstream = FALSE,      # Keep FALSE for DE analysis
alra_k = NULL,                        # Auto-detect SVD rank

# Enable CLTS re-normalization (Module 07b)
run_clts_renormalization = TRUE,
clts_clustering_source = "scice",
```

Ensure required environments are installed:
- afMF: ~/.conda/envs/afMF_SCImputation_env/
- ALRA: install.packages("ALRA") in R
- Julia/scICE: For Module 06 (scICE)
- R package IDclust: For Module 06 (IDclust)

------------------------------------------------------------------------------
Step 6: Validate and Run
------------------------------------------------------------------------------

```bash
cd /path/to/new/project

# Validate setup
bash validate_pipeline.sh

# List available groups
./submit_all_jobs.sh --list-groups

# Dry run
./submit_all_jobs.sh --preprocess --dry-run

# Full run
./submit_all_jobs.sh --preprocess

# Downstream (after preprocessing completes)
./submit_all_jobs.sh --downstream --1    # First group
./submit_all_jobs.sh --downstream --2    # Second group

# OR using legacy syntax (still works)
./submit_all_jobs.sh --downstream --ventricle LV
./submit_all_jobs.sh --downstream --ventricle 4V

# With optional modules
./submit_all_jobs.sh --downstream --1 --modules 2b,3,3b,4,5,6,7,7b,8,9,10,11
```

================================================================================
10. TROUBLESHOOTING
================================================================================

------------------------------------------------------------------------------
10.1 Common Issues
------------------------------------------------------------------------------

ISSUE: "Array task X exceeds sample count N"
CAUSE: Array size doesn't match samplesheet
FIX:   Let submit_all_jobs.sh set array size automatically

ISSUE: "Samplesheet not found"
CAUSE: PROJECT_ROOT or SAMPLESHEET not set
FIX:   Run via submit_all_jobs.sh or export variables manually

ISSUE: "Column 'X' MISSING" in validation
CAUSE: Samplesheet missing required column
FIX:   Add column to samplesheet.csv

ISSUE: Cell Ranger fails with filename error
CAUSE: FASTQ files don't match Cell Ranger naming convention
FIX:   Verify fastq_samplesheet.csv read1_renamed/read2_renamed columns

ISSUE: "No samples found for ventricle X"
CAUSE: Ventricle column missing or values don't match
FIX:   Check samplesheet.csv ventricle column values

ISSUE: "No samples found for group X"
CAUSE: Group ID doesn't exist in samplesheet
FIX:
  - Run ./submit_all_jobs.sh --list-groups to see available groups
  - Check samplesheet.csv group_id column values
  - Ensure group_id column exists (column 13)

ISSUE: "GROUP_LABEL not set" or empty output directory name
CAUSE: group_label column missing or empty
FIX:
  - Add group_label column (column 14) to samplesheet.csv
  - Ensure all rows have non-empty group_label values

ISSUE: "PREPROCESS_DIR not found" or "Input file not found"
CAUSE: Environment variable not set or wrong path
FIX:
  - Run via submit_all_jobs.sh (sets automatically), OR
  - Check DATASET_NAME in samplesheet column 12
  - Verify Output_dir_<dataset_name>/ exists

ISSUE: "Output going to wrong directory"
CAUSE: DATASET_NAME or GROUP_LABEL not set correctly
FIX:
  - Check samplesheet column 12 has correct dataset_name value
  - Check samplesheet column 14 has correct group_label value

ISSUE: "--ventricle LV not working after upgrade"
CAUSE: Backward compatibility issue (should not happen)
FIX:
  - Ensure ventricle column still exists in samplesheet
  - The --ventricle syntax should automatically map to groups

ISSUE: "afMF imputation failed" or "has_afmf is FALSE"
CAUSE: afMF Python environment not available
FIX:
  - Check afmf_python path in params.R
  - Verify environment exists: ls ~/.conda/envs/afMF_SCImputation_env/
  - Install afMF:
    ```bash
    conda activate afMF_SCImputation_env
    cd path/to/SCImputation/afMF
    pip install .
    ```
  - Module 00 checks availability and sets has_afmf flag

ISSUE: "Module 02b skipped but imputation_method includes afmf"
CAUSE: afMF environment check failed in Module 00
FIX:
  - Review Module 00 output for afMF status
  - Check has_afmf value in pipeline_environment.RData
  - Ensure afMF Python environment is properly installed

ISSUE: "ALRA imputation failed" or "has_alra is FALSE"
CAUSE: ALRA R package not installed
FIX:
  - Install ALRA: install.packages("ALRA") or
    devtools::install_github("KlugerLab/ALRA")
  - Module 00 checks availability and sets has_alra flag

ISSUE: "Module 03b skipped but imputation_method includes alra"
CAUSE: Either ALRA not installed OR SCTransform was selected
FIX:
  - If ALRA not installed: install.packages("ALRA")
  - If SCTransform selected: This is expected! ALRA is incompatible with
    SCTransform. To use ALRA, disable SCTransform:
    ```R
    run_sctransform = FALSE
    ```
  - Or force a specific normalization method:
    ```R
    integration_normalization_method = "scran"  # or "LogNormalize"
    ```

ISSUE: "ALRA incompatible with SCTransform" warning
CAUSE: SCTransform won normalization benchmarking
NOTE: This is correct behavior - ALRA cannot work with SCTransform's
      Pearson residuals. Module 03b is automatically skipped.
FIX:  If you need ALRA, use scran or LogNormalize normalization

ISSUE: "CLTS failed - no cluster column found"
CAUSE: Clustering source object missing expected column
FIX:
  - Check clts_clustering_source setting in params.R
  - Verify clustering was run: check for scice_subcluster or idclust_subcluster column
  - Try different source: clts_clustering_source = "leiden"

ISSUE: "DE results differ between scice and clts sources"
CAUSE: This is expected - CLTS removes composition bias
NOTE: This is not an error; CLTS-normalized results may show
      fewer false positives from composition effects

ISSUE: "Biological conservation metrics: NOT AVAILABLE"
CAUSE: No cell type annotations found in metadata
FIX:
  - Set celltype_column in params.R to the correct metadata column name
  - Auto-detection searches for common names (cell_type, celltype,
    CellType, annotation, etc.)
  - The pipeline falls back to batch-only scoring if no annotations found

ISSUE: "celltype_asw/nmi/ari all NA"
CAUSE: Cell type column exists but has <2 unique types or too few cells
FIX:
  - Check the annotation column has meaningful cell type labels
  - Ensure at least 2 distinct cell types with >=10 cells each
  - Try a different annotation granularity

------------------------------------------------------------------------------
10.2 Checking Job Status
------------------------------------------------------------------------------

```bash
# All jobs
squeue -u $USER

# Specific job
sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed

# Array job tasks
sacct -j JOB_ID --format=JobID%20,State,ExitCode
```

------------------------------------------------------------------------------
10.3 Log Files
------------------------------------------------------------------------------

SLURM logs: logs/{step}/{step}_%a_%A.out/err
  - %a = array task ID (1, 2, 3...)
  - %A = job ID

Per-step README: {OUTPUT_DIR}/README.txt
  - Processing status per sample
  - Summary statistics

------------------------------------------------------------------------------
10.4 Verifying Environment Variables
------------------------------------------------------------------------------

```bash
# In bash script
echo "PROJECT_ROOT: $PROJECT_ROOT"
echo "PREPROCESS_DIR: $PREPROCESS_DIR"
echo "DOWNSTREAM_DIR: $DOWNSTREAM_DIR"
echo "DATASET_NAME: $DATASET_NAME"
echo "GROUP_ID: $GROUP_ID"
echo "GROUP_LABEL: $GROUP_LABEL"
echo "VENTRICLE_FILTER: $VENTRICLE_FILTER"

# In R
Sys.getenv("PROJECT_ROOT")
Sys.getenv("PREPROCESS_DIR")
Sys.getenv("DOWNSTREAM_DIR")
Sys.getenv("DATASET_NAME")
Sys.getenv("GROUP_ID")
Sys.getenv("GROUP_LABEL")
Sys.getenv("VENTRICLE_FILTER")
```

------------------------------------------------------------------------------
10.5 Verifying Optional Module Status
------------------------------------------------------------------------------

```R
# In R, after Module 00
load("objects/pipeline_environment.RData")

# Check afMF availability
cat("afMF available:", has_afmf, "\n")
cat("afMF Python:", afmf_python, "\n")

# Check ALRA availability
cat("ALRA available:", has_alra, "\n")

# Check module settings
cat("imputation_method:", params$imputation_method, "\n")
cat("run_scice_subclustering:", params$run_scice_subclustering, "\n")
cat("run_idclust_subclustering:", params$run_idclust_subclustering, "\n")
cat("subclustering_methods:", paste(params$subclustering_methods, collapse = ", "), "\n")
cat("subclustering_source:", params$subclustering_source, "\n")
cat("use_afmf_for_normalization:", params$use_afmf_for_normalization, "\n")
cat("use_alra_for_downstream:", params$use_alra_for_downstream, "\n")
cat("run_clts_renormalization:", params$run_clts_renormalization, "\n")
cat("clts_clustering_source:", params$clts_clustering_source, "\n")
```

NORMALIZATION METHODS:
-----------------------------
Module 03 performs FOUR normalization methods for benchmarking:

1. SCTransform - Variance stabilization (Pearson residuals, no log)
2. scran - Deconvolution size factors (log2 transformation)
3. LogNormalize - Simple library size normalization (log1p transformation)
4. scKWARN - Kernel Weighted Adjusted Regularized Normalization (log1p transformation)

Controlled by: run_sctransform, run_scran, run_lognorm, run_sckwarn
Best method selected by: scIB-style composite score (40% bio + 60% batch) via
  majority vote across integration methods, or batch-only composite if no cell
  type annotations are available (see Benchmarking Metrics below)
Force specific method: integration_normalization_method = "scKWARN"

scKWARN Details:
- Uses LocASN (Local Adaptive Size Normalization) algorithm
- Applies log1p() as recommended by scKWARN authors
- Output includes 'scKWARN_linear' assay with linear (non-log) data
- Install: devtools::install_github('cyhsuTN/scKWARN')

BENCHMARKING METRICS:
-----------------------------
Modules 03 and 04 use scIB-style composite scoring to select the best
normalization and integration methods. When cell type annotations are
available (e.g., from label transfer), both batch correction and biological
conservation are evaluated.

Batch correction metrics (always computed):
  - batch_variance: Variance of batch proportions across neighborhoods
    (lower = better mixing)
  - batch_asw: Average Silhouette Width on batch labels
    (closer to 0 = better mixing)
  - LISI: Local Inverse Simpson Index on batch labels
    (higher = better mixing)

Biological conservation metrics (computed when celltype_column is available):
  - celltype_asw: Silhouette Width on cell type labels
    (higher = better separation of cell types)
  - celltype_lisi: LISI on cell type labels
    (lower = better local purity of cell types)
  - NMI: Normalized Mutual Information between Louvain clustering and
    cell type annotations (higher = better agreement)
  - ARI: Adjusted Rand Index between Louvain clustering and cell type
    annotations (higher = better agreement)

Composite scoring:
  - With bio metrics: 40% bio + 60% batch (scIB paper default)
  - Without bio metrics: equal-weight batch metrics only
  - Weights configurable via bio_weight/batch_weight in params.R

Method selection (Module 03):
  - Each integration method independently scores all normalizations
  - Best normalization per integration method determined by composite score
  - Majority vote across integration methods selects the winner
  - Tie-breaking by average composite score

Method selection (Module 04):
  - Existing scIB framework activated by setting celltype_column in params.R
  - Uses same composite scoring (bio_weight/batch_weight)

Cell type annotation detection priority:
  1. Explicit: params$celltype_column (e.g., "cell_type")
  2. Auto-detect common names: cell_type, celltype, CellType, annotation, etc.
  3. Fallback: batch-only scoring (fully backward compatible)

params.R settings for benchmarking:
```R
  celltype_column = "cell_type",           # Cell type annotation column
  bio_weight = 0.4,                        # Biological conservation weight
  batch_weight = 0.6,                      # Batch correction weight
```

Shared bio metric functions (functions.R):
  - compute_celltype_asw(embeddings, celltype_labels)
  - compute_celltype_lisi(embeddings, celltype_labels, perplexity)
  - compute_nmi_score(embeddings, celltype_labels, resolution)
  - compute_ari_score(embeddings, celltype_labels, resolution)
  All functions subsample to 5000 cells if needed, handle NA labels, and
  require >=2 cell types with >=10 cells.

================================================================================
11. QUICK REFERENCE
================================================================================

------------------------------------------------------------------------------
Essential Commands
------------------------------------------------------------------------------

# Validate setup
bash validate_pipeline.sh

# List available groups
./submit_all_jobs.sh --list-groups

# Full preprocessing (Steps 1-10)
./submit_all_jobs.sh --preprocess

# Downstream analysis
./submit_all_jobs.sh --downstream --1              # Group 1
./submit_all_jobs.sh --downstream --2              # Group 2
./submit_all_jobs.sh --downstream --group 1        # Explicit

# Legacy syntax (still works)
./submit_all_jobs.sh --downstream --ventricle LV
./submit_all_jobs.sh --downstream --ventricle 4V

# Downstream with optional modules
./submit_all_jobs.sh --downstream --1 --modules 2b,3       # afMF + norm
./submit_all_jobs.sh --downstream --1 --modules 3,3b,4     # Norm + ALRA + integration
./submit_all_jobs.sh --downstream --1 --modules 7b,8       # CLTS + DE

# Check job status
squeue -u $USER

------------------------------------------------------------------------------
Key Files to Check
------------------------------------------------------------------------------

BEFORE RUNNING:
  - samplesheet.csv (correct samples? dataset_name? group_id? group_label?)
  - fastq_samplesheet.csv (correct FASTQ paths?)
  - Complete_FASTQ/ (files exist?)
  - params.R (imputation_method? afmf_python path? alra settings?
              celltype_column? bio_weight/batch_weight?)

AFTER EACH STEP:
  - Output_dir_<dataset_name>/Single_cell_preprocessed/{STEP}_output/README.txt
  - logs/{step}/*.out (errors?)

AFTER MODULE 00:
  - objects/pipeline_environment.RData (has_afmf? has_alra? has_scice?)

FINAL OUTPUTS:
  - Output_dir_<dataset_name>/Single_cell_preprocessed/7_scCDC_correction/
  - Output_dir_<dataset_name>/Single_cell_clustering/11_Downstream_Analysis_*/
  - .../11_Downstream_Analysis_*/reports/Analysis_Report_MultiSample.html

------------------------------------------------------------------------------
Environment Variables
------------------------------------------------------------------------------

PROJECT_ROOT      : Base directory for project
SAMPLESHEET       : Path to samplesheet.csv
DATASET_NAME      : Dataset identifier (from samplesheet column 12)
PREPROCESS_DIR    : Output_dir_<dataset_name>/Single_cell_preprocessed
DOWNSTREAM_DIR    : Output_dir_<dataset_name>/Single_cell_clustering
PIPELINE_DIR      : Scripts/scrnaseq_pipeline
GROUP_ID          : Numeric group ID for filtering
GROUP_LABEL       : Human-readable group label
VENTRICLE_FILTER  : LV, 4V, or ALL (backward compatible)

------------------------------------------------------------------------------
Directory Structure Summary
------------------------------------------------------------------------------

Output_dir_<dataset_name>/
├── Single_cell_preprocessed/     # PREPROCESS_DIR (Steps 1-10)
│   ├── 1_CellRanger_output/
│   ├── 2_DropletQC_output/
│   ├── 3_qClus_empty_droplets/
│   ├── 4_Vaeda_doublet_detection/
│   ├── 5_Seurat_conversion/
│   ├── 6_Doublet_consensus/
│   ├── 7_scCDC_correction/       # ← R pipeline input
│   ├── 8_DecontX_correction/
│   ├── 9_CHOIR_integration/
│   │   ├── Individual/
│   │   ├── {group_label}_only/
│   │   └── All_combined/
│   └── 10_QC_visualization/
│
└── Single_cell_clustering/       # DOWNSTREAM_DIR (Step 11+)
    ├── 11_Downstream_Analysis_<group_label_A>/   # Group 1 (--1)
    │   ├── objects/
    │   │   ├── merged_object.rds
    │   │   ├── merged_object_imputed.rds           # [if imputation afmf]
    │   │   ├── merged_normalized/
    │   │   │   ├── merged_LogNormalize_unintegrated.rds
    │   │   │   ├── merged_LogNormalize_alra.rds    # [if imputation alra]
    │   │   │   └── ...
    │   │   ├── scice_subclustered_object.rds
    │   │   ├── 06_scice_data.RData              # scICE + IDclust results
    │   │   └── scice_clts_renormalized_object.rds  # [if run_clts]
    │   ├── 02b_Imputation_afMF/                    # [if imputation afmf]
    │   ├── 03b_Imputation_ALRA/                    # [if imputation alra]
    │   ├── 07b_CLTS_Normalization/                 # [if run_clts]
    │   └── ...
    └── 11_Downstream_Analysis_<group_label_B>/     # Group 2 (--2)
        └── [same structure as above]

------------------------------------------------------------------------------
Samplesheet Column Reference
------------------------------------------------------------------------------

Column  | Name          | Required | Description
--------|---------------|----------|------------------------------------------
1       | sample_id     | YES      | Unique sample identifier
2       | sex           | YES      | M or F
3       | age           | No       | Sample age (e.g., P60, 8w)
4       | batch         | YES      | Batch identifier
5       | ventricle     | No       | Anatomical region - backward compat
6       | condition     | No       | Treatment condition
7       | SRA_accession | No       | SRA accession number
8       | ftp_path      | No       | FTP download path
9       | local_fastq_dir| No      | Local FASTQ directory
10      | genome        | No       | Cell Ranger reference path
11      | include       | No       | TRUE/FALSE to include sample
12      | dataset_name  | Recommended| Dataset identifier for output naming
13      | group_id      | Recommended| Numeric group ID
14      | group_label   | Recommended| Human-readable group name

------------------------------------------------------------------------------
Module Reference
------------------------------------------------------------------------------

Module | Name                    | Optional | Controlled by
-------|-------------------------|----------|------------------------------
0      | Environment Setup       | No       | Always runs
1      | Load Data               | No       | Always runs
2      | QC Validation           | No       | Always runs
2b     | Imputation (afMF)       | YES      | imputation_method ∈ {afmf, both}
3      | Normalization           | No       | Always runs
3b     | Imputation (ALRA)       | YES      | imputation_method ∈ {alra, both}
       |                         |          | AND normalization ∈ {LogNorm, scran}
4      | Integration             | No       | Always runs
5      | CHOIR Clustering        | No       | run_choir_clustering
6      | Subclustering (scICE + IDclust) | No | subclustering_methods / run_scice_subclustering / run_idclust_subclustering
7      | Leiden Clustering       | No       | run_leiden_clustering
7b     | CLTS Re-normalization   | YES      | run_clts_renormalization = TRUE
8      | Differential Expression | No       | Always runs
9      | Gene Visualization      | No       | Always runs
10     | Final Summary           | No       | Always runs
11     | HTML Report             | No       | Always runs

------------------------------------------------------------------------------
Imputation Method Summary
------------------------------------------------------------------------------

Method | Module | Stage            | Input              | Output Assay
-------|--------|------------------|--------------------|--------------
afMF   | 02b    | BEFORE norm      | Raw counts         | "imputed"
ALRA   | 03b    | AFTER norm       | Log-normalized     | "ALRA"

ALRA Compatibility:
  LogNormalize : YES
  scran        : YES
  SCTransform  : NO (auto-skipped)

Recommended Settings for DE Analysis:
  imputation_method = "both"           # Run both methods
  use_afmf_for_normalization = FALSE   # Don't use afMF counts downstream
  use_alra_for_downstream = FALSE      # Don't use ALRA data downstream

------------------------------------------------------------------------------
Conda Environments
------------------------------------------------------------------------------

R_4_5           : R-based steps (2, 5, 6, 7, 8, 10, 11-21)
qclus           : Step 3 (QClus)
vaeda_env       : Step 4 (VAEDA)
afMF_SCImputation_env : Module 02b (afMF imputation)

================================================================================
END OF DOCUMENTATION
================================================================================
