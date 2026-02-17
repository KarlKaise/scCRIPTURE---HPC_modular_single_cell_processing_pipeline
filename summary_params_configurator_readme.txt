summary_params_configurator_readme.txt
========================================

scRNA-seq Pipeline Configurator — Functional Summary
-----------------------------------------------------

PURPOSE:
A self-contained, browser-based GUI for configuring all parameters of a
modular single-cell RNA-seq downstream analysis pipeline. It replaces
manual editing of params.R files with a guided, step-by-step wizard that
validates inputs and exports a ready-to-run configuration.

WORKFLOW (two phases, 15 steps total):

  Phase A — Sample Sheet (steps A1–A5)
    - Load a sample_sheet.csv (drag-and-drop or file picker)
    - Assign column roles (sample_name, sex, batch, condition, etc.)
    - Inspect experimental design via cross-tabulations
    - Edit, validate, and export the finalised sample sheet

  Phase B — params.R Configuration (steps B0–B9)
    - Optionally load an existing params.R or start from a template
    - Configure parameters across eight thematic sections:
        B1  Project metadata & input paths
        B2  QC thresholds & cell filtering
        B3  Normalization method selection
        B4  Imputation (afMF / ALRA / both / none)
        B5  Batch integration & benchmarking strategy
        B6  Clustering (CHOIR, Leiden, scICE, scCobra)
        B7  Differential expression settings
        B8  Visualization options & output paths
    - Review all parameters with a live R-code preview
    - Export the full params.R (template-injected) or a flat params list

KEY FEATURES:
  - Tag-based template injection: user-configured values are inserted into
    the full params.R template at <<CFG:key>> tags, preserving all
    pipeline infrastructure code, validation logic, and inline docs.
  - Per-parameter approve/modify workflow with visual status tracking.
  - Input validation and contextual warnings (e.g., SCTransform ↔ ALRA
    incompatibility, batch variable recommendations).
  - Session save/restore (JSON) for interrupted work.
  - Progress bar and sidebar navigation with step-locking to enforce
    sequential completion.
  - Fully client-side — no server required; runs from a single HTML file.

PIPELINE MODULES COVERED:
  00  Environment Setup
  01  Load & QC
  02  Merge Samples
  02b Imputation – afMF
  03  Normalization (SCTransform / scran / LogNormalize)
  03b Imputation – ALRA
  04  Integration & Benchmarking
  05  CHOIR Clustering
  05b scCobra
  06  scICE Subclustering
  07  Leiden Clustering
  07b CLTS Re-normalization
  08  Differential Expression
  09  Visualization

OUTPUT FILES:
  - params.R          Full pipeline configuration (template-injected)
  - params_flat.R     Minimal params list (reference only)
  - sample_sheet.csv  Validated sample metadata
