# Changelog

All notable changes to **scCRIPTURE** will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- MapMyCells cell type annotation transfer scripts (Python + R)
- SLURM wrapper for MapMyCells label transfer
- Analytical block framework for post-pipeline analyses
- Browser-based pipeline configurator GUI (`params_scrna_pipeline_configurator.html`)

### Changed
- Group-based filtering system replacing hardcoded ventricle logic
- Backward-compatible `--ventricle` syntax mapped to generic groups

## [1.0.0] — 2025-XX-XX

### Added
- Full preprocessing pipeline (Steps 1–10) with SLURM array job support
- Downstream analysis pipeline (Steps 11–21) with modular R framework
- Multi-method doublet detection with consensus voting (VAEDA + scDblFinder + DoubletFinder)
- Dual contamination correction (scCDC + DecontX)
- Optional imputation modules (afMF + ALRA) with automatic compatibility checks
- Optional CLTS re-normalisation for improved DE accuracy
- Multi-method integration benchmarking (scVI, Harmony, BBKNN, Scanorama)
- Subclustering via scICE (Julia) and IDclust (R)
- Normalisation benchmarking (SCTransform vs scran vs LogNormalise)
- Master job submission script (`submit_all_jobs.sh`) with dependency tracking
- Pipeline validation script (`validate_pipeline.sh`)
- Comprehensive QC visualisation and HTML report generation
- Session save/restore in configurator GUI
