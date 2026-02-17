#!/usr/bin/env python
# ==============================================================================
# run_afmf_imputation.py
# ==============================================================================
# Standalone afMF imputation script - called via system2() from R
# Part of the scrnaseq_pipeline (Module 02b)
#
# Usage:
#   python run_afmf_imputation.py input.csv output.csv [options]
#
# Options:
#   --transpose       Input is cells × genes (will transpose for afMF)
#   --min-genes N     Filter genes expressed in fewer than N cells (default: 10)
#   --log-file PATH   Write log to file instead of stdout
#   --quiet           Suppress progress messages
#
# Input format:
#   CSV with row names (genes) and column names (cells)
#   Format: genes × cells (default) or cells × genes (with --transpose)
#
# Output format:
#   CSV with same dimensions as input, imputed values
#
# Author: scrnaseq_pipeline
# Date: 2026-01-09
# ==============================================================================

import sys
import os
import argparse
import time
from datetime import datetime

import pandas as pd
import numpy as np

# Suppress warnings during import
import warnings
warnings.filterwarnings('ignore')

def log_message(msg, quiet=False, log_file=None):
    """Print message with timestamp."""
    if quiet:
        return
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    formatted = f"[{timestamp}] {msg}"
    print(formatted, flush=True)
    if log_file:
        with open(log_file, 'a') as f:
            f.write(formatted + "\n")

def filter_genes(df, min_cells=10):
    """Filter genes expressed in fewer than min_cells."""
    # Count non-zero cells per gene
    nonzero_counts = (df > 0).sum(axis=1)
    genes_to_keep = nonzero_counts >= min_cells
    n_removed = (~genes_to_keep).sum()
    return df.loc[genes_to_keep], n_removed

def run_afmf_imputation(input_file, output_file, transpose=False, 
                        min_genes=10, quiet=False, log_file=None):
    """
    Run afMF imputation on count matrix.
    
    Parameters:
    -----------
    input_file : str
        Path to input CSV (genes × cells or cells × genes)
    output_file : str
        Path to output CSV
    transpose : bool
        If True, input is cells × genes (will transpose)
    min_genes : int
        Filter genes expressed in fewer than this many cells
    quiet : bool
        Suppress progress messages
    log_file : str
        Path to log file (optional)
    
    Returns:
    --------
    bool : Success status
    """
    
    start_time = time.time()
    
    log_message("=" * 60, quiet, log_file)
    log_message("afMF IMPUTATION", quiet, log_file)
    log_message("=" * 60, quiet, log_file)
    
    # Import afMF
    log_message("Loading afMF module...", quiet, log_file)
    try:
        from afMF.runafMF import afMF
        log_message("  afMF loaded successfully", quiet, log_file)
    except ImportError as e:
        log_message(f"ERROR: Failed to import afMF: {e}", False, log_file)
        log_message("  Make sure afMF is installed: pip install .", False, log_file)
        return False
    
    # Load input data
    log_message(f"Loading input data: {input_file}", quiet, log_file)
    try:
        dat = pd.read_csv(input_file, index_col=0)
        log_message(f"  Loaded matrix: {dat.shape[0]} rows × {dat.shape[1]} columns", quiet, log_file)
    except Exception as e:
        log_message(f"ERROR: Failed to load input file: {e}", False, log_file)
        return False
    
    # Store original orientation info
    original_shape = dat.shape
    original_index = dat.index.tolist()
    original_columns = dat.columns.tolist()
    
    # Transpose if needed (afMF expects genes × cells)
    if transpose:
        log_message("Transposing input (cells × genes → genes × cells)...", quiet, log_file)
        dat = dat.T
        log_message(f"  After transpose: {dat.shape[0]} genes × {dat.shape[1]} cells", quiet, log_file)
    else:
        log_message(f"Input format: {dat.shape[0]} genes × {dat.shape[1]} cells", quiet, log_file)
    
    # Filter genes
    if min_genes > 0:
        log_message(f"Filtering genes expressed in < {min_genes} cells...", quiet, log_file)
        dat_filtered, n_removed = filter_genes(dat, min_genes)
        log_message(f"  Removed {n_removed} genes, kept {dat_filtered.shape[0]} genes", quiet, log_file)
        
        # Store removed genes for later
        removed_genes = dat.index.difference(dat_filtered.index)
        removed_data = dat.loc[removed_genes]
        dat = dat_filtered
    else:
        removed_genes = pd.Index([])
        removed_data = pd.DataFrame()
    
    # Check minimum size
    if dat.shape[0] < 100:
        log_message(f"WARNING: Only {dat.shape[0]} genes after filtering (< 100)", quiet, log_file)
        log_message("  afMF may produce suboptimal results with few features", quiet, log_file)
    
    if dat.shape[1] < 50:
        log_message(f"WARNING: Only {dat.shape[1]} cells (< 50)", quiet, log_file)
        log_message("  afMF may produce suboptimal results with few cells", quiet, log_file)
    
    # Run afMF imputation
    log_message("Running afMF imputation...", quiet, log_file)
    log_message(f"  Input: {dat.shape[0]} genes × {dat.shape[1]} cells", quiet, log_file)
    
    imputation_start = time.time()
    try:
        # Suppress runtime warnings from afMF
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            imputed_dat = afMF(dat)
        imputation_time = time.time() - imputation_start
        log_message(f"  Imputation completed in {imputation_time:.1f} seconds", quiet, log_file)
    except Exception as e:
        log_message(f"ERROR: afMF imputation failed: {e}", False, log_file)
        return False
    
    # Verify output
    if imputed_dat.shape != dat.shape:
        log_message(f"WARNING: Output shape {imputed_dat.shape} differs from input {dat.shape}", quiet, log_file)
    
    # Add back filtered genes (with original values, not imputed)
    if len(removed_genes) > 0:
        log_message(f"Adding back {len(removed_genes)} filtered genes (unimputed)...", quiet, log_file)
        imputed_dat = pd.concat([imputed_dat, removed_data])
        # Reorder to match original gene order
        original_gene_order = dat.index.tolist() + removed_genes.tolist()
        imputed_dat = imputed_dat.reindex(original_gene_order)
    
    # Transpose back if needed
    if transpose:
        log_message("Transposing output back (genes × cells → cells × genes)...", quiet, log_file)
        imputed_dat = imputed_dat.T
        log_message(f"  Final output: {imputed_dat.shape[0]} cells × {imputed_dat.shape[1]} genes", quiet, log_file)
    else:
        log_message(f"Final output: {imputed_dat.shape[0]} genes × {imputed_dat.shape[1]} cells", quiet, log_file)
    
    # Ensure non-negative values
    n_negative = (imputed_dat < 0).sum().sum()
    if n_negative > 0:
        log_message(f"  Clipping {n_negative} negative values to 0", quiet, log_file)
        imputed_dat = imputed_dat.clip(lower=0)
    
    # Save output
    log_message(f"Saving imputed matrix to: {output_file}", quiet, log_file)
    try:
        # Create output directory if needed
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        imputed_dat.to_csv(output_file)
        
        # Verify file was created
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file) / (1024 * 1024)  # MB
            log_message(f"  Saved successfully ({file_size:.1f} MB)", quiet, log_file)
        else:
            log_message("ERROR: Output file was not created", False, log_file)
            return False
            
    except Exception as e:
        log_message(f"ERROR: Failed to save output: {e}", False, log_file)
        return False
    
    # Summary statistics
    log_message("", quiet, log_file)
    log_message("--- IMPUTATION SUMMARY ---", quiet, log_file)
    
    # Compare sparsity
    if transpose:
        # Compare in cells × genes orientation
        pass  # Already in correct orientation
    
    original_zeros = (dat == 0).sum().sum()
    original_total = dat.size
    original_sparsity = original_zeros / original_total * 100
    
    imputed_zeros = (imputed_dat == 0).sum().sum()
    imputed_total = imputed_dat.size
    imputed_sparsity = imputed_zeros / imputed_total * 100
    
    log_message(f"  Original sparsity: {original_sparsity:.1f}% zeros", quiet, log_file)
    log_message(f"  Imputed sparsity:  {imputed_sparsity:.1f}% zeros", quiet, log_file)
    log_message(f"  Dropout recovery:  {original_sparsity - imputed_sparsity:.1f}% points", quiet, log_file)
    
    # Value ranges
    log_message(f"  Original value range: [{dat.values.min():.2f}, {dat.values.max():.2f}]", quiet, log_file)
    log_message(f"  Imputed value range:  [{imputed_dat.values.min():.2f}, {imputed_dat.values.max():.2f}]", quiet, log_file)
    
    total_time = time.time() - start_time
    log_message("", quiet, log_file)
    log_message(f"Total runtime: {total_time:.1f} seconds ({total_time/60:.2f} minutes)", quiet, log_file)
    log_message("=" * 60, quiet, log_file)
    log_message("afMF IMPUTATION COMPLETE", quiet, log_file)
    log_message("=" * 60, quiet, log_file)
    
    return True


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Run afMF imputation on scRNA-seq count matrix",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Standard usage (genes × cells input):
  python run_afmf_imputation.py counts.csv imputed.csv
  
  # If input is cells × genes:
  python run_afmf_imputation.py counts.csv imputed.csv --transpose
  
  # With custom gene filter:
  python run_afmf_imputation.py counts.csv imputed.csv --min-genes 5
  
  # Save log to file:
  python run_afmf_imputation.py counts.csv imputed.csv --log-file imputation.log
"""
    )
    
    parser.add_argument("input_file", help="Input CSV file (genes × cells)")
    parser.add_argument("output_file", help="Output CSV file")
    parser.add_argument("--transpose", action="store_true",
                        help="Input is cells × genes (will transpose for afMF)")
    parser.add_argument("--min-genes", type=int, default=10,
                        help="Filter genes expressed in < N cells (default: 10)")
    parser.add_argument("--log-file", type=str, default=None,
                        help="Write log to file")
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress progress messages")
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(args.input_file):
        print(f"ERROR: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Run imputation
    success = run_afmf_imputation(
        input_file=args.input_file,
        output_file=args.output_file,
        transpose=args.transpose,
        min_genes=args.min_genes,
        quiet=args.quiet,
        log_file=args.log_file
    )
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
