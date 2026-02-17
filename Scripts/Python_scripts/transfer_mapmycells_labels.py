#!/usr/bin/env python3
# ==============================================================================
# Transfer MapMyCells annotations to integrated CHOIR AnnData object (.h5ad)
# ==============================================================================
# This script transfers cell-type labels from MapMyCells hierarchical mapping
# output to the integrated CHOIR AnnData object.
#
# MapMyCells cell_id format: BARCODE_SAMPLEID_all
# ==============================================================================

import pandas as pd
import anndata as ad
import numpy as np
import sys
import os
from datetime import datetime

print("=" * 78)
print("  MapMyCells Label Transfer to AnnData Object (Python)")
print("=" * 78)
print()

# --- Configuration -----------------------------------------------------------
base_dir = "/scicore/home/doetsch/kaiser0001/Single_cell_paper/Datasets/Human_Covid_LV_ChP_PMID_34153974"

mapmycells_csv = os.path.join(base_dir,
    "MapMyCells_Human_Covid/All_combined_CHOIR_integrated_10xWholeHumanBrain(CCN202210140)_HierarchicalMapping_UTC_1770459613213.csv")

h5ad_path = os.path.join(base_dir,
    "Output_dir_Human_Covid_LV_ChP/Single_cell_preprocessed/9_CHOIR_integration/All_combined/All_combined_CHOIR_integrated.h5ad")

sample_ids = ["c15_12", "c15_17", "c16_18", "CP_90", "CP_91"]

annotation_cols = [
    "supercluster_label", "supercluster_name", "supercluster_bootstrapping_probability",
    "cluster_label", "cluster_name", "cluster_bootstrapping_probability",
    "subcluster_label", "subcluster_name", "subcluster_alias",
    "subcluster_bootstrapping_probability"
]

# --- Load MapMyCells annotations ----------------------------------------------
print("[1] Loading MapMyCells annotations...")
mmc = pd.read_csv(mapmycells_csv)
print(f"    Total annotated cells: {len(mmc)}")
print(f"    Columns: {', '.join(mmc.columns.tolist())}")

# Parse cell_id: strip trailing _all
mmc["cell_id_clean"] = mmc["cell_id"].str.replace(r"_all$", "", regex=True)

# Parse sample and barcode
mmc["parsed_sample"] = None
mmc["parsed_barcode"] = None

for sid in sample_ids:
    pattern = f"_{sid}$"
    mask = mmc["cell_id_clean"].str.endswith(f"_{sid}")
    mmc.loc[mask, "parsed_sample"] = sid
    mmc.loc[mask, "parsed_barcode"] = mmc.loc[mask, "cell_id_clean"].str.replace(
        pattern, "", regex=True)

print("\n    Parsed sample distribution:")
sample_counts = mmc["parsed_sample"].value_counts(dropna=False)
for s, c in sample_counts.items():
    label = "UNMATCHED" if pd.isna(s) else s
    print(f"      {label:<10} : {c} cells")

n_unmatched = mmc["parsed_sample"].isna().sum()
if n_unmatched > 0:
    print(f"\n    WARNING: {n_unmatched} cells could not be parsed!")
    examples = mmc.loc[mmc["parsed_sample"].isna(), "cell_id"].head(5).tolist()
    for ex in examples:
        print(f"      {ex}")

# --- Load h5ad ----------------------------------------------------------------
print(f"\n[2] Loading AnnData object: {os.path.basename(h5ad_path)}")
if not os.path.exists(h5ad_path):
    print(f"    ERROR: File not found: {h5ad_path}")
    sys.exit(1)

adata = ad.read_h5ad(h5ad_path)
print(f"    Shape: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"    Existing obs columns: {', '.join(adata.obs.columns.tolist()[:20])}...")

int_barcodes = adata.obs_names.tolist()
print(f"\n    AnnData barcode examples:")
for i, bc in enumerate(int_barcodes[:5]):
    print(f"      [{i+1}] {bc}")

print(f"\n    MapMyCells cell_id_clean examples:")
for i, cid in enumerate(mmc["cell_id_clean"].head(5).tolist()):
    print(f"      [{i+1}] {cid}")

# --- Matching strategies ------------------------------------------------------
print("\n[3] Testing matching strategies...")

# Strategy 1: Direct match cell_id_clean vs obs_names
int_set = set(int_barcodes)
direct_match = sum(1 for c in mmc["cell_id_clean"] if c in int_set)
print(f"    Strategy 1 - Direct (cell_id_clean): {direct_match} / {len(mmc)}")

# Strategy 2: Full cell_id vs obs_names
full_match = sum(1 for c in mmc["cell_id"] if c in int_set)
print(f"    Strategy 2 - Full cell_id:           {full_match} / {len(mmc)}")

# Strategy 3: Strip -1 from adata barcodes
int_stripped = [bc.replace("-1", "", 1) if bc.endswith("-1") else bc for bc in int_barcodes]
int_stripped_set = set(int_stripped)
stripped_match = sum(1 for c in mmc["cell_id_clean"] if c in int_stripped_set)
print(f"    Strategy 3 - Stripped '-1':           {stripped_match} / {len(mmc)}")

# Choose best
strategies = {1: direct_match, 2: full_match, 3: stripped_match}
best = max(strategies, key=strategies.get)
print(f"\n    Best strategy: {best} (matched {strategies[best]} cells)")

if best == 1:
    mmc_key_col = "cell_id_clean"
    int_keys = int_barcodes
elif best == 2:
    mmc_key_col = "cell_id"
    int_keys = int_barcodes
else:
    mmc_key_col = "cell_id_clean"
    int_keys = int_stripped

# Build mapping: int_key -> original barcode
key_to_original = dict(zip(int_keys, int_barcodes))

# Build mmc lookup
mmc_lookup = mmc.set_index(mmc_key_col)

matched_keys = set(int_keys) & set(mmc_lookup.index)
n_matched = len(matched_keys)
match_pct = n_matched / len(int_barcodes) * 100

print(f"\n    MATCHING SUMMARY:")
print(f"      Cells in AnnData:      {len(int_barcodes)}")
print(f"      Cells in MapMyCells:   {len(mmc)}")
print(f"      Successfully matched:  {n_matched} ({match_pct:.2f}% of AnnData)")
print(f"      Unmatched in AnnData:  {len(int_barcodes) - n_matched}")
print(f"      Unmatched in MMC:      {len(mmc) - n_matched}")

# --- Transfer annotations ----------------------------------------------------
print("\n[4] Transferring annotations...")

for col in annotation_cols:
    new_col = f"{col}_MapMyCells"
    # Initialize with NaN/None
    adata.obs[new_col] = None

    if n_matched > 0:
        # Build a series for matched cells
        for int_key, orig_bc in key_to_original.items():
            if int_key in mmc_lookup.index:
                val = mmc_lookup.at[int_key, col]
                adata.obs.at[orig_bc, new_col] = val

    n_filled = adata.obs[new_col].notna().sum()
    print(f"    {new_col}: {n_filled} cells annotated")

# --- Distribution summaries ---------------------------------------------------
print("\n    Supercluster distribution (matched cells):")
sc_dist = adata.obs["supercluster_name_MapMyCells"].value_counts(dropna=False)
for name, count in sc_dist.items():
    label = "NA (unmatched)" if pd.isna(name) else name
    print(f"      {label:<30} : {count}")

print("\n    Cluster distribution (top 15):")
cl_dist = adata.obs["cluster_name_MapMyCells"].value_counts(dropna=True).head(15)
for name, count in cl_dist.items():
    print(f"      {name:<30} : {count}")

# --- Save ---------------------------------------------------------------------
print(f"\n[5] Saving to: {h5ad_path}")

# Ensure string columns are properly typed for h5ad
for col in annotation_cols:
    new_col = f"{col}_MapMyCells"
    if "probability" in col:
        adata.obs[new_col] = pd.to_numeric(adata.obs[new_col], errors="coerce")
    else:
        adata.obs[new_col] = adata.obs[new_col].astype(str).replace("None", np.nan).replace("nan", np.nan)

adata.write_h5ad(h5ad_path)
print("    Saved successfully.")

# --- Reload and verify -------------------------------------------------------
print("\n[6] Reloading to verify...")
adata_verify = ad.read_h5ad(h5ad_path)

mmc_cols = [c for c in adata_verify.obs.columns if c.endswith("_MapMyCells")]
print(f"    MapMyCells columns found: {len(mmc_cols)}")
for c in mmc_cols:
    print(f"      {c}")

n_annotated = adata_verify.obs["supercluster_name_MapMyCells"].notna().sum()
print(f"\n    Cells with annotations after reload: {n_annotated} / {adata_verify.shape[0]}")

print(f"\n    First 5 barcodes (reloaded):")
for i in range(min(5, adata_verify.shape[0])):
    bc = adata_verify.obs_names[i]
    sc = adata_verify.obs.at[bc, "supercluster_name_MapMyCells"]
    cl = adata_verify.obs.at[bc, "cluster_name_MapMyCells"]
    sub = adata_verify.obs.at[bc, "subcluster_name_MapMyCells"]
    has_label = "YES" if pd.notna(sc) and str(sc) != "nan" else "NO"
    print(f"      [{i+1}] {bc}")
    print(f"           supercluster={sc} | cluster={cl} | subcluster={sub} | has_label={has_label}")

# --- Per-sample verification ---------------------------------------------------
print("\n    Per-sample annotation verification:")
if "parsed_sample" not in adata_verify.obs.columns:
    # Reconstruct sample info from barcodes
    for sid in sample_ids:
        mask = adata_verify.obs_names.str.contains(f"_{sid}")
        n_total = mask.sum()
        n_ann = (adata_verify.obs.loc[mask, "supercluster_name_MapMyCells"].notna() &
                 (adata_verify.obs.loc[mask, "supercluster_name_MapMyCells"].astype(str) != "nan")).sum()
        pct = n_ann / n_total * 100 if n_total > 0 else 0
        print(f"      {sid:<10}: {n_ann}/{n_total} annotated ({pct:.1f}%)")

print("\n" + "=" * 78)
print(f"[DONE] Python label transfer complete.")
print(f"Timestamp: {datetime.now()}")
print("=" * 78)
