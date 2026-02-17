#!/bin/bash
# ============================================================================
# Universal Modular scRNA-seq Pipeline — Master Submission Script
# ============================================================================
#
# This script submits all pipeline steps to SLURM with proper dependencies.
# It reads samplesheet.csv to determine samples, groups, and output paths.
#
# All paths are resolved relative to the script's location (PROJECT ROOT).
# Just copy the project folder anywhere and run - no path changes needed.
#
# ============================================================================
# QUICK COMMAND EXAMPLES
# ============================================================================
#
# -------------------------
# VALIDATION & INFO
# -------------------------
#
#   ./submit_all_jobs.sh --validate
#
#   # List available groups
#   ./submit_all_jobs.sh --list-groups
#
#   # Dry run (see what would be submitted)
#   ./submit_all_jobs.sh --preprocess --dry-run
#
# -------------------------
# PREPROCESSING (Steps 1-10)
# -------------------------
#
#   # Run full preprocessing pipeline (all steps 1-10)
#   ./submit_all_jobs.sh --preprocess
#
#   # Run only specific preprocessing step
#   ./submit_all_jobs.sh --only 1                    # Only Cell Ranger
#   ./submit_all_jobs.sh --only 3                    # Only step 3 (QClus)
#   ./submit_all_jobs.sh --only 7                    # Only step 7 (scCDC)
#   ./submit_all_jobs.sh --only 8                    # Only step 8 (DecontX)
#
#   # Run preprocessing steps in a range
#   ./submit_all_jobs.sh --from 1 --to 5             # Steps 1,2,3,4,5
#   ./submit_all_jobs.sh --from 6 --to 10            # Steps 6,7,8,9,10
#
# -------------------------
# DOWNSTREAM (Modules 0-11)
# -------------------------
#
#   # Run full downstream pipeline
#   ./submit_all_jobs.sh --downstream --1            # Group 1
#   ./submit_all_jobs.sh --downstream --2            # Group 2
#   ./submit_all_jobs.sh --downstream --3            # Group 3
#
#   # Alternative syntax
#   ./submit_all_jobs.sh --downstream --group 1
#   ./submit_all_jobs.sh --downstream --group 2
#   ./submit_all_jobs.sh --downstream --group 3
#
#   # Backward compatible (if ventricle column matches)
#   ./submit_all_jobs.sh --downstream --ventricle LV
#
#   # Dry run downstream
#   ./submit_all_jobs.sh --downstream --1 --dry-run
#   ./submit_all_jobs.sh --downstream --2 --dry-run
#   ./submit_all_jobs.sh --downstream --3 --dry-run
#
#   # Run specific downstream modules only
#   ./submit_all_jobs.sh --downstream --1 --modules 2        # Only QC validation
#   ./submit_all_jobs.sh --downstream --2 --modules 5,6,7    # Clustering only
#   ./submit_all_jobs.sh --downstream --3 --modules 8,9,10   # DE + viz + summary
#
#   # Run downstream modules in a range
#   ./submit_all_jobs.sh --downstream --1 --start 0 --stop 4   # Setup through integration
#   ./submit_all_jobs.sh --downstream --2 --start 5            # Clustering to end
#   ./submit_all_jobs.sh --downstream --3 --start 8 --stop 11  # DE + reports
#
# -------------------------
# COMMON WORKFLOWS
# -------------------------
#
#   # Full pipeline from scratch (all samples together)
#   ./submit_all_jobs.sh --preprocess
#   # Wait for completion, then:
#   ./submit_all_jobs.sh --downstream --1
#
#   # Full pipeline with 3-group analysis
#   ./submit_all_jobs.sh --preprocess
#   # Wait for completion, then:
#   ./submit_all_jobs.sh --downstream --1
#   ./submit_all_jobs.sh --downstream --2
#   ./submit_all_jobs.sh --downstream --3
#
#   # Re-run clustering for all groups
#   ./submit_all_jobs.sh --downstream --1 --start 5
#   ./submit_all_jobs.sh --downstream --2 --start 5
#   ./submit_all_jobs.sh --downstream --3 --start 5
#
#   # Re-run only differential expression for group 2
#   ./submit_all_jobs.sh --downstream --2 --modules 8
#
#   # Regenerate final report for group 3
#   ./submit_all_jobs.sh --downstream --3 --modules 10,11
#
#   # Monitor jobs
#   squeue -u $USER
#
#   # Cancel specific job
#   scancel <job_id>
#
# ============================================================================
# PIPELINE STRUCTURE
# ============================================================================
#
# PREPROCESSING STEPS (1-10):
#   1  - Cell Ranger Count       (Array: N jobs, ~6-12h each)
#   2  - DropletQC Nuclear Fraction (Array: N jobs, ~1-2h each)
#   3  - QClus Empty Droplets    (Array: N jobs, ~1-2h each)
#   4  - VAEDA Doublet Detection (Array: N jobs, ~2-4h each)
#   5  - AnnData to Seurat       (Array: N jobs, ~30min each)
#   6  - scDblFinder/DoubletFinder (Array: N jobs, ~1-2h each)
#   7  - scCDC Correction        (Array: N jobs, ~1-2h each)
#   8  - DecontX Correction      (Array: N jobs, ~1-2h each)
#   9  - CHOIR Integration       (Single job, ~2-4h)
#   10 - QC Visualization        (Single job, ~1-2h)
#
# DOWNSTREAM MODULES (0-11):
#   0  - Environment Setup
#   1  - Load Data (reads corrected .rds files for selected group)
#   2  - QC Validation
#   3  - Normalization (SCTransform)
#   4  - Integration (Harmony with batch correction)
#   5  - CHOIR Clustering
#   6  - scICE Subclustering
#   7  - Leiden Clustering
#   8  - Differential Expression
#   9  - Gene Visualization
#   10 - Final Summary
#   11 - HTML Report
#
# ============================================================================
# SAMPLESHEET STRUCTURE
# ============================================================================
#
# Main samplesheet.csv (14 columns):
#   sample_id,sex,age,batch,fastq_path,cellranger_ref,species,tissue,
#   ventricle,condition,notes,dataset_name,group_id,group_label
#
# Column reference:
#   Col 1:  sample_id      - Unique sample name (e.g., sample_A, sample_B)
#   Col 2:  sex            - Male/Female
#   Col 3:  age            - Sample age (e.g., P60)
#   Col 4:  batch          - Library prep batch (e.g., batch_A, batch_B, batch_C)
#   Col 5:  fastq_path     - Path to FASTQ directory
#   Col 6:  cellranger_ref - Path to Cell Ranger reference
#   Col 7:  species        - Species identifier
#   Col 8:  tissue         - Tissue type
#   Col 9:  ventricle      - Region identifier
#   Col 10: condition      - Experimental condition
#   Col 11: notes          - Additional sample metadata
#   Col 12: dataset_name   - Dataset identifier (determines output dir)
#   Col 13: group_id       - Numeric group for downstream filtering (1, 2, 3, ...)
#   Col 14: group_label    - Human-readable group name for output dir
#
# Downstream sample_sheet.csv (9 columns):
#   sample_name,sex,age,batch,ventricle,condition,include,group_id,group_label
#
# IMPORTANT: Both samplesheets must have matching group_id and group_label values!
#
# ============================================================================

set -e

# ---------------------------------------------------------------------------
# CONFIGURATION - AUTOMATIC PATH DETECTION
# ---------------------------------------------------------------------------

# Get the directory where this script is located (PROJECT ROOT)
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SCRIPT_DIR="${BASE_DIR}/Scripts/Slurm_scripts"
CONFIG_DIR="${BASE_DIR}/config"

# SAMPLESHEET - Always in project root with standard name
SAMPLESHEET="${BASE_DIR}/samplesheet.csv"
FASTQ_SAMPLESHEET="${BASE_DIR}/fastq_samplesheet.csv"

# Preprocessing Scripts (Steps 1-10)
PREPROCESS_SCRIPTS=(
    "${SCRIPT_DIR}/01_cellranger_count.sh"
    "${SCRIPT_DIR}/02_DropletQC_nuclear_fraction.sh"
    "${SCRIPT_DIR}/03_qClus_empty_droplets.sh"
    "${SCRIPT_DIR}/04_Vaeda_doublet_analysis.sh"
    "${SCRIPT_DIR}/05_Vaeda_anndata_to_seurat.sh"
    "${SCRIPT_DIR}/06_scDblFinder_DoubletFinder_analysis.sh"
    "${SCRIPT_DIR}/07_scCDC_correction.sh"
    "${SCRIPT_DIR}/08_DecontX_correction.sh"
    "${SCRIPT_DIR}/09_CHOIR_integration.sh"
    "${SCRIPT_DIR}/10_QC_visualization.sh"
)

PREPROCESS_NAMES=(
    "Cell Ranger Count"
    "DropletQC Nuclear Fraction & Empty Droplet ID"
    "QClus Contaminated Droplet Filtering"
    "VAEDA Doublet Detection"
    "AnnData to Seurat Conversion"
    "scDblFinder + DoubletFinder Consensus"
    "scCDC Contamination Correction"
    "DecontX Ambient RNA Correction"
    "CHOIR Integration (Individual, by Group, All)"
    "QC Visualization"
)

# Downstream Analysis Script
DOWNSTREAM_SCRIPT="${SCRIPT_DIR}/11_downstream_analysis.sh"

# All step names
ALL_STEP_NAMES=(
    "Cell Ranger Count"
    "DropletQC Nuclear Fraction"
    "QClus Empty Droplets"
    "VAEDA Doublet Detection"
    "AnnData to Seurat"
    "scDblFinder/DoubletFinder"
    "scCDC Correction"
    "DecontX Correction"
    "CHOIR Integration"
    "QC Visualization"
    "Environment Setup"
    "Load Data"
    "QC Validation"
    "Normalization"
    "Integration"
    "CHOIR Clustering"
    "scICE Subclustering"
    "Leiden Clustering"
    "Differential Expression"
    "Gene Visualization"
    "Final Summary & HTML Report"
)

# Downstream module names (for display)
DOWNSTREAM_MODULE_NAMES=(
    "Environment Setup"
    "Load Data"
    "QC Validation"
    "Normalization"
    "Integration"
    "CHOIR Clustering"
    "scICE Subclustering"
    "Leiden Clustering"
    "Differential Expression"
    "Gene Visualization"
    "Final Summary"
    "HTML Report"
)

# ---------------------------------------------------------------------------
# SAMPLESHEET PARSING FUNCTIONS
# ---------------------------------------------------------------------------

validate_samplesheet() {
    if [[ ! -f "$SAMPLESHEET" ]]; then
        echo "ERROR: samplesheet.csv not found at: $SAMPLESHEET"
        echo ""
        echo "Please create samplesheet.csv with the following format:"
        echo "sample_id,sex,age,batch,fastq_path,cellranger_ref,species,tissue,ventricle,condition,notes,dataset_name,group_id,group_label"
        return 1
    fi

    HEADER=$(head -1 "$SAMPLESHEET")
    for col in sample_id sex batch dataset_name; do
        if [[ ! "$HEADER" == *"$col"* ]]; then
            echo "ERROR: samplesheet.csv missing '$col' column"
            return 1
        fi
    done

    N_SAMPLES=$(tail -n +2 "$SAMPLESHEET" | grep -v '^$' | wc -l)
    if [[ $N_SAMPLES -eq 0 ]]; then
        echo "ERROR: samplesheet.csv contains no samples"
        return 1
    fi

    return 0
}

validate_fastq_samplesheet() {
    if [[ ! -f "$FASTQ_SAMPLESHEET" ]]; then
        echo "WARNING: fastq_samplesheet.csv not found at: $FASTQ_SAMPLESHEET"
        echo "         Cell Ranger step will fail without this file."
        return 1
    fi

    HEADER=$(grep -v '^#' "$FASTQ_SAMPLESHEET" | head -1)
    for col in sample_id fastq_dir read1_file read2_file; do
        if [[ ! "$HEADER" == *"$col"* ]]; then
            echo "ERROR: fastq_samplesheet.csv missing '$col' column"
            return 1
        fi
    done

    return 0
}

get_samples() {
    tail -n +2 "$SAMPLESHEET" | cut -d',' -f1 | grep -v '^$'
}

get_samples_by_ventricle() {
    local vent=$1
    tail -n +2 "$SAMPLESHEET" | awk -F',' -v v="$vent" '$9 == v {print $1}'
}

get_ventricles() {
    tail -n +2 "$SAMPLESHEET" | cut -d',' -f9 | grep -v '^$' | sort -u
}

get_dataset_name() {
    # Get dataset_name from first data row (column 12)
    tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f12
}

# ---------------------------------------------------------------------------
# GROUP-BASED FILTERING FUNCTIONS
# ---------------------------------------------------------------------------

get_column_index() {
    # Get 1-based column index for a given column name
    local col_name=$1
    head -1 "$SAMPLESHEET" | tr ',' '\n' | grep -n "^${col_name}$" | cut -d: -f1
}

get_group_label_by_id() {
    # Get group_label for a given group_id
    local gid=$1
    local gid_col=$(get_column_index "group_id")
    local glabel_col=$(get_column_index "group_label")

    if [[ -z "$gid_col" || -z "$glabel_col" ]]; then
        echo ""
        return
    fi

    tail -n +2 "$SAMPLESHEET" | awk -F',' -v g="$gid" -v gc="$gid_col" -v lc="$glabel_col" '$gc == g {print $lc; exit}'
}

get_group_id_by_ventricle() {
    # Get group_id for a given ventricle value (for backward compatibility)
    local vent=$1
    local vent_col=$(get_column_index "ventricle")
    local gid_col=$(get_column_index "group_id")

    if [[ -z "$vent_col" || -z "$gid_col" ]]; then
        echo ""
        return
    fi

    tail -n +2 "$SAMPLESHEET" | awk -F',' -v v="$vent" -v vc="$vent_col" -v gc="$gid_col" '$vc == v {print $gc; exit}'
}

get_samples_by_group() {
    # Get sample IDs for a given group_id
    local gid=$1
    local gid_col=$(get_column_index "group_id")

    if [[ -z "$gid_col" ]]; then
        get_samples
        return
    fi

    tail -n +2 "$SAMPLESHEET" | awk -F',' -v g="$gid" -v gc="$gid_col" '$gc == g {print $1}'
}

get_ventricle_by_group() {
    # Get ventricle value for a given group_id (for backward compatibility)
    local gid=$1
    local gid_col=$(get_column_index "group_id")
    local vent_col=$(get_column_index "ventricle")

    if [[ -z "$gid_col" || -z "$vent_col" ]]; then
        echo ""
        return
    fi

    tail -n +2 "$SAMPLESHEET" | awk -F',' -v g="$gid" -v gc="$gid_col" -v vc="$vent_col" '$gc == g {print $vc; exit}'
}

list_available_groups() {
    # List all available groups from samplesheet
    local gid_col=$(get_column_index "group_id")
    local glabel_col=$(get_column_index "group_label")

    if [[ -z "$gid_col" || -z "$glabel_col" ]]; then
        echo "  WARNING: group_id and/or group_label columns not found in samplesheet"
        echo "  Falling back to ventricle-based filtering:"
        echo ""
        for vent in $(get_ventricles); do
            local count=$(get_samples_by_ventricle "$vent" | wc -l)
            local samples=$(get_samples_by_ventricle "$vent" | tr '\n' ' ')
            echo "    --ventricle ${vent}  →  ${count} samples: ${samples}"
        done
        return
    fi

    echo "  Available groups (use --<group_id> or --group <group_id>):"
    echo ""

    # Get unique group_id,group_label pairs
    tail -n +2 "$SAMPLESHEET" | awk -F',' -v gc="$gid_col" -v lc="$glabel_col" '{print $gc","$lc}' | sort -t',' -k1,1n -u | while IFS=',' read -r gid glabel; do
        if [[ -n "$gid" && -n "$glabel" ]]; then
            local count=$(get_samples_by_group "$gid" | wc -l)
            local samples=$(get_samples_by_group "$gid" | tr '\n' ' ')
            printf "    --%-2s  →  %-20s (%d samples: %s)\n" "$gid" "$glabel" "$count" "$samples"
        fi
    done

    echo ""
    echo "  Backward compatible (--ventricle):"
    for vent in $(get_ventricles); do
        local gid=$(get_group_id_by_ventricle "$vent")
        if [[ -n "$gid" ]]; then
            echo "    --ventricle ${vent}  →  maps to group ${gid}"
        else
            echo "    --ventricle ${vent}  →  (no group mapping)"
        fi
    done
}

print_samplesheet_summary() {
    echo ""
    echo "============================================================================"
    echo "SAMPLESHEET SUMMARY"
    echo "============================================================================"
    echo "Location: $SAMPLESHEET"
    echo ""

    local n_samples=$(get_samples | wc -l)
    local dataset_name=$(get_dataset_name)
    echo "Dataset name: $dataset_name"
    echo "Total samples: $n_samples"
    echo ""

    local gid_col=$(get_column_index "group_id")
    local glabel_col=$(get_column_index "group_label")

    echo "Samples:"
    echo "  Sample_ID    Sex      Age    Ventricle  Batch     Group"
    echo "  ----------   ------   ----   ---------  ------    -----"
    while IFS=',' read -r line; do
        [[ "$line" == sample_id* ]] && continue
        [[ -z "$line" ]] && continue

        local sample_id=$(echo "$line" | cut -d',' -f1)
        local sex=$(echo "$line" | cut -d',' -f2)
        local age=$(echo "$line" | cut -d',' -f3)
        local batch=$(echo "$line" | cut -d',' -f4)
        local ventricle=$(echo "$line" | cut -d',' -f9)
        local group_id=""
        local group_label=""

        if [[ -n "$gid_col" ]]; then
            group_id=$(echo "$line" | cut -d',' -f"$gid_col")
        fi
        if [[ -n "$glabel_col" ]]; then
            group_label=$(echo "$line" | cut -d',' -f"$glabel_col")
        fi

        if [[ -n "$group_id" ]]; then
            printf "  %-12s %-8s %-6s %-10s %-8s %s (%s)\n" "$sample_id" "$sex" "$age" "$ventricle" "$batch" "$group_id" "$group_label"
        else
            printf "  %-12s %-8s %-6s %-10s %s\n" "$sample_id" "$sex" "$age" "$ventricle" "$batch"
        fi
    done < "$SAMPLESHEET"
    echo ""

    echo "By Group:"
    if [[ -n "$gid_col" && -n "$glabel_col" ]]; then
        tail -n +2 "$SAMPLESHEET" | awk -F',' -v gc="$gid_col" -v lc="$glabel_col" '{print $gc","$lc}' | sort -t',' -k1,1n -u | while IFS=',' read -r gid glabel; do
            if [[ -n "$gid" ]]; then
                local samples=$(get_samples_by_group "$gid" | tr '\n' ' ')
                echo "  Group $gid ($glabel): $samples"
            fi
        done
    else
        echo "  (group_id/group_label columns not found - using ventricle)"
        for vent in $(get_ventricles); do
            local samples=$(get_samples_by_ventricle "$vent" | tr '\n' ' ')
            echo "  $vent: $samples"
        done
    fi
    echo ""

    echo "Sex distribution:"
    local n_male=$(tail -n +2 "$SAMPLESHEET" | cut -d',' -f2 | grep -ci "^male$" || echo 0)
    local n_female=$(tail -n +2 "$SAMPLESHEET" | cut -d',' -f2 | grep -ci "^female$" || echo 0)
    echo "  Male: $n_male"
    echo "  Female: $n_female"
    echo ""
    echo "============================================================================"
}

# ---------------------------------------------------------------------------
# PARSE ARGUMENTS
# ---------------------------------------------------------------------------

PIPELINE_START=1
PIPELINE_END=21
DRY_RUN=false
PREPROCESS_ONLY=false
DOWNSTREAM_ONLY=false
VALIDATE_ONLY=false
LIST_GROUPS=false
VENTRICLE=""
GROUP_ID=""
GROUP_LABEL=""
DOWNSTREAM_START=""
DOWNSTREAM_STOP=""
DOWNSTREAM_MODULES=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --from)
            PIPELINE_START="$2"
            shift 2
            ;;
        --to)
            PIPELINE_END="$2"
            shift 2
            ;;
        --only)
            PIPELINE_START="$2"
            PIPELINE_END="$2"
            shift 2
            ;;
        --preprocess)
            PREPROCESS_ONLY=true
            PIPELINE_START=1
            PIPELINE_END=10
            shift
            ;;
        --downstream)
            DOWNSTREAM_ONLY=true
            PIPELINE_START=11
            PIPELINE_END=21
            shift
            ;;
        --ventricle)
            # Backward compatibility: map ventricle to group_id
            VENTRICLE="$2"
            shift 2
            ;;
        --group)
            # Explicit group ID specification
            GROUP_ID="$2"
            shift 2
            ;;
        --[0-9]|--[0-9][0-9]|--[0-9][0-9][0-9])
            # Shorthand: --1, --2, --10, etc.
            GROUP_ID="${1#--}"
            shift
            ;;
        --list-groups)
            LIST_GROUPS=true
            shift
            ;;
        --start)
            DOWNSTREAM_START="$2"
            shift 2
            ;;
        --stop)
            DOWNSTREAM_STOP="$2"
            shift 2
            ;;
        --modules)
            DOWNSTREAM_MODULES="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --validate)
            VALIDATE_ONLY=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "PREPROCESSING OPTIONS:"
            echo "  --preprocess            Run only preprocessing (steps 1-10)"
            echo "  --from N                Start from preprocessing step N (1-10)"
            echo "  --to N                  Stop at preprocessing step N (1-10)"
            echo "  --only N                Run only preprocessing step N"
            echo ""
            echo "DOWNSTREAM OPTIONS:"
            echo "  --downstream            Run downstream analysis (GPU: A100)"
            echo "  --<N>                   Filter by group ID shorthand (e.g., --1, --2)"
            echo "  --group <N>             Filter by group ID (alternative syntax)"
            echo "  --ventricle LV|4V       Filter by ventricle (backward compatible)"
            echo "  --list-groups           Show available groups from samplesheet"
            echo "  --start N               Start from downstream module N (0-11)"
            echo "  --stop N                Stop at downstream module N (0-11)"
            echo "  --modules N,M,O         Run only specific modules (comma-separated)"
            echo ""
            echo "GENERAL OPTIONS:"
            echo "  --dry-run               Show what would be submitted without submitting"
            echo "  --validate              Validate samplesheet only"
            echo "  -h, --help              Show this help message"
            echo ""
            echo "PREPROCESSING STEPS:"
            echo "   1  - Cell Ranger Count          6 - scDblFinder/DoubletFinder"
            echo "   2  - DropletQC Nuclear Fraction  7 - scCDC Correction"
            echo "   3  - QClus Empty Droplets        8 - DecontX Correction"
            echo "   4  - VAEDA Doublet Detection     9 - CHOIR Integration"
            echo "   5  - AnnData to Seurat          10 - QC Visualization"
            echo ""
            echo "DOWNSTREAM MODULES:"
            echo "   0 - Environment Setup      6 - scICE Subclustering"
            echo "   1 - Load Data              7 - Leiden Clustering"
            echo "   2 - QC Validation          8 - Differential Expression"
            echo "   3 - Normalization          9 - Gene Visualization"
            echo "   4 - Integration           10 - Final Summary"
            echo "   5 - CHOIR Clustering      11 - HTML Report"
            echo ""
            echo "EXAMPLES:"
            echo "  $0 --preprocess"
            echo "  $0 --only 8"
            echo "  $0 --from 3 --to 8"
            echo "  $0 --list-groups"
            echo "  $0 --downstream --1"
            echo "  $0 --downstream --group 2"
            echo "  $0 --downstream --ventricle LV"
            echo "  $0 --downstream --1 --start 2 --stop 4"
            echo "  $0 --downstream --1 --modules 2,8,11"
            echo "  $0 --downstream --1 --dry-run"
            echo ""
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# VALIDATE SAMPLESHEET (needed for --list-groups and further processing)
# ---------------------------------------------------------------------------

if [[ "$LIST_GROUPS" == "true" ]]; then
    if ! validate_samplesheet; then
        exit 1
    fi
    echo ""
    echo "============================================================================"
    echo "AVAILABLE GROUPS FOR DOWNSTREAM ANALYSIS"
    echo "============================================================================"
    echo ""
    list_available_groups
    echo ""
    echo "============================================================================"
    exit 0
fi

# ---------------------------------------------------------------------------
# RESOLVE GROUP FROM VENTRICLE (backward compatibility)
# ---------------------------------------------------------------------------

if [[ -n "$VENTRICLE" && -z "$GROUP_ID" ]]; then
    # Map ventricle to group_id
    GROUP_ID=$(get_group_id_by_ventricle "$VENTRICLE")
    if [[ -z "$GROUP_ID" ]]; then
        # No group mapping, use ventricle directly for filtering
        echo "Note: Using ventricle '$VENTRICLE' directly (no group_id mapping found)"
    else
        echo "Note: Mapped --ventricle $VENTRICLE to --group $GROUP_ID"
    fi
fi

# ---------------------------------------------------------------------------
# RESOLVE GROUP LABEL
# ---------------------------------------------------------------------------

if [[ -n "$GROUP_ID" ]]; then
    GROUP_LABEL=$(get_group_label_by_id "$GROUP_ID")
    if [[ -z "$GROUP_LABEL" ]]; then
        echo "ERROR: Group ID '$GROUP_ID' not found in samplesheet"
        echo ""
        list_available_groups
        exit 1
    fi
    # Also get the corresponding ventricle for backward compatibility
    if [[ -z "$VENTRICLE" ]]; then
        VENTRICLE=$(get_ventricle_by_group "$GROUP_ID")
    fi
fi

# ---------------------------------------------------------------------------
# VALIDATE DOWNSTREAM GROUP REQUIREMENT
# ---------------------------------------------------------------------------

if [[ "$DOWNSTREAM_ONLY" == "true" && -z "$GROUP_ID" && -z "$VENTRICLE" ]]; then
    echo "ERROR: Group filter required for downstream analysis"
    echo ""
    echo "Use one of the following:"
    echo "  --<N>            Group ID shorthand (e.g., --1, --2)"
    echo "  --group <N>      Group ID"
    echo "  --ventricle X    Ventricle (backward compatible)"
    echo ""
    list_available_groups
    exit 1
fi

# Validate ventricle if using legacy mode without group
if [[ -n "$VENTRICLE" && -z "$GROUP_ID" ]]; then
    valid_ventricles=$(get_ventricles | tr '\n' '|' | sed 's/|$//')
    if [[ ! "$VENTRICLE" =~ ^($valid_ventricles)$ ]]; then
        echo "ERROR: --ventricle must be one of: $(get_ventricles | tr '\n' ', ' | sed 's/, $//')"
        exit 1
    fi
fi

# Set defaults for downstream module range
if [[ "$DOWNSTREAM_ONLY" == "true" && -z "$DOWNSTREAM_MODULES" ]]; then
    if [[ -n "$DOWNSTREAM_START" && -z "$DOWNSTREAM_STOP" ]]; then
        DOWNSTREAM_STOP=11
    fi
    if [[ -z "$DOWNSTREAM_START" && -n "$DOWNSTREAM_STOP" ]]; then
        DOWNSTREAM_START=0
    fi
    if [[ -z "$DOWNSTREAM_START" && -z "$DOWNSTREAM_STOP" ]]; then
        DOWNSTREAM_START=0
        DOWNSTREAM_STOP=11
    fi
fi

# ---------------------------------------------------------------------------
# VALIDATE SAMPLESHEET
# ---------------------------------------------------------------------------

echo "============================================================================"
echo "scRNA-seq Pipeline"
echo "============================================================================"
echo ""
echo "Project root: ${BASE_DIR}"
echo "Samplesheet: ${SAMPLESHEET}"
echo ""

if ! validate_samplesheet; then
    exit 1
fi

if [[ $PIPELINE_START -le 1 && "$DOWNSTREAM_ONLY" != "true" ]]; then
    if ! validate_fastq_samplesheet; then
        echo "ERROR: FASTQ samplesheet required for Cell Ranger step"
        exit 1
    fi
else
    validate_fastq_samplesheet || true
fi

print_samplesheet_summary

if [[ "$VALIDATE_ONLY" == "true" ]]; then
    echo "Samplesheet validation complete. Exiting."
    exit 0
fi

# ---------------------------------------------------------------------------
# SETUP OUTPUT DIRECTORIES BASED ON DATASET NAME
# ---------------------------------------------------------------------------

DATASET_NAME=$(get_dataset_name)

if [[ -z "$DATASET_NAME" ]]; then
    echo "ERROR: dataset_name column is empty in samplesheet.csv"
    exit 1
fi

# Define output directory structure
OUTPUT_BASE="${BASE_DIR}/Output_dir_${DATASET_NAME}"
PREPROCESS_DIR="${OUTPUT_BASE}/Single_cell_preprocessed"
DOWNSTREAM_DIR="${OUTPUT_BASE}/Single_cell_clustering"

echo ""
echo "Output configuration:"
echo "  Dataset name:    ${DATASET_NAME}"
echo "  Output base:     ${OUTPUT_BASE}"
echo "  Preprocessing:   ${PREPROCESS_DIR}"
echo "  Downstream:      ${DOWNSTREAM_DIR}"
echo ""

# Export for child scripts
export BASE_DIR
export SAMPLESHEET
export FASTQ_SAMPLESHEET
export PROJECT_ROOT="${BASE_DIR}"
export DATASET_NAME
export OUTPUT_BASE
export PREPROCESS_DIR
export DOWNSTREAM_DIR

# Get samples based on group or ventricle filter
if [[ -n "$GROUP_ID" ]]; then
    SAMPLES=($(get_samples_by_group "$GROUP_ID"))
    echo "Filtered to group $GROUP_ID ($GROUP_LABEL) samples: ${SAMPLES[*]}"
elif [[ -n "$VENTRICLE" ]]; then
    SAMPLES=($(get_samples_by_ventricle "$VENTRICLE"))
    echo "Filtered to $VENTRICLE samples: ${SAMPLES[*]}"
else
    SAMPLES=($(get_samples))
fi

N_SAMPLES=${#SAMPLES[@]}
echo "Processing ${N_SAMPLES} samples."
echo ""

# Export group and ventricle filter for downstream scripts
export GROUP_ID="${GROUP_ID}"
export GROUP_LABEL="${GROUP_LABEL}"
export VENTRICLE_FILTER="${VENTRICLE}"

# ---------------------------------------------------------------------------
# CREATE DIRECTORIES
# ---------------------------------------------------------------------------

echo "Creating directory structure..."

# Create preprocessing directories
mkdir -p "${PREPROCESS_DIR}/1_CellRanger_output"
mkdir -p "${PREPROCESS_DIR}/2_DropletQC_output"
mkdir -p "${PREPROCESS_DIR}/3_qClus_empty_droplets"
mkdir -p "${PREPROCESS_DIR}/4_Vaeda_doublet_detection"
mkdir -p "${PREPROCESS_DIR}/5_Seurat_conversion"
mkdir -p "${PREPROCESS_DIR}/6_Doublet_consensus"
mkdir -p "${PREPROCESS_DIR}/7_scCDC_correction"
mkdir -p "${PREPROCESS_DIR}/8_DecontX_correction"
mkdir -p "${PREPROCESS_DIR}/9_CHOIR_integration/Individual"
mkdir -p "${PREPROCESS_DIR}/9_CHOIR_integration/All_combined"
mkdir -p "${PREPROCESS_DIR}/10_QC_visualization"

# Create group-specific CHOIR integration directories from samplesheet
if [[ -n "$(get_column_index 'group_label')" ]]; then
    tail -n +2 "$SAMPLESHEET" | awk -F',' -v lc="$(get_column_index 'group_label')" '{print $lc}' | sort -u | while read -r glabel; do
        if [[ -n "$glabel" ]]; then
            mkdir -p "${PREPROCESS_DIR}/9_CHOIR_integration/${glabel}_only"
        fi
    done
else
    # Fallback to ventricle-based directories
    for vent in $(get_ventricles); do
        mkdir -p "${PREPROCESS_DIR}/9_CHOIR_integration/${vent}_only"
    done
fi

# Create downstream directories based on available groups
if [[ -n "$(get_column_index 'group_id')" && -n "$(get_column_index 'group_label')" ]]; then
    tail -n +2 "$SAMPLESHEET" | awk -F',' -v gc="$(get_column_index 'group_id')" -v lc="$(get_column_index 'group_label')" '{print $lc}' | sort -u | while read -r glabel; do
        if [[ -n "$glabel" ]]; then
            mkdir -p "${DOWNSTREAM_DIR}/11_Downstream_Analysis_${glabel}"
        fi
    done
else
    # Fallback to ventricle-based directories
    for vent in $(get_ventricles); do
        mkdir -p "${DOWNSTREAM_DIR}/11_Downstream_Analysis_${vent}"
    done
fi

# Create log directories (keep in project root for easy access)
mkdir -p "${BASE_DIR}/logs/cellranger"
mkdir -p "${BASE_DIR}/logs/dropletqc"
mkdir -p "${BASE_DIR}/logs/qclus"
mkdir -p "${BASE_DIR}/logs/vaeda"
mkdir -p "${BASE_DIR}/logs/h5ad2rds"
mkdir -p "${BASE_DIR}/logs/doublets"
mkdir -p "${BASE_DIR}/logs/scCDC"
mkdir -p "${BASE_DIR}/logs/decontx"
mkdir -p "${BASE_DIR}/logs/CHOIR"
mkdir -p "${BASE_DIR}/logs/qc_viz"
mkdir -p "${BASE_DIR}/logs/downstream"

# Create sample-specific directories
for SAMPLE in "${SAMPLES[@]}"; do
    mkdir -p "${PREPROCESS_DIR}/2_DropletQC_output/${SAMPLE}"
    mkdir -p "${PREPROCESS_DIR}/3_qClus_empty_droplets/${SAMPLE}"
    mkdir -p "${PREPROCESS_DIR}/4_Vaeda_doublet_detection/${SAMPLE}"
    mkdir -p "${PREPROCESS_DIR}/5_Seurat_conversion/${SAMPLE}"
    mkdir -p "${PREPROCESS_DIR}/6_Doublet_consensus/${SAMPLE}"
    mkdir -p "${PREPROCESS_DIR}/7_scCDC_correction/${SAMPLE}"
    mkdir -p "${PREPROCESS_DIR}/8_DecontX_correction/${SAMPLE}"
    mkdir -p "${PREPROCESS_DIR}/9_CHOIR_integration/Individual/${SAMPLE}"
done

# Create config directory
mkdir -p "${CONFIG_DIR}"

echo "Created directory structure for ${N_SAMPLES} samples."
echo ""

# ---------------------------------------------------------------------------
# PRINT EXECUTION PLAN
# ---------------------------------------------------------------------------

echo "============================================================================"
echo "Execution Plan"
echo "============================================================================"
echo ""
echo "Samples: ${SAMPLES[*]}"
if [[ -n "$GROUP_ID" ]]; then
    echo "Group filter: $GROUP_ID ($GROUP_LABEL)"
elif [[ -n "$VENTRICLE" ]]; then
    echo "Ventricle filter: $VENTRICLE"
fi

if [[ "$DOWNSTREAM_ONLY" == "true" ]]; then
    echo ""
    echo "Downstream Analysis Configuration:"
    if [[ -n "$GROUP_LABEL" ]]; then
        echo "  Group ID:    ${GROUP_ID}"
        echo "  Group Label: ${GROUP_LABEL}"
        echo "  Output: ${DOWNSTREAM_DIR}/11_Downstream_Analysis_${GROUP_LABEL}/"
    else
        echo "  Ventricle: ${VENTRICLE}"
        echo "  Output: ${DOWNSTREAM_DIR}/11_Downstream_Analysis_${VENTRICLE}/"
    fi
    echo "  GPU: A100 partition (gpu1day QoS)"

    if [[ -n "$DOWNSTREAM_MODULES" ]]; then
        echo "  Mode: Specific modules"
        echo "  Modules: $DOWNSTREAM_MODULES"
        echo ""
        echo "  Modules to run:"
        IFS=',' read -ra MOD_ARRAY <<< "$DOWNSTREAM_MODULES"
        for mod in "${MOD_ARRAY[@]}"; do
            if [[ $mod -ge 0 && $mod -lt ${#DOWNSTREAM_MODULE_NAMES[@]} ]]; then
                printf "    %2d - %s\n" "$mod" "${DOWNSTREAM_MODULE_NAMES[$mod]}"
            fi
        done
    else
        echo "  Mode: Sequential range"
        echo "  Start module: $DOWNSTREAM_START"
        echo "  Stop module: $DOWNSTREAM_STOP"
        echo ""
        echo "  Modules to run:"
        for mod in $(seq $DOWNSTREAM_START $DOWNSTREAM_STOP); do
            if [[ $mod -ge 0 && $mod -lt ${#DOWNSTREAM_MODULE_NAMES[@]} ]]; then
                printf "    %2d - %s\n" "$mod" "${DOWNSTREAM_MODULE_NAMES[$mod]}"
            fi
        done
    fi
else
    echo "Steps: ${PIPELINE_START} to ${PIPELINE_END}"
    echo ""
    for step_num in $(seq $PIPELINE_START $PIPELINE_END); do
        step_idx=$((step_num - 1))
        if [[ $step_idx -lt ${#ALL_STEP_NAMES[@]} ]]; then
            echo "  Step ${step_num}: ${ALL_STEP_NAMES[$step_idx]}"
        fi
    done
fi
echo ""

# ---------------------------------------------------------------------------
# SUBMIT JOBS
# ---------------------------------------------------------------------------

echo "============================================================================"
echo "Submitting Jobs"
echo "============================================================================"
echo ""

PREV_JOB_ID=""
SUBMITTED_JOBS=()

# Submit preprocessing steps (1-10)
if [[ "$DOWNSTREAM_ONLY" != "true" ]]; then
    for step_num in $(seq $PIPELINE_START $((PIPELINE_END < 10 ? PIPELINE_END : 10))); do
        if [[ $step_num -ge 1 && $step_num -le 10 ]]; then
            step_idx=$((step_num - 1))
            script="${PREPROCESS_SCRIPTS[$step_idx]}"
            step_name="${PREPROCESS_NAMES[$step_idx]}"

            echo "Step ${step_num}: ${step_name}"
            echo "  Script: $(basename "$script")"

            SBATCH_CMD="sbatch"

            if [[ -n "$PREV_JOB_ID" ]]; then
                SBATCH_CMD="$SBATCH_CMD --dependency=afterok:${PREV_JOB_ID}"
                echo "  Dependency: afterok:${PREV_JOB_ID}"
            fi

            # Steps 1-8 are per-sample array jobs; steps 9-10 are single jobs
            if [[ $step_num -ge 1 && $step_num -le 8 ]]; then
                SBATCH_CMD="$SBATCH_CMD --array=1-${N_SAMPLES}"
                echo "  Array job: 1-${N_SAMPLES} samples"
            fi

            SBATCH_CMD="$SBATCH_CMD --export=ALL,SAMPLESHEET=${SAMPLESHEET},PROJECT_ROOT=${BASE_DIR},PREPROCESS_DIR=${PREPROCESS_DIR},DOWNSTREAM_DIR=${DOWNSTREAM_DIR},DATASET_NAME=${DATASET_NAME}"
            SBATCH_CMD="$SBATCH_CMD $script"

            if [[ "$DRY_RUN" == "true" ]]; then
                echo "  [DRY-RUN] Would submit: $SBATCH_CMD"
                PREV_JOB_ID="DRYRUN_${step_num}"
            else
                if [[ -f "$script" ]]; then
                    JOB_OUTPUT=$($SBATCH_CMD)
                    JOB_ID=$(echo "$JOB_OUTPUT" | grep -oP '\d+')
                    PREV_JOB_ID="$JOB_ID"
                    SUBMITTED_JOBS+=("$JOB_ID")
                    echo "  Submitted: Job ID ${JOB_ID}"
                else
                    echo "  SKIPPED: Script not found"
                fi
            fi
            echo ""
        fi
    done
fi

# Submit downstream analysis (GPU)
if [[ "$DOWNSTREAM_ONLY" == "true" || ($PIPELINE_END -ge 11 && $PIPELINE_START -le 21) ]]; then
    echo "Downstream Analysis Pipeline (GPU: A100)"

    if [[ -z "$GROUP_ID" && -z "$VENTRICLE" ]]; then
        echo "  WARNING: No group or ventricle specified. Skipping downstream."
        echo "  Use: $0 --downstream --1"
        echo "  Or:  $0 --downstream --ventricle LV"
    else
        echo "  Script: $(basename "$DOWNSTREAM_SCRIPT")"

        if [[ -n "$GROUP_LABEL" ]]; then
            echo "  Group: $GROUP_ID ($GROUP_LABEL)"
            echo "  Output: ${DOWNSTREAM_DIR}/11_Downstream_Analysis_${GROUP_LABEL}/"
        else
            echo "  Ventricle: $VENTRICLE"
            echo "  Output: ${DOWNSTREAM_DIR}/11_Downstream_Analysis_${VENTRICLE}/"
        fi
        echo "  Partition: a100 (GPU)"

        SBATCH_CMD="sbatch"

        if [[ -n "$PREV_JOB_ID" && "$PREV_JOB_ID" != DRYRUN_* ]]; then
            SBATCH_CMD="$SBATCH_CMD --dependency=afterok:${PREV_JOB_ID}"
            echo "  Dependency: afterok:${PREV_JOB_ID}"
        fi

        # Export all group-related variables
        SBATCH_CMD="$SBATCH_CMD --export=ALL,SAMPLESHEET=${SAMPLESHEET},PROJECT_ROOT=${BASE_DIR},PREPROCESS_DIR=${PREPROCESS_DIR},DOWNSTREAM_DIR=${DOWNSTREAM_DIR},DATASET_NAME=${DATASET_NAME},GROUP_ID=${GROUP_ID},GROUP_LABEL=${GROUP_LABEL},VENTRICLE_FILTER=${VENTRICLE}"

        # Determine output label for downstream script
        if [[ -n "$GROUP_LABEL" ]]; then
            OUTPUT_LABEL="${GROUP_LABEL}"
        else
            OUTPUT_LABEL="${VENTRICLE}"
        fi

        DOWNSTREAM_ARGS="--group-label ${OUTPUT_LABEL}"

        if [[ -n "$DOWNSTREAM_MODULES" ]]; then
            DOWNSTREAM_ARGS="$DOWNSTREAM_ARGS --modules ${DOWNSTREAM_MODULES}"
            echo "  Modules: $DOWNSTREAM_MODULES"
        else
            DOWNSTREAM_ARGS="$DOWNSTREAM_ARGS --start ${DOWNSTREAM_START} --stop ${DOWNSTREAM_STOP}"
            echo "  Module range: ${DOWNSTREAM_START} to ${DOWNSTREAM_STOP}"
        fi

        SBATCH_CMD="$SBATCH_CMD $DOWNSTREAM_SCRIPT $DOWNSTREAM_ARGS"

        if [[ "$DRY_RUN" == "true" ]]; then
            echo "  [DRY-RUN] Would submit: $SBATCH_CMD"
        else
            if [[ -f "$DOWNSTREAM_SCRIPT" ]]; then
                JOB_OUTPUT=$($SBATCH_CMD)
                JOB_ID=$(echo "$JOB_OUTPUT" | grep -oP '\d+')
                SUBMITTED_JOBS+=("$JOB_ID")
                echo "  Submitted: Job ID ${JOB_ID}"
            else
                echo "  SKIPPED: Script not found at $DOWNSTREAM_SCRIPT"
            fi
        fi
    fi
    echo ""
fi

# ---------------------------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------------------------

echo "============================================================================"
echo "Submission Summary"
echo "============================================================================"
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "DRY RUN - No jobs were actually submitted."
    echo ""
    echo "To submit for real, remove the --dry-run flag."
else
    if [[ ${#SUBMITTED_JOBS[@]} -gt 0 ]]; then
        echo "Submitted ${#SUBMITTED_JOBS[@]} job(s):"
        for job in "${SUBMITTED_JOBS[@]}"; do
            echo "  - ${job}"
        done
        echo ""
        echo "Monitor with: squeue -u \$USER"
        echo "Cancel with:  scancel <job_id>"
    else
        echo "No jobs were submitted."
    fi
fi

echo ""
echo "Output directories:"
echo "  Project root:   ${BASE_DIR}"
echo "  Output base:    ${OUTPUT_BASE}"
echo "  Preprocessing:  ${PREPROCESS_DIR}/"
if [[ -n "$GROUP_LABEL" ]]; then
    echo "  Downstream:     ${DOWNSTREAM_DIR}/11_Downstream_Analysis_${GROUP_LABEL}/"
elif [[ -n "$VENTRICLE" ]]; then
    echo "  Downstream:     ${DOWNSTREAM_DIR}/11_Downstream_Analysis_${VENTRICLE}/"
else
    echo "  Downstream:     ${DOWNSTREAM_DIR}/"
fi
echo ""
echo "============================================================================"
echo "Pipeline submission complete!"
echo "============================================================================"
