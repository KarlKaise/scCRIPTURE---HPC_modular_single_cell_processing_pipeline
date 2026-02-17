#!/usr/bin/env bash
# ==============================================================================
# install_tools.sh
# ==============================================================================
# Post-install script for:
#   - uv (fast pip replacement) + Python pip dependencies
#   - scAURA (GitHub clone)
#   - Julia/scICE (GitHub clone + Julia packages)
#
# Run AFTER creating the conda environment and R post-install:
#   conda activate kaiser_test_py3.11
#   bash install_tools.sh
#
# This script can be re-run safely (idempotent).
# ==============================================================================

set -euo pipefail

# --- Configuration -----------------------------------------------------------
TOOLS_DIR="${TOOLS_DIR:-$(pwd)/tools}"
SCAURA_DIR="${TOOLS_DIR}/scAURA"
SCICE_DIR="${TOOLS_DIR}/scICE"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================================="
echo "  External Tools Installer                                    "
echo "============================================================="
echo "Tools directory: ${TOOLS_DIR}"
echo ""

mkdir -p "${TOOLS_DIR}"

# ==============================================================================
# 0. Install uv + Python pip dependencies
# ==============================================================================
# uv is a fast drop-in replacement for pip. We install it via pip, then use
# it for all subsequent Python package installations.
# ==============================================================================

echo "--- [0/3] uv + Python pip dependencies ---"

if ! command -v uv &> /dev/null; then
    echo "Installing uv..."
    pip install uv
else
    echo "uv already installed: $(uv --version)"
fi

# Install Python pip dependencies from requirements.txt
REQUIREMENTS_FILE="${SCRIPT_DIR}/requirements.txt"
if [ -f "${REQUIREMENTS_FILE}" ]; then
    echo "Installing Python pip packages via uv..."
    echo "  Source: ${REQUIREMENTS_FILE}"
    uv pip install -r "${REQUIREMENTS_FILE}"
    echo "Python pip packages installed."
else
    echo "Warning: ${REQUIREMENTS_FILE} not found."
    echo "  Looked in: ${SCRIPT_DIR}"
    echo "  If requirements.txt is elsewhere, run manually:"
    echo "    uv pip install -r /path/to/requirements.txt"
fi

echo ""

# ==============================================================================
# 1. scAURA - Graph Contrastive Learning for scRNA-seq Clustering
# ==============================================================================
# Repository: https://github.com/bozdaglab/scAURA
# Uses: PyTorch, DGL, scanpy (all already installed above)
# GPU version (scAURA_gpu.py) for datasets >= 2,500 cells
# CPU version (scAURA.py) for smaller datasets
# ==============================================================================

echo "--- [1/3] scAURA ---"

if [ -d "${SCAURA_DIR}" ]; then
    echo "scAURA already cloned at ${SCAURA_DIR}"
    echo "Updating..."
    cd "${SCAURA_DIR}"
    git pull --ff-only || echo "Warning: git pull failed, using existing version"
else
    echo "Cloning scAURA..."
    git clone https://github.com/bozdaglab/scAURA.git "${SCAURA_DIR}"
    cd "${SCAURA_DIR}"
fi

# Install scAURA's Python dependencies.
# Most should already be satisfied by the main requirements.txt, but this
# ensures any additional packages specific to scAURA are picked up.
if [ -f "requirements.txt" ]; then
    echo "Installing scAURA pip requirements via uv (skipping already-installed)..."
    uv pip install --no-deps -r requirements.txt 2>/dev/null || \
    uv pip install -r requirements.txt
    echo "scAURA requirements installed."
else
    echo "Warning: No requirements.txt found in scAURA repo."
fi

echo ""
echo "scAURA installed at: ${SCAURA_DIR}"
echo "  CPU version: python ${SCAURA_DIR}/scAURA.py"
echo "  GPU version: python ${SCAURA_DIR}/scAURA_gpu.py"
echo ""

# ==============================================================================
# 2. scICE / scLENS - Julia-based clustering (with CUDA GPU support)
# ==============================================================================
# Repository: https://github.com/Mathbiomed/scICE
# Requires: Julia <1.12.3 (installed via conda), CUDA drivers for GPU
# Key Julia packages: CUDA, CSV, scLENS, CairoMakie
# These are auto-installed via Pkg.instantiate() from the project environment.
# ==============================================================================

echo "--- [2/3] Julia scICE / scLENS ---"

# Verify Julia is available
if ! command -v julia &> /dev/null; then
    echo "ERROR: Julia not found in PATH."
    echo "  Julia should be installed via conda (see kaiser_test_py3.11.yml)."
    echo "  If conda solver failed for Julia, install manually:"
    echo "    curl -fsSL https://install.julialang.org | sh"
    echo "    # Choose version 1.11.x (avoid 1.12.3 - PyCall crash)"
    exit 1
fi

JULIA_VERSION=$(julia --version 2>/dev/null | grep -oP '\d+\.\d+\.\d+')
echo "Julia version: ${JULIA_VERSION}"

# Check Julia version compatibility
JULIA_MAJOR=$(echo "${JULIA_VERSION}" | cut -d. -f1)
JULIA_MINOR=$(echo "${JULIA_VERSION}" | cut -d. -f2)
JULIA_PATCH=$(echo "${JULIA_VERSION}" | cut -d. -f3)

if [ "${JULIA_MAJOR}" -eq 1 ] && [ "${JULIA_MINOR}" -eq 12 ] && [ "${JULIA_PATCH}" -ge 3 ]; then
    echo "WARNING: Julia ${JULIA_VERSION} detected. PyCall has a known crash"
    echo "  issue with Julia >=1.12.3. Consider downgrading to 1.11.x."
fi

# Clone scICE
if [ -d "${SCICE_DIR}" ]; then
    echo "scICE already cloned at ${SCICE_DIR}"
    echo "Updating..."
    cd "${SCICE_DIR}"
    git pull --ff-only || echo "Warning: git pull failed, using existing version"
else
    echo "Cloning scICE..."
    git clone https://github.com/Mathbiomed/scICE.git "${SCICE_DIR}"
    cd "${SCICE_DIR}"
fi

# Install Julia packages from the scICE project environment
echo "Installing Julia packages (CUDA, CSV, scLENS, CairoMakie, etc.)..."
echo "This may take several minutes on first install..."

julia --project="${SCICE_DIR}" -e '
    import Pkg
    Pkg.activate(".")
    Pkg.instantiate()
    println("\n--- Installed Julia packages ---")
    for (uuid, pkg) in Pkg.dependencies()
        if pkg.is_direct_dep
            println("  $(pkg.name) $(pkg.version)")
        end
    end
    println("\n--- CUDA status ---")
    try
        using CUDA
        if CUDA.functional()
            println("  CUDA available: $(CUDA.device())")
            println("  CUDA version: $(CUDA.runtime_version())")
        else
            println("  CUDA not functional (CPU-only mode)")
        end
    catch e
        println("  CUDA check failed: $(e)")
    end
'

echo ""
echo "scICE installed at: ${SCICE_DIR}"
echo "  Run: julia --project=${SCICE_DIR} your_script.jl"
echo ""

# ==============================================================================
# 3. Configure Julia <-> Python interop (PyCall)
# ==============================================================================

echo "--- [3/3] Configuring PyCall for current Python ---"

# Set PYTHON environment variable so PyCall uses the conda Python
CONDA_PYTHON=$(which python)
echo "Setting PyCall to use: ${CONDA_PYTHON}"

julia --project="${SCICE_DIR}" -e "
    ENV[\"PYTHON\"] = \"${CONDA_PYTHON}\"
    import Pkg
    Pkg.add(\"PyCall\")
    Pkg.build(\"PyCall\")
    using PyCall
    println(\"PyCall Python: \", PyCall.python)
" 2>/dev/null || echo "Warning: PyCall configuration optional - skip if not needed."

# ==============================================================================
# Summary
# ==============================================================================

echo ""
echo "============================================================="
echo "  Installation Summary                                        "
echo "============================================================="
echo ""
echo "  uv:        $(uv --version 2>/dev/null || echo 'not found')"
echo "  Python:    $(python --version 2>/dev/null || echo 'not found')"
echo "  R:         $(R --version 2>/dev/null | head -1 || echo 'not found')"
echo "  Julia:     $(julia --version 2>/dev/null || echo 'not found')"
echo ""
echo "  scAURA:    ${SCAURA_DIR}"
echo "    CPU:     python ${SCAURA_DIR}/scAURA.py"
echo "    GPU:     python ${SCAURA_DIR}/scAURA_gpu.py"
echo ""
echo "  scICE:     ${SCICE_DIR}"
echo "    Run:     julia --project=${SCICE_DIR} your_script.jl"
echo ""
echo "=== External tools installation complete ==="
