================================================================================
kaiser_test_py3.11 — Installation Guide
================================================================================

Environment: Python 3.11 + R 4.4+ + Julia 1.10–1.11.x
GPU support: CUDA 12.1 (PyTorch, TensorFlow, RAPIDS, DGL)

--------------------------------------------------------------------------------
Quick Start (4 steps)
--------------------------------------------------------------------------------

  1. CONDA — base environment (Python, R, Julia, conda-only packages):

       conda env create -f kaiser_test_py3.11.yml
       conda activate kaiser_test_py3.11

  2. R — GitHub-only R packages (DropletQC, DoubletFinder, CHOIR, etc.):

       Rscript install_r_packages.R

  3. TOOLS — uv, Python pip packages, scAURA, Julia/scICE:

       bash install_tools.sh

     This script handles everything:
       - Installs uv (fast pip replacement)
       - Runs uv pip install -r requirements.txt
       - Clones & configures scAURA (GitHub)
       - Clones & configures scICE/scLENS (Julia packages)
       - Sets up PyCall for Julia <-> Python interop

  4. VERIFY:

       python -c "import scanpy, scvi, torch; print('Python OK')"
       Rscript -e "library(Seurat); message('R OK')"
       julia -e "println(\"Julia OK\")"

--------------------------------------------------------------------------------
Files
--------------------------------------------------------------------------------

  kaiser_test_py3.11.yml   Conda env definition (Python + R + Julia + conda pkgs)
  requirements.txt         Python pip packages (installed via uv)
  install_r_packages.R     R packages from CRAN, Bioconductor, GitHub
  install_tools.sh         uv + pip deps + scAURA + Julia/scICE setup

--------------------------------------------------------------------------------
No GPU?
--------------------------------------------------------------------------------

Remove these sections from requirements.txt before step 3:
  - PyTorch + CUDA 12.1 → replace with CPU-only torch
  - RAPIDS (GPU-accelerated analytics)
  - DGL → replace with CPU-only dgl

--------------------------------------------------------------------------------
Manual uv install (alternative to install_tools.sh for pip packages only)
--------------------------------------------------------------------------------

  pip install uv
  uv pip install -r requirements.txt
