#!/bin/bash
# install_analyzer.sh
# Run on the headnode (not via sbatch — no GPU needed).
# Installs complex-analyzer into the existing /shared/miniconda3/envs/boltzgen env.

set -e

echo "=========================================="
echo "  Installing complex-analyzer"
echo "=========================================="

# Check we're on the headnode
if [ ! -d "/shared" ]; then
    echo "Error: /shared not found. Are you on the ParallelCluster headnode?"
    exit 1
fi

# Activate conda
source /shared/miniconda3/etc/profile.d/conda.sh
conda activate boltzgen

echo "[1/3] Installing dependencies..."
pip install biopython pandas numpy freesasa prodigy-prot 2>/dev/null || \
pip install biopython pandas numpy prodigy-prot  # freesasa optional

# Install DSSP for secondary structure
echo "[2/3] Installing DSSP..."
sudo apt-get update -qq && sudo apt-get install -y dssp 2>/dev/null || \
    echo "Warning: DSSP not installed (secondary structure will use CA heuristic)"

echo "[3/3] Installing complex-analyzer..."

# If repo is already in /shared, install from there; otherwise install from tarball
if [ -d "/shared/complex-analyzer" ]; then
    pip install -e /shared/complex-analyzer
elif [ -f "/shared/complex-analyzer.tar.gz" ]; then
    cd /shared
    tar xzf complex-analyzer.tar.gz
    pip install -e /shared/complex-analyzer
else
    echo "Error: Place the repo at /shared/complex-analyzer/ or tarball at /shared/complex-analyzer.tar.gz"
    exit 1
fi

echo ""
echo "=========================================="
echo "  Installation complete!"
echo "=========================================="
echo ""
echo "Test it:"
echo "  source /shared/miniconda3/etc/profile.d/conda.sh && conda activate boltzgen"
echo "  complex-analyzer --help"
echo ""
echo "Submit analysis job:"
echo "  /shared/complex-analyzer/scripts/submit_analysis.sh"
echo "  /shared/complex-analyzer/scripts/submit_analysis.sh /path/to/structures/"
