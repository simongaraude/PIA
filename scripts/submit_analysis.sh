#!/bin/bash
# submit_analysis.sh — Submit complex-analyzer as a Slurm job
#
# Usage:
#   /shared/complex-analyzer/scripts/submit_analysis.sh [INPUT_DIR] [OPTIONS]
#
# Examples:
#   # Default: analyse /shared/all_structures, output CSV there
#   /shared/complex-analyzer/scripts/submit_analysis.sh
#
#   # Custom input directory
#   /shared/complex-analyzer/scripts/submit_analysis.sh /shared/my_results/structures
#
#   # Custom input + output + parental reference
#   /shared/complex-analyzer/scripts/submit_analysis.sh /shared/my_results/structures \
#       --parental_pdb /shared/references/parental.pdb
#
# The CSV is written next to the input structures and optionally synced to Google Drive.

set -e

# =====================================================================
# Configuration
# =====================================================================
CONDA_SH="/shared/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="boltzgen"
RCLONE_BIN="/shared/rclone"
RCLONE_CONF="/shared/config/rclone.conf"
GDRIVE_REMOTE="alceus_gdrive"
GDRIVE_FOLDER="analysis_results_v2"

DEFAULT_INPUT="/shared/all_structures"
DEFAULT_BINDER_CHAIN="B"
DEFAULT_RECEPTOR_CHAIN="A"

# =====================================================================
# Parse arguments
# =====================================================================
INPUT_DIR="${1:-$DEFAULT_INPUT}"
shift 2>/dev/null || true   # remaining args passed through to complex-analyzer
EXTRA_ARGS="$@"

# Resolve output path
OUTPUT_CSV="${INPUT_DIR}/complex_metrics.csv"

echo "=========================================="
echo "  Complex Analyzer — Job Submission"
echo "=========================================="
echo ""

# Validate input directory
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR"
    echo ""
    echo "Usage: $0 [INPUT_DIR] [extra complex-analyzer flags]"
    echo ""
    echo "Examples:"
    echo "  $0                                          # uses $DEFAULT_INPUT"
    echo "  $0 /shared/my_structures"
    echo "  $0 /shared/my_structures --parental_pdb /shared/ref/parental.pdb"
    exit 1
fi

# Count structures
N_FILES=$(ls "$INPUT_DIR"/*.pdb "$INPUT_DIR"/*.cif "$INPUT_DIR"/*.mmcif "$INPUT_DIR"/*.ent 2>/dev/null | wc -l)
if [ "$N_FILES" -eq 0 ]; then
    echo "Error: No PDB/CIF files found in $INPUT_DIR"
    exit 1
fi

# Estimate time (~30-60s per structure on CPU)
EST_MINUTES=$(( (N_FILES * 45 + 59) / 60 ))
EST_HOURS=$(( (EST_MINUTES + 59) / 60 ))

echo "Configuration:"
echo "  Input dir:     $INPUT_DIR"
echo "  Output CSV:    $OUTPUT_CSV"
echo "  Structures:    $N_FILES"
echo "  Binder chain:  $DEFAULT_BINDER_CHAIN"
echo "  Receptor chain:$DEFAULT_RECEPTOR_CHAIN"
echo "  Extra args:    ${EXTRA_ARGS:-none}"
echo "  Est. runtime:  ~${EST_MINUTES} min (~${EST_HOURS}h)"
echo ""

# =====================================================================
# Create job directory & script
# =====================================================================
JOB_DIR="${INPUT_DIR}/analysis_job"
mkdir -p "$JOB_DIR"

JOB_SCRIPT="${JOB_DIR}/analyze_job.sh"

cat > "$JOB_SCRIPT" << JOBEOF
#!/bin/bash
#SBATCH --job-name=complex-analyzer
#SBATCH --output=${JOB_DIR}/slurm_%j.out
#SBATCH --error=${JOB_DIR}/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=${EST_HOURS}:30:00

echo "=========================================="
echo "  Complex Analyzer — Running"
echo "=========================================="
echo "Node:       \$(hostname)"
echo "Date:       \$(date)"
echo "Input:      ${INPUT_DIR}"
echo "Structures: ${N_FILES}"
echo ""

# Activate conda environment
source ${CONDA_SH}
conda activate ${CONDA_ENV}

# Run analysis
echo "Starting analysis..."
START=\$(date +%s)

complex-analyzer \\
    --input_dir "${INPUT_DIR}" \\
    --output "${OUTPUT_CSV}" \\
    --binder_chain ${DEFAULT_BINDER_CHAIN} \\
    --receptor_chain ${DEFAULT_RECEPTOR_CHAIN} \\
    ${EXTRA_ARGS}

EXIT_CODE=\$?
END=\$(date +%s)
ELAPSED=\$(( END - START ))

echo ""
echo "Finished in \${ELAPSED}s (exit code: \${EXIT_CODE})"

# Sync results to Google Drive
if [ \$EXIT_CODE -eq 0 ] && [ -f "${RCLONE_BIN}" ] && [ -f "${RCLONE_CONF}" ]; then
    echo ""
    echo "Syncing results to Google Drive..."
    ${RCLONE_BIN} --config ${RCLONE_CONF} copy \\
        "${OUTPUT_CSV}" \\
        ${GDRIVE_REMOTE}:${GDRIVE_FOLDER}/ \\
        2>&1 | grep -v "Failed to save config" || true
    echo "Synced to ${GDRIVE_REMOTE}:${GDRIVE_FOLDER}/complex_metrics.csv"
fi

echo ""
echo "=========================================="
echo "  Done!"
echo "=========================================="
echo "Results: ${OUTPUT_CSV}"
JOBEOF

chmod +x "$JOB_SCRIPT"

# =====================================================================
# Submit
# =====================================================================
echo "Submitting job to Slurm..."
JOB_ID=$(sbatch "$JOB_SCRIPT" | awk '{print $4}')

echo "Job submitted!"
echo ""
echo "  Job ID:      $JOB_ID"
echo "  Output CSV:  $OUTPUT_CSV"
echo "  Job script:  $JOB_SCRIPT"
echo "  Logs:        ${JOB_DIR}/slurm_${JOB_ID}.out"
echo ""
echo "Monitoring commands:"
echo "  Check status:  squeue"
echo "  View logs:     tail -f ${JOB_DIR}/slurm_${JOB_ID}.out"
echo "  Cancel job:    scancel ${JOB_ID}"
echo ""
echo "When done, results at: $OUTPUT_CSV"
