#!/usr/bin/env bash

# Submit with runner/runner.py, not by hand.
# Required env vars:
#   DIAZODB_JOB_ID   - API job UUID
#   DIAZODB_INPUT    - input FASTA path in the job workspace
#   DIAZODB_OUTDIR   - output directory for intermediate files
#   DIAZODB_OUTPUT   - final CSV/HTML/PDF path to post back to the API

#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH -J diazodb_classify
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err

set -euo pipefail

: "${DIAZODB_JOB_ID:?DIAZODB_JOB_ID is required}"
: "${DIAZODB_INPUT:?DIAZODB_INPUT is required}"
: "${DIAZODB_OUTDIR:?DIAZODB_OUTDIR is required}"
: "${DIAZODB_OUTPUT:?DIAZODB_OUTPUT is required}"

mkdir -p "$DIAZODB_OUTDIR"

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID      : $DIAZODB_JOB_ID"
echo "Input       : $DIAZODB_INPUT"
echo "Output      : $DIAZODB_OUTPUT"
echo "Slurm ID    : ${SLURM_JOBID:-local}"
echo "====================================================="

if [[ -n "${DIAZODB_CONDA_ENV:-}" ]]; then
  eval "$(conda shell.bash hook)"
  conda activate "$DIAZODB_CONDA_ENV"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_SCRIPT="${DIAZODB_PIPELINE_SCRIPT:-$SCRIPT_DIR/diazodb_classify_pipeline.sh}"

bash "$PIPELINE_SCRIPT" "$DIAZODB_INPUT" "$DIAZODB_OUTDIR" "$DIAZODB_OUTPUT"

test -s "$DIAZODB_OUTPUT"

echo "====================================================="
echo "End Time    : $(date)"
echo "====================================================="
