#!/usr/bin/env bash

# Submit with runner/runner.py, not by hand.
# Required env vars:
#   DIAZODB_JOB_ID   - API job UUID
#   DIAZODB_INPUT    - input FASTA path in the job workspace
#   DIAZODB_OUTDIR   - output directory for intermediate files
#   DIAZODB_OUTPUT   - final CSV/HTML/PDF path to post back to the API

#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=4GB
#SBATCH -J diazodb_classify
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Prefer env vars injected by runner; fall back to positional args if provided.
DIAZODB_JOB_ID="${DIAZODB_JOB_ID:-${1:-}}"
INPUT_FASTA="${DIAZODB_INPUT:-${2:-}}"
OUTDIR="${DIAZODB_OUTDIR:-${3:-}}"
FINAL_OUTPUT="${DIAZODB_OUTPUT:-${4:-}}"

if [[ -z "$DIAZODB_JOB_ID" || -z "$INPUT_FASTA" || -z "$OUTDIR" || -z "$FINAL_OUTPUT" ]]; then
  echo "Missing required inputs. Provide DIAZODB_JOB_ID, DIAZODB_INPUT, DIAZODB_OUTDIR, DIAZODB_OUTPUT (or positional args)." >&2
  exit 2
fi

mkdir -p "$OUTDIR"

HMMSEARCH_BIN="${DIAZODB_HMMSEARCH_BIN:-hmmsearch}"
PRODIGAL_BIN="${DIAZODB_PRODIGAL_BIN:-prodigal}"
USE_PRODIGAL="${DIAZODB_USE_PRODIGAL:-false}"

QUERY_FASTA="$INPUT_FASTA"
HMM_DB="$REPO_ROOT/HMMs/combined_nif_03192026.hmm"

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID      : $DIAZODB_JOB_ID"
echo "Input       : $INPUT_FASTA"
echo "Output      : $FINAL_OUTPUT"
echo "Slurm ID    : ${SLURM_JOBID:-local}"
echo "====================================================="

if [[ -n "${DIAZODB_CONDA_ENV:-}" ]]; then
  eval "$(conda shell.bash hook)"
  conda activate "$DIAZODB_CONDA_ENV"
fi

module load mafft/7.505-gcc-13.2.0-nklkvtc

# if [[ "$USE_PRODIGAL" == "true" ]]; then
#   QUERY_FASTA="$OUTDIR/predicted_proteins.faa"
#   "$PRODIGAL_BIN" \
#     -i "$INPUT_FASTA" \
#     -a "$QUERY_FASTA" \
#     -p meta \
#     -q
# fi

# HMM search
"$HMMSEARCH_BIN" "$HMM_DB" "$QUERY_FASTA" > "$OUTDIR/hmmsearch.out"
#chgrp hpc_enviromics ../results/archaea/hmmsearch_results/*

# Parse HMM
python "$SCRIPT_DIR/Parse_hmm_results.py" --hits "$OUTDIR/hmmsearch.out" --outdir "$OUTDIR/hmmsearch_results"
python "$SCRIPT_DIR/Parse_tophits.py" --hits "$OUTDIR/hmmsearch_results/hits.feather" --outdir "$OUTDIR/hmmsearch_results"

# Emit final result for the runner
cp "$OUTDIR/hmmsearch_results/tophits.csv" "$FINAL_OUTPUT"
test -s "$FINAL_OUTPUT"


#SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#PIPELINE_SCRIPT="${DIAZODB_PIPELINE_SCRIPT:-$SCRIPT_DIR/diazodb_classify_pipeline.sh}"

#bash "$PIPELINE_SCRIPT" "$DIAZODB_INPUT" "$DIAZODB_OUTDIR" "$DIAZODB_OUTPUT"

#test -s "$DIAZODB_OUTPUT"

echo "====================================================="
echo "End Time    : $(date)"
echo "====================================================="
