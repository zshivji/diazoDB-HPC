#!/usr/bin/env bash

# Submit with runner/runner.py, not by hand.
# Required env vars:
#   DIAZODB_JOB_ID   - API job UUID
#   DIAZODB_INPUT    - input FASTA path in the job workspace
#   DIAZODB_OUTDIR   - output directory for intermediate files
#   DIAZODB_OUTPUT   - final CSV/HTML/PDF path to post back to the API

#SBATCH --time=08:04:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40GB
#SBATCH -J diazodb_classify
#SBATCH -o /resnick/scratch/zshivji/diazoDB-HPC/slurm/slurm-%x-%j.out
#SBATCH -e /resnick/scratch/zshivji/diazoDB-HPC/slurm/slurm-%x-%j.err

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

mkdir -p "$OUTDIR" "$(dirname "$FINAL_OUTPUT")"

HMMSEARCH_BIN="${DIAZODB_HMMSEARCH_BIN:-hmmsearch}"
PRODIGAL_BIN="${DIAZODB_PRODIGAL_BIN:-prodigal}"
USE_PRODIGAL="${DIAZODB_USE_PRODIGAL:-false}"
DIAZODB_CONDA_ENV="${DIAZODB_CONDA_ENV:-/resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm}"

QUERY_FASTA="$INPUT_FASTA"
HMM_DB="${DIAZODB_HMM_PROFILE:-$REPO_ROOT/HMMs/combined_nif_07202026.hmm}"
JOB_HMM_DIR="$OUTDIR/hmmsearch"
JOB_PARSE_DIR="$OUTDIR/parse_hmm"
SHARED_HMM_DIR="$REPO_ROOT/results/hmmsearch_results/bacteria"
SHARED_ARCHAEA_DIR="$REPO_ROOT/results/hmmsearch_results/archaea"
SHARED_FINAL_RESULTS="$REPO_ROOT/results/final/nif_clusters.csv"
PIPELINE_LOCK="$REPO_ROOT/results/.runner-pipeline.lock"

if [[ ! -f "$INPUT_FASTA" ]]; then
  echo "Input FASTA does not exist: $INPUT_FASTA" >&2
  exit 2
fi
if [[ ! -f "$HMM_DB" ]]; then
  echo "HMM profile does not exist: $HMM_DB" >&2
  exit 2
fi

mkdir -p "$JOB_HMM_DIR" "$JOB_PARSE_DIR" "$SHARED_HMM_DIR" "$SHARED_ARCHAEA_DIR"

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID      : $DIAZODB_JOB_ID"
echo "Input       : $INPUT_FASTA"
echo "Output      : $FINAL_OUTPUT"
echo "Slurm ID    : ${SLURM_JOBID:-local}"
echo "====================================================="

eval "$(conda shell.bash hook)"
conda activate "$DIAZODB_CONDA_ENV"

module load mafft/7.505-gcc-13.2.0-nklkvtc

if [[ "$USE_PRODIGAL" == "true" ]]; then
  QUERY_FASTA="$OUTDIR/predicted_proteins.faa"
  "$PRODIGAL_BIN" \
    -i "$INPUT_FASTA" \
    -a "$QUERY_FASTA" \
    -p meta \
    -q
fi

INPUT_STEM="$(basename "$INPUT_FASTA")"
INPUT_STEM="${INPUT_STEM%.*}"
SAFE_INPUT_STEM="$(printf '%s' "$INPUT_STEM" | tr -c 'A-Za-z0-9._-' '_')"
SAFE_JOB_ID="$(printf '%s' "$DIAZODB_JOB_ID" | tr -c 'A-Za-z0-9._-' '_')"
SAFE_INPUT_STEM="${SAFE_INPUT_STEM:-input}"
HMM_PREFIX="${SAFE_INPUT_STEM}__job_${SAFE_JOB_ID}_nif"
JOB_DOMTBLOUT="$JOB_HMM_DIR/${HMM_PREFIX}.domtblout"
JOB_HMM_OUT="$JOB_HMM_DIR/${HMM_PREFIX}.out"

# Match bin/2_hmmsearch.sh: parse_hmm.py consumes HMMER domtblout files.
"$HMMSEARCH_BIN" \
  --domtblout "$JOB_DOMTBLOUT" \
  -o "$JOB_HMM_OUT" \
  "$HMM_DB" \
  "$QUERY_FASTA"

# Preserve a parser result containing only this job's HMM output.
cd "$SCRIPT_DIR"
python parse_hmm.py \
  --hits "$JOB_HMM_DIR" \
  --outdir "$JOB_PARSE_DIR/bacteria" \
  --min_genes 3 \
  --gene_range 15

# The current parser and conserved-residue code use shared repository paths.
# Serialize that section, expose only this job's domtblout there, and retain
# the job-specific parser output above.
exec 9>"$PIPELINE_LOCK"
flock 9

SHARED_DOMTBLOUT="$SHARED_HMM_DIR/${HMM_PREFIX}.domtblout"
cleanup_shared_input() {
  if [[ -L "$SHARED_DOMTBLOUT" ]]; then
    rm -f "$SHARED_DOMTBLOUT"
  fi
}
trap cleanup_shared_input EXIT
ln -sfn "$JOB_DOMTBLOUT" "$SHARED_DOMTBLOUT"

# Match bin/3_parse_hmm.sh. conserved-res.py loads both generated hit tables.
python parse_hmm.py \
  --hits "$SHARED_ARCHAEA_DIR" \
  --outdir "${SHARED_ARCHAEA_DIR%/}" \
  --min_genes 3 \
  --gene_range 15
python parse_hmm.py \
  --hits "$SHARED_HMM_DIR" \
  --outdir "${SHARED_HMM_DIR%/}" \
  --min_genes 3 \
  --gene_range 15

# Match bin/4_conserved_res.sh. This stage creates results/final/nif_clusters.csv.
python conserved-res.py --align_fasta

# Emit final result for the runner
if [[ ! -s "$SHARED_FINAL_RESULTS" ]]; then
  echo "Expected result was not created: $SHARED_FINAL_RESULTS" >&2
  exit 1
fi
cp "$SHARED_FINAL_RESULTS" "$FINAL_OUTPUT"
test -s "$FINAL_OUTPUT"
cleanup_shared_input
trap - EXIT

echo "====================================================="
echo "End Time    : $(date)"
echo "====================================================="
