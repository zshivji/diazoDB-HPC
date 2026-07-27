#!/usr/bin/env bash

# Submit with runner/runner.py, not by hand.
# Required env vars:
#   DIAZODB_JOB_ID   - API job UUID
#   DIAZODB_INPUT    - input FASTA path in the job workspace
#   DIAZODB_OUTDIR   - output directory for intermediate files
#   DIAZODB_OUTPUT   - final CSV/HTML/PDF path to post back to the API

#SBATCH --time=02:04:00
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
JOB_HMM_DIR="$OUTDIR/hmmsearch_results/hmm_out"
JOB_PARSE_DIR="$OUTDIR/hmmsearch_results"
JOB_PARSED_HITS="$JOB_PARSE_DIR/hits.csv"
JOB_CONSERVED_DIR="$OUTDIR/conserved_res"
JOB_PROTEINS_DIR="$OUTDIR/proteins"
FINAL_DIR="$(dirname "$FINAL_OUTPUT")"
REFERENCE_PROTEINS="$REPO_ROOT/runner/references/RS_GCF_001571105.1_protein.faa"
REFERENCE_IDS="$REPO_ROOT/runner/references/reference_genome.ids"

if [[ ! -f "$INPUT_FASTA" ]]; then
  echo "Input FASTA does not exist: $INPUT_FASTA" >&2
  exit 2
fi
if [[ ! -f "$HMM_DB" ]]; then
  echo "HMM profile does not exist: $HMM_DB" >&2
  exit 2
fi
if [[ ! -f "$REFERENCE_PROTEINS" || ! -f "$REFERENCE_IDS" ]]; then
  echo "Packaged conserved-residue references are missing under runner/references." >&2
  exit 2
fi

mkdir -p \
  "$JOB_HMM_DIR" \
  "$JOB_PARSE_DIR" \
  "$JOB_CONSERVED_DIR" \
  "$JOB_PROTEINS_DIR/user" \
  "$JOB_PROTEINS_DIR/reference"

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID      : $DIAZODB_JOB_ID"
echo "Input       : $INPUT_FASTA"
echo "Output      : $FINAL_OUTPUT"
echo "Slurm ID    : ${SLURM_JOBID:-local}"
echo "====================================================="

eval "$(conda shell.bash hook)"
conda activate "$DIAZODB_CONDA_ENV"

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
GENOME_ID="$SAFE_INPUT_STEM"
if [[ "$SAFE_INPUT_STEM" =~ ([[:alnum:]_]+_GC[AF]_[0-9]+\.[0-9]+) ]]; then
  GENOME_ID="${BASH_REMATCH[1]}"
fi

# Present only the uploaded proteins and small packaged conserved-residue
# references through the directory layout understood by helper.get_seq().
ln -sfn "$QUERY_FASTA" "$JOB_PROTEINS_DIR/user/${GENOME_ID}_protein.faa"
ln -sfn \
  "$REFERENCE_PROTEINS" \
  "$JOB_PROTEINS_DIR/reference/RS_GCF_001571105.1_protein.faa"

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
  --outdir "$JOB_PARSE_DIR" \
  --output_file "$JOB_PARSED_HITS" \
  --min_genes 3 \
  --gene_range 15 \
  --skip_taxonomy

# Run conserved-residue classification entirely within this job's directories.
python conserved-res.py \
  --reload_fasta \
  --hits_file "$JOB_PARSED_HITS" \
  --results_dir "$JOB_CONSERVED_DIR" \
  --final_dir "$FINAL_DIR" \
  --proteins_dir "$JOB_PROTEINS_DIR" \
  --config_file "$SCRIPT_DIR/nif-config.json" \
  --ref_seq_ids "$REFERENCE_IDS" \
  --skip_metadata

# Emit final result for the runner
if [[ ! -s "$FINAL_OUTPUT" ]]; then
  echo "Expected result was not created: $FINAL_OUTPUT" >&2
  exit 1
fi

# to do --> add figure creation etc.

echo "====================================================="
echo "End Time    : $(date)"
echo "====================================================="
