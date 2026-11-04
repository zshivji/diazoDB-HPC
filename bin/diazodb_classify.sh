#!/usr/bin/env bash

# Submit with runner/runner.py, not by hand.
# Required env vars:
#   DIAZODB_JOB_ID   - API job UUID
#   DIAZODB_INPUT    - input FASTA path in the job workspace
#   DIAZODB_INTERMEDIATE   - output directory for intermediate files
#   DIAZODB_OUTPUT   - path to results to post back to the API

# Operon annotation is substantially slower than HMM classification.
#SBATCH --time=0:10:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -J diazodb_classify
#SBATCH -o /resnick/scratch/zshivji/diazoDB-HPC/logs/%x-%j.out

set -euo pipefail

REPO_ROOT="/resnick/groups/enviromics/zahra/diazoDB-HPC"
SCRIPT_DIR="/resnick/groups/enviromics/zahra/diazoDB-HPC/bin"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Prefer env vars injected by runner; fall back to positional args if provided.
DIAZODB_JOB_ID="${DIAZODB_JOB_ID:-${1:-}}"
INPUT_FASTA="${DIAZODB_INPUT:-${2:-}}"
INTERMEDIATE="${DIAZODB_INTERMEDIATE:-${3:-}}"
FINAL_DIR="${DIAZODB_OUTPUT:-${4:-}}"

if [[ -z "$DIAZODB_JOB_ID" || -z "$INPUT_FASTA" || -z "$INTERMEDIATE" || -z "$FINAL_DIR" ]]; then
  echo "Missing required inputs. Provide DIAZODB_JOB_ID, DIAZODB_INPUT, DIAZODB_INTERMEDIATE, DIAZODB_OUTPUT (or positional args)." >&2
  exit 2
fi

mkdir -p "$INTERMEDIATE" "$FINAL_DIR"

HMMSEARCH_BIN="${DIAZODB_HMMSEARCH_BIN:-hmmsearch}"
PRODIGAL_BIN="${DIAZODB_PRODIGAL_BIN:-prodigal}"
USE_PRODIGAL="${DIAZODB_USE_PRODIGAL:-false}"
DIAZODB_CONDA_ENV="${DIAZODB_CONDA_ENV:-/resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm}"
DIAZODB_PRODIGAL_ENV="${DIAZODB_PRODIGAL_ENV:-/resnick/groups/enviromics/zahra/miniconda3/envs/prodigal}"
CONDA_BIN="${DIAZODB_CONDA_BIN:-/resnick/groups/enviromics/zahra/miniconda3/bin/conda}"
MICROBE_ENV="${DIAZODB_MICROBE_ENV:-/resnick/groups/enviromics/zahra/miniconda3/envs/microbeannotator}"
MICROBE_DB="${DIAZODB_MICROBE_DB:-/resnick/groups/enviromics/databases/microbeannotator-db}"
MICROBE_THREADS="${DIAZODB_MICROBE_THREADS:-1}"

QUERY_FASTA="$INPUT_FASTA"
# Allow deployments to pin a profile explicitly; keep the repository's current
# profile as the default for the shared runner environment.
HMM_DB="${DIAZODB_HMM_PROFILE:-${REPO_ROOT}/HMMs/combined_nif_07292026.hmm}"
JOB_HMM_DIR="$INTERMEDIATE/hmmsearch_results"
JOB_PARSED_HITS="$JOB_HMM_DIR/hits.csv"
TMP="$INTERMEDIATE/tmp"

if [[ ! -f "$INPUT_FASTA" ]]; then
  echo "Input FASTA does not exist: $INPUT_FASTA" >&2
  exit 2
fi
if [[ ! -f "$HMM_DB" ]]; then
  echo "HMM profile does not exist: $HMM_DB" >&2
  exit 2
fi

mkdir -p \
  "$JOB_HMM_DIR" \
  "$TMP"

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID      : $DIAZODB_JOB_ID"
echo "Input       : $INPUT_FASTA"
echo "Results Dir : $FINAL_DIR"
echo "Slurm ID    : ${SLURM_JOBID:-local}"
echo "====================================================="

if [[ ! -x "$CONDA_BIN" ]]; then
  echo "Conda executable not found or not executable: $CONDA_BIN" >&2
  exit 127
fi
eval "$("$CONDA_BIN" shell.bash hook)"
conda activate "$DIAZODB_CONDA_ENV"

if [[ "$USE_PRODIGAL" == "true" ]]; then
  QUERY_FASTA="$INTERMEDIATE/predicted_proteins.faa"
  conda run --no-capture-output -p "$DIAZODB_PRODIGAL_ENV" "$PRODIGAL_BIN" \
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
# GENOME_ID="$SAFE_INPUT_STEM"
# if [[ "$SAFE_INPUT_STEM" =~ ([[:alnum:]_]+_GC[AF]_[0-9]+\.[0-9]+) ]]; then
#   GENOME_ID="${BASH_REMATCH[1]}"
# fi

# Match bin/2_hmmsearch.sh: parse_hmm.py consumes HMMER domtblout files.
"$HMMSEARCH_BIN" \
  --domtblout "$JOB_DOMTBLOUT" \
  -o "$JOB_HMM_OUT" \
  "$HMM_DB" \
  "$QUERY_FASTA"

# run jobs from diazoDB-HPC/bin to ensure relative paths to config files are correct
cd "$SCRIPT_DIR"

# Preserve a parser result containing only this job's HMM output.
python parse_hmm.py \
  --hits "$JOB_HMM_DIR" \
  --outdir "$JOB_HMM_DIR" \
  --output_file "$JOB_PARSED_HITS" \
  --min_genes 3 \
  --gene_range 15 \
  --skip_taxonomy

# Run conserved-residue classification entirely within this job's directories.
python conserved_res.py \
  --reload_fasta \
  --hits_file "$JOB_PARSED_HITS" \
  --proteins_dir "$(dirname "$QUERY_FASTA")" \
  --tmp "$TMP" \
  --final_dir "$FINAL_DIR" \
  --skip_metadata \
  --external

# Create operon organization diagrams in this job's isolated workspace.
# The classifier has already produced nif_clusters.csv and nif_final.csv in
# FINAL_DIR.  The metadata helper uses those files to select neighborhoods,
# while MicrobeAnnotator supplies annotations for surrounding genes.
OPERON_DIR="$INTERMEDIATE/operon-org"
OPERON_INPUT_DIR="$OPERON_DIR/input-fastas"
OPERON_ANNOT_DIR="$OPERON_DIR/microbeannotator"
OPERON_METADATA="$FINAL_DIR/operon_metadata.json"
OPERON_CLUSTERS="$FINAL_DIR/nif_clusters.csv"
OPERON_NIF_FINAL="$FINAL_DIR/nif_final.csv"

log "Preparing operon FASTA inputs in $OPERON_INPUT_DIR"
mkdir -p "$OPERON_INPUT_DIR" "$OPERON_ANNOT_DIR"
if [[ -s "$OPERON_CLUSTERS" && -s "$OPERON_NIF_FINAL" ]]; then
  (
    cd "$SCRIPT_DIR"
    conda run -p "$DIAZODB_CONDA_ENV" python diazoDB-metadata.py \
      --prepare \
      --proteins_dir "$(dirname "$QUERY_FASTA")" \
      --clusters_file "$OPERON_CLUSTERS" \
      --operon_dir "$OPERON_DIR"
  )

  mapfile -t OPERON_INPUTS < <(find "$OPERON_INPUT_DIR" -maxdepth 1 -type f -name '*.fasta' | sort)
  if ((${#OPERON_INPUTS[@]})); then
    if command -v module >/dev/null 2>&1; then
      module load diamond/2.1.7-gcc-13.2.0-cfkl5pd
    fi
    conda activate "$MICROBE_ENV"
    if [[ ! -d "$MICROBE_DB" ]]; then
      echo "MicrobeAnnotator database does not exist: $MICROBE_DB" >&2
      exit 2
    fi

    log "Annotating ${#OPERON_INPUTS[@]} operon neighborhoods"
    microbeannotator \
      --input "${OPERON_INPUTS[@]}" \
      --outdir "$OPERON_ANNOT_DIR" \
      --method diamond \
      --database "$MICROBE_DB" \
      -t "$MICROBE_THREADS" \
      --refine \
      --light \
      --no_plot

    (
      cd "$SCRIPT_DIR"
      conda run -p "$DIAZODB_CONDA_ENV" python diazoDB-metadata.py \
        --data \
        --clusters_file "$OPERON_CLUSTERS" \
        --nif_final_file "$OPERON_NIF_FINAL" \
        --operon_dir "$OPERON_DIR" \
        --metadata_file "$OPERON_METADATA"
    )
  else
    log "No candidate operons found; skipping operon diagram"
  fi
else
  log "No nif result tables found; skipping operon diagram"
fi

# show identified nif genes on nifH reference tree
#module load blast/2.15.0
#makeblastdb 
#blastp

# Emit final result for the runner
FINAL_CSV="$FINAL_DIR/nif_clusters.csv"
FINAL_DETAIL_CSV="$FINAL_DIR/nif_final.csv"
if [[ ! -s "$FINAL_CSV" ]]; then
  echo "Expected result was not created: $FINAL_CSV" >&2
  exit 1
fi
if [[ ! -s "$FINAL_DETAIL_CSV" ]]; then
  echo "Expected result was not created: $FINAL_DETAIL_CSV" >&2
  exit 1
fi

echo "====================================================="
echo "End Time    : $(date)"
echo "====================================================="
