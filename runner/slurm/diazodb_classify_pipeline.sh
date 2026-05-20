#!/usr/bin/env bash

# Single-submission DiazoDB pipeline wrapper.
# This intentionally does not modify the older database-build scripts in bin/.
# Replace the placeholder commands below with the final per-upload classifier.

set -euo pipefail

INPUT_FASTA="${1:?input FASTA is required}"
OUTDIR="${2:?output directory is required}"
FINAL_OUTPUT="${3:?final output path is required}"

mkdir -p "$OUTDIR"

HMM_PROFILE="${DIAZODB_HMM_PROFILE:-}"
HMMSEARCH_BIN="${DIAZODB_HMMSEARCH_BIN:-hmmsearch}"
PRODIGAL_BIN="${DIAZODB_PRODIGAL_BIN:-prodigal}"
USE_PRODIGAL="${DIAZODB_USE_PRODIGAL:-false}"

QUERY_FASTA="$INPUT_FASTA"

if [[ "$USE_PRODIGAL" == "true" ]]; then
  QUERY_FASTA="$OUTDIR/predicted_proteins.faa"
  "$PRODIGAL_BIN" \
    -i "$INPUT_FASTA" \
    -a "$QUERY_FASTA" \
    -p meta \
    -q
fi

if [[ -n "$HMM_PROFILE" ]]; then
  "$HMMSEARCH_BIN" \
    --tblout "$OUTDIR/hmmsearch.tbl" \
    "$HMM_PROFILE" \
    "$QUERY_FASTA" \
    > "$OUTDIR/hmmsearch.out"
else
  echo "DIAZODB_HMM_PROFILE is not set; writing placeholder HMM outputs." >&2
  : > "$OUTDIR/hmmsearch.tbl"
  : > "$OUTDIR/hmmsearch.out"
fi

if [[ -n "${DIAZODB_PARSE_COMMAND:-}" ]]; then
  # Example:
  # DIAZODB_PARSE_COMMAND='python /path/to/parse_uploaded_job.py --input "$QUERY_FASTA" --tblout "$OUTDIR/hmmsearch.tbl" --output "$FINAL_OUTPUT"'
  eval "$DIAZODB_PARSE_COMMAND"
else
  {
    echo "sequence_id,prediction,confidence,source"
    awk '
      /^>/ {
        id=$1
        sub(/^>/, "", id)
        print id ",pending_parser,,placeholder"
      }
    ' "$QUERY_FASTA"
  } > "$FINAL_OUTPUT"
fi
