#!/bin/bash
# Run `bash bin/8_pangenome.sh` to batch inputs and submit the array.

#SBATCH --time=36:10:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120GB
#SBATCH --job-name=pangenome
#SBATCH --output=logs/%x-%A_%a.out

set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT="${ROOT}/bin/8_pangenome.sh"
DIR="${ROOT}/pangenome"
INPUT_DIR="${DIR}/input-fastas"
BATCH_DIR="${DIR}/pangenome-batch"
OUTPUT_DIR="${DIR}/microbeannotator"
DB="/resnick/groups/enviromics/databases/microbeannotator-db"
MICROBE_ENV="/resnick/groups/enviromics/zahra/miniconda3/envs/microbeannotator"
BATCH_SIZE=50

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

prepare_batches() {
    mkdir -p "${BATCH_DIR}" "${OUTPUT_DIR}" "${ROOT}/logs"
    find "${BATCH_DIR}" -maxdepth 1 -type f -name 'batch_*' -delete
    find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.faa' | sort \
        > "${DIR}/genomes.txt"
    [[ -s "${DIR}/genomes.txt" ]] || {
        log "No pangenome FASTA files found in ${INPUT_DIR}"
        exit 1
    }
    split -l "${BATCH_SIZE}" -d -a 3 \
        "${DIR}/genomes.txt" "${BATCH_DIR}/batch_"
}

run_batch() {
    local batch output
    batch=$(printf '%s/batch_%03d' "${BATCH_DIR}" "${SLURM_ARRAY_TASK_ID}")
    output=$(printf '%s/batch_%03d' "${OUTPUT_DIR}" "${SLURM_ARRAY_TASK_ID}")
    [[ -s "${batch}" ]] || {
        log "Missing or empty batch file: ${batch}"
        exit 1
    }

    mapfile -t inputs < "${batch}"
    mkdir -p "${output}"
    eval "$(conda shell.bash hook)"
    conda activate "${MICROBE_ENV}"
    module load diamond/2.1.7-gcc-13.2.0-cfkl5pd

    log "Annotating ${#inputs[@]} pangenome FASTA files"
    microbeannotator \
        --input "${inputs[@]}" \
        --outdir "${output}" \
        --method diamond \
        --database "${DB}" \
        -p 8 \
        -t 4 \
        --refine \
        --no_plot
}

submit_workflow() {
    local batch_count array_job
    prepare_batches
    batch_count=$(find "${BATCH_DIR}" -maxdepth 1 \
        -type f -name 'batch_[0-9][0-9][0-9]' | wc -l)

    array_job=$(sbatch --parsable \
        --array="0-$((batch_count - 1))" "${SCRIPT}" run)

    log "Submitted MicrobeAnnotator array ${array_job} (${batch_count} batches)"
    log "No downstream job was submitted"
}

case "${1:-submit}" in
    submit) submit_workflow ;;
    run) run_batch ;;
    *)
        echo "Usage: $0 [submit|run]" >&2
        exit 2
        ;;
esac
