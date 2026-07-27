#!/bin/bash
# Run `bash bin/6_operon_org.sh` to batch inputs and submit the workflow.

#SBATCH --time=3:10:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120GB
#SBATCH --job-name=operon-org
#SBATCH -o logs/%x-%j.out # STDOUT

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

set -euo pipefail

ROOT="/resnick/groups/enviromics/zahra/diazoDB-HPC"
SCRIPT="${ROOT}/bin/6_operon_org.sh"
DIR="${ROOT}/operon-org"
INPUT_DIR="${DIR}/input-fastas"
BATCH_DIR="${DIR}/microbeannotator-batches"
OUTPUT_DIR="${DIR}/microbeannotator"
DB="/resnick/groups/enviromics/databases/microbeannotator-db"
MICROBE_ENV="/resnick/groups/enviromics/zahra/miniconda3/envs/microbeannotator"
PARSE_ENV="/resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm"
BATCH_SIZE=50

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

prepare_batches() {
    log "Creating operon FASTA inputs"
    mkdir -p "$BATCH_DIR"
#    (
#        cd "${ROOT}/bin"
#        conda run -p "${PARSE_ENV}" python diazoDB-metadata.py --prepare
#    )
    find "${BATCH_DIR}" -maxdepth 1 -type f -name 'batch_*' -delete
    find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.fasta' | sort \
        > "${BATCH_DIR}/inputs.txt"
    [[ -s "${BATCH_DIR}/inputs.txt" ]] || {
        log "No operon FASTA files found in ${INPUT_DIR}"
        exit 1
    }
    split -l "${BATCH_SIZE}" -d -a 3 \
        "${BATCH_DIR}/inputs.txt" "${BATCH_DIR}/batch_"
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
    eval "$(conda shell.bash hook)"
    conda activate "${MICROBE_ENV}"
    module load diamond/2.1.7-gcc-13.2.0-cfkl5pd

    log "Annotating ${#inputs[@]} operon FASTA files"
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

postprocess() {
    local annotation annotation_dir="${OUTPUT_DIR}/annotation_results"
    find "${annotation_dir}" -maxdepth 1 -type l -delete

    while IFS= read -r annotation; do
        ln -sf "${annotation}" "${annotation_dir}/$(basename "${annotation}")"
    done < <(
        find "${OUTPUT_DIR}"/batch_* \
            -type f -path '*/annotation_results/*.annot' | sort
    )

    log "Calling get_plot_data() after all array tasks succeeded"
    (
        cd "${ROOT}/bin"
        conda run -p "${PARSE_ENV}" python diazoDB-metadata.py --data
    )
}

submit_workflow() {
    local batch_count array_job post_job
    prepare_batches
    batch_count=$(find "${BATCH_DIR}" -maxdepth 1 \
        -type f -name 'batch_[0-9][0-9][0-9]' | wc -l)

    array_job=$(sbatch --parsable \
        --output="logs/%x-%A_%a.out" \
        --array="0-$((batch_count - 1))" \
        "${SCRIPT}" run)
    post_job=$(sbatch --parsable \
        --output="logs/%x-%j.out" \
        --dependency="afterok:${array_job}" \
        "${SCRIPT}" postprocess)

    log "Submitted MicrobeAnnotator array ${array_job} (${batch_count} batches)"
    log "Submitted get_plot_data job ${post_job} after ${array_job}"
}

case "${1:-submit}" in
    submit) submit_workflow ;;
    run) run_batch ;;
    postprocess) postprocess ;;
    *)
        echo "Usage: $0 [submit|run|postprocess]" >&2
        exit 2
        ;;
esac

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
