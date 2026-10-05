#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=40GB
#SBATCH --job-name=clusterDB
#SBATCH -o logs/%x-%j.out

set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm

ROOT="/resnick/groups/enviromics/zahra/diazoDB-HPC"
DIR="${ROOT}/diazoDB-comparison"

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 <file-gene-table.csv> <identity-threshold>" >&2
    echo "CSV must have headers: file,gene" >&2
    exit 2
fi

TABLE="$1"
THRESHOLD="$2"

if [[ "${TABLE}" = /* ]]; then
    TABLE_PATH="${TABLE}"
elif [[ -e "${TABLE}" ]]; then
    TABLE_PATH="${TABLE}"
elif [[ -e "${ROOT}/${TABLE}" ]]; then
    TABLE_PATH="${ROOT}/${TABLE}"
elif [[ -e "${DIR}/${TABLE}" ]]; then
    TABLE_PATH="${DIR}/${TABLE}"
else
    echo "CSV table not found: ${TABLE}" >&2
    exit 1
fi

if ${THRESHOLD} == "1"; then
    COV_MODE=1
    C=1
else
    COV_MODE=0
    C=0.8
fi

OUTPUT_DIR="${DIR}/clusterDB_${THRESHOLD}"
mkdir -p "${OUTPUT_DIR}" "${OUTPUT_DIR}/tmp"

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : ${SLURM_JOBID:-NA} / ${SLURM_JOB_NAME:-clusterDB}"
echo "======================================================"
echo ""

for gene in $(tail -n +2 "${TABLE_PATH}" | cut -d, -f2 | sort -u); do
    input="${OUTPUT_DIR}/${gene}_allDB.fasta"
    cluster="${OUTPUT_DIR}/${gene}_allDB_clustered"

    echo "Combining ${gene} sequences"
    : > "${input}"
    while IFS=, read -r file file_gene; do
        [[ "${file_gene}" == "${gene}" ]] || continue

        file="${DIR}/${file}"
        database="${file#${DIR}/}"
        database="${database%%/*}"
        sed "s/^>/>${database}_/" "${file}" >> "${input}"
    done < <(tail -n +2 "${TABLE_PATH}")

    echo "Clustering ${gene} at ${THRESHOLD} identity"
    mmseqs easy-cluster "${input}" "${cluster}" "${OUTPUT_DIR}/tmp/${gene}" \
        --min-seq-id "${THRESHOLD}" -c "${C}" --cov-mode "${COV_MODE}" --threads 4
done

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
