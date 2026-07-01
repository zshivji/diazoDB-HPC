#!/bin/bash
#===============================================================================
# SLURM Job: Diazotroph pangenome
#===============================================================================
# Usage:  sbatch gtdbtk.sh
#
# Description:
#   Runs microbeannotator on diazotrophs (GTDB genomes)
#
# Prerequisites:
#   - List of genomes + their protein files
#   - Microbeannotator env with database configured
#
# Inputs:
#   - genome list from /resnick/groups/enviromics/zahra/diazoDB-HPC/results/final/nif_genomes.csv
#   - GTDB genomes from /resnick/groups/enviromics/zahra/diazoDB-HPC/proteins_faa_reps
#
# Outputs:
#   - DIR/ : full microbeannotator output (annotation results, etc.)
#===============================================================================

#--- SLURM directives ---------------------------------------------------------
#SBATCH --time=36:10:00   # walltime --> 24 hrs for both Archaea, Bacteria for de novo
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=32
#SBATCH --mem 120GB   # memory per CPU core
#SBATCH --job-name=pangenome   # job name
#SBATCH --array=0-103

# Notify at the beginning, end of job and on failure.
##SBATCH --mail-user=zshivji@caltech.edu   # email address
##SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH -o logs/%x-%j.out # STDOUT

#--- Safety flags -------------------------------------------------------------
#set -uo pipefail

#--- Configuration ------------------------------------------------------------
THREADS=${SLURM_CPUS_PER_TASK:-32}
CONDA_ENV="/resnick/groups/enviromics/zahra/miniconda3/envs/microbeannotator"

# ===== EDIT THESE PATHS =====
DIR="/resnick/groups/enviromics/zahra/diazoDB-HPC/pangenome"
DB="/resnick/groups/enviromics/databases/microbeannotator-db"

#--- Functions ----------------------------------------------------------------

log() {
echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

activate_env() {
    local env_name="$1"
    eval "$(conda shell.bash hook)"
    conda activate "$env_name"
    log "Activated conda env: ${env_name}"
}

load_files() {
    local GTDB_dir="$1"
    mkdir -p "${DIR}"

    find "$GTDB_dir" -type f > all_files.txt
    : > missing_files.txt

    while IFS=, read -r genome _
    do
        file=$(grep -E "/${genome}.*\.faa$" all_files.txt | head -n1)
        if [[ -n "$file" ]]; then
            ln -sf "$file" "${DIR}/input-fastas"
            echo "$genome"
        else
            echo "$genome" >> missing_files.txt
        fi
    done < ../results/final/nif_genomes.csv

    rm -f all_files.txt
}

batch_files() {
    find "${DIR}/input-fastas" -name "*.faa" > "${DIR}/genomes.txt"
    split -l 50 -d -a 3 "${DIR}/genomes.txt" batch_
}

run_microbeannotator() {
    mkdir -p "${DIR}/microbeannotator"
    module load diamond/2.1.7-gcc-13.2.0-cfkl5pd

    batch=$(printf "batch_%03d" "$SLURM_ARRAY_TASK_ID")
    inputs=$(tr '\n' ' ' < "$batch")

    microbeannotator \
    --input $inputs  \
    --outdir "${DIR}/batch_${SLURM_ARRAY_TASK_ID}" \
    --method diamond \
    --database "${DB}" \
    -p 8 \
    -t 4  \
    --refine \
    --no_plot

}

#--- Main workflow ------------------------------------------------------------
# should be split into file which runs load_files, batch_files then call to slurm array job for microbeannotator then back to main job for merging files

main() {
    echo "======== Microbeannotator for Pangenome Analysis started ======"
    echo "Start Time  : $(date)"
    echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
    echo "Node: $(hostname), CPUs: ${THREADS}, Mem: ${SLURM_MEM_PER_NODE:-unknown}MB"
    echo "==============================================================="
    echo ""

    activate_env "$CONDA_ENV"

    # pull all nif genomes in to a flat dir (still_missing.txt contains genomes no longer in either GTDB folder)
#    log "Linking input genomes to a single directory"
#    load_files /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps/

    # batch files in first array job
    if [[ "$SLURM_ARRAY_TASK_ID" -eq 0 ]]; then
        log "Batching input files"
        batch_files
    fi

    # after batch is done, run microbeannotator
    while [[ ! -f batch_000 ]]; do
        sleep 5
    done
   
    log "Running Microbeannotator"
    run_microbeannotator


#    log "========== GTDB-Tk classification finished =========="
#    log "Results: ${GTDBTK_OUTDIR}/de_novo_wf"
#    log "  Bacterial summary: ${GTDBTK_OUTDIR}/de_novo_wf/bacteria/gtdbtk.bac120.summary.tsv"
#    log "  Archaeal summary:  ${GTDBTK_OUTDIR}/de_novo_wf/archaea/gtdbtk.ar53.summary.tsv"
}

main


