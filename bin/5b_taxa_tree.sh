#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=18:15:00   # walltime # about 3hrs for ~300 seqs, 20+ hrs for 7000+
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 150GB   # memory per node
#SBATCH --job-name=tree   # job name
#SBATCH -o logs/%x-%j.out # STDOUT

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

# load env, software
eval "$(conda shell.bash hook)"
CONDA_ENV="/resnick/groups/enviromics/gsolini/miniconda3/envs/gtdbtk"
conda activate "$CONDA_ENV"

BASE="/resnick/groups/enviromics/zahra/diazoDB-HPC"
DIR="${BASE}/regulon"

# download GTDB reference trees
wget -P "$DIR" https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_r232.tree
wget -P "$DIR" https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_r232.tree

# make sure GTDB taxonomy file is up to date (same release as trees)

gtdb_itol_decorate 
    "${DIR}/ar53_r232.tree" \
    ../bin/GTDB_taxonomy_r232.tsv \
    "${DIR}/ar53_r232.itol.tree" 

gtdb_itol_decorate
    "${DIR}/bac120_r232.tree" \
    ../bin/GTDB_taxonomy_r232.tsv \
    "${DIR}/bac120_r232.itol.tree"




echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
