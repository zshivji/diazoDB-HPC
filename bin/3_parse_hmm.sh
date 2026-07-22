#!/bin/bash

#SBATCH --time=2:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per CPU core
#SBATCH --job-name=parse-hmm   # job name
#SBATCH -o logs/%x-%j.out # STDOUT

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm/

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

python parse_hmm.py \
	--hits /resnick/groups/enviromics/zahra/diazoDB-HPC/results/hmmsearch_results/archaea/ \
	--outdir /resnick/groups/enviromics/zahra/diazoDB-HPC/results/hmmsearch_results/archaea \
    --min_genes 3 \
    --gene_range 15

echo "Archaea hits found!"

python parse_hmm.py \
	--hits /resnick/groups/enviromics/zahra/diazoDB-HPC/results/hmmsearch_results/bacteria/ \
	--outdir /resnick/groups/enviromics/zahra/diazoDB-HPC/results/hmmsearch_results/bacteria \
    --min_genes 3 \
    --gene_range 15

echo "Bacteria hits found!"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""

