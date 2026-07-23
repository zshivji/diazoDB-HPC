#!/bin/bash

#SBATCH --time=06:04:00   # walltime #8hrs?
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 40GB   # memory per CPU core
#SBATCH --job-name=conserved-res   # job name
#SBATCH -o logs/%x-%j.out # STDOUT
##SBATCH --dependency=afterok:48728583

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

python conserved-res.py #--align_fasta #--reload_fasta

#python final-fasta-export.py

#python diazoDB-check.py

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""

