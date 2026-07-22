#!/bin/bash

#SBATCH --time=18:10:00   # walltime --> 15 hrs for both Archaea, Bacteria
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per CPU core
#SBATCH --job-name=hmmsearch   # job name
#SBATCH -o logs/%x-%j.out # STDOUT

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

#HMM search on GTDB archaea
for file in /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps_232/archaea/*; do
        f="${file##*/}"
        hmmsearch \
            --domtblout "/resnick/groups/enviromics/zahra/diazoDB-HPC/results/hmmsearch_results/archaea/${f%.faa}_nif.domtblout" \
            -o "/resnick/groups/enviromics/zahra/diazoDB-HPC/results/hmmsearch_results/archaea/${f%.faa}_nif.out" \
            /resnick/groups/enviromics/zahra/diazoDB-HPC/HMMs/combined_nif_07202026.hmm "$file"
done

chgrp hpc_enviromics /resnick/groups/enviromics/zahra/diazoDB-HPC/results/hmmsearch_results/archaea/*

echo "Archaea completed"

#HMM search on GTDB bacteria
for file in /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps_232/bacteria/*; do
        f="${file##*/}"
        hmmsearch \
            --domtblout "/resnick/groups/enviromics/zahra/diazoDB-HPC/results/hmmsearch_results/bacteria/${f%.faa}_nif.domtblout" \
            -o "/resnick/groups/enviromics/zahra/diazoDB-HPC/results/hmmsearch_results/bacteria/${f%.faa}_nif.out" \
            /resnick/groups/enviromics/zahra/diazoDB-HPC/HMMs/combined_nif_07202026.hmm "$file" 
done

for file in /resnick/groups/enviromics/zahra/diazoDB-HPC/results/bacteria/hmmsearch_results/*; do
        chgrp hpc_enviromics "$file"
done

echo "Bacteria completed"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""

