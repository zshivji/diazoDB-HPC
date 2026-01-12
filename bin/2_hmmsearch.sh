#!/bin/bash

#SBATCH --time=18:10:00   # walltime --> 15 hrs for both Archaea, Bacteria
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per CPU core
#SBATCH -J hmmsearch.e%j   # job name

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=zshivji@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

## /SBATCH -p general # partition (queue)
## /SBATCH -o slurm.%N.%j.out # STDOUT
## /SBATCH -e slurm.%N.%j.err # STDERR

eval "$(conda shell.bash hook)"
conda activate /central/groups/enviromics/miniconda3/envs/hmmer

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

#HMM search on GTDB archaea
#for file in ../all_rep_proteins_aa/archaea/*; do
#        f="${file##*/}"
#        hmmsearch ../HMMs/nifV_07242025.hmm "$file" >> "../results/archaea/hmmsearch_results/${f%.faa}_nif.out"
#done

#chgrp hpc_enviromics ../results/archaea/hmmsearch_results/*

#echo "Archaea completed"

#HMM search on GTDB bacteria
#for file in ../all_rep_proteins_aa/bacteria/*; do
#        f="${file##*/}"
#        hmmsearch ../HMMs/nifV_07242025.hmm "$file" >> "../results/bacteria/hmmsearch_results/${f%.faa}_nif.out"
#done

#for file in ../results/bacteria/hmmsearch_results/*; do
#        chgrp hpc_enviromics ../results/bacteria/hmmsearch_results/"$file"
#done

#echo "Bacteria completed"

for file in ../Zostera_marina_isolates/*.faa; do # takes 2 min
        f="${file##*/}"
        hmmsearch ../HMMs/combined_nif_04012025.hmm "$file" >> "../results/Zostera/hmmsearch_results/${f%.faa}_nif.out"
done

chgrp hpc_enviromics ../results/Zostera/hmmsearch_results/*

echo "Zostera completed"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""

