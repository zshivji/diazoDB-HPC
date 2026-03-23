#!/bin/bash

#SBATCH --time=2:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per CPU core
#SBATCH -J parse_hmm.e%j   # job name

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=zshivji@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

## /SBATCH -p general # partition (queue)
## /SBATCH -o slurm.%N.%j.out # STDOUT
## /SBATCH -e slurm.%N.%j.err # STDERR

eval "$(conda shell.bash hook)"
conda activate /home/zshivji/miniconda3/envs/parse_hmm

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

python Parse_hmm_results.py archaea

echo "Found hits!"

python Parse_tophits.py archaea

echo "Found nif!"

echo "Archaea done!"

python Parse_hmm_results.py bacteria

echo "Found hits!"

python Parse_tophits.py bacteria

echo "Found nif!"

echo "Bacteria done!"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""

