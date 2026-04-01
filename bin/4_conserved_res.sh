#!/bin/bash

#SBATCH --time=04:04:00   # walltime #8hrs?
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 40GB   # memory per CPU core
#SBATCH -J aln_hits   # job name
##SBATCH --dependency=afterok:48728583

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=zshivji@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

## /SBATCH -p general # partition (queue)
## /SBATCH -o slurm.%N.%j.out # STDOUT
## /SBATCH -e slurm.%N.%j.err # STDERR

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

python aln_nif_hits.py #--reload_fasta

python conserved-res.py

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""

