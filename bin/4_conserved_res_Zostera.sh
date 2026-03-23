#!/bin/bash

#SBATCH --time=2:14:00   # walltime (2hrs?)
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 8GB   # memory per CPU core
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
conda activate /home/zshivji/miniconda3/envs/parse_hmm

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

#mkdir ../results/fasta_splits/
#python aln_ars_hits.py WW
#python conserved-res.py WW 09042025
#python final-fasta-export.py WW 09042025
#rm -R ../results/fasta_splits/
#echo "WW done!"


mkdir ../results/fasta_splits/
python aln_ars_hits.py GTDB
python conserved-res.py GTDB 08272025
python final-fasta-export.py GTDB 08272025
rm -R ../results/fasta_splits/
echo "GTDB done!"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""

