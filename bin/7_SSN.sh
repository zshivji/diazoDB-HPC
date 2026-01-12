#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=3:10:00   # walltime >6 hrs
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks) 32
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 50GB   # memory per node 150 GB
#SBATCH -J SSN   # job name

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=zshivji@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

## /SBATCH -p general # partition (queue)
## /SBATCH -o slurm.%N.%j.out # STDOUT
## /SBATCH -e slurm.%N.%j.err # STDERR

echo "====================================================="
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

# load env, software
eval "$(conda shell.bash hook)"
conda activate /home/zshivji/miniconda3/envs/make_trees/

module load mafft/7.505-gcc-13.2.0-nklkvtc

# blastn (all-against-all)
mmseqs createdb ../SSN/PF00148-unreviewed.fasta ../SSN/queryDB
mmseqs createdb ../SSN/PF00148-unreviewed.fasta ../SSN/targetDB
mmseqs createindex ../SSN/targetDB ../SSN/tmp
mmseqs search ../SSN/queryDB ../SSN/targetDB ../SSN/resultDB tmp
mmseqs convertalis ../SSN/queryDB ../SSN/targetDB ../SSN/resultDB ../SSN/resultDB.m8

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
