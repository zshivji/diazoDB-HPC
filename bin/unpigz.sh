#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=04:10:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 8GB   # memory per CPU core
#SBATCH -J unpigz.e%j   # job name

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

# unzip bacteria files --> need to split up bc too many files
unpigz ../all_rep_proteins_aa/bacteria/GB_GCA_00*.gz

unpigz ../all_rep_proteins_aa/bacteria/GB_GCA_01*.gz
unpigz ../all_rep_proteins_aa/bacteria/GB_GCA_02*.gz
unpigz ../all_rep_proteins_aa/bacteria/GB_GCA_09*.gz
unpigz ../all_rep_proteins_aa/bacteria/GB_GCA_08*.gz
unpigz ../all_rep_proteins_aa/bacteria/GB_GCA_05*.gz
unpigz ../all_rep_proteins_aa/bacteria/GB_GCA_03*.gz

unpigz ../all_rep_proteins_aa/bacteria/GB_GCA_0*.gz

unpigz ../all_rep_proteins_aa/bacteria/GB_GCA*.gz

unpigz ../all_rep_proteins_aa/bacteria/RS_GCF_9*.gz
unpigz ../all_rep_proteins_aa/bacteria/RS_GCF_8*.gz
unpigz ../all_rep_proteins_aa/bacteria/RS_GCF_7*.gz
unpigz ../all_rep_proteins_aa/bacteria/RS_GCF_6*.gz
unpigz ../all_rep_proteins_aa/bacteria/RS_GCF_5*.gz

unpigz ../all_rep_proteins_aa/bacteria/RS_GCF_*.gz
unpigz ../all_rep_proteins_aa/bacteria/*.gz

# change all files to hpc_enviromics
chgrp -R hpc_enviromics ../all_rep_proteins_aa/bacteria/*

echo ""
echo "===================================================="
echo "End Time   : $(date)"
echo "======================================================"
echo ""
