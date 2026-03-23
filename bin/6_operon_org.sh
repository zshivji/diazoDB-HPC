#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=1:10:00   # walltime >18 hrs
#SBATCH --ntasks=32   # number of processor cores (i.e. tasks) 32
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 40GB   # memory per node 100 GB
#SBATCH -J operon-org   # job name

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
conda activate /home/zshivji/miniconda3/envs/parse_hmm/

echo "getting operon fastas"
#python get_operon.py

echo "running microbeannotator"
#conda deactivate
conda activate /groups/enviromics/miniconda3/envs/microbeannotator

module load diamond/2.1.7-gcc-13.2.0-cfkl5pd

microbeannotator \
--input $(ls ../Zostera/operon-org/input-fastas/*11062025.fasta)  \
--outdir ../Zostera/operon-org/microbeannotator/nifEN \
--method diamond \
--database /groups/enviromics/db \
-p 8 \
-t 4  \
--refine \
--no_plot

# --continue_run

#echo "plotting"
#python operon-org-plot.py 0

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
